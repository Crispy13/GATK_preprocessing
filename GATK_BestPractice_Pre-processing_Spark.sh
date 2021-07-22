#! /bin/bash

# Help message
usage='
Necessary arguments:
-f1	fastq 1
-f2	fastq 2
-r      reference file [Default="/data/eck/Data/GATK/bundle/hg19/ucsc.hg19.fasta"]
-t      number of threads [Default=32]
-m      Including process that makes reference index for bwa. (true/false) [Default=false]\n'

# Get Arguments
while :;
do
case "$1" in
	--help|-h)
		echo -e "${usage}"
		echo -e "\nThis script must be executed at the directory where the fastq files are."
		exit 0
		;;
	-f1)
		if [ "$2" ]; then
		f1=$2
		shift
		fi
		;;
	-f2)
		if [ "$2" ]; then
		f2=$2
		shift
		fi
		;;
	-?*)
		echo "Unknown options."
		;;
	*)
		break
	esac

	shift
done

# Check Arguments.
mi_s="false"
t=32
r="/data/eck/Data/GATK/bundle/hg19/ucsc.hg19.fasta"

echo -e "${f1}\n${f2}"

if [ -z ${f1} ] || [ -z ${f2} ] || [ -z "$r" ] || [ -z "$mi_s" ] || [ -z "$t" ]
then
        echo "Required arguments are not given."
        exit 1
fi

if ! { [ "$mi_s" = "true" ] || [ "$mi_s" = "false" ] ; } ; then
        echo "-m option should receive only 'true' or 'false'."
        exit 1
fi

# Activate ES Virtual Environment
# source /data/eck/software/anaconda3/bin/activate ES

# Pre-process for separating fastq files.
echo -e "Preprocess for separting fastq files. $(date +%y%m%d_%H%M)"
f=`echo ${f1} | sed -r "s/_1\.fastq/rgfasta/"`

f1_rg_ar=() # fastq1 read group array
f2_rg_ar=() # fastq2 read group array

if [ -d "./${f}.temp" ] # Make rgfasta.temp directory
then
	mv ./${f}.temp ./${f}.temp_backup_on_$(date +%y%m%d)
	mkdir ${f}.temp
else
	mkdir ${f}.temp
fi

for i in `echo -e "${f1}\n${f2}"`
do
	# Grep Headers.
        echo "awk 'NR%4==1 {print}' ${i} > ${i}.header"
        awk 'NR%4==1 {print}' ${i} > ${i}.header &
done
wait

for i in `echo -e "${f1}\n${f2}"`
do
	# Readgroup details
	echo 'cut -f 1,2,3,4 -d ":" ${i}.header | uniq > ${i}.rgd_uniq_list'
	cut -f 1,2,3,4 -d ":" ${i}.header | uniq > ${i}.rgd_uniq_list &
done
wait

for i in `echo -e "${f1}\n${f2}"`
do
##	((f_index+=1))

	# Separate fastq file.
	echo "Separating the fastq file."
	
	for i2 in `cat ${i}.rgd_uniq_list`
	do
		f_index=`echo $i | sed -r "s/.*_//" | sed -r "s/\.fastq//"`
		rg=`echo $i2 | sed -r "s/:/_/g" | sed -r "s/(.*)/\1_${f_index}/" | sed -r "s/^@//"`

		if [ ${f_index} -eq 1 ]; then # append rg var to the array.
			f1_rg_ar+=(${rg})
		elif [ ${f_index} -eq 2 ]; then
			f2_rg_ar+=(${rg})
		fi

		grep -A3 ${i2} ${i} > ${f}.temp/${rg}.fastq &
	done

	rm ${i}.rgd_uniq_list
	rm ${i}.header

	echo "Separating complete."
	echo "${i} was done."
done
wait

echo -e "f1_rg_ar=${f1_rg_ar[@]}\nf2_rg_ar=${f2_rg_ar[@]}"

# Make uBam file.
readarray -t f1_rg_ars < <(for a in "${f1_rg_ar[@]}"; do echo "$a"; done | sort)
readarray -t f2_rg_ars < <(for a in "${f2_rg_ar[@]}"; do echo "$a"; done | sort)
f1_rg_ars=("${f1_rg_ars[@]/%_1/}")
f2_rg_ars=("${f2_rg_ars[@]/%_2/}")

echo "${f1_rg_ars[@]}"
echo "${f2_rg_ars[@]}"

A="${f1_rg_ars[@]}" # To fix "Too many arguments."
B="${f2_rg_ars[@]}"

if [ "$A" = "$B" ] # Check the number of array elements.
	then
		:
	else
		echo "The arrays are different!"
		exit 1
fi	

SM=`echo ${f1/#\.\//} | sed -r "s/(.*)(_[12].fastq)/\1/"`
cd ./${f}.temp

for ((i=0; i<${#f1_rg_ars[@]}; i++))
do
	echo -e "\n> Make a uBam file for readgroup. $(date +%y%m%d_%H%M)\n"
	gatk FastqToSam -F1 "${f1_rg_ars[i]}_1.fastq" -F2 "${f2_rg_ars[i]}_2.fastq" -O "${f1_rg_ars[i]}.bam" -RG "${f1_rg_ars[i]}" -SM ${SM} -LB ${SM} -PL "illumina" &
done
wait

# MarkIlluminaAdapters.
if [ -d "./MIA_metrics" ]; then # Make a folder for metrics files.
        :
        else
                mkdir ./MIA_metrics
fi

for i in `ls *.bam`
do
        echo -e "\n> Marking Illumina Adapters in $i file... $(date +%y%m%d_%H%M)\n"
        gatk MarkIlluminaAdapters -I $i -O ${i/\.bam/_MIA.bam} -M "./MIA_metrics/${i/\.bam/_MIA_metrics.txt}" &
done
wait

# SamToFastq
for i in `find -maxdepth 1 -name "*_MIA.bam"`
do
        echo -e "\n> Starts SamToFastq process using $i file. $(date +%y%m%d_%H%M)\n"
        gatk SamToFastq -I $i -F ${i/\.bam/_SamToFastq.fastq} -CLIP_ATTR XT -CLIP_ACT 2 -INTER true -NON_PF true &
done
wait

# Make reference index files.
if [ "$mi_s" = "true" ]; then
        bwa index -a bwtsw $r
fi

# Align sequences to reference.
for i in `find -maxdepth 1 -name "*_SamToFastq.fastq"`
do
        echo -e "\n> Aligning $i to referece... $(date +%y%m%d_%H%M)\n"
        bwa mem -M -t $t -p $r $i > ${i/.fastq/_aligned.sam}
done

# MergeBamAlignment
for i in `find -maxdepth 1 -name "*_aligned.sam"`
do
        echo -e "\n> Do MergeBamAlignment process for $i file. $(date +%y%m%d_%H%M)\n"
        gatk MergeBamAlignment -R $r --UNMAPPED_BAM ${i/_MIA_SamToFastq_aligned.sam/.bam} --ALIGNED_BAM $i -O ${i/.sam/_merged.bam} --CREATE_INDEX true --ADD_MATE_CIGAR true --CLIP_ADAPTERS false --CLIP_OVERLAPPING_READS true --INCLUDE_SECONDARY_ALIGNMENTS true --MAX_INSERTIONS_OR_DELETIONS -1 --PRIMARY_ALIGNMENT_STRATEGY MostDistant --ATTRIBUTES_TO_RETAIN XS &
done
wait

# SortSam (for Spark only)
for i in $(find -maxdepth 1 -name "*_merged.bam")
do
	echo -e "\n> Do SortSam process for $i file. $(date +%y%m%d_%H%M)\n"
	gatk SortSamSpark --spark-master local[$t] -I $i -O ${i/.bam/_sorted.bam} --sort-order queryname
done

# MarkDuplicates
MD_input=$(find -maxdepth 1 -name "*_merged_sorted.bam" | sed -r "s/(\.\/)/--input \1/g" | sed -r "s/\.\///")
MD_output=$(echo ${MD_input/*.bam[[:space:]]/} | sed -r "s/(.*)(_MIA.*)(\.bam)/${SM}\2_MD\3/") # parameter expansion matches maximally

if [ -d "./MD_metrics" ]; then # Make a folder for MarkDuplicates metrics files.
        :
        else
                mkdir ./MD_metrics
fi

echo -e "\n> Starts MarkDuplicates process.\n MD_input=\n${MD_input}\nMD_output=\n${MD_output} $(date +%y%m%d_%H%M)\n"
gatk MarkDuplicatesSpark --spark-master local[$t] ${MD_input} --output ${MD_output} --metrics-file "./MD_metrics/${MD_output/.bam/_MD_metrics.txt}" --optical-duplicate-pixel-distance 2500 --create-output-bam-index true # --conf 'spark.local.dir=tmp'

# BaseRecalibrator
echo -e "\n> Starts to analyze patterns of covariation in the sequence dataset for ${MD_output} file. $(date +%y%m%d_%H%M)\n"
gatk BaseRecalibratorSpark --spark-master local[$t] -R $r -I ${MD_output} --known-sites /data/eck/Data/GATK/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf --known-sites /data/eck/Data/GATK/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf --known-sites /data/eck/Data/GATK/bundle/hg19/dbsnp_138.hg19.vcf -O ${MD_output/.bam/_recal_data.table}

echo -e "\n> ApplyBQSR for ${MD_output} file. $(date +%y%m%d_%H%M)\n"
gatk ApplyBQSRSpark --spark-master local[$t] -R $r -I ${MD_output} -bqsr ${MD_output/.bam/_recal_data.table} -O ${MD_output/.bam/_recal_reads.bam}

echo -e "\n> Makes post_recal_data.table for $MD_output/.bam/_recal_reads.bam} file. $(date +%y%m%d_%H%M)\n"
gatk BaseRecalibratorSpark --spark-master local[$t] -R $r -I ${MD_output/.bam/_recal_reads.bam} --known-sites /data/eck/Data/GATK/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf --known-sites /data/eck/Data/GATK/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf --known-sites /data/eck/Data/GATK/bundle/hg19/dbsnp_138.hg19.vcf -O ${MD_output/.bam/_post_recal_data.table}

echo -e "\n> Starts to generate before/after plots for ${MD_output} file. $(date +%y%m%d_%H%M)\n"
gatk AnalyzeCovariates -before ${MD_output/.bam/_recal_data.table} -after ${MD_output/.bam/_post_recal_data.table} -plots ${MD_output/.bam/_AnalyzeCovariates.pdf}

# Leave the files needed at the next step except the others.
echo -e "\n`find ! -name "*_recal_reads.bam" \( -name "*.bam" -o -name "*sam" -o -name "*fastq" \)`\nwill be removed.\n"
find ! -name "*_recal_reads.bam" \( -name "*.bam" -o -name "*sam" -o -name "*fastq" \) -exec rm {} \;

echo "End of $0. $(date +%y%m%d_%H%M)"
