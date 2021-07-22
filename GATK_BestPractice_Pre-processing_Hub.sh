#! /bin/bash

# list to do
: <<'END'
1. Input folder as full path.
END


# Help message
usage="
Necessary arguments:
-fl	input folder list
"
caveat="
** You should notice the followings:

- This script must be executed at the parent directory of the subs that include fastq files.
- If *rgfasta.temp folder exists, it will be removed and be made newly.
"

while :;
do
case "$1" in
	--help|-h)
		echo -e "${caveat}"
                echo -e "${usage}"
                exit 0
                ;;
        -fl)
                if [ "$2" ]; then
                fl=$2
                shift
                fi
                ;;
	--wes)
		if [ "$2" ]; then
                wes=$2
                shift
                fi
                ;;
	-s)
                s=1
                ;;
        -?*)
                echo "Unknown options."
                ;;
	*)
                break
esac
shift
done

# Default Arguments
if [ -z "${fl}" ]; then
fl=$(find -maxdepth 1 -type d | tr "\n" " " | sed -r "s/\.[[:space:]]//")
fi

if [ -z "${wes}" ]; then
wes="/data/eck/Workspace/SMK_Colon_WES"
fi

# Convert relative path to absolute path
for i in ${fl};
do
        i=$(echo $i | sed -r "s/\.\///")
        fl_arr+=($(find "${wes}" -maxdepth 3 -name "${i}"))
done

# Caveat
echo -e "${caveat}\nStarting in 10 secs...\n"
echo -e "fl=$fl"
sleep 10

:<<'END'
Are you sure you want to launch this script? [yes|no]
"

echo "${caveat}"
read -r ans

while [ "$ans" != "yes" ] && [ "$ans" != "Yes" ] && [ "$ans" != "YES" ] && \
          [ "$ans" != "no" ]  && [ "$ans" != "No" ]  && [ "$ans" != "NO" ]
do
        printf "Please answer 'yes' or 'no':'\\n"
        printf ">>> "
        read -r ans
done

if [ "$ans" != "yes" ] && [ "$ans" != "Yes" ] && [ "$ans" != "YES" ]
        then
                "Abort running."
                exit 2
fi
END

# Activate ES venv
source /data/eck/software/anaconda3/bin/activate ES

# Parameters
cs="/data/eck/CustomScriptsArchive"

# notice Spark
if [ "$s" == "1" ]; then
echo -e "\nSpark will be used."; else
echo -e "\nSpark will not be used."
fi

# Feed f1, f2 for the main script.
for i in ${fl};
do
	cd $i
	echo -e "\n>> Current working folder:\n$(pwd)\n"
	f1=`find -maxdepth 1 -name "*_1.fastq"`
	f2=`find -maxdepth 1 -name "*_2.fastq"`
#	echo -e "${f1}\n${f2}"
	if [ "$s" == "1" ]
	then
		${cs}/GATK_BestPractice_Pre-processing_Spark.sh -f1 ${f1} -f2 ${f2};
	else
		${cs}/GATK_BestPractice_Pre-processing.sh -f1 ${f1} -f2 ${f2}
	fi
	cd ../
done
