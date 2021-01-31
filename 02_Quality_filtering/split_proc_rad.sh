#!/bin/bash
# define variables from input
INSEQ=$1
PNO=$2

# check input and print usage
if [[ -z $INSEQ || -z $PNO ]]
	then
	echo "Usage: script_name input.fastq cpu_no"; 
	echo "Make sure your adaptors file is present in the same directory as reads.";
	exit
fi

# preset variables and thresholds. modify if you understand what youre changing
TSTART=1	# beginning of region to include during truncation
TEND=36		# end of region to include during truncation
QT=20		# quality threshold
NLQ=18		# number of low quality bases allowed
MLEN=18		# minimum length of adaptor match (cross_match)
MSCORE=12	# minimum score of adaptor match (cross_match)

BS=3000000

NOLINES=`wc -l $INSEQ | perl -pi -e "s/\s.+//"`
NOSEQS=`expr $NOLINES / 4`
split -d -a 3 -l $BS $INSEQ
ls x* >list
perl -pi -e "s/x/\nx/g" list
perl -pi -e "s/^\n//g" list
cat list

exec < "list"
while read LINE
do
	mv $LINE $LINE.fastq
        mkdir $LINE
        mv $LINE.fastq ./$LINE
	ls $LINE
done

exec < "list"
j=0
while read LINE
do
	if [ $j == $PNO ] 
	then
	j=0
	fi
	j=$[ j + 1 ]
	echo "$LINE  $j"
	cd $LINE
	../TruncateFastq.pl x*.fastq $TSTART $TEND trunc.fastq ;
	../QualFilterFastq.pl trunc.fastq $QT $NLQ hq.fastq ; 
	sh ~/Programs/bbmap/bbduk.sh in=hq.fastq ref=../adaptors.fasta k=$MSCORE stats=stats.txt out=clean.fastq ;
	cd ..
	if [ $j == $PNO ]
	then
	  wait
	fi
done
wait;

cat x*/clean.fastq > clean.fastq;

exec < "list"
while read LINE
do
	rm -rf $LINE
done

#mv clean.fastq $file_clean.fastq
