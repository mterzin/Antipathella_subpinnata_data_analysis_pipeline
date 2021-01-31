#!/bin/bash
top=$(readlink -f $(dirname $0)/..)

# STEP 16: Go back to the top directory.
cd $top

# STEP 17-B-i: Align every sample's reads.
bwa_db=genome/bwa/gac
for sample in $(cut -f1 info/popmap.tsv) ;do
	fq_file=cleaned/$sample.fq.gz
	bam_file=alignments/$sample.bam
	log_file=alignments/$sample.oe
	{
		bwa mem -M $bwa_db $fq_file | samtools view -b > $bam_file
	} &> $log_file
done

# STEP 17-B-ii: Check that all alignments have completed.
echo "Checking that all BWA runs have completed..."
ls alignments/*.oe | wc
wc -l alignments/*.oe
grep -iE '\b(err|e:|warn|w:|fail|abort)' alignments/*.oe
tail -n5 alignments/*.oe | head -n100

# STEP 17-B-iii: Change directory.
cd stacks.ref

# STEP 17-B-iv: Run pstacks for every sample.
popmap=../info/popmap.tsv
index=1
for sample in $(cut -f1 $popmap) ;do
	bam_file=../alignments/$sample.bam
	log_file=$sample.pstacks.oe
	pstacks -f $bam_file -i $index -o ./ &> $log_file
	index=$(( $index + 1 ))
done

# STEP 17-B-v: Check that the pstacks runs have completed.
echo "Checking that all pstacks runs have completed..."
ls *.pstacks.oe | wc
wc -l *.pstacks.oe
grep -iE '\b(err|e:|warn|w:|fail|abort)' *.pstacks.oe
grep -L 'pstacks is done' *.pstacks.oe

# STEP 17-B-vi: Run cstacks.
cstacks --aligned -P ./ -M $popmap &> cstacks.oe

# STEP 17-B-vii: Run sstacks on every sample.
for sample in $(cut -f1 $popmap) ;do
	log_file=$sample.sstacks.oe
	sstacks --aligned -c ./batch_1 -s ./$sample -o ./ &> $log_file
done

# STEP 18: We use rxstacks to improve calls. We create a new subdirectory.
mkdir rxstacks
cd rxstacks
popmap=../$popmap # (After cd'ing we must update the popmap path.)

# STEP 19: Run rxstacks.
rxstacks -P ../ -o ./ --prune_haplo --conf_lim=0.10 &> rxstacks.oe

# STEP 20: Re-run cstacks and sstacks
# (Same commands as above, with updated paths.)
cstacks --aligned -P ./ -M $popmap &> cstacks.oe
for sample in $(cut -f1 $popmap) ;do
	log_file=$sample.sstacks.oe
	sstacks --aligned -c ./batch_1 -s ./$sample -o ./ &> $log_file
done

# STEP 21: Filter genotypes with populations.
min_samples=0.80
min_maf=0.05
max_obs_het=0.70
populations -P ./ -r $min_samples --min_maf $min_maf --max_obs_het $max_obs_het &> populations.oe
