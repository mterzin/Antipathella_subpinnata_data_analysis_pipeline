#!/bin/bash
top=$(readlink -f $(dirname $0)/..)

# STEP 5:
echo "We expect:"
echo "<6-nt barcode>TGCAGG<unique 89-nt sequence>\n"

# STEP 6: Inspect the first read of each lane.
cd $top/raw
for lane in lane* ;do
	echo "$lane:"
	zcat $(ls $lane/*.fq.gz | awk NR==2) | head -n4
done
