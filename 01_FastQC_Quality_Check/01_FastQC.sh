#!/bin/bash

#Used this tutorial for quality check: https://metagenomics-workshop.readthedocs.io/en/latest/reads-qc/index.html

#This is done for all reads, and they will all be placed in one directory output directory [/store/usr7eco1/Fede_store/Marko/01_FastQC]

#This is the output directory where I want to place all the files after the quality check
  #It only needs to be done once, so I will comment it out

  # mkdir /store/usr7eco1/Fede_store/Marko/01_FastQC

#Making the PREFIX variable to extract sample name
for DATASET in $(ls /store/usr7eco1/Fede_store/Marko/02_Quality_adapter_filtering_truncation/As*.fastq); 
do
  PREFIX=$(basename $DATASET | cut -d "." -f 1)
  echo "for the file $DATASET, the prefix is $PREFIX"
  echo "and the current file is $(basename $DATASET)"

  #Running FastQC
  #I needed to add the directory in the specified path
  export PATH=/usr/arb/bin:/home/usr7eco1/Programs/SHRiMP/bin:/home/usr7eco1/Programs/bin:/home/usr7eco1/Paolo/Bacteria_exp_field_copy/Analyses/ncbi-blast-2.9.0+/bin:/home/usr7eco1/Programs/Stacks/bin:/home/usr7eco1/perl5/bin:/usr/arb/bin:/home/usr7eco1/Programs/SHRiMP/bin:/home/usr7eco1/Programs/bin:/home/usr7eco1/Paolo/Bacteria_exp_field_copy/Analyses/ncbi-blast-2.9.0+/bin:/home/usr7eco1/Programs/Stacks/bin:/home/usr7eco1/perl5/bin:/usr/lib64/qt-3.3/bin:/usr/lib64/ccache:/usr/local/bin:/usr/local/sbin:/usr/bin:/usr/sbin:/bin:/sbin:/home/usr7eco1/.local/bin:/home/usr7eco1/bin${PATH} && /home/usr7eco1/Programs/FastQC/fastqc /store/usr7eco1/Fede_store/Marko/02_Quality_adapter_filtering_truncation/$PREFIX.fastq

  #FastQC will generate 2 files for each input file: an HTML report (unzipped) and a zipped file. The file of interest is the unzipped one
  
done

#Then move the files to the output directory: /store/usr7eco1/Fede_store/Marko/01_FastQC
  mv /store/usr7eco1/Fede_store/Marko/02_Quality_adapter_filtering_truncation/*.html /store/usr7eco1/Fede_store/Marko/01_FastQC/
  mv /store/usr7eco1/Fede_store/Marko/02_Quality_adapter_filtering_truncation/*.zip /store/usr7eco1/Fede_store/Marko/01_FastQC/ # be careful not to have other zipped files/filders in the directory!!! You only want to remove the output from this step

#Finally, remove the part in the name that is not needed: _R1_001

for file in /store/usr7eco1/Fede_store/Marko/01_FastQC/*;
do
  mv "${file}" "${file/_R1_001/}";
done

    #FastQC will generate two files for each input file, one zipped and one unzipped (html report)
    #The results can be seen from the HTML file