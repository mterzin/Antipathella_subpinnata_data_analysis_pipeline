#!/bin/bash

      ###################################################################################################
      #                                                                                                 #
#######                                      IMPORTANT                                                  #######
      #                                                                                                 #
      ###################################################################################################
      # The script must be ran from the directory where you have your demultiplexed raw reads stored (if using relative paths)
      # [/store/usr7eco1/Fede_store/Marko/02_Quality_adapter_filtering_truncation] in my case
      #
      # I placed the scripts in the same folder, although they can be put elsewhere as well (use absolute paths in this case)

######################################################################################################################
#Concatenate all input fastq files to one file to make the reference
######################################################################################################################

for file in As*.fastq # This will select all the raw files of Antipathella subpinnata
do
	bash split_proc_rad.sh $file 16
	mv clean.fastq clean.$file
done



  

	cat clean.* >> all.fastq 





######################################################################################################################
#Fastq to fasta
######################################################################################################################

	echo -e '\n> Running FastqToFasta...'
	perl FastqToFasta.pl all.fastq all.fasta all.qual

######################################################################################################################
#Make the reference
######################################################################################################################

	echo -e '\n> Making the reference...'
	perl CDR_CsCpI.pl all.fasta all.qual

######################################################################################################################
#Fastq to fasta, but SEPARATELY, not all the files together. This is needed for the 2bRADmodElisa script
######################################################################################################################

for file in clean*.fastq
do
	perl FastqToFasta.pl $file $file.fasta $file.qual
done

######################################################################################################################
#Now change .fastq.fasta extension to .fasta
######################################################################################################################

for file in clean*fastq.fasta;
	do mv "${file}" "${file/.fastq/}";
done

######################################################################################################################
#Now run the 2bRADmodElisa script
#This is needed to remove the restriction enzyme recognition sites
######################################################################################################################

# First we move the "clean" fasta files into the elaborazione folder

	mkdir elaborazione
	mv clean*fasta ./elaborazione

#Now we run the script and make the output folder Trimmed
#This is where we store files that will be used in STACKS

	mkdir Trimmed
	  bash 2bRADmodElisa_v1.sh
	
#	mv clean.fastq clean.$file
#done

./FastaStats.pl $PWD/CDR >cdr_stats.txt
