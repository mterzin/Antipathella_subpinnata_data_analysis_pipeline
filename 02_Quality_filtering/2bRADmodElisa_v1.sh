#Script per preparare files fastq da analizzare con Stacks. Lo script è pensato per l'enzima CsCpI,senza il passaggio del mapping su un genoma di riferimento.
#!/bin/bash
a=elaborazione

b=Trimmed
#prima parte corretta
#test -d $a || mkdir -p $a
#test -d $b || mkdir -p $b
#for file in *fastq;
#do
#Taglio le reads dalla posizione 1 alla posizione 36.
#cat $file | perl TruncateFastq.pl $file 1 36 ./elaborazione/$file
#done

#seconda parte corretta
#for file in ./elaborazione/*fastq 
#do
#cat $file | perl FastqToFasta.pl $file $file.fasta $file.qual
#done

#da correggere Voglio cambiare l'estensione dei file da fastq.fasta to .fasta
#for file in ./elaborazione/*.fastq.fasta
#do
#mv fastq.fasta $file.fasta
#done

#rm ./elaborazione/*fasta
#comando centOS
#rename fastq.fasta fasta 
#./elaborazione/*fastq.fasta
#comando debian
#rename 's/\.fastq//' ./elaborazione/*fastq.fasta


#for file in ./elaborazione/*fasta
#do

#Utlizzo lo script 2b_extract per selezionare le sequenze contenenti il sito di riconoscimneto di CsCpI(For e Rew ).Salvo i file così trimmati nella cartella Trimmed,pronti per essere usati come input per Stacks.
#cat $file | perl CspCI_Extract.pl $file ".{10}CAA.{5}GTGG.{10}" $file.forw

#cat $file | perl CspCI_Extract.pl $file $file.forw
#done

#questo è il giusto...solo controllare il percorso per Trimmed
for file in ./elaborazione/*fasta
do
cat $file | perl CspCI_Extract.pl $file $file.forw && perl CspCI_Extract.pl $file $file.rew && perl ReverseComplement.pl $file.rew $file.rew2 && cat $file.forw $file.rew2 > ../Trimmed/T_$file

done



#perl CspCI_Extract.pl AsSAV_8r.fasta TRIM_AsSAV_8r.forw
#&& perl CspCI_Extract.pl $file ".{10}CCAC.{5}TTG.{10}" $file.rew && perl ReverseComplement.pl $file.rew $file.rew2 && cat $file.forw $file.rew2 > ../Trimmed/T_$file
#rm *.rew
#rm *.forw
#rm *.rew2
done
#rm -R ../elaborazione 

#exit
