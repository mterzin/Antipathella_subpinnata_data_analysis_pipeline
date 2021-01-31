#!/usr/bin/perl
$scriptname=$0; $scriptname =~ s/.+\///g;
use Bio::Seq;
use Bio::SeqIO;
use warnings;

# -- program name
print "-"x60, "\n";
print "$scriptname v 1.0\n";
print "Created 04 Feb 2013 E Meyer\n";
print "Last modified 04 Feb 2013\n";
print "-"x60, "\n";

# -- program description and required arguments
unless ($#ARGV == 1)
        {print "Builds a cluster derived reference (CDR) from Illumina reads\n";
        print "Output:\t fasta-formatted nucleotide sequence file for use as a reference\n";
        print "Usage:\t $scriptname seq.file qual.file\n";
        print "Arguments:\n";
        print "\t seq.file\t nucleotide reads, fasta format\n";
        print "\t qual.file\t quality scores for those reads\n";
        print "\n"; exit;
        }

# -- user input --ho aggiunto "$site2" per trovare il sito di riconoscimento in reverse
my $csfile = $ARGV[0];
my $qualfile = $ARGV[1];

my $qthd = 10;					# scores lower than this are LQ
my $maxlq = 5;					# this many LQ scores allowed
my $site = "\\w{10}CAA\\w{5}GTGG\\w{10}";	# recognition site in Perl regexp
my $site2 = "\\w{10}CCAC\\w{5}TTG\\w{10}";
my $minobs = 2;					# min number of perfect matches for PAs
my $cpct = 0.90;				# percent identity for clustering PAs
my $word = 8;					# selected as specified below (from manual)
						# -n 8,9,10 for thresholds 0.90 ~ 1.0
						# -n 7      for thresholds 0.88 ~ 0.9
						# -n 6      for thresholds 0.85 ~ 0.88
						# -n 5      for thresholds 0.80 ~ 0.85
						# -n 4      for thresholds 0.75 ~ 0.8 
my $endtrunc = 32;				# end of region to be retained (base 1)
my $begtrunc = 1;				# beginning of region to be retained (base 1)

print "\nBeginning process...\n";

# -- stringently filter and truncate raw reads and scores
my $nqual = 0;
open(QF, $qualfile);
open(OQ, ">VHQ.qual");
while(<QF>)
	{chomp;
	if ($_ =~ />/) {$name = $_; $name =~ s/>//; $name =~ s/\s+.*//; $inqual++; next;}
	else	{@qa = split(" ", $_);
		$lqn = 0; $nn = 0; $poscount = 0;
		foreach $q (@qa) 
			{if ($poscount>$endtrunc) {next;}
			$poscount++; 
			if ($q < $qthd) {$lqn++;} 
			if ($q == -1) {$nn++;}
			}
		if ($nn>0) {$nqual++; next;}
		if ($lqn>$maxlq) {$nlq++; next;}
		elsif ($lqn <= $maxlq)
			{print OQ ">", $name, "\n";
			print OQ join(" ", @qa[$begtrunc-1..$endtrunc-1]), "\n";
			$gh{$name}++;
			$outqual++;
			}
		}
	}
close(QF);

open(SF, $csfile);
open(OS, ">VHQ.fasta");
while(<SF>)
	{chomp;
	if ($_ =~ />/) {$name = $_; $name =~ s/>//; $name =~ s/\s+.*//; $inseq++; next;}
	else	{
		if (exists($gh{$name})) 
			{print OS ">", $name, "\n";
			print OS substr($_, $begtrunc-1, $endtrunc+1-$begtrunc+1), "\n";
			$outseq++;
			}
		}
	}
close(SF);
close(OS);
print $inseq, " sequences in.\n";
print $inqual, " scores in.\n";
print $nqual, " records with Ns excluded.\n";
print $nlq, " records with low quality excluded.\n";
print $outseq, " sequences out.\n";
print $outqual, " scores out.\n";
print "Done filtering raw reads.\n";

# -- Filter for perfect match to restriction site --crea 2 files con i siti in For e Rev e poi li unisce in unico file matches.fasta
system("grep -P \"$site\" VHQ.fasta -B 1 > matches1.fasta");
system("grep -P \"$site2\" VHQ.fasta -B 1 > matches2.fasta");
system("cat matches1.fasta matches2.fasta > matches.fasta");
system("perl -pi -e \"s/^--\n//g\" matches.fasta");
print "Done filtering for site matches. Matches found:\n";
system("grep \">\" matches.fasta -c");

# -- cluster to identify repeatedly observed sequences (pseudoalleles, PA)
system("cd-hit-est -i matches.fasta -d 0 -M 0 -o 100pct -c 1 >pseudoalleles.log");
open(IN, "100pct.clstr");
while(<IN>)
	{chomp;
	if ($_ =~ /^>/) {$name = $_; $name =~ s/ /_/; $name =~ s/>//; next;}
	else	{
		$ch{$name}++;
		if ($_ =~ /\*/)
			{$rsn = $_; $rsn =~ s/.+>//; $rsn =~ s/\s.+//; $rsn =~ s/\.+//;
			$rh{$name} = $rsn;
			}
		}
	}
close(IN);
if (-e "PA.fasta") {system("rm PA.fasta");}
foreach $n (sort(keys(%ch)))
	{
	if ($ch{$n}>=$minobs)
		{system("grep $rh{$n} 100pct -A 1 >>PA.fasta");
		}
	}
print "Done clustering site matches and identifying pseudoalleles. Number found:\n";
system("grep \">\" PA.fasta -c");

# -- cluster PA to collapse alleles
system("cd-hit-est -d 0 -i PA.fasta -o CDR -c $cpct -n $word  >cdr.log");
print "Done clustering pseudoalleles into pseudosites. Number found:\n";
system("grep \">\" CDR -c");

print "Process completed.\n\n";
