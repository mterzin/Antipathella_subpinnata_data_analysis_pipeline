#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions
use Bio::SeqIO;

$scriptname=$0; $scriptname =~ s/.+\///g;
# -- program name
print "-"x60, "\n";
print "FastqToFasta.pl v 1.11 E Meyer\n";
print "Last modified 01 Oct 2013\n";
print "-"x60, "\n";

# -- program description and required arguments
unless ($#ARGV == 2)
        {print "Converts a fastq file from Illumina into fasta sequence and quality score files.\n";
        print "Output:\t fasta files of sequence and quality scores.\n";
        print "Usage:\t $scriptname fastq out_fasta out_qual\n";
        print "Arguments:\n";
        print "\t fastq\t name of fastq input file \n";
        print "\t out_fasta\t name of fasta output file \n";
        print "\t out_qual\t name of qual output file \n";
        print "\n"; exit;
        }

use warnings;
use Bio::SeqIO;
use Bio::Seq::Quality;

my $wrap = 1000;		# edit to be longer than maximum read length
my $fastqfile = $ARGV[0];
my $inseqs = new Bio::SeqIO(-file=>$fastqfile, -format=>"fastq");
my $osfile = $ARGV[1];
my $oqfile = $ARGV[2];
my $outseqs = new Bio::SeqIO(-file=>">$osfile", -format=>"fasta", -width=>$wrap);
my $outquals = new Bio::SeqIO(-file=>">$oqfile", -format=>"qual", -width=>$wrap);

my %sh; my $scount = 0;
while ($seq = $inseqs->next_seq) 
	{$sh{$seq->display_id} = $seq->seq;
	$outseqs->write_seq($seq);
	$outquals->write_seq($seq);
	$scount++;
	}

print "Converted ", $scount, " reads from FASTQ to FASTA and QUAL.\n";

