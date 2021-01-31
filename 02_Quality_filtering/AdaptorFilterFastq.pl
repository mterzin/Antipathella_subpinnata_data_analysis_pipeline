#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions
use Bio::Perl;
use Bio::SeqIO;
$scriptname=$0; $scriptname =~ s/.+\///g;

# -- program description and required arguments
unless ($#ARGV == 4)
        {print "\nFilters a set of short reads in FASTQ format, excluding any matching\n";
	print "the specified adaptor sequences\n";
        print "Usage:\t $scriptname sequences adaptors min_bp min_score output\n";
        print "Arguments:\n";
        print "\t sequences\t file of short reads to be filtered, fastq format\n";
        print "\t adaptors\t file of adaptor sequences to screen for, fasta format\n";
        print "\t min_bp\t\t length threshold; alignments this long are removed\n";
        print "\t min_score\t score threshold; alignments scoring this high are removed\n";
        print "\t output\t\t a name for the output file (fastq format)\n";
        print "\n"; exit;
        }

my $seqfile = $ARGV[0];		# raw reads, fastq format
my $adfile = $ARGV[1];		# adaptors, fasta format
my $maxn = 0;			# max number of Xs allowed in screened file
my $minmatch = $ARGV[2];	# min length of alignment between adaptor and read
my $minscore = $ARGV[3];	# min score (matching bases) of alignment
my $outfile = $ARGV[4];		# name for output file, fastq format

# convert fastq to fasta
open (IN, $seqfile);
open (TMP, ">tmp.fasta");
my $switch = 0;
my $count = 0;
while(<IN>)
	{
	chomp;
	$count++;
	if ($count==1) {$ss = substr($_, 0, 4);}
	if ($_ =~ /^$ss/) {$_ =~ s/\@//; print TMP ">", $_, "\n"; $switch = 1; next;}
	if ($_ =~ /^\+/) {$switch = 0; next;}
	if ($switch == 1) {print TMP $_, "\n";}
	}
close(IN);
close(TMP);

# adaptor searching
system("cross_match.manyreads tmp.fasta $adfile -minmatch $minmatch -minscore $minscore -screen >cross_match.log");

# sequence filtering -- first count number of Xs
my $seqs = new Bio::SeqIO(-file=>"tmp.fasta.screen", -format=>"fasta");
my $incount = 0; my $goodcount = 0; my $bcount = 0;
while ($seq = $seqs->next_seq)
        {
        $incount++;
        my $ss = $seq->seq;
        my $sid = $seq->display_id;
        $ncount = ($ss =~ tr/X//);
        if ($ncount > $maxn)    {$bcount++; next;}
	else
		{$goodcount++;
		$gh{$sid}++;
	        next;
                }
        }
print "Output from ", $scriptname, "\n";
print $incount, " reads input\n";
print $bcount, " reads failed\n";
print $goodcount, " reads passed\n";

# finally loop through fastq file and write out only the passing sequences
open (IN, $seqfile);
open (OUT, ">$outfile");
my $switch = 0;
while(<IN>)
	{
	chomp;
	if ($_ =~ /^$ss/) 
		{$_ =~ s/\@//; 
		$_ =~ s/\s.+//;
		if (exists($gh{$_}))
			{
			print OUT "@", $_, "\n";
			$switch = 1; 
			next;
			}
		else {$switch = 0;}
		}
	if ($switch == 1) {print OUT $_, "\n";}
	}
close(IN);
system("date");
system("rm tmp.fasta*; rm cross_match.log");

