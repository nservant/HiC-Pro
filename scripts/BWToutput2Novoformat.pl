#!/usr/bin/perl -w

#The script is used to convert BWT output to Nonoalign-like format, i.e. put information of unmapped or max-mapped reads to BWT ouput to keep the same order of reads as given sequence file


use Getopt::Std;
use utf8;
use strict;
use POSIX;

## get options from command line
my %opts= ('a'=>'','b'=>'','c'=>'','o'=>'output.data');

sub usage{
print STDERR <<EOF;
    usage:  $0 -a raw_sequence_file -b bwt_output -c maxmapped_sequence_file -o output [-h]

     -h   : help message;
     -a   : the raw sequence file;
     -b   : the default bowtie1 output;
     -c   : all reads with a number of valid alignments exceeding the limit set in bowtie1;
     -o   : the name of output, the default name is output.data;
EOF
    exit;
}

getopts('a:b:c:o:h', \%opts) || usage();
usage() if ($opts{h});

my ($raw_seqfile,$bwt_output,$maxmapped_seqfile,$output)=($opts{a},$opts{b},$opts{c},$opts{o});

#load and store bwt output
my %bwt_res;
open(BWT,"$bwt_output") || die;
while(<BWT>){
    chomp $_;
    my @array=split(/\t/,$_);
    my $cur_name=$array[0];
    my $cur_strand="F";
    $cur_strand="R" if ($array[1] eq "-");
    my $cur_seq=$array[4];
    $cur_seq=revcompl($cur_seq) if ($cur_strand eq "R");
    my $cur_quality=$array[5];
    $cur_quality=reverse($cur_quality) if ($cur_strand eq "R");
    if($array[7]){
	$bwt_res{$cur_name}="@" . $cur_name . "\tS\t" . $cur_seq . "\t" . $cur_quality . "\tU\t0\t150\t>" . $array[2] . "\t" . $array[3] . "\t" . $cur_strand . "\t.\t.\t.\t" . $array[7] . "\n";
    }
    else{
	$bwt_res{$cur_name}="@" . $cur_name . "\tS\t" . $cur_seq . "\t" . $cur_quality . "\tU\t0\t150\t>" . $array[2] . "\t" . $array[3] . "\t" . $cur_strand . "\t.\t.\t.\n";
    }
}
close(BWT);


#load and store maxmapped reads
my %maxmapped;
open(MAX,"$maxmapped_seqfile") || die;
while(<MAX>){
    if(/\@(.*)/){
	$maxmapped{$1}=1;
	<MAX>;
	<MAX>;
	<MAX>;
    }
}
close(MAX);

#read raw sequence file and output into novoalign-like result format.
open(RAW,"$raw_seqfile") || die;
open(OUTPUT,">$output");
while(<RAW>){
    if(/\@(.*)/){
	my $cur_name=$1;
	#get information of current read
	my $cur_seq=<RAW>;
	<RAW>;
	my $cur_quality=<RAW>;

	if($bwt_res{$cur_name}){
	    print OUTPUT $bwt_res{$cur_name};
	    next;
	}
	
	chomp $cur_seq;
	chomp $cur_quality;
	if ($maxmapped{$cur_name}){
	    print OUTPUT "@",$cur_name,"\tS\t",$cur_seq,"\t",$cur_quality,"\tR\t2\n";
	}
	else{
	    print OUTPUT "@",$cur_name,"\tS\t",$cur_seq,"\t",$cur_quality,"\tNM\n";
	}
    }    
}
close(RAW);
close(OUTPUT);


## return reverse complementary sequence
sub revcompl{
    my ($s)=@_;
    $s=reverse($s);
    $s=~ tr/ACGTacgt/TGCAtgca/;
    return $s;
}
