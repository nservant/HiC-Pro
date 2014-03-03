#!/usr/bin/perl -w

###################################
## The script is used to convert BWT2 output to Nonoalign-like format. 
## Only output unique valid alignment 
###################################

use Getopt::Std;
use utf8;
use strict;
use POSIX;

## get options from command line
my %opts= ('a'=>'','o'=>'output.data');

sub usage{
print STDERR <<EOF;
 Description:

     Process Bowtie2 output file in SAM/BAM format. 
     Only unmapped/repeated/unique valid alignment are reported.
 
 Usage:  $0 -a <bwt2_output.sam\/.bam> -o <output> [-h]
     
     -h   : help message;
     -a   : the default bowtie2 output (SAM or BAM files);
     -o   : the name of output, the default name is output.data;
EOF
    exit;
}


getopts('a:o:h', \%opts) || usage();
usage() if ($opts{h});
#usage() if (@ARGV == 0 && -t STDIN);

my ($bwt2_output,$output)=($opts{a},$opts{o});

#load and store bwt2 output
my $bam=($bwt2_output =~ /.bam$/)? 1:0;
if($bam){
    open(BWT2, "samtools view $bwt2_output |") or die "$0: can't open ".$bwt2_output.":$!\n";
}else{
    open BWT2, "<".$bwt2_output or die "$0: can't open ".$bwt2_output.":$!\n";
}

#open(BWT2,"$bwt2_output") || die;
open(OUT,">$output");

while(<BWT2>){
    next if (/^\@/);
    my $aln_type;
    ## Unaligned reads
    if (/^(\S+)\s+4\s+\S+\s+\S+\s+\d+\s+\S+\s+\S+\s+\d+\s+\d+\s+(\S+)\s+(\S+)/){
	print OUT "@",$1,"\tS\t",$2,"\t",$3,"\tUN\n";
	next;
    }
    ## Multiple hits based on XS Flag (i.e. secondary aligment)
    if(/^(\S+)\s+\d+\s+\S+\s+\S+\s+\d+\s+\S+\s+\S+\s+\d+\s+\d+\s+(\S+)\s+(\S+).*\s+AS\:i\:(\d+)\s+XS\:i\:(\d+)/){
	my ($cur_seq_name,$cur_seq,$cur_quality,$best_score,$secondary_score)=($1,$2,$3,$4,$5);
	if ($best_score==$secondary_score){
	    print OUT "@",$cur_seq_name,"\tS\t",$cur_seq,"\t",$cur_quality,"\tR\t2\n";
	    next;
	}
    }
    ## Unique reads
    if(/^(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s+\d+\s+(\S+)\s+\S+\s+\d+\s+\d+\s+(\S+)\s+(\S+)/){
	my ($cur_seq_name,$cur_strand,$cur_chrm,$cur_coord,$cur_cigar,$cur_seq,$cur_quality)=($1,$2,$3,$4,$5,$6,$7);
	my $match_len;
	if($cur_cigar=~/^(\d+)S(\d+)M/){
	    my $soft_clip_len=$1;
	    $match_len=$2;
	    $cur_seq=substr($cur_seq,$soft_clip_len-1,$match_len);
	    $cur_quality=substr($cur_quality,$soft_clip_len-1,$match_len);
	}
	elsif($cur_cigar=~/^(\d+)M/){
	    $match_len=$1;
	    $cur_seq=substr($cur_seq,0,$match_len);
	    $cur_quality=substr($cur_quality,0,$match_len);
	}
	## Get strand information
	if($cur_strand==16){
	    $cur_seq=revcompl($cur_seq);
	    $cur_quality=reverse($cur_quality);
	    $cur_strand="R";
	}
	else{
	    $cur_strand="F";
	}
	print OUT "@",$cur_seq_name,"\tS\t",$cur_seq,"\t",$cur_quality,"\tU\t0\t150\t>",$cur_chrm,"\t",$cur_coord,"\t",$cur_strand,"\t.\t.\t.\n";
    }
}
close(BWT2);
close(OUT);


## return reverse complementary sequence
sub revcompl{
    my ($s)=@_;
    $s=reverse($s);
    $s=~ tr/ACGTacgt/TGCAtgca/;
    return $s;
}


