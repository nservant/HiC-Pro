#!/usr/bin/perl -w

#The script is used to merge alignment from bowtie1 and bowtie2

use strict;
use Getopt::Std;
use IO::Handle;

## get options from command line
my %opts= ('a'=>'','b'=>'','o'=>'output.data');

sub usage{
    print STDERR <<EOF;
    usage:  $0 -a aln_file_bowtie1 -b aln_file_bowtie2 -o output_name [-h]

    -h   : help message;
    -a   : alignment file from bowtie1;
    -b   : alignment file from bowtie2;
    -o   : output file;
EOF
    exit;
}


getopts('a:b:o:h', \%opts) || usage();
usage() if ($opts{h});

my ($aln_file1,$aln_file2,$output_name)=($opts{a},$opts{b},$opts{o});

open(g2_aln,"$aln_file2") || die;
my $cur_g2_line=g2_aln->getline();
my %g2_res;
while($cur_g2_line) {
    my @cur_g2_arr=split(/\s/,$cur_g2_line);
    my $cur_g2_name=$cur_g2_arr[0];
    $g2_res{$cur_g2_name}=$cur_g2_line;
    $cur_g2_line=g2_aln->getline();
}
close(g2_aln);


open(g1_aln,"$aln_file1") || die;
my $cur_g1_line=g1_aln->getline();
open(OUT, ">$output_name");
while($cur_g1_line) {
    my @cur_g1_arr=split(/\s/,$cur_g1_line);
    my $cur_g1_name=$cur_g1_arr[0];
    if($g2_res{$cur_g1_name}){
	print OUT $g2_res{$cur_g1_name};
    }
    else{
	print OUT $cur_g1_line;
    }
    $cur_g1_line=g1_aln->getline();
}
close(g1_aln);
close(OUT);

