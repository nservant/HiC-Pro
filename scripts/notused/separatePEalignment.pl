#!/usr/bin/perl -w

#The script is used to separate strand paired end alignments. 

use strict;
use Getopt::Std;
use IO::Handle;

## get options from command line
my %opts= ('a'=>'','b'=>'','g'=>'','o'=>'output');

sub usage{
print STDERR <<EOF;
    usage:  $0 -a aln_file1 -b aln_file2 -g genome -o output_prefix [-h]

     -h   : help message;
     -a   : alignment file from left end of paired end read;
     -b   : alignment file from right end of paired end read;
     -g   : name of genome;
     -o   : the prefix of output name, default is output;
EOF
    exit;
}

getopts('a:b:g:o:h', \%opts) || usage();
usage() if ($opts{h});

my ($aln_file1,$aln_file2,$g_name,$output_name)=($opts{a},$opts{b},$opts{g},$opts{o});

open(g1_aln,"$aln_file1") || die;
open(g2_aln,"$aln_file2") || die;
my $cur_g1_line=g1_aln->getline();
my $cur_g2_line=g2_aln->getline();

open (OUT1,">$output_name.$g_name.1.out");
open (OUT2,">$output_name.$g_name.2.out");

while($cur_g1_line || $cur_g2_line) {
    if($cur_g1_line=~/^\#/ || $cur_g1_line=~/^\s*$/) {
	$cur_g1_line=g1_aln->getline();
	$cur_g2_line=g2_aln->getline();
	next;
    }
    
    chomp $cur_g1_line;
    chomp $cur_g2_line;
    
    my @cur_g1_arr=split(/\t/,$cur_g1_line);
    my @cur_g2_arr=split(/\t/,$cur_g2_line);
    
    my $cur_g1_alnstat=$cur_g1_arr[4];
    my $cur_g2_alnstat=$cur_g2_arr[4];
    
    if(($cur_g1_alnstat eq "NM") || ($cur_g2_alnstat eq "NM") || ($cur_g1_alnstat eq "UN") || ($cur_g2_alnstat eq "UN") || ($cur_g1_alnstat eq "?") || ($cur_g2_alnstat eq "?") || ($cur_g1_alnstat eq "R") || ($cur_g2_alnstat eq "R") ){
	$cur_g1_line=g1_aln->getline();
	$cur_g2_line=g2_aln->getline();
	next;
    }
    else{
	$cur_g1_line=~s/\_/\-/g;
	$cur_g2_line=~s/\_/\-/g;	
	print OUT1 $cur_g1_line,"\n";
	print OUT2 $cur_g2_line,"\n";
    }
    $cur_g1_line=g1_aln->getline();
    $cur_g2_line=g2_aln->getline();
}

close(g1_aln);
close(g2_aln);

close(OUT1);
close(OUT2);



