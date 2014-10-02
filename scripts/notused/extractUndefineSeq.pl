#!/usr/bin/perl -w

#The script is used to extract sequences whose alignment type is not yet defined. 

use strict;
use Getopt::Std;
use IO::Handle;

## get options from command line
my %opts= ('a'=>'','o'=>'output.data');

sub usage{
print STDERR <<EOF;
    usage:  $0 -a aln_file -o output_name [-h]

     -h   : help message;
     -a   : alignment file from merge step;
     -o   : the output name, default is output.data;

EOF
    exit;
}

getopts('a:o:h', \%opts) || usage();
usage() if ($opts{h});

my ($aln_file,$output_name)=($opts{a},$opts{o});

open(g_aln,"$aln_file") || die;

my $cur_g_line=g_aln->getline();


open (OUT,">$output_name");
while($cur_g_line) {
    #store read mapping information in array
    chomp $cur_g_line;
    my @cur_g_arr=split(/\t/,$cur_g_line);

    #algnment state: U/R/NM/? 
    my $cur_g_alnstat=$cur_g_arr[4];

    if(($cur_g_alnstat eq "UN") || ($cur_g_alnstat eq "?")){
	print OUT $cur_g_arr[0],"\n",$cur_g_arr[2],"\n+\n",$cur_g_arr[3],"\n";
    }
    $cur_g_line=g_aln->getline();
}

close(g_aln);
close(OUT);
