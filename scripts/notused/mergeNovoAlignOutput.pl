#!/usr/bin/perl -w

#The script is used to merge alignment results (Novoalign raw format) coming from aligning single-end reads to two similar (allele-specific) genomes which differ in snps.Two alignment files will be merged into one file. 

use strict;
use Getopt::Std;
use IO::Handle;

## get options from command line
my %opts= ('a'=>'','b'=>'','c'=>'','d'=>'','e'=>'','o'=>'output.data');

sub usage{
print STDERR <<EOF;
    usage:  $0 -a aln_file1 -b aln_file2 -c genome1 -d genome2 -e common_genome -o output_dir [-h]

     -h   : help message;
     -a   : alignment file from genome1 (novoalign raw output);
     -b   : alignment file from genome2 (novoalign raw output);
     -c   : name of genome1, eg. mm9-129S1;
     -d   : name of genome2, eg. mm9-CAST;
     -e   : name of common genome, eg. mm9-common;
     -o   : the output name, default is output.data;

EOF
    exit;
}

getopts('a:b:c:d:e:o:h', \%opts) || usage();
usage() if ($opts{h});

my ($aln_file1,$aln_file2,$g1_name,$g2_name,$g_name,$output_name)=($opts{a},$opts{b},$opts{c},$opts{d},$opts{e},$opts{o});

open(g1_aln,"$aln_file1") || die;
open(g2_aln,"$aln_file2") || die;
my $cur_g1_line=g1_aln->getline();
my $cur_g2_line=g2_aln->getline();

open (OUT,">$output_name");
while($cur_g1_line || $cur_g2_line) {
    if($cur_g1_line=~/^\#/ || $cur_g1_line=~/^\s*$/) {
		$cur_g1_line=g1_aln->getline();
		$cur_g2_line=g2_aln->getline();
		next;
	}
    
    chomp $cur_g1_line;
    chomp $cur_g2_line;

    #store read mapping information in array
    my @cur_g1_arr=split(/\t/,$cur_g1_line);
    my @cur_g2_arr=split(/\t/,$cur_g2_line);

    #algnment state: U/R/QC/NM/QL for novoalgn raw output
    my $cur_g1_alnstat=$cur_g1_arr[4];
    my $cur_g2_alnstat=$cur_g2_arr[4];
    
    if($cur_g1_alnstat eq "R") {
	print OUT $cur_g1_line,"\n";
    } elsif($cur_g2_alnstat eq "R") {
	print OUT $cur_g2_line,"\n";
    } else {
	
	my $cur_g1_len=$#cur_g1_arr;
	my $cur_g2_len=$#cur_g2_arr;
	
	my $cur_g1_seq=$cur_g1_arr[2];
	my $cur_g2_seq=$cur_g2_arr[2];

	my @output_arr;
	if((length($cur_g1_seq)>length($cur_g2_seq)) && ($cur_g1_len>4)){
	    @output_arr=@cur_g1_arr;
	    $output_arr[7].="_" .  $g1_name;		
	}
	elsif((length($cur_g1_seq)<length($cur_g2_seq)) && ($cur_g2_len>4)){
	    @output_arr=@cur_g2_arr;
	    $output_arr[7].="_" .  $g2_name;		
	}
	else{
	    if($cur_g1_len==$cur_g2_len) {
		@output_arr=@cur_g1_arr;
		if(($cur_g1_len==12) || (($cur_g1_len==13) && ($cur_g1_arr[13] eq $cur_g2_arr[13]))){
		    if($cur_g1_arr[8]==$cur_g2_arr[8]){
			$output_arr[7].="_" . $g_name;
		    }
		    else{
			$output_arr[7].="_undefine";
			$output_arr[4]="R";	    
		    }
		}
		elsif($cur_g1_len==13){
		    $output_arr[7].="_undefine";
		    $output_arr[4]="?";
		}
	    } elsif((($cur_g1_len<$cur_g2_len) && ($cur_g1_len>4)) || ($cur_g2_len==4)) {
		@output_arr=@cur_g1_arr;
		$output_arr[7].="_" .  $g1_name;	    
	    } else {
		@output_arr=@cur_g2_arr;
		$output_arr[7].="_" .  $g2_name;	
	    }
	}	
	my $output_str=join("\t",@output_arr);
	print OUT $output_str,"\n";
    }
    
    $cur_g1_line=g1_aln->getline();
    $cur_g2_line=g2_aln->getline();
	
}

close(g1_aln);
close(g2_aln);
close(OUT);
