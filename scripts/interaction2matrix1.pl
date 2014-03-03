#!/usr/bin/perl -w

#The script is used to tranform interaction txt to matrix.

use strict;
use Getopt::Std;

## get options from command line
my %opts= ('a'=>'');

sub usage{
print STDERR <<EOF;
    usage:  $0 -a interaction_txt [-h]

     -h   : help message;
     -a   : interaction file from scripts/overlapMapped2HiCFragments.pl;

EOF
    exit;
}

getopts('a:h', \%opts) || usage();
usage() if ($opts{h});

my ($txt_file)=($opts{a});

my %interaction;
my %chr_idx;
my %max_chr_idx;
open (IN,"$txt_file") || die;
while(<IN>){
    if(/(\S+?\_chr(\S+?)\_(\d+?)\|\S+?\|\S+?)\@\d+\s+(\S+?\_chr(\S+?)\_(\d+?)\|\S+?\|\S+?)\@\d+/){
	my ($lstr,$lchr,$lidx,$rstr,$rchr,$ridx)=($1,$2,$3,$4,$5,$6);
	#interaction between i and j, ie. C(i,j) is the same as C(j,i), only count once
	if($interaction{$lstr}{$rstr}){
	    $interaction{$lstr}{$rstr}++;
}
	elsif($interaction{$rstr}{$lstr}){
	    $interaction{$rstr}{$lstr}++;
}
	else{
	    $interaction{$lstr}{$rstr}=1;
}
	#full information of index in each chromosome
	$chr_idx{$lchr}{$lidx}=$lstr;
	$chr_idx{$rchr}{$ridx}=$rstr;
	#the maximum index in each each chromosome
	$max_chr_idx{$lchr}=$lidx if (!$max_chr_idx{$lchr} || ($max_chr_idx{$lchr}<$lidx));
	$max_chr_idx{$rchr}=$ridx if (!$max_chr_idx{$rchr} || ($max_chr_idx{$rchr}<$ridx));
}
}
close(IN);


##create matrix
my %count_matrix;
foreach my $chrm1 (sort {$a <=> $b} keys %chr_idx){
    next if ($chrm1 ne "10");
    foreach my $chrm2 (sort {$a <=> $b} keys %chr_idx){
	next if ($chrm2 ne "10");
	warn $chrm1,"\t",$chrm2,"\n";
	my $max_idx1=$max_chr_idx{$chrm1};
	my $max_idx2=$max_chr_idx{$chrm2};
	warn $max_idx1,"\t",$max_idx2,"\n";
	for(my $id1=1;$id1<=$max_idx1;$id1++){
	    next if (!$chr_idx{$chrm1}{$id1});
	    warn $id1,"\n" if ($id1%1000==0);
	    my $cur_str1=$chr_idx{$chrm1}{$id1};
	    for(my $id2=1;$id2<=$max_idx2;$id2++){
		next if(!$chr_idx{$chrm2}{$id2});
		my $cur_str2=$chr_idx{$chrm2}{$id2};
		my $cur_count=0;
		if(!$interaction{$cur_str1}{$cur_str2} && !$interaction{$cur_str2}{$cur_str1}){
		    $cur_count=0;
}
		else{
		    $cur_count=exists($interaction{$cur_str1}{$cur_str2}) ? $interaction{$cur_str1}{$cur_str2} : $interaction{$cur_str2}{$cur_str1};
}
		$count_matrix{$chrm1}{$id1}.="\t" . $cur_count;
}
}
}
}


## print matrix

foreach my $chrm (keys %chr_idx){
    next if ($chrm ne "10");
    my $max_idx=$max_chr_idx{$chrm};
    for(my $i=1;$i<=$max_idx;$i++){
	next if (!$chr_idx{$chrm}{$i});
	print "\t",$chr_idx{$chrm}{$i};
}
}
print "\n";

foreach my $chrm (keys %chr_idx){
    my $max_idx=$max_chr_idx{$chrm};
    for(my $i=1;$i<=$max_idx;$i++){
	next if (!$chr_idx{$chrm}{$i});
	print $chr_idx{$chrm}{$i},$count_matrix{$chrm}{$i},"\n";
}
}

