#!/usr/bin/perl -w

#The script is used to tranform interaction txt to matrix.

use strict;
use Getopt::Std;

## get options from command line
my %opts= ('a'=>'','b'=>'chrX_mm9.chrX_mm9');

sub usage{
print STDERR <<EOF;
    usage:  $0 -a interaction_txt -b chrN_genome.chrN_genome [-h]

     -h   : help message;
     -a   : interaction file from scripts/overlapMapped2HiCFragments.pl;
     -b   : the two chromosomes chosen to create heatmap,format: chrN_genome.chrN_genome;

EOF
    exit;
}

getopts('a:b:h', \%opts) || usage();
usage() if ($opts{h});

my ($txt_file,$chrchr)=($opts{a},$opts{b});

my ($chrm_left,$chrm_right,$genome_left,$genome_right)=("","","","");
if($chrchr=~/chr(\S+?)\_(\S+)\.chr(\S+?)\_(\S+)/){
    ($chrm_left,$genome_left,$chrm_right,$genome_right)=($1,$2,$3,$4);
}
else{
    exit("wrong chromosome format\n");
}

my %interaction;
my %chr_idx;
my %max_chr_idx;
open (IN,"$txt_file") || die;
while(<IN>){
    if(/(\S+?\_chr(\S+?)\_(\d+?)\|(\S+?)\|\S+?\:\d+\-\d+).*\s+(\S+?\_chr(\S+?)\_(\d+?)\|(\S+?)\|\S+?\:\d+\-\d+).*/){
	my ($lstr,$lchr,$lidx,$lgenome,$rstr,$rchr,$ridx,$rgenome)=($1,$2,$3,$4,$5,$6,$7,$8);
	next if ((($lchr ne $chrm_left) || ($rchr ne $chrm_right)) && (($lchr ne $chrm_right) || ($rchr ne $chrm_left)));
	next if ((($lgenome ne $genome_left) || ($rgenome ne $genome_right)) && (($lgenome ne $genome_right) || ($rgenome ne $genome_left)));
	$lchr.="-" . $lgenome;
	$rchr.="-" . $rgenome;
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

my $file_prefix="";
if($txt_file=~/(\S+)\./){
    $file_prefix=$1;
}

##print matrix
my $genome_chrm_left=$chrm_left . "-" . $genome_left;
my $genome_chrm_right=$chrm_right .  "-" . $genome_right;
open(MATRIX,">$file_prefix.$chrchr.matrix");
open(XBED,">$file_prefix.$chrchr.xbed");
open(YBED,">$file_prefix.$chrchr.ybed");
foreach my $chrm1 (sort {$a cmp $b} keys %chr_idx){
    next if ($chrm1 ne $genome_chrm_right);
    foreach my $chrm2 (sort {$a cmp $b} keys %chr_idx){
	next if ($chrm2 ne $genome_chrm_left);
	warn $chrm1,"\t",$chrm2,"\n";
	my $max_idx1=$max_chr_idx{$chrm1};
	my $max_idx2=$max_chr_idx{$chrm2};
	warn $max_idx1,"\t",$max_idx2,"\n";
	for(my $i=1;$i<=$max_idx2;$i++){
	    next if (!$chr_idx{$chrm2}{$i});
	    print MATRIX "\t",$chr_idx{$chrm2}{$i};
	    if($chr_idx{$chrm2}{$i}=~/(\S+)\|\S+\|(\S+)\:(\d+)\-(\d+)/){
		print XBED $2,"\t",$3,"\t",$4,"\t",$1,"\n";
	    }
	}
	print MATRIX "\n";	
	for(my $id1=1;$id1<=$max_idx1;$id1++){
	    next if (!$chr_idx{$chrm1}{$id1});
	    warn $id1,"\n" if ($id1%1000==0);
	    my $cur_str1=$chr_idx{$chrm1}{$id1};
	    if($cur_str1=~/(\S+)\|\S+\|(\S+)\:(\d+)\-(\d+)/){
		print YBED $2,"\t",$3,"\t",$4,"\t",$1,"\n";
	    }
	    print MATRIX $cur_str1;
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
		print MATRIX "\t" . $cur_count;
}
	    print MATRIX "\n";
}
}
}
close(MATRIX);
close(XBED);
close(YBED);



