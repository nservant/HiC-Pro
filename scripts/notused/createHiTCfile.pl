#!/usr/bin/perl -w
#The script is used to create files for HiTC package

use strict;
use Getopt::Std;

## get options from command line
my %opts= ('a'=>'');

sub usage{
    print STDERR <<EOF;
    usage:  $0 -a interaction_file [-h]

	-h   : help message;
    -a   : interaction file

EOF
    exit;
}

getopts('a:h', \%opts) || usage();
usage() if ($opts{h});

my ($input_file)=($opts{a});

open(IN,"$input_file") || die;
my %lhash;
my %rhash;
open(MAPOUT,">$input_file.map");
while(<IN>){
    if(/(\S+)\|\S+\|(\S+)\@\d+\s+(\S+)\|\S+\|(\S+)\@\d+/){
	my ($lsite,$lcoord,$rsite,$rcoord)=($1,$2,$3,$4);
	print MAPOUT $lsite,"\t",$rsite,"\t1\n";
	$lhash{$lcoord}=$lsite;
	$rhash{$rcoord}=$rsite;
#	warn $rcoord,"\t",$rsite,"\n";
#	sleep 2;
}
}
close(IN);
close(MAPOUT);

open (XBED,">$input_file.xbed");
foreach my $item (sort {$a cmp $b} keys %lhash){
    if($item=~/(\S+)\:(\d+)\-(\d+)/){
	print XBED $1,"\t",$2,"\t",$3,"\t",$lhash{$item},"\n";
}
}
close(XBED);

open (YBED,">$input_file.ybed");
foreach my $item1 (sort {$a cmp $b} keys %rhash){
    if($item1=~/(\S+)\:(\d+)\-(\d+)/){
        print YBED $1,"\t",$2,"\t",$3,"\t",$rhash{$item1},"\n";
    }
}
close(YBED);

