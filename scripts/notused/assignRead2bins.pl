#!/usr/bin/perl -w
# chongjian chen 23/07/2012
# The script is used to assign reads to bins, output one matrix file, and two bed files (xbed and ybed)

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use POSIX qw(ceil floor);


##----
#options
##----

sub check_options {
    my $opts = shift;
    my ($map_txt,$gsize_txt,$bin_size,$bin_adjust,$step,$size_norm,$chrchr,$output,$obed)=("","",0,0,0,0,"","output",0);

    $bin_adjust=1 if ( exists($opts->{'binAdjust'}) );
    $size_norm=1 if ( exists($opts->{'sizeNorm'}) );
    $obed=1 if ( exists($opts->{'obed'}) );

    $chrchr=$opts->{'chrchr'} if ( exists($opts->{'chrchr'}) );
    $output=$opts->{'output'} if ( exists($opts->{'output'}) );
   
    if( exists($opts->{'mapTxt'}) ) {
	$map_txt = $opts->{'mapTxt'};
    } else {
	print "Option mapTxt is required.\n";
	exit;
    }
	
    if( exists($opts->{'binSize'}) ) {
	$bin_size = $opts->{'binSize'};
    } else {
	print "Option binSize is required.\n";
	exit;
    }
	
    if( exists($opts->{'genomeDesc'}) ) {
	$gsize_txt = $opts->{'genomeDesc'};
    } else {
	print "Option genomeDesc is required.\n";
	exit;
    }
    
    if( exists($opts->{'step'}) ) {
	$step = $opts->{'step'};
    } else {
	print "Option step is required.\n";
	exit;
    }
    
    return ($map_txt,$gsize_txt,$bin_size,$bin_adjust,$step,$size_norm,$chrchr,$output, $obed);
}	

my %options;
GetOptions( \%options,'help|?','man','mapTxt=s','binSize=i','genomeDesc=s','binAdjust','step=i','sizeNorm','chrchr=s','output=s','obed');

##usage
pod2usage(1) if ($options{"help"});
pod2usage(-verbose => 2) if ($options{"man"});

##check options
my ($map_txt,$gsize_txt,$bin_size,$bin_adjust,$step,$size_norm,$chrchr,$output, $obed)=check_options(\%options);


##----
#main
##----


##check if only select two chromosomes
my ($chrm_left,$chrm_right,$genome_left,$genome_right)=("","","","");
if($chrchr){
    if($chrchr=~/(\S+?)\_(\S+)\.(\S+?)\_(\S+)/){
	($chrm_left,$genome_left,$chrm_right,$genome_right)=($1,$2,$3,$4);
    }
    else{
	exit("wrong chromosome format\n");
    }
}


##load size of each chromosome into hash
open(DESC,"$gsize_txt") || die;
my %chrm_size;
while(<DESC>){
    if(/(\S+)\s+(\d+)/){
	$chrm_size{$1}=$2;
}
}
close(DESC);


##get bin/step information
my %chrm_bin_def;

foreach my $chrm (keys %chrm_size){
    ##in case binAdjust option is used
    if($bin_adjust){
#	warn $chrm,"\t",$chrm_size{$chrm},"\n";
	if ($chrm_size{$chrm}<$bin_size){
	    $chrm_bin_def{$chrm}{"binsize"}=$chrm_size{$chrm};
	    $chrm_bin_def{$chrm}{"stepsize"}=$chrm_size{$chrm};
	    $chrm_bin_def{$chrm}{"nbin"}=1;
	    next;
	}
	my ($cur_binsize,$cur_stepsize,$cur_nbin)=adjust_bin($chrm_size{$chrm},$bin_size,$step);
	$chrm_bin_def{$chrm}{"binsize"}=$cur_binsize;
	$chrm_bin_def{$chrm}{"stepsize"}=$cur_stepsize;
	$chrm_bin_def{$chrm}{"nbin"}=$cur_nbin;
##	warn $chrm,"\t",$cur_binsize,"\t",$cur_stepsize,"\t",$cur_nbin,"\n";
    }
    ##default
    else{
	if ($chrm_size{$chrm}<$bin_size){
	    $chrm_bin_def{$chrm}{"binsize"}=$chrm_size{$chrm};
	    $chrm_bin_def{$chrm}{"stepsize"}=$chrm_size{$chrm};
	    $chrm_bin_def{$chrm}{"nbin"}=1;
	    next;
	}
	$chrm_bin_def{$chrm}{"binsize"}=$bin_size;
	$chrm_bin_def{$chrm}{"stepsize"}=floor($bin_size/$step);
	my $cur_remainder=($chrm_size{$chrm}-$bin_size)%$chrm_bin_def{$chrm}{"stepsize"};
	my $cur_nbin=1+floor(($chrm_size{$chrm}-$bin_size)/$chrm_bin_def{$chrm}{"stepsize"});
	$chrm_bin_def{$chrm}{"nbin"}=$cur_remainder>0 ? $cur_nbin+1 : $cur_nbin;
##	warn $chrm,"\t",$chrm_bin_def{$chrm}{"binsize"},"\t",$chrm_bin_def{$chrm}{"stepsize"},"\t",$chrm_bin_def{$chrm}{"nbin"},"\n";
    }
}


##read interaction txt file
my %bin_sites;
my %checked_marks;
my %interaction;
my %chr_idx;
my %genomes; #all names of genomes
open (IN,"$map_txt") || die ("$map_txt can not be located!!!");
while(<IN>){
    if(/^\S+?\_(\S+?\_\d+?)\|(\S+?)\|(\S+?)\:(\d+)\-\d+\@(\d+)\s+\S+?\_(\S+?\_\d+?)\|(\S+?)\|(\S+?)\:(\d+)\-\d+\@(\d+)/){
	my ($lmark,$lgenome,$lchr,$lst,$ldist,$rmark,$rgenome,$rchr,$rst,$rdist)=($1,$2,$3,$4,$5,$6,$7,$8,$9,$10);
	$genomes{$lgenome}=1;
	$genomes{$rgenome}=1;
#	warn $_;
#	warn $chrm_left,"\t",$lchr,"\t",$chrm_right,"\t",$rchr,"\n";
#	warn $genome_left,"\t",$lgenome,"\t",$genome_right,"\t",$rgenome,"\n";
#	sleep 2;

	#filter based on chr if specified
	next if ($chrm_left && (($lchr ne $chrm_left) || ($rchr ne $chrm_right)) && (($lchr ne $chrm_right) || ($rchr ne $chrm_left)));
	next if ($genome_left && (($lgenome ne $genome_left) || ($rgenome ne $genome_right)) && (($lgenome ne $genome_right) || ($rgenome ne $genome_left)));
#	warn $_;
#	warn $chrm_left,"\t",$lchr,"\t",$chrm_right,"\t",$rchr,"\n";
#	warn $genome_left,"\t",$lgenome,"\t",$genome_right,"\t",$rgenome,"\n";
#	sleep 2;

	#get bin information for each read
	my ($lhit,$lbin_idx)=assign_bin($lchr,$lst,$ldist,$lgenome);
	my ($rhit,$rbin_idx)=assign_bin($rchr,$rst,$rdist,$rgenome);


	#count number of restriction sites in each bin
	my $lbin_idx1=$lchr . "_" . $lbin_idx . "|" . $lgenome;
	my $rbin_idx1=$rchr . "_" . $rbin_idx . "|" . $rgenome;
	
	$lmark.="_" . $lbin_idx1;
	$rmark.="_" . $rbin_idx1;

	$bin_sites{$lbin_idx1}++ if (!$checked_marks{$lmark});
	$bin_sites{$rbin_idx1}++ if (!$checked_marks{$rmark});

	$checked_marks{$lmark}=1;
	$checked_marks{$rmark}=1;

	##store information of interactions

	#add genome name to be able to handle allelic interactions
	$lchr.="@" . $lgenome;
	$rchr.="@" . $rgenome;

	#interaction between i and j, ie. C(i,j) is the same as C(j,i), only count once
	if($interaction{$lhit}{$rhit}){
	    $interaction{$lhit}{$rhit}++;
}
	elsif($interaction{$rhit}{$lhit}){
	    $interaction{$rhit}{$lhit}++;
}
	else{
	    $interaction{$lhit}{$rhit}=1;
}

	#full information of index in each chromosome
	$chr_idx{$lchr}{$lbin_idx}=$lhit;
	$chr_idx{$rchr}{$rbin_idx}=$rhit;

#	warn $_;
#	warn "-------------\n";
#	warn $lmark,"\t",$lchr,"\t",$lbin_idx,"\t",,$lhit,"\n";
#	warn $rmark,"\t",$rchr,"\t",$rbin_idx,"\t",,$rhit,"\n";
}
}
close(IN);


##output normalized or raw interaction map
my ($genome_lchr,$genome_rchr)=("","");
if($chrm_left){
    $genome_lchr=$chrm_left . "@" . $genome_left;
    $genome_rchr=$chrm_right .  "@" . $genome_right;
}

$chrchr="All" if (!$chrchr);

open(MATRIX,">$output.$chrchr.matrix") || die "$0: can't open output file\n";;
if ($obed == 1){
    open(XBED,">$output.$chrchr.xbed");
    open(YBED,">$output.$chrchr.ybed");
}
my %count_matrix;
my $chrm_ct=0;

my @genome_ids=keys %genomes;
if($#genome_ids>1){
    warn "unable to handle more than 2 different genomes\n";
    exit(0);
}


foreach my $chrm1 (sort {$a cmp $b} keys %chr_idx){

    #only focus on given chrm
    next if ($genome_rchr && ($chrm1 ne $genome_rchr));
    #to avoid to output both allelic information to the same axis
    next if (!$genome_rchr && ($chrm1!~/$genome_ids[0]/));

    $chrm_ct++;

    #get real chrm id and name of genome
    my ($realchrm1,$cur_genome1)=get_realchrm($chrm1);    
    #get the maximum id of bin in each chromosome
    my $max_idx1=$chrm_bin_def{$realchrm1}{"nbin"};
   
    foreach my $chrm2 (sort {$a cmp $b} keys %chr_idx){
	#only focus on given chrm
	next if ($genome_lchr && ($chrm2 ne $genome_lchr));

	next if (!$genome_lchr && (($#genome_ids==0 && ($chrm2!~/$genome_ids[0]/)) || ($#genome_ids==1 && ($chrm2!~/$genome_ids[1]/))));
	
#	warn $chrm1,"\t",$chrm2,"\n";

	#get real chrm id and name of genome
	my ($realchrm2,$cur_genome2)=get_realchrm($chrm2);
	#get the maximum id of bin in each chromosome
	my $max_idx2=$chrm_bin_def{$realchrm2}{"nbin"};

#	warn $max_idx1,"\t",$max_idx2,"\n";

	#During the first chrm1 loop, print head of matrix (X, colums in matrix) to matrix file, and X coordinates to bed file
	if($chrm_ct==1){
	    for(my $i=1;$i<=$max_idx2;$i++){		
		if($chr_idx{$chrm2}{$i}){
		    print MATRIX "\t",$chr_idx{$chrm2}{$i};
		    if($obed == 1){
			if($chr_idx{$chrm2}{$i}=~/(\S+)\|\S+\|(\S+)\:(\d+)\-(\d+)/){
			    print XBED $2,"\t",$3,"\t",$4,"\t",$1,"\n";
			}
		    }
		}
		else{
		    my $cur_binst=$chrm_bin_def{$realchrm2}{"stepsize"}*($i-1)+1;
		    my $cur_binend=$cur_binst+$chrm_bin_def{$realchrm2}{"binsize"}-1;
		    $cur_binend=$chrm_size{$realchrm2} if ($cur_binend>$chrm_size{$realchrm2});
		    my $cur_binname="HIC_" . $realchrm2 . "_" . $i . "|" . $cur_genome2 . "|" . $realchrm2 . ":" . $cur_binst . "-" . $cur_binend;
		    print MATRIX "\t",$cur_binname;
		    if ($obed == 1){
			print XBED $realchrm2,"\t",$cur_binst,"\t",$cur_binend,"\tHIC_",$realchrm2,"_",$i,"\n";
		    }
		}
	    }
	}
	
	for(my $id1=1;$id1<=$max_idx1;$id1++){
	    warn $id1,"\n" if ($id1%1000==0);

	    my $cur_str1="NNNNN";
	    #get information of current bin in chrm1, if this is a new bin, print head of Y coordinates (rows in matrix) to bed file
	    if($chr_idx{$chrm1}{$id1}){
		$cur_str1=$chr_idx{$chrm1}{$id1};
		if(!$count_matrix{$chrm1}{$id1}){
		    if($chr_idx{$chrm1}{$id1}=~/(\S+)\|\S+\|(\S+)\:(\d+)\-(\d+)/ && $obed == 1){
			print YBED $2,"\t",$3,"\t",$4,"\t",$1,"\n";
		    }
		}
	    }
	    else{
		my $cur_binst=$chrm_bin_def{$realchrm1}{"stepsize"}*($id1-1)+1;
		my $cur_binend=$cur_binst+$chrm_bin_def{$realchrm1}{"binsize"}-1;
		$cur_binend=$chrm_size{$realchrm1} if ($cur_binend>$chrm_size{$realchrm1});
		$cur_str1="HIC_" . $realchrm1 . "_" . $id1 . "|" . $cur_genome1 . "|" . $realchrm1 . ":" . $cur_binst . "-" . $cur_binend;
		if ($obed == 1){
		    print YBED $realchrm1,"\t",$cur_binst,"\t",$cur_binend,"\tHIC_",$realchrm1,"_",$id1,"\n" if (!$count_matrix{$chrm1}{$id1});
		}
	    }

	    #put rowname(name of bin in $chrm1)
	    $count_matrix{$chrm1}{$id1}=$cur_str1 if (!$count_matrix{$chrm1}{$id1});

	    for(my $id2=1;$id2<=$max_idx2;$id2++){
		my $cur_str2="XXXXX";
		$cur_str2=$chr_idx{$chrm2}{$id2} if ($chr_idx{$chrm2}{$id2});
		my $cur_count=0;
		if(!$interaction{$cur_str1}{$cur_str2} && !$interaction{$cur_str2}{$cur_str1}){
		    $cur_count=0;
}
		else{
		    $cur_count=exists($interaction{$cur_str1}{$cur_str2}) ? $interaction{$cur_str1}{$cur_str2} : $interaction{$cur_str2}{$cur_str1};
#		    warn "map: ",$cur_str1,"\t",$cur_str2,"\n";

		    #normalize read counts by the number of restriction sites in interacting bins
		    if($size_norm){
			my $cur_chrm1_bin_sites=get_site_count($cur_str1);
			my $cur_chrm2_bin_sites=get_site_count($cur_str2);
			$cur_count=$cur_count/($cur_chrm1_bin_sites*$cur_chrm2_bin_sites)
		    }
		}
		$count_matrix{$chrm1}{$id1}.= "\t" . $cur_count;
	    }
	}
    }
}


if ($obed == 1){
    close(XBED);
    close(YBED);
}
print MATRIX "\n";

foreach my $chrm (sort {$a cmp $b} keys %chr_idx){

    next if ($genome_rchr && ($chrm ne $genome_rchr));
    
#to avoid to output both allelic information to the same axis
    next if (!$genome_rchr && ($chrm!~/$genome_ids[0]/));

#    warn "output: ",$chrm,"\n";

    #get real chrm id and name of genome                                                                                                                  
    my ($realchrm,$cur_genome)=get_realchrm($chrm);

    #get the maximum id of bin in each chromosome                                                                                                         
    my $max_idx=$chrm_bin_def{$realchrm}{"nbin"};
   
    for(my $i=1;$i<=$max_idx;$i++){
	print MATRIX $count_matrix{$chrm}{$i},"\n";
    }
}

close(MATRIX);


##----
#functions
##----

#return adjusted bin size and bin step
sub adjust_bin{
    my ($gsize,$bin_size,$step)=@_;
    my $step_size=floor($bin_size/$step);
    my $nbin=1+floor(($gsize-$bin_size)/$step_size);
    my $cur_remainder=($gsize-$bin_size)%$step_size;
    my $bin_size_mod=$bin_size+floor($cur_remainder/$nbin);
    my $step_size_mod=floor($bin_size_mod/$step);
    if($bin_size==$bin_size_mod){
	return($bin_size_mod,$step_size_mod,$nbin);
}
    else{
	adjust_bin($gsize,$bin_size_mod,$step);
}
}

#return the index of bin given genomic coordinate
sub assign_bin{
    my ($chr,$st,$dist,$genome)=@_;
    my $loc=$st+$dist-1;
    my $binsize=$chrm_bin_def{$chr}{"binsize"};
    my $stepsize=$chrm_bin_def{$chr}{"stepsize"};
    my $cur_binidx=1+ceil(($loc-$binsize)/$stepsize);
    my $cur_binst=$stepsize*($cur_binidx-1)+1;
    my $cur_binend=$cur_binst+$binsize-1;
    $cur_binend=$chrm_size{$chr} if ($cur_binend>$chrm_size{$chr});
    my $cur_hit="HIC_" . $chr . "_" . $cur_binidx . "|" . $genome . "|" . $chr . ":" . $cur_binst . "-" . $cur_binend;
    return ($cur_hit,$cur_binidx);
}


#return real chrm and genome
sub get_realchrm{
    my ($chrm)=@_;
    my $genome="";
    if($chrm=~/(\S+)\@(\S+)/){
	($chrm,$genome)=($1,$2);
}
    return($chrm,$genome)
}

#return the number of restriction sites in given bin
sub get_site_count{
    my ($hit)=@_;
    my $binct=0;
    if($hit=~/HIC\_(\S+?\|\S+?)\|/){
	my $cur_id=$1;
	$binct=$bin_sites{$cur_id};
    }
    return ($binct);
}

__END__

=head1 NAME

- Using assignRead2bins.pl program

=head1 SYNOPSIS

perl assignRead2bins.pl [options]

  Options:
    --help               brief help message
    --man                full documentation
    --mapTxt=s           output file from overlapMapped2HiCFragments.pl (*.interaction)
    --genomeDesc=s       genome description file, format: chrm length;  can be generated using "samtools faidx ref.fa"
    --binSize=i          size of bin
    --step=i             number of windows in each bin used to compute the real size of step: binSize/step
    --sizeNorm           normalize read count in each bin by the number of restriction sites
    --binAdjust          slightly adjust the size of bin and step for each chromosome to avoid smaller size of the last bin
    --chrchr=s           choose only two chromosomes, format: chrN_genome.chrN_genome;         
    --output=s           prefix of output filename
    --obed                generate BED outputs of genome intervals

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-[other options]>

Check information in Usage

=back

=head1 DESCRIPTION

B<This program> is used to  assign reads to user-defined bins for each chromosome
    
=cut
