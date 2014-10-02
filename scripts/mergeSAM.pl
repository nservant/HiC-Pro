#!/usr/bin/perl -w

## Nicolas Servant
## Merge two SAM files (SE) into one BAM file
## Clean de SAM files by MAPQ, uniqueness, mapped

use Getopt::Std;
use strict;

## get options from command line
my $arg_num = scalar @ARGV;
my %opts= ('f'=>'','r'=>'', 'o'=>'', 'u'=>'', 'm'=>'', 'q'=>'', 'v'=>'');

sub usage{
print STDERR <<EOF;
    Merge two SE SAM files into one PE SAM file.
    Reads can be filtered according to MAPQ, mappability and uniqueness.
    A pair is reported only if both mates pass the different filters.

    usage:  $0 [-u -m -q -h] -f forward_mapped_file -r reverse_mapped_file -o output.sam 

     -h   : help message;
     -f   : forward read mapped file
     -r   : reverse read mapped file
     -u   : Do NOT report unmapped read pairs
     -m   : Do NOT report multiple mapped reads
     -q   : Do NOT report reads with a mapping quality lower than 'q'
     -v   : verbose mode
     -o   : output file name

    example: $0 -f mapping_R1.sam -r mapping_R2.sam
EOF
    exit;
}

getopts('r:f:o:q:umhv', \%opts) || usage();
usage() if ($opts{h});
usage() if $arg_num==0;

my($forward_file, $reverse_file, $output_file, $rm_unmapped, $rm_multi, $mapq, $quiet)=($opts{f},$opts{r}, ${opts{o}}, $opts{u}, $opts{m}, $opts{q}, $opts{v});

## Default - report all
$rm_unmapped=0 unless $rm_unmapped;
$rm_multi=0 unless $rm_multi;
$mapq=0 unless $mapq;
$quiet=0 unless $quiet;

paire_sam($forward_file, $reverse_file, $output_file, $rm_unmapped, $rm_multi, $mapq, $quiet);


###################
## paire_sam
## Merge two SE files in one PE file

sub paire_sam{
    my ($fileForward, $fileReverse, $filePair, $rm_unmap, $rm_mult, $rm_mapq, $verb) = @_;
    if ($verb == 1){
	print "##Pairing '$fileForward' and '$fileReverse' in '$filePair'\n";
	print "##Discard unmapped reads: $rm_unmap\n";
	print "##Discard multiple reads: $rm_mult\n";
	print "##Report reads with MAPQ >= $rm_mapq\n";
    }
    
    my ($f_read, $r_read, $isMulti);
    my $r_read_id="";
    my $f_read_id="";
    my @f_split_read;
    my @r_split_read;

    ##stats
    my ($tot_pairs_counter, $multi_pairs_counter, $uniq_pairs_counter, $unmapped_pairs_counter, $lowq_pairs_counter, $paired_reads_counter)=(0, 0, 0, 0, 0, 0);

    open (PAIRED, ">$filePair") or die "Can't create \'$filePair\' : $!";
    open (FORWARD, $fileForward) or die "Can't read \'$fileForward\' : $!";
    open (REVERSE, $fileReverse) or die "Can't read \'$fileReverse\' : $!";
    
    while(<FORWARD>){
	$f_read = $_;
	if(/^\s*$/){
	    next;
	}elsif(/^@/){    #Print SAM header lines so conversion to BAM is possible
	    print PAIRED $_;
	    next;
	}else{
	    @f_split_read=split /\t/, $f_read;
	    $f_read_id = $f_split_read[0];

	    while ($r_read_id ne $f_read_id){    
		$r_read = scalar <REVERSE> or last;
		if($r_read =~ /^\s*$/){
		    next;
		}elsif($r_read =~ /^@/){  #Ignore SAM format header lines  
		    next;
		}else{
		    @r_split_read=split /\t/, $r_read;
		    $r_read_id = $r_split_read[0];
		}
		
		if($f_read_id eq $r_read_id){
		    $isMulti=0;
		    $tot_pairs_counter++;
		    @f_split_read = split /\t/, $f_read;
		    @r_split_read = split /\t/, $r_read;

		    ## Unmapped reads
		    if(  (($f_split_read[1] & 0x4) or ($r_split_read[1] & 0x4)) and $rm_unmap eq 1){
			$unmapped_pairs_counter++;
			next;
		    }

		    ## Unique versus Multiple hits (bowtie 2 only) - note that should give the same results than filter out MAPQ=0 or 1
		    ## http://biofinysics.blogspot.fr/2014/05/how-does-bowtie2-assign-mapq-scores.html	  		    
		    if ($f_read =~ /\s+AS\:i\:(-?\d+)\s+XS\:i\:(-?\d+)/){
			my ($best_score,$secondary_score)=($1,$2);
			if ($best_score==$secondary_score){
			    $isMulti=1;
			}
		    }
		    if($isMulti==0 and $r_read =~ /\s+AS\:i\:(-?\d+)\s+XS\:i\:(-?\d+)/){
			my ($best_score,$secondary_score)=($1,$2);
			if ($best_score==$secondary_score){
			    $isMulti=1;
			}
		    }
		    ## Multihits
		    if ($isMulti == 1){
			$multi_pairs_counter++;
		    }else{
			$uniq_pairs_counter++;
		    }
		    if ($isMulti == 1 and $rm_mult == 1){
			next;
		    }

		    ## MAQ
		    if( ($f_split_read[4] < $rm_mapq) or ($r_split_read[4] < $rm_mapq)){
			$lowq_pairs_counter++;
			next;
		    }
		    
		    ## Reporting
		    if($f_read and $r_read){
		 	## Build pairs
		 	($f_read, $r_read) = sam_flag($f_read, $r_read);
			
			print PAIRED $f_read;
			print PAIRED $r_read;
		 	$paired_reads_counter++;
		    }
		}
	    }
	}    
    }
    close(PAIRED);
    close(FORWARD);
    close(REVERSE);
    if ($verb == 1){
	printf ("File\t%s\nTotal_pairs_processed\t%d\t(100)\nUnmapped_pairs\t%d\t(%.3f)\nLow_qual_pairs\t%d\t(%.3f)\nUnique_paired_alignments\t%d\t(%.3f)\nMultiple_pairs_alignments\t%d\t(%.3f)\nReported pairs\t%d\t(%.3f)\n",$filePair, $tot_pairs_counter, $unmapped_pairs_counter, $unmapped_pairs_counter/$tot_pairs_counter*100, $lowq_pairs_counter, $lowq_pairs_counter/$tot_pairs_counter*100, $uniq_pairs_counter, $uniq_pairs_counter/$tot_pairs_counter*100, $multi_pairs_counter, $multi_pairs_counter/$tot_pairs_counter*100, $paired_reads_counter, $paired_reads_counter/$tot_pairs_counter*100);
    }
}


############################
##
## Receives 2 single-end reads in SAM format and converts to paired-end read SAM format
##
sub sam_flag{
  my @f_read = split(/\t/, $_[0]);
  my @r_read = split(/\t/, $_[1]);
  ##my $report_unmapped=$_[2];

  #Relevant bitwise flags (flag in an 11-bit binary number)
  #1 The read is one of a pair
  #2 The alignment is one end of a proper paired-end alignment
  #4 The read has no reported alignments
  #8 The read is one of a pair and has no reported alignments
  #16 The alignment is to the reverse reference strand
  #32 The other mate in the paired-end alignment is aligned to the reverse reference strand
  #64 The read is the first (#1) mate in a pair
  #128 The read is the second (#2) mate in a pair
  
  #The reads were mapped as single-end data, so should expect flags of 
  #0 (map to the '+' strand) or 16 (map to the '-' strand)
  #Output example: a paired-end read that aligns to the reverse strand 
  #and is the first mate in the pair will have flag 83 (= 64 + 16 + 2 + 1)
  
  my $f_bitwise = $f_read[1];   
  my $r_bitwise = $r_read[1];

  #Ignore non-alignments
  ##if(  (($f_bitwise & 0x4) or ($r_bitwise & 0x4)) and $report_unmapped eq 0){
  ##    return(0,0);
  ##}else{
  if ($f_bitwise & 0x4){
      $r_bitwise = $r_bitwise | 0x8;
  }
  if($r_bitwise & 0x4){
      $f_bitwise = $f_bitwise | 0x8;
  }
  ##}
  
  if(  !($f_bitwise & 0x4) and !($r_bitwise & 0x4)){
      #The flag should now indicate this is paired-end data
      $f_bitwise = $f_bitwise | 0x1;
      $f_bitwise = $f_bitwise | 0x2;
      $r_bitwise = $r_bitwise | 0x1;
      $r_bitwise = $r_bitwise | 0x2;
  }
    
  #Indicate if the pair is on the reverse strand
  if($f_bitwise & 0x10){
      $r_bitwise = $r_bitwise | 0x20;
  }
  if($r_bitwise & 0x10){
      $f_bitwise = $f_bitwise | 0x20;
  }
  
  #Is this first or the second pair?
  $f_bitwise = $f_bitwise | 0x40;
  $r_bitwise = $r_bitwise | 0x80;
  
  #Insert the modified bitwise flags into the reads
  $f_read[1] = $f_bitwise;
  $r_read[1] = $r_bitwise;

  #Determine the RNEXT and PNEXT values (i.e. the positional values of a read's pair)
  #RNEXT
  if($f_read[2] eq $r_read[2]){
    $f_read[6] = '=';
    $r_read[6] = '=';
  }else{
    $f_read[6] = $r_read[2];
    $r_read[6] = $f_read[2];
  }
  #PNEXT
  $f_read[7] = $r_read[3];
  $r_read[7] = $f_read[3];

  my $f_read_string = join("\t", @f_read);
  my $r_read_string = join("\t", @r_read);

  return($f_read_string, $r_read_string);
}
