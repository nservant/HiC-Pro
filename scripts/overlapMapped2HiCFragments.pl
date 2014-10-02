use strict;
use POSIX qw(ceil floor);
use List::Util qw[min max];
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Cwd;

sub check_options {
    my $opts = shift;
	
    my ($fragmentFile_1,$fragmentFile_2,$mappedReadsFile_1,$mappedReadsFile_2);
 
	if( exists($opts->{'fragmentFile_1'}) ) {
		$fragmentFile_1 = $opts->{'fragmentFile_1'};
	} else {
		print "Option fragmentFile_1|f1 is required.\n";
		exit;
	}
	
	if( exists($opts->{'fragmentFile_2'}) ) {
		$fragmentFile_2 = $opts->{'fragmentFile_2'};
	} else {
		print "Option fragmentFile_2|f2 is required.\n";
		exit;
	}
	
	if( exists($opts->{'mappedReadsFile_1'}) ) {
		$mappedReadsFile_1 = $opts->{'mappedReadsFile_1'};
	} else {
		print "Option mappedReadsFile_1|m is required.\n";
		exit;
	}
	
	if( exists($opts->{'mappedReadsFile_2'}) ) {
		$mappedReadsFile_2 = $opts->{'mappedReadsFile_2'};
	} else {
		print "Option mappedReadsFile_2|m is required.\n";
		exit;
	}
	
	return($fragmentFile_1,$fragmentFile_2,$mappedReadsFile_1,$mappedReadsFile_2);
}

sub badFormat($$$) {
	my $line=shift;
	my $lineNum=shift;
	my $errorType=shift;
	
	print "bad format @ line # $lineNum ($errorType) | $line\n";
	exit;
}

sub round($$) {
	my $num=shift;#the number to work on
	my $digs_to_cut=shift;# the number of digits after 
  
	if ($num=~/\d+\.(\d){$digs_to_cut,}/) {
		$num=sprintf("%.".($digs_to_cut-1)."f", $num);
	}		
	return $num;
}

## Read the fragment gff file and store the fragment position in %fragments and 
## the number of fragment per chromosome in %fragmentSizes
## Note : this is not the restriction site themselves, but the expected fragments

sub getFragments($) {
	my $dataFile=shift;

	my %fragments=();
	my %fragmentSizes=();
	
	open(IN,$dataFile);
	while(my $line = <IN>) {
		
		chomp($line);
		$line =~ s/>//;
		
		my ($fragmentName,$assembly,$coordinates)=split(/\|/,$line);
		my ($chr,$position)=split(/:/,$coordinates);
		my ($startPosition,$endPosition)=split(/-/,$position);
		
		my $index=0;
		$index=$fragmentSizes{$chr} if(exists($fragmentSizes{$chr}));
		$fragments{$chr}{$index}{'name'}=$line;
		$fragments{$chr}{$index}{'start'}=$startPosition;
		$fragments{$chr}{$index}{'end'}=$endPosition;
		$fragmentSizes{$chr}++;		
	}
	close(IN);
	
	return(\%fragments,\%fragmentSizes);
}


sub strand2direction($) {
	my $mappingStrand=shift;
	
	my $direction="";
	$direction="->" if($mappingStrand eq "F");
	$direction="<-" if($mappingStrand eq "R");
	
	return($direction);
}

sub commify {
   (my $num = shift) =~ s/\G(\d{1,3})(?=(?:\d\d\d)+(?:\.|$))/$1,/g; 
   return $num; 
}

sub annotate($$$$$$) {
	my $mappingStrand1=shift;
	my $mappingStrand2=shift;
	my $match1=shift;
	my $match2=shift;
	my $offset1=shift;
	my $offset2=shift;
	
	my @tmp1_1=split(/\|/,$match1);
	my @tmp1_2=split(/:/,$tmp1_1[2]);
	my $chr1=$tmp1_2[0];
	my @tmp1_3=split(/-/,$tmp1_2[1]);
	my $start1=$tmp1_3[0];
	my $end1=$tmp1_3[1];
	my $size1=$end1-$start1;
	
	my @tmp2_1=split(/\|/,$match2);
	my @tmp2_2=split(/:/,$tmp2_1[2]);
	my $chr2=$tmp2_2[0];
	my @tmp2_3=split(/-/,$tmp2_2[1]);
	my $start2=$tmp2_3[0];
	my $end2=$tmp2_3[1];
	my $size2=$end2-$start2;
	
	
	my ($direction1,$direction2,$direction,$interaction,$interactionOffset);
	$direction1=$direction2=$direction=$interaction=$interactionOffset="";
	
	($direction1)=strand2direction($mappingStrand1);
	($direction2)=strand2direction($mappingStrand2);
	
	if($match1 eq $match2) {
		if($offset1 < $offset2) {
			$direction=$direction1.".".$direction2;
			$interaction=$match1."\t".$match2;
			$interactionOffset=$match1."@".$offset1."\t".$match2."@".$offset2;
		} else {
			$direction=$direction2.".".$direction1;
			$interaction=$match2."\t".$match1;
			$interactionOffset=$match2."@".$offset2."\t".$match1."@".$offset1;
		}
	} else {
		if($chr1 eq $chr2) {
			if($start1 < $start2) {
				$direction=$direction1.".".$direction2;
				$interaction=$match1."\t".$match2;
				$interactionOffset=$match1."@".$offset1."\t".$match2."@".$offset2;
			} else {
				$direction=$direction2.".".$direction1;
				$interaction=$match2."\t".$match1;
				$interactionOffset=$match2."@".$offset2."\t".$match1."@".$offset1;
			}
		} else {
			if($chr1 lt $chr2) {
				$direction=$direction1.".".$direction2;
				$interaction=$match1."\t".$match2;
				$interactionOffset=$match1."@".$offset1."\t".$match2."@".$offset2;
			} else {
				$direction=$direction2.".".$direction1;
				$interaction=$match2."\t".$match1;
				$interactionOffset=$match2."@".$offset2."\t".$match1."@".$offset1;
			}
		}
	}
	
	return($interaction,$interactionOffset,$direction1,$direction2,$direction);
}

sub distanceDirection($$$$$$$$$$$) {
	my $junctionType = shift;
	my $direction1 = shift;
	my $match1 = shift;
	my $offset1 = shift;
	my $seqlen1 = shift;
	my $junction1 = shift;
	my $direction2 = shift;
	my $match2 = shift;
	my $offset2 = shift;
	my $seqlen2 = shift;
	my $junction2 = shift;
	
	my @tmp1_1=split(/\|/,$match1);
	my @tmp1_2=split(/:/,$tmp1_1[2]);
	my $chr1=$tmp1_2[0];
	my @tmp1_3=split(/-/,$tmp1_2[1]);
	my $start1=$tmp1_3[0];
	my $end1=$tmp1_3[1];
	my $size1=$end1-$start1;
	
	my @tmp2_1=split(/\|/,$match2);
	my @tmp2_2=split(/:/,$tmp2_1[2]);
	my $chr2=$tmp2_2[0];
	my @tmp2_3=split(/-/,$tmp2_2[1]);
	my $start2=$tmp2_3[0];
	my $end2=$tmp2_3[1];
	my $size2=$end2-$start2;
	
	my ($distanceNORM,$distanceSNP,$moleculeSize);
	$distanceNORM=$distanceSNP=$moleculeSize=0;

	if($junctionType eq "danglingEnd") {	
		$moleculeSize=($offset2-$offset1) if(($direction1 eq "->") and ($direction2 eq "<-"));
		$moleculeSize=($offset1-$offset2) if(($direction1 eq "<-") and ($direction2 eq "->"));
	} elsif($junctionType eq "selfCircle") {
		$moleculeSize=($offset1-$offset2) if(($direction1 eq "->") and ($direction2 eq "<-"));
		$moleculeSize=($offset1+($size2-$offset2)) if(($direction1 eq "<-") and ($direction2 eq "->"));
	} elsif($junctionType eq "interaction") {
		
		my ($moleculeSize1,$moleculeSize1NORM,$moleculeSize1SNP);
		$moleculeSize1=$moleculeSize1NORM=$moleculeSize1SNP=0;
		$moleculeSize1 = $moleculeSize1NORM = ($size1-$offset1-1) if($direction1 eq "->");
		$moleculeSize1 = $moleculeSize1SNP = $seqlen1 if(($junction1 != -1) and ($direction1 eq "->"));
		$moleculeSize1 = $moleculeSize1NORM = ($offset1-1) + $seqlen1 if($direction1 eq "<-");
		$moleculeSize1 = $moleculeSize1SNP = $seqlen1 if(($junction1 != -1) and ($direction1 eq "<-"));

		my ($moleculeSize2,$moleculeSize2NORM,$moleculeSize2SNP);
		$moleculeSize2=$moleculeSize2NORM=$moleculeSize2SNP=0;
		$moleculeSize2 = $moleculeSize2NORM = ($size2-$offset2-1) if($direction2 eq "->");
		$moleculeSize2 = $moleculeSize2SNP = $seqlen2 if(($junction2 != -1) and ($direction2 eq "->"));
		$moleculeSize2 = $moleculeSize2NORM = ($offset2-1) + $seqlen2 if($direction2 eq "<-");
		$moleculeSize2 = $moleculeSize2SNP = $seqlen2 if(($junction2 != -1) and ($direction2 eq "<-"));
		
		my ($inFrontDistance1,$behindDistance1,$inFrontDistance2,$behindDistance2);
		$inFrontDistance1=$behindDistance1=$inFrontDistance2=$behindDistance2=0;
		
		if($direction1 eq "->") {			
			$inFrontDistance1=($size1-$offset1-1);
			$behindDistance1=($offset1-1);
			$inFrontDistance1=$behindDistance1=$seqlen1 if($junction1 != -1);
		} elsif($direction1 eq "<-") {
			$inFrontDistance1=($offset1-1)+$seqlen1;
			$behindDistance1=($size1-$offset1-1);
			$inFrontDistance1=$behindDistance1=$seqlen1 if($junction1 != -1);
		}
		
		if($direction2 eq "->") {			
			$inFrontDistance2=($size2-$offset2-1);
			$behindDistance2=($offset2-1);
			$inFrontDistance2=$behindDistance2=$seqlen2 if($junction2 != -1);
		} elsif($direction2 eq "<-") {
			$inFrontDistance2=($offset2-1)+$seqlen2;
			$behindDistance2=($size2-$offset2-1);
			$inFrontDistance2=$behindDistance2=$seqlen2 if($junction2 != -1);
		}
		
		$moleculeSize=($moleculeSize1+$moleculeSize2);
		
		$moleculeSize1=min($inFrontDistance1,$behindDistance1);
		$moleculeSize2=min($inFrontDistance2,$behindDistance2);
		
		$moleculeSize=($moleculeSize1+$moleculeSize2);		
		
	}		
	
	return($moleculeSize);
	
}

my %dna2complement=();
$dna2complement{'A'}='T';
$dna2complement{'C'}='C';
$dna2complement{'G'}='G';
$dna2complement{'T'}='A';

sub normalRS2Super($) {
	my $restrictionSite=shift;
	
	my $antiRestrictionSite = reverse $restrictionSite;
	$antiRestrictionSite =~ tr/ACGTacgt/TGCAtgca/;
	$antiRestrictionSite = reverse $restrictionSite;
	
	my ($leftCut,$rightCut)=split(/\^/,$restrictionSite);
	my ($leftAntiCut,$rightAntiCut)=split(/\^/,$antiRestrictionSite);
	
	my $leftSuperRestrictionSite = $leftAntiCut;
	$leftSuperRestrictionSite =~ tr/ACGTacgt/TGCAtgca/;
	my $rightSuperRestrictionSite = $rightCut;
	my $superRestrictionSite=$leftSuperRestrictionSite.$rightSuperRestrictionSite;
	
	
	my $antiSuperRestrictionSite = reverse $superRestrictionSite;
	$antiSuperRestrictionSite =~ tr/ACGTacgt/TGCAtgca/;
	
	my $dnaSequence=$restrictionSite;
	$dnaSequence =~ s/\^//;
	return($dnaSequence,$superRestrictionSite,$antiSuperRestrictionSite);
}

# @HWI-ST570:35:D0ELHACXX:6:1101:7138:2249 1:N:0: S       NTAAGCCACTCTCGCCTGGCCCAGCCACACCTCAGGTGCTGAGAGCCCAG      #1=DDFFFHHHHHJJJJJJJJJJJJJJIGJJJJJJJGIJJJJJIJIJJJG    U       6       86      >chr3-cast      79656566        R       .       .       .
# @HWI-ST570:35:D0ELHACXX:6:1101:11664:2228 1:N:0:        S       NGTGTTGTGGATTGGCTGCTGTGCTCTGCTGATGCTTTACAGTTTGAAGA      #1:DDDDDHHDHHGIIIIIGGAFEDG>FEEHGIGGGIIE<DHDHIIDE>?    U       6       150     >chr1-cast      190661730       F       .       .       .

sub overlapMappedReads($$$$$$) {
	my $mappedReadsFile_1=shift;
	my $mappedReadsFile_2=shift;
	my $fragments_1={};
	$fragments_1=shift;
	my $fragmentSizes_1={};
	$fragmentSizes_1=shift;
	my $fragments_2={};
	$fragments_2=shift;
	my $fragmentSizes_2={};
	$fragmentSizes_2=shift;
	
	my %lastIndex=();
	
	my @tmp_1=split(/\//,$mappedReadsFile_1);
	my $mappedReadsFileName_1=$tmp_1[@tmp_1-1];
	my $name_1=$mappedReadsFileName_1;
	$name_1 =~ s/\.out//;
	$name_1 =~ s/\_perfect//;
	$name_1 =~ s/\.1$//;
	
	my @tmp_2=split(/\//,$mappedReadsFile_2);
	my $mappedReadsFileName_2=$tmp_2[@tmp_2-1];
	my $name_2=$mappedReadsFileName_2;
	$name_2 =~ s/\.out//;
	$name_2 =~ s/\_perfect//;
	$name_2 =~ s/\.2$//;
	
	die("error with files...($name_1 vs $name_2)\n") if($name_1 ne $name_2);
	
	open(OUT,">",$name_1.".1");
	%lastIndex=();
	
	#####################################################
        ## Firt step
        ## Reading input file and output the pairs information
        ## chr1	181511146	R	chr1	181510959	F
        ## chr2	78673481	R	chr16	98206868	R
        ## chr12	33698600	R	chr8	113172753	F
	## ...

	print "standardizing file [".$name_1.".1]...\n";
	open (SIDE1, $mappedReadsFile_1) or die $!;
	open (SIDE2, $mappedReadsFile_2) or die $!;
	while((!eof(SIDE1)) || (!eof(SIDE2))) { 
		
		# Reading both file line by line but only one line at a time.
		my $line1 = <SIDE1>;
		my $line2 = <SIDE2>;
		chomp ($line1);
		chomp ($line2);
		
		my ($chr1,$origPos1,$strand1)=(split(/\t/,$line1))[7,8,9];
		my ($chr2,$origPos2,$strand2)=(split(/\t/,$line2))[7,8,9];
		$chr1 =~ s/>//;
		$chr1=(split(/-/,$chr1))[0];
		$chr2 =~ s/>//;
		$chr2=(split(/-/,$chr2))[0];
		my $pos1=$origPos1+1;
		my $pos2=$origPos2+1;
		
		print OUT "$chr1\t$origPos1\t$strand1\t$chr2\t$origPos2\t$strand2\n";				
	}
	close(IN);
	close(OUT);
	
	print "\tsorting file ...\n";
	system("sort -k1,1 -k2,2n ".$name_1.".1 > ".$name_1.".1.sorted");
	##system("rm ".$name_1.".1");
	print "\tdone\n";
	
	## Assign the first read to a HindIII fragment
	print "processing first interval [".$name_1.".2]...\n";
	open(OUT,">",$name_1.".2");
	open(IN,$name_1.".1.sorted");
	while(my $line = <IN>) {
		
		chomp($line);
		
		my ($chr1,$origPos1,$strand1,$chr2,$origPos2,$strand2,$overlap1)=split(/\t/,$line);		
		$chr1=(split(/-/,$chr1))[0];
		$chr2=(split(/-/,$chr2))[0];
		my $pos1=$origPos1+1;
		my $pos2=$origPos2+1;
		
		$pos1+=(50-5) if($strand1 eq "R"); ## TO CHECK
		
		my $index=0;
		$index=$lastIndex{$chr1} if(exists($lastIndex{$chr1}));
		
		my $found=0;
		my $chrSize=$fragmentSizes_1->{$chr1};
		## For each reads pairs, loop on chromosome HindIII sites and look for a site 
		for(my $i=$index;$i<$chrSize;$i++) {
			my $tmpName=$fragments_1->{$chr1}->{$i}->{'name'};
			my $tmpStart=$fragments_1->{$chr1}->{$i}->{'start'};
			my $tmpEnd=$fragments_1->{$chr1}->{$i}->{'end'};
			
			$lastIndex{$chr1}=$i if($origPos1 > $tmpEnd);
			next if($origPos1 > $tmpEnd);
			
			if(($pos1 >= $tmpStart) and ($pos1 <= $tmpEnd)) {
				$found=1;
				#print "$chr1\t$origPos1\t$pos1\t$strand1\t$tmpName\t$chr2\t$origPos2\t$pos2\t$strand2\n";
				print OUT "$chr1\t$origPos1\t$strand1\t$chr2\t$origPos2\t$strand2\t$tmpName\n";				
				last;
			}
		}
		print "WARNING - could not find match for $chr1 | $pos1 | $strand1\n$line\n" if($found == 0);
	}
	close(IN);
	
	##system("rm ".$name_1.".1.sorted");
	
	close(OUT);
	
	print "\tsorting file...\n";
	# now sort by second fragment
	system("sort -k4,4 -k5,5n ".$name_1.".2 > ".$name_1.".2.sorted");
	##system("rm ".$name_1.".2");
	# now sort by second fragment
	print "\tdone\n";

	## Assign the second read to a HindIII fragment
	print "processing second interval...\n";
	open(OUT,">",$name_1.".3");
	%lastIndex=();
	
	open(IN,$name_1.".2.sorted");
	while(my $line = <IN>) {
		
		chomp($line);
		
		my ($chr1,$origPos1,$strand1,$chr2,$origPos2,$strand2,$overlap1)=split(/\t/,$line);		
		$chr1=(split(/-/,$chr1))[0];
		$chr2=(split(/-/,$chr2))[0];
		my $pos1=$origPos1+1;
		my $pos2=$origPos2+1;
		
		$pos2+=(50-5) if($strand2 eq "R");
		
		my $index=0;
		$index=$lastIndex{$chr2} if(exists($lastIndex{$chr2}));
		
		my $chrSize=$fragmentSizes_2->{$chr2};
		for(my $i=$index;$i<$chrSize;$i++) {
			my $tmpName=$fragments_2->{$chr2}->{$i}->{'name'};
			my $tmpStart=$fragments_2->{$chr2}->{$i}->{'start'};
			my $tmpEnd=$fragments_2->{$chr2}->{$i}->{'end'};
			
			$lastIndex{$chr2}=$i if($origPos2 > $tmpEnd);
			next if($origPos2 > $tmpEnd);
			
			if(($pos2 >= $tmpStart) and ($pos2 <= $tmpEnd)) {
				print OUT "$chr1\t$origPos1\t$strand1\t$chr2\t$origPos2\t$strand2\t$overlap1\t$tmpName\n";
				last;
			}
		}
		
	}
	close(IN);
	
	##system("rm ".$name_1.".2.sorted");
	
	close(OUT);
	print "\tdone\n";
	
	# now parse all reads
	
	my ($stats);
	$stats->{'noMap'} = 0;
	$stats->{'single'} = 0;
	$stats->{'peMapped'} = 0;
	$stats->{'side1_match'} = 0;
	$stats->{'side2_match'} = 0;
	$stats->{'valid'} = 0;
	$stats->{'invalid'} = 0;

	$stats->{'same'}{'->.->'} = 0;
	$stats->{'same'}{'->.<-'} = 0;
	$stats->{'same'}{'<-.<-'} = 0;
	$stats->{'same'}{'<-.->'} = 0;
	$stats->{'different'}{'->.->'} = 0;
	$stats->{'different'}{'->.<-'} = 0;
	$stats->{'different'}{'<-.->'} = 0;
	$stats->{'different'}{'<-.<-'} = 0;

	my %dist=();
	
	my $statFile=$name_1;

	#valid interaction
	my $interactionFile=$statFile.".interaction";
	open(INTERACTION,">$interactionFile") || die "SC file error $!\n";
	#only 1 end mapped
	my $singleFile=$statFile.".single";
	open(SINGLE,">$singleFile") || die "SC file error $!\n";	
	#self circle interaction
	my $selfCircleFile=$statFile.".selfCircle";
	open(SELFCIRCLE,">$selfCircleFile") || die "SC file error $!\n";
	#dangling end interaction
	my $danglingEndFile=$statFile.".danglingEnd";
	open(DANGLINGEND,">$danglingEndFile") || die "DE file error $!\n";
	#impossible paired read
	my $errorFile=$statFile.".error";
	open(ERROR,">$errorFile") || die "ERROR file error $!\n";

	#selfCirlce molecule sizes
	my $selfCircleSizesFile=$statFile.".selfCircle.sizes";
	open(SELFCIRCLESIZES,">$selfCircleSizesFile") || die "SC file error $!\n";
	#danglingEnd molecule sizes
	my $danglingEndSizesFile=$statFile.".danglingEnd.sizes";
	open(DANGLINGENDSIZES,">$danglingEndSizesFile") || die "DE file error $!\n";
	#valid interaction molecule sizez
	my $interactionSizesFile=$statFile.".interaction.sizes";
	open(INTERACTIONSIZES,">$interactionSizesFile") || die "SC file error $!\n";
	
	print "determing HiC junction...\n";
	open(IN,$name_1.".3");
	while(my $line = <IN>) {
		
		chomp($line);
		
		my ($chr1,$side1Offset,$side1Strand,$chr2,$side2Offset,$side2Strand,$side1Align,$side2Align)=split(/\t/,$line);		
		my $read1pos=$side1Offset;
		my $read2pos=$side2Offset;

		$side1Offset++;
		$side2Offset++;
		
		my ($fragmentName1,$assembly1,$coordinates1)=split(/\|/,$side1Align);
		my ($chr1,$position1)=split(/:/,$coordinates1);
		my ($startPosition1,$endPosition1)=split(/-/,$position1);
		$side1Offset=$side1Offset-$startPosition1;
		$chr1=(split(/-/,$chr1))[0];
		
		my ($fragmentName2,$assembly2,$coordinates2)=split(/\|/,$side2Align);
		my ($chr2,$position2)=split(/:/,$coordinates2);
		my ($startPosition2,$endPosition2)=split(/-/,$position2);
		$side2Offset=$side2Offset-$startPosition2;
		$chr2=(split(/-/,$chr2))[0];
		
		#override missing data
		my $side1Junction=-1;
		my $side2Junction=-1;
		my $side1Length=50;
		my $side2Length=50;
		
		# Count number of mached and non matched reads.
		$stats->{'side1_match'}++ if($side1Align ne "");
		$stats->{'side2_match'}++ if($side2Align ne "");
		$stats->{'peReads'}++;
		
		my ($moleculeSize);
		$moleculeSize=-1;
		
		## Both reads are assigned to a fragment
		if(($side1Align ne "") and ($side2Align ne "")) {

		    my ($interaction,$interactionOffset,$direction1,$direction2,$direction)=annotate($side1Strand,$side2Strand,$side1Align,$side2Align,$side1Offset,$side2Offset);
		    print "$chr1\t$read1pos\t$side1Strand\t$side1Align\t$chr2\t$read2pos\t$side2Strand\t$side2Align\t";
		    #print "==> $interaction - $interactionOffset - $direction1 - $direction2 - $direction \n";
		    
		    ## to the same fragment
		    if($side1Align eq $side2Align) {  
			$stats->{'same'}{$direction}++;
			if($direction eq "<-.->") { #selfCircle
			    print SELFCIRCLE "$interactionOffset\n";
			    ($moleculeSize)=distanceDirection('selfCircle',$direction1,$side1Align,$side1Offset,$side1Length,$side1Junction,$direction2,$side2Align,$side2Offset,$side2Length,$side2Junction);
			    print SELFCIRCLESIZES "$moleculeSize\n";

			    print "SELFCYCLE\n";## -$moleculeSize\n";

			} elsif($direction eq "->.<-") { #danglingEnd
			    print DANGLINGEND "$interactionOffset\n";
			    ($moleculeSize)=distanceDirection('danglingEnd',$direction1,$side1Align,$side1Offset,$side1Length,$side1Junction,$direction2,$side2Align,$side2Offset,$side2Length,$side2Junction);
			    print DANGLINGENDSIZES "$moleculeSize\n";

			    print "DANGLINGEND\n";## -$moleculeSize\n";

			} else {
			    print ERROR "$interactionOffset\n";
			    print "ERROR\n";
			}
		    }else{		  
			## To different fragments
			$stats->{'different'}{$direction}++;
			$stats->{'valid'}++;				
			print INTERACTION "$interactionOffset\n";
			($moleculeSize)=distanceDirection('interaction',$direction1,$side1Align,$side1Offset,$side1Length,$side1Junction,$direction2,$side2Align,$side2Offset,$side2Length,$side2Junction);
			print INTERACTIONSIZES "$moleculeSize\n";
			print "INTERACTION\n";## $moleculeSize\n";
		    }
		    $stats->{'peMapped'}++;	
		} elsif(($side1Align ne "") and ($side2Align eq "")) {
		    my $singleOffset=$side1Align."@".$side1Offset;
		    $stats->{'single'}++;
		    print SINGLE "$singleOffset\n";
		    print "SINGLE\n";

		} elsif(($side1Align eq "") and ($side2Align ne "")) {
		    my $singleOffset=$side2Align."@".$side2Offset;
		    $stats->{'single'}++;
		    print SINGLE "$singleOffset\n";
		    print "SINGLE\n";

		} else {
		    print "ERROR\n";
		    exit;
		}
		
	}
	close(IN);

	##system("rm ".$name_1.".3");
	
	my $peReads = $stats->{'peReads'};

	#both sides mapped - valid interaction.
	my $peMapped = $stats->{'peMapped'};
	my $peMapped_pc = 0;
	$peMapped_pc = (($peMapped / $peReads)*100) if($peReads != 0);
	$peMapped_pc=sprintf "%.1f", $peMapped_pc;

	#side1 or side2 mapped.
	my $side1_good =  $stats->{'side1_match'};
	my $side1_mappc = ($side1_good / $peReads)*100;
	$side1_mappc=sprintf "%.1f", $side1_mappc;
	my $side2_good =  $stats->{'side2_match'};
	my $side2_mappc = ($side2_good / $peReads)*100;
	$side2_mappc=sprintf "%.1f", $side2_mappc;

	#neither side1 nor side1 mapped.
	my $noMap = $stats->{'noMap'};
	my $noMap_pc = ($noMap / $peReads)*100;
	$noMap_pc=sprintf "%.1f", $noMap_pc;

	#only 1 side mapped.
	my $single = $stats->{'single'};
	my $single_pc = ($single / $peReads)*100;
	$single_pc=sprintf "%.1f", $single_pc;

	#invalid 'C' pairs.
	my $invalid = $stats->{'same'}{'->.<-'}+$stats->{'same'}{'<-.->'};
	my $invalid_pc = ($invalid/$peReads)*100;
	$invalid_pc=sprintf "%.1f", $invalid_pc;

	#valid 'C' pairs.
	my $valid = $stats->{'valid'};
	my $valid_pc = ($valid/$peReads)*100;
	$valid_pc=sprintf "%.1f", $valid_pc;

	#impossible (error) 'C' pairs.
	my $error=$stats->{'same'}{'<-.<-'}+$stats->{'same'}{'->.->'};
	my $error_pc = ($error/$peReads)*100;
	$error_pc=sprintf "%.1f", $error_pc;

	#self circle ligation -> same fragment ligated to itself. 
	my $selfCircle = $stats->{'same'}{'<-.->'};
	my $selfCircle_pc = ($selfCircle/($peMapped+1))*100;
	$selfCircle_pc=sprintf "%.1f", $selfCircle_pc;
	#dangling end ligation -> an unligated fragment that was pulled down.
	my $danglingEnd = $stats->{'same'}{'->.<-'};
	my $danglingEnd_pc = ($danglingEnd/($peMapped+1))*100;
	$danglingEnd_pc=sprintf "%.1f", $danglingEnd_pc;

	close(INTERACTION);
	close(SINGLE);
	close(NOMAP);
	close(SELFCIRCLE);
	close(DANGLINGEND);
	close(ERROR);
	
	# EV: 2014-10-30 added -f option
	system("gzip -f $interactionFile");
	system("gzip -f $singleFile");
	system("gzip -f $selfCircleFile");
	system("gzip -f $danglingEndFile");
	system("gzip -f $errorFile");

	open(HEAD,">".$statFile.".head");
	print HEAD "## Mapping\n";
	print HEAD "side1Mapped\t".commify($side1_good)."\t($side1_mappc%)\n";
	print HEAD "side2Mapped\t".commify($side2_good)."\t($side2_mappc%)\n";
	print HEAD "noSideMapped\t".commify($noMap)."\t($noMap_pc%)\n";
	print HEAD "oneSideMapped\t".commify($single)."\t($single_pc%)\n";
	print HEAD "bothSideMapped\t".commify($peMapped)."\t($peMapped_pc%)\n";
	print HEAD "errorPairs\t".commify($error)."\t($error_pc%)\n";
	print HEAD "invalidPairs\t".commify($invalid)."\t($invalid_pc%)\n";
	print HEAD "validPairs\t".commify($valid)."\t($valid_pc%)\n";
	print HEAD "## Artifacts\n";
	print HEAD "selfCircle\t".commify($selfCircle)."\t($selfCircle_pc%)\n";
	print HEAD "danglingEnd\t".commify($danglingEnd)."\t($danglingEnd_pc%)\n";
	print HEAD "## Advanced\n";

	my $stat_pc=sprintf "%.1f", ($stats->{'same'}{'->.->'}/($error+$invalid))*100;
	print HEAD "same|->.->\t".commify($stats->{'same'}{'->.->'})."\t($stat_pc%)\n";
	$stat_pc=sprintf "%.1f", ($stats->{'same'}{'->.<-'}/($error+$invalid))*100;
	print HEAD "same|->.<-\t".commify($stats->{'same'}{'->.<-'})."\t($stat_pc%)\n";
	$stat_pc=sprintf "%.1f", ($stats->{'same'}{'<-.->'}/($error+$invalid))*100;
	print HEAD "same|<-.->\t".commify($stats->{'same'}{'<-.->'})."\t($stat_pc%)\n";
	$stat_pc=sprintf "%.1f", ($stats->{'same'}{'<-.<-'}/($error+$invalid))*100;
	print HEAD "same|<-.<-\t".commify($stats->{'same'}{'<-.<-'})."\t($stat_pc%)\n";

	$stat_pc=sprintf "%.1f", ($stats->{'different'}{'->.->'}/($valid))*100;
	print HEAD "different|->.->\t".commify($stats->{'different'}{'->.->'})."\t($stat_pc%)\n";
	$stat_pc=sprintf "%.1f", ($stats->{'different'}{'->.<-'}/($valid))*100;
	print HEAD "different|->.<-\t".commify($stats->{'different'}{'->.<-'})."\t($stat_pc%)\n";
	$stat_pc=sprintf "%.1f", ($stats->{'different'}{'<-.->'}/($valid))*100;
	print HEAD "different|<-.->\t".commify($stats->{'different'}{'<-.->'})."\t($stat_pc%)\n";
	$stat_pc=sprintf "%.1f", ($stats->{'different'}{'<-.<-'}/($valid))*100;
	print HEAD "different|<-.<-\t".commify($stats->{'different'}{'<-.<-'})."\t($stat_pc%)\n";
	close(HEAD);
}

my %options;
my $results = GetOptions( \%options,'fragmentFile_1|f1=s','fragmentFile_2|f2=s','mappedReadsFile_1|m1=s','mappedReadsFile_2|m2=s');

my ($fragmentFile_1,$fragmentFile_2,$mappedReadsFile_1,$mappedReadsFile_2);
($fragmentFile_1,$fragmentFile_2,$mappedReadsFile_1,$mappedReadsFile_2)=check_options( \%options );

print "\n";
print "fragmentFile_1\t$fragmentFile_1\n";
print "fragmentFile_2\t$fragmentFile_2\n";
print "mappedReadsFile_1\t$mappedReadsFile_1\n";
print "mappedReadsFile_2\t$mappedReadsFile_2\n";
print "\n";

my ($dnaSequence,$superRestrictionSite,$antiSuperRestrictionSite)=normalRS2Super("A^AGCTT");
print "$dnaSequence\t$superRestrictionSite\t$antiSuperRestrictionSite\n";

print "parsing fragment file for side1 ($fragmentFile_1)...\n";
my $fragments_1={};
my $fragmentSizes_1={};
($fragments_1,$fragmentSizes_1)=getFragments($fragmentFile_1);

print "parsing fragment file for side2 ($fragmentFile_2)...\n";
my $fragments_2={};
my $fragmentSizes_2={};
($fragments_2,$fragmentSizes_2)=getFragments($fragmentFile_2);

#print "$_ ${$fragmentSizes_2}{$_}\n" for (keys %{$fragmentSizes_2});

print "overlapping mapped reads ($mappedReadsFile_1 - $mappedReadsFile_2)...\n";
overlapMappedReads($mappedReadsFile_1,$mappedReadsFile_2,$fragments_1,$fragmentSizes_1,$fragments_2,$fragmentSizes_2);
