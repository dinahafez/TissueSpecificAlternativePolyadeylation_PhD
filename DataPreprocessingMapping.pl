use strict;
use Switch;
use List::Util qw(sum);
use Statistics::R ;
use List::Util qw(shuffle);
use List::Util qw(sum);
use Math::Random qw(random_multinomial);
use Math::Random qw(random_permutation);
use List::MoreUtils qw(firstidx);
 use List::MoreUtils qw(indexes);
 use Statistics::Zscore qw();

 use POSIX;
 use Storable qw(dclone);


	
#files ending with filter1 are already filtered in this criterea: 
	#a. annotate 5' tags and considers those falling in 3'UTR or exon
	#b. have >=2   5' different tags
#After preprocessing, all filtered files are mapped back to the genome. After trying Bowtie and BWA, I decided to stick with BWA
#This step is done one the cluster with another multi-mapping script
#Mapped files are then used in here 

#1. Generate TSS for each of the cell lines
#2. Run Fseq for each cell, to get the cutoff, or consider the general one??? now using 12
#3. Cluster PAS into clusters, and separate different sites according to overlap and median representation
#3. Annotate each cell to get the corresponding gene. If it got annotated as intergenic, I can always check in filter1  to get the corresponding gene
#4. create a hash: Gene->CellLine->pos3,count  (Counts are normalized)
#This is the output of this script

my  %used_chr = (
        chr1 => 1,
        chr2 => 1,
        chr3 => 1,
        chr4 => 1,
        chr5 => 1,
        chr6 => 1,
        chr7 => 1,
        chr8 => 1,
        chr9 => 1,
        chr10 => 1,
        chr11 => 1,
        chr12 => 1,
        chr13 => 1,
        chr14 => 1,
        chr15 => 1,
        chr16 => 1,
        chr17 => 1,
        chr18 => 1,
        chr19 => 1,
        chr20 => 1,
        chr21 => 1,
        chr22 => 1,
        chrM => 1,
        chrX=>1,
        chrY=>1
    );
    
my %cell_Types = (
	brain => 1,
	kidney => 1,
	liver => 1,
	
);

my %dataset;
$dataset{'145'}->{'Kidney'} = '../../analysis/145_Kidney_hg19.uniq';
$dataset{'145'}->{'Liver'} = '../../analysis/145_Liver_hg19.uniq';
$dataset{'145'}->{'Brain'} = '../../analysis/Brain_PA_hg19.uniq';

my ( $brain_genes_sorted_ref,$brain_resid_sorted_ref ) = ();
my ( $kidney_genes_sorted_ref,$kidney_resid_sorted_ref ) = ();
my ( $liver_genes_sorted_ref,$liver_resid_sorted_ref ); 


my %cellTypes_hash;

$cellTypes_hash{'Liver'}->{'resid'} = \$liver_resid_sorted_ref;
$cellTypes_hash{'Liver'}->{'genes'} = \$liver_genes_sorted_ref;

$cellTypes_hash{'Kidney'}->{'resid'} = \$kidney_resid_sorted_ref;
$cellTypes_hash{'Kidney'}->{'genes'} = \$kidney_genes_sorted_ref;

$cellTypes_hash{'Brain'}->{'resid'} = \$brain_resid_sorted_ref;
$cellTypes_hash{'Brain'}->{'genes'} = \$brain_genes_sorted_ref;



my %PAS_modes =();

#functions to filter the data, and preprocessing
filterHapRandom();
FilterInternalPriming();
filter2Uniq5Read();

#Map to genome hg19
annotate();
CountUniq5prime();
count5primeTagsPerGene();

#Peak calling
fseq();

#Sampling 
sample_wig("145_Kidney_hg19.uniq");

PAS();
run_generate_tss();

#Annotation based on refSeq
annotatePAS();

#Cluster polyA sites
groupPASByGene();
getPASMedians();
FilterOnlyGenesInThree ();
AdjustedCount_review();


sub filterHapRandom
{
	print ("in hap\n");
	while (my ($lane, $cellHash) = each %dataset)
	{
		while (my ($cell, $cell_file) = each %$cellHash)
		{
			my $cmd = "grep -v random $cell_file > 1";
			system ($cmd);
			$cmd = "grep -v hap 1 > 2";
			system ($cmd);
			$cmd = "grep -v Un 2 > $cell_file";
			system($cmd);
			
		}
	}
}

sub FilterInternalPriming
{
	print ("filtereing ");
	while (my ($lane, $cellHash) = each %dataset)
	{
		while (my ($cell, $cell_file) = each %$cellHash)
		{
			my $seqIn = "$cell_file.positionIn.fa";
			my $seqOut = "$cell_file.positionOut.fa";
	
	print("$seqIn");
			open IN, "<$cell_file" or die "Can not open file :$cell_file";
			open OUT, ">$seqIn" or die "Can not open file: $seqIn";
	
			while (my $line = <IN>){
				chomp ($line);
				$line =~ s/^\s+//; #trim leading space
				
				my ($uniq_count,$chr, $strand,$pos_5, $pos_3) = split(/\s+/, $line);
				
				my $pos_25 = 0;
				
				if($strand eq "+")
				{	$pos_25 = $pos_3 + 26;
					if(($chr eq "chrM") && ($pos_25 > 16571)) {$pos_25 = 16571;}
					if(($chr eq "chrM") && ($pos_3 > 16570)) {$pos_3 = 16570;}
					print OUT ("$chr:$pos_3-$pos_25\n");
				}
				else
				{	$pos_25 = $pos_3-26;
					if($pos_25<0){$pos_25=0;}
					if(($chr eq "chrM") && ($pos_25 > 16570)) {$pos_25 = 16570;}
					if(($chr eq "chrM") && ($pos_3 > 16571)) {$pos_3 = 16571;}
					print OUT ("$chr:$pos_25-$pos_3\n");
				}		
				
			}	
			my $cmd = "/nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/twoBitToFa /nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/hg19.2bit -seqList=$seqIn $seqOut";
			#print ("$cmd\n");
			system ($cmd);
			
			close IN;
			close OUT;
			
			RemoveSeq($cell_file);
		}
	}
}

sub RemoveSeq
{
	
	my ($cell_file) = @_;	
	my $seqOut = "$cell_file.positionOut.fa";
	
	#print $cell_file;
	open INseq, "<$seqOut" or die "Can not open file :$seqOut";
	open INtag, "<$cell_file" or die "Can not open file :$cell_file";
	
	my $count13=0;
	my $count=0;
	my $count14=0;
	my $count16=0;
	my $filtered_InternalPriming = $cell_file;
	my $tissuePath = "analysis/tissueSpecific";
	$filtered_InternalPriming =~ s/analysis/$tissuePath/g;
	$filtered_InternalPriming =~ s/uniq/filtered.InternalPriming/g;
	open OUT, ">$filtered_InternalPriming" or die "Can not open file: $filtered_InternalPriming";
	
	my $count=0;
	
	while (my $line = <INseq> ){
		my $seq = <INseq>;
		my $tag = <INtag>;
		#print ("Tag : $tag\t");
		chomp ($seq);	
	
		if($seq =~ m/AAAAAAAAAAAAA/i)
		{
			$count ++;
		}
		else
		{print OUT ($tag);}
	}
	print ("$filtered_InternalPriming : internal priming filtered out $count\n");
}

sub filter2Uniq5Read
{
	
	while (my ($lane, $cellHash) = each %dataset)
	{
		while (my ($cell, $cell_file) = each %$cellHash)
		{
			my %reads=();
			my $filtered_InternalPriming = $cell_file;
			my $tissuePath = "analysis/tissueSpecific";
			$filtered_InternalPriming =~ s/analysis/$tissuePath/g;
			$filtered_InternalPriming =~ s/uniq/filtered.InternalPriming/g;
			open IN, "<$filtered_InternalPriming" or die "Can not open file :$filtered_InternalPriming";
			
			my $filteredUniq5 = $filtered_InternalPriming;
			$filteredUniq5 =~ s/filtered.InternalPriming/filtered.PrimingUniq5/;
			open OUT, ">$filteredUniq5" or die "Can not open file: $filteredUniq5";
			my $line;
			while ( $line = <IN>)
			{
					chomp ($line);
					$line =~ s/^\s+//; #trim leading space
				
					my ($uniq_count,$chr, $strand,$pos_5, $pos_3) = split(/\s+/, $line);
			 		$reads{$chr}->{$strand}->{$pos_3}->{"count5"} ++;
			 		$reads{$chr}->{$strand}->{$pos_3}->{"noOfReads"} += $uniq_count;
			 		push(@{$reads{$chr}->{$strand}->{$pos_3}->{'pos_5'}}, $pos_5=>$uniq_count);
			}
			close IN;
			
			while (my ($chr, $strand_hash) = each %reads) 
			{	while (my ($strand, $pos_hash) = each %$strand_hash)
				{
					while(my ($pos,$hash) = each %$pos_hash)
					{
						if($reads{$chr}->{$strand}->{$pos}->{"count5"} >=2)
						{
							#chr,strand,pos_5,pos_3,uniq_count
							my @array = @{$reads{$chr}->{$strand}->{$pos}->{'pos_5'}};
							 for (my $i=0; $i<@array; $i+=2)
							 {
							 	print OUT ("$chr,$strand,$array[$i],$pos,$array[$i+1]\n");
							 
							 }
							
						}
					}
				}
			}
			close OUT;
			
	    	
		}
	}
}
sub annotate
{
	#Folder /nfs/labs/ohlerlab/sata/data/dina/PA_all_2/analysis/tissueSpecific/
	my %all_cells;
	
	while (my ($lane, $cellHash) = each %dataset)
	{
		while (my ($cell, $cell_file) = each %$cellHash)
		{
			my $filteredUniq5 = $cell_file;
			my $tissuePath = "analysis/tissueSpecific";
			$filteredUniq5 =~ s/analysis/$tissuePath/g;
			$filteredUniq5 =~ s/uniq/filtered.PrimingUniq5/g;
			
			#prepare the input file for annotator chr strand start end (take care for the -ve strand)
			my $annotateInput = $filteredUniq5;
			$annotateInput =~ s/filtered.PrimingUniq5/annotateInput/;
			 
			
			open IN, "<$filteredUniq5" or die "Can not open file :$filteredUniq5";
       		open OUT, ">$annotateInput" or die "Can not open file: $annotateInput";
			my $line;
			while ( $line = <IN>)
			{
					chomp ($line);
					$line =~ s/^\s+//; #trim leading space
				
					my ($chr, $strand,$pos_5, $pos_3,$uniq_count) = split(/\,/, $line);
			 		
		 			#if($strand eq "+")
			 		#{
			 			print OUT ("$chr,$strand,$pos_5,$pos_5,$pos_3,$uniq_count\n");
			 		#}
			 		#else
			 		#{
			 		#	print OUT ("$chr,$strand,$pos_3,$pos_5,$uniq_count\n"); #coz pos_3 is smaller
			 		#}
			 		
			}
			close IN;
			close OUT;
			
			my $sortedFile = $annotateInput."sorted";
			#sort according to the first position
			my $cmd = "sort -t, -k1r,1 -k 2r,2 -k 3n,3 $annotateInput  > $sortedFile";
    		system($cmd);
			
			#call the annotator script	
			my $annotateOutput = $annotateInput;
			$annotateOutput =~ s/annotateInput/annotateOutput/;
    		my $profile = "my_profile.yaml";
    		
    		$cmd = "annotator -s $sortedFile -b hg19 -p $profile -o $annotateOutput";
    		print($cmd."\n");
    		system($cmd);
		}
	}
}

sub CountUniq5prime
{

	my $bad = "../../analysis/tissueSpecific/badENST";
	open OUT_BAD, ">$bad" or die "Can not open file: $bad";
	while (my ($lane, $cellHash) = each %dataset)
	{
		while (my ($cell, $cell_file) = each %$cellHash)
		{
			my $annotatorFile = $cell_file;
			my $tissuePath = "analysis/tissueSpecific";
			$annotatorFile =~ s/analysis/$tissuePath/g;
			$annotatorFile =~ s/uniq/annotateOut/g;
	
			print ("$annotatorFile\n");
			open IN_2, "<$annotatorFile" or die "Can not open file :$annotatorFile";
   
		    my $line = <IN_2>;
		    my %Genes;
			my $countOk=0;
			my %GenesTotalCount;
			
    		while($line = <IN_2>)
		    {
		    	chomp ($line);	
		    	#we only consider the line that has a transcriptID
		    	my @values = split(/\s+/, $line);  	#cut on space or tab	
		    	if(($line =~ /ENST/)&&($values[4] !~ /5\'utr/)) #check that the annotation is not 5'utr
		    	{	
					my $chr = $values[0];
					my $strand = $values[1];
					my $count_5_3 = $values[@values-1];
					my $pos_5 = $values[2];
					my $pos_3 = $values[3];
					my $annotation = $values[4];
					my $ENST = $values[5];
					my $geneID = $values[@values-3];
					my $geneName = $values[@values-2];
				
					my $curr_5;
					$curr_5 = $values [2];
					
						
					#neglect reads that fall upstream of 3'utr (according to strand)
					
					my $transcriptID = $values[5];
					if(($annotation =~ /3\'utr/) or ($annotation =~ /coding/))
					{
						if (grep {$_ eq $curr_5} @{$Genes{$geneID}{'5pos'}}){}
						else
						{
							push (@{$Genes{$geneID}{'5pos'}}, $curr_5);
						}
						$GenesTotalCount{$geneID}{'count'} += $count_5_3;
						
					}
					else
					{	
						if((($strand eq '+') && ($transcriptID =~ /\[\d*\]/)) || (($strand eq '-')&& ($transcriptID =~ /\[-\d*\]/)))
						{
							print OUT_BAD ("$chr, $strand, $annotation, $pos_5, $pos_3, $transcriptID	\n");
						}
						else
						{
							if (grep {$_ eq $curr_5} @{$Genes{$geneID}{'5pos'}}){}
							else
							{
								push (@{$Genes{$geneID}{'5pos'}}, $curr_5);
							}
							$GenesTotalCount{$geneID}{'count'} += $count_5_3;
						}
					}
		    	}
		    }#end while
		  
		    my $Count5file = $annotatorFile;
    		$Count5file =~ s/annotateOut/count5prime/g;
    		open OUT, ">$Count5file" or die "Can not open file: $Count5file";
    		print OUT ("Gene,NoOfUniq5primeTags\n");
    		foreach my $gene (keys %Genes)
    		{
    			my $size = scalar (@{$Genes{$gene}{'5pos'}});
    			print OUT ("$gene,$size\n");
    		}
    		my $Totalfile = $annotatorFile;
    		$Totalfile =~ s/annotateOut/totalCount/g;
    		open OUT_2, ">$Totalfile" or die "Can not open file: $Totalfile";
    		print OUT_2 ("Gene,TotalNoOfReads\n");
    		foreach my $gene (keys %GenesTotalCount)
    		{
    			my $size = scalar (@{$Genes{$gene}{'5pos'}});
    			print OUT_2 ("$gene,$GenesTotalCount{$gene}{'count'}\n");
    		}
    		close OUT_2; 
    		close OUT; 
    		close IN_2;
    		close OUT_BAD;
		}
	}

    
}	

sub count5primeTagsPerGene
{
		while (my ($lane, $cellHash) = each %dataset)
		{
			while (my ($cell, $cell_file) = each %$cellHash)
			{
				my $Count5file = $cell_file;
				my $tissuePath = "analysis/tissueSpecific";
				$Count5file =~ s/analysis/$tissuePath/g;
				$Count5file =~ s/uniq/totalCount/g;
				open IN, "<$Count5file" or die "Can not open file :$Count5file";
				my %geneCount;
				#my $line = <IN>;
				while(my $line = <IN>)
			    {
			    	chomp ($line);					
					my ($gene, $count) = split(/\,/, $line); 
					 		
					$geneCount{$count}{'number'} ++;
			    }
			    my $Count5perGene = $cell_file;
				my $tissuePath = "analysis/tissueSpecific";
				$Count5perGene =~ s/analysis/$tissuePath/g;
				$Count5perGene =~ s/uniq/TotalCountPerGene/g;
			    open OUT, ">$Count5perGene" or die "Can not open file: $Count5perGene";
	    		print OUT ("noOfUniq5primeReads,noOfGenes\n");
	    		
	    		foreach (sort { $a <=> $b } keys(%geneCount) )
	      		
	    		{
	    			print OUT ("$_,$geneCount{$_}{'number'}\n");
	    		}
	    		close OUT;
	    		close IN;
	    		
	    		###
	    		my $Count5file = $cell_file;
				my $tissuePath = "analysis/tissueSpecific";
				$Count5file =~ s/analysis/$tissuePath/g;
				$Count5file =~ s/uniq/count5prime2/g;
				open IN, "<$Count5file" or die "Can not open file :$Count5file";
				my %geneCount;
				#my $line = <IN>;
				while(my $line = <IN>)
			    {
			    	chomp ($line);					
					my ($gene, $count) = split(/\,/, $line); 
					 		
					$geneCount{$count}{'number'} ++;
			    }
			    my $Count5perGene = $cell_file;
				my $tissuePath = "analysis/tissueSpecific";
				$Count5perGene =~ s/analysis/$tissuePath/g;
				$Count5perGene =~ s/uniq/Uniq5CountPerGene/g;
			    open OUT, ">$Count5perGene" or die "Can not open file: $Count5perGene";
	    		print OUT ("noOfUniq5primeReads,noOfGenes\n");
	    		
	    		foreach (sort { $a <=> $b } keys(%geneCount) )
	      		
	    		{
	    			print OUT ("$_,$geneCount{$_}{'number'}\n");
	    		}
	    		close OUT;
	    		close IN;
			}
		}
}

sub fseq
{
	#THE FOR LOOP SHOULD BE REMOVED!!!!
	my %all_cells;
	my $FseqCutoff = "../../analysis/tissueSpecific/fseq_cutoff";
	open FSEQ, ">$FseqCutoff" or die "Can not open file: $FseqCutoff";
	
	my $badfile = "../../analysis/bad";
	open BAD , ">$badfile" or die "Can not open file: $badfile";
	
	while (my ($lane, $cellHash) = each %dataset)
	{
		while (my ($cell, $cell_file) = each %$cellHash)
		{
			my $annotatorFile = $cell_file;
			my $tissuePath = "analysis/tissueSpecific";
			$annotatorFile =~ s/analysis/$tissuePath/g;
			$annotatorFile =~ s/uniq/annotateOutput/g;
	
			print ("$annotatorFile\n");
			open IN_2, "<$annotatorFile" or die "Can not open file :$annotatorFile";
   
   			my $fseqInput = $annotatorFile;
			$fseqInput =~ s/annotateOutput/fseq/;
			open OUT_P, ">$fseqInput" or die "Can not open file: $fseqInput";
		    
		    my $line = <IN_2>;
		    my %Genes;
			my $countOk=0;
			my %GenesTotalCount;
			
    		while($line = <IN_2>)
		    {
		    	chomp ($line);	
		    	#we only consider the line that has a geneID
		    	my @values = split(/\s+/, $line);  	#cut on space or tab	
		    		my $chr = $values[0];
					my $strand = $values[1];
					
					my $pos_5 = $values[2];
					#my $pos_3 = $values[3];
					my $annotation = $values[4];
					my $trailing = $values[@values-1];
					my($pos_3,$count_5_3) = split(/\,/, $trailing);
					
					switch ($annotation) 
					{
						case "intergenic" 
						{
							# if 7 elements:tes only
							#if 8 elements: utr and tes, consider tes
							if((@values ==7)&&($values[5] =~ /^ENST/))
							{
								if(($strand eq "+")&&($values[5]=~/\[\-/))
								{
									#polya site
							#		for(my $i=0;$i<$count_5_3;$i++)
									{
									print OUT_P ("chr$chr+\t$pos_3\t$pos_3\t3_utr\t100\t$strand\n");
									}
								}
								elsif(($strand eq "-")&&($values[5]=~/\[\d+\]/))
								{
									#polya site
									
								#	for(my $i=0;$i<$count_5_3;$i++)
									{
									print OUT_P ("chr$chr-\t$pos_3\t$pos_3\t3_utr\t100\t$strand\n");
									}
								}
								
							}
							elsif((@values ==8)&&($values[6] =~ /^ENST/))
							{
							#	for(my $i=0;$i<$count_5_3;$i++)
									{
								if(($strand eq "+")&&($values[6]=~/\[\-/))
								{
									#polya site
									print OUT_P ("chr$chr+\t$pos_3\t$pos_3\t3_utr\t100\t$strand\n");
								}
								elsif(($strand eq "-")&&($values[6]=~/\[\d+\]/))
								{
									#polya site
									print OUT_P ("chr$chr-\t$pos_3\t$pos_3\t3_utr\t100\t$strand\n");
								}
									}
							}
							else
							{
								print BAD ("$line\n");
							}
								

						}
						case "coding"
						{
							if(($values[5] =~ /^ENST/)|| ($values[5] =~/^\+ENST/))
							{
								#polya site
						#		for(my $i=0;$i<$count_5_3;$i++)
									{
								if($strand eq "+")
								{
									print OUT_P ("chr$chr+\t$pos_3\t$pos_3\t3_utr\t100\t$strand\n");
								}
								elsif($strand eq "-")
								{
									print OUT_P ("chr$chr-\t$pos_3\t$pos_3\t3_utr\t100\t$strand\n");
								}
									}
							}
							else
								{
									print BAD ("$line\n");
								}
						}
						case "3\'utr"
						{	
							if(($values[5] =~ /^ENST/)|| ($values[5] =~/^\+ENST/))
							{
								#polya site
				#	for(my $i=0;$i<$count_5_3;$i++)
									{
								if($strand eq "+")
								{
									print OUT_P ("chr$chr+\t$pos_3\t$pos_3\t3_utr\t100\t$strand\n");
								}
								elsif($strand eq "-")
								{
									print OUT_P ("chr$chr-\t$pos_3\t$pos_3\t3_utr\t100\t$strand\n");
								}
									
							}}
							else
							{
								print BAD ("$line\n");
							}
						}
						case "tes"
						{	
					#		for(my $i=0;$i<$count_5_3;$i++)
									{
							if($strand eq "+")
							{
								print OUT_P ("chr$chr+\t$pos_3\t$pos_3\t3_utr\t100\t$strand\n");
							}
							elsif($strand eq "-")
							{
								print OUT_P ("chr$chr-\t$pos_3\t$pos_3\t3_utr\t100\t$strand\n");
							}
									}
						}
						
						else
						{
							print BAD ("$line\n");
						}
					}
		    }#end while
					
			close IN;
			close OUT_P;
		
			
		my $cellType = $cell_file;
		my $path = '../../analysis/';
		$cellType =~ s/$path//g;
		
		my $cmd = "mkdir ../../analysis/tissueSpecific/$cellType.fseq";
		system($cmd);
		
		my $fseq = "../../analysis/tissueSpecific/$cellType.fseq";
		
		$cmd = "/home/dc97/biotools/fseq/bin/fseq -f 0 -l 30 -of wig -d ../../analysis/tissueSpecific -o $fseq $fseqInput";
	    print($cmd."\n");
	    system($cmd);
	   
#	    my $cutoff = sample_wig($cellType);
#	    print FSEQ  ("$cellType\t$cutoff\n");
	}
	}
	close FSEQ;
}

#This function computes the cutoff value for the fseq generated files
sub sample_wig {
	my ($cellType) = @_;
	print ("in sample wig $cellType\n");
	my $total =0;
	my $score=0;

	foreach my $chr (keys(%used_chr)) 
	{
		my $file;
		
		for( my $j = 0; $j <= 1; $j++ ) 
		{
			if( $j == 0 ) 
			{ $file = "../../analysis/tissueSpecific/$cellType.fseq/".$chr."+.wig"; 
		
			}
			else 
			{ $file = "../../analysis/tissueSpecific/$cellType.fseq/".$chr."-.wig";
			 }
			
			my $cmd = "touch $file";
			open IN, "<$file" ;
			print ("$file\t");
			while (<IN>) {
				chomp $_;
				if( /start=(\w+)/ ) {
					my $location = $1;
					my $start_location = $1;
					
				}
				elsif ((/\d/)&&($_>0.0)) { $score +=$_; $total++;} #$val{$location} = $_; $location++; print("$location\t");}	
			}
			close IN;
			
		}	
	}
	my $cut = $score / $total;

	print ("\nscore / total No of wigs: $cut\n");
	return $cut;
}

sub run_generate_tss #davis's java code, has to be run from inside preProcess
{#Input file chromosome,strand,3'readLocation,readCount
        while (my ($lane, $cellHash) = each %dataset)
		{
			while (my ($cell, $cell_file) = each %$cellHash)
			{
	        	my $annotatorFile = $cell_file;
				my $tissuePath = "analysis/tissueSpecific";
				$annotatorFile =~ s/analysis/$tissuePath/g;
				$annotatorFile =~ s/uniq/annotateOutput/g;
		
				print ("$annotatorFile\n");
				open IN_2, "<$annotatorFile" or die "Can not open file :$annotatorFile";
	   
	   			 my $PAS_in = $cell_file;
	             $PAS_in =~ s/analysis/$tissuePath/g;
	             $PAS_in =~ s/uniq/PASinput/g;
	             open OUT, ">$PAS_in" or die "Can not open file: $PAS_in";
			    
			    my $line = <IN_2>;
			    my %Genes;
				my $countOk=0;
				my %GenesTotalCount;
				
				my %reads;
	    		while($line = <IN_2>)
			    {
			    	chomp ($line);	
			    	#we only consider the line that has a geneID
			    	my @values = split(/\s+/, $line);  	#cut on space or tab	
			    	if(($line =~ /ENSG/)&&($values[4] !~ /5\'utr/)) #check that the annotation is not 5'utr
			    	{	
						my $chr = $values[0];
						my $strand = $values[1];
						my $trailing = $values[@values-1];
						my $pos_5 = $values[2];
						#my $pos_3 = $values[3];
						my $annotation = $values[4];
						my $ENST = $values[5];
						my $geneID = $values[@values-3];
						my $geneName = $values[@values-2];
						my ($pos_3, $count_5_3) = split(/\,/, $trailing);
					
						#neglect reads that fall upstream of 3'utr (according to strand)
						
						my $transcriptID = $values[5];
						my $utr_pos;
						if(($annotation =~ /3\'utr/) or ($annotation =~ /coding/))
						{
							#if($strand eq "+")
							#{
								$utr_pos = $pos_3;
							#}
							#else
							#{
							#	$utr_pos = $pos_5;
							#}
	                               
	                        $reads{$chr}{$strand}{$utr_pos} ++;
							
						}
						else
						{	
							if((($strand eq '+') && ($transcriptID =~ /\[\d*\]/)) || (($strand eq '-')&& ($transcriptID =~ /\[-\d*\]/)))
							{
								#print OUT_BAD ("$chr, $strand, $annotation, $pos_5, $pos_3, $transcriptID	\n");
							}
							else
							{
								if($strand eq "+")
								{	$utr_pos = $pos_3;
								}
								else
								{
									$utr_pos = $pos_3;
								}
		                               
		                        $reads{$chr}{$strand}{$utr_pos} ++;
							}
						}
			    	}
			    }#end while
						
	          
                        while (my ($chr, $strandHash) = each %reads)
                        {
                                while (my ($strand, $posHash) = each %$strandHash)
                                {
                                        while (my ($position, $count) = each %$posHash)
                                        {
                                                print OUT ("chr$chr,$strand,$position,$count\n");
                                        }
                                }
                        }
                        close IN;
						close OUT;
						
                        my $output = $PAS_in.".sorted";
                        my $cmd = "sort -t, -k1r,1 -k 2r,2 -k 3n,3 $PAS_in  > $output";
                		system($cmd);
                        $cmd = "mv $output $PAS_in";
                    	system($cmd);
                        
                        my $cellType = $cell_file;
						my $path = '../../analysis/';
						$cellType =~ s/$path//g;
		
                        
                        my $PAS_out =  $PAS_in;
                        $PAS_out =~ s/PASinput/PASoutput/g;
                        my $cmd = "java GenerateTSSs -s 0 -n 5 -w ../../../analysis/tissueSpecific/$cellType.fseq/ -c $PAS_in -o $PAS_out";
                        print($cmd."\n");
                		system($cmd);
                }
        }
}

sub annotatePAS
{   
        while (my ($lane, $cellHash) = each %dataset)
        {
                while (my ($cell, $cell_file) = each %$cellHash)
                {
      				my $PAS_file=  $cell_file;
      				my $tissuePath = "analysis/tissueSpecific";
      				$PAS_file =~ s/analysis/$tissuePath/g;
                    $PAS_file =~ s/uniq/PASoutput/g;
                  
                    open IN, "<$PAS_file" or die "Can not open file: $PAS_file";
                    
                    my $PAS_annotate_in = $PAS_file;
                    $PAS_annotate_in =~ s/PASoutput/PAS.annotateIn/g;
                 	
                    open OUT, ">$PAS_annotate_in" or die "Can not open file: $PAS_annotate_in";
                    
				 	while (my $line = <IN>)
				 	{
			          	chomp ($line);
			          	my ($chr, $strand, $PAS_start, $PAS_end, $count, $mode, $mode_count, $type)=  split(/\,/, $line);
			          	
			          	print OUT ("$chr,$strand,$PAS_start,$PAS_end,$count,$mode,$mode_count,$type\n");
				 	}
				 	my $PAS_annotate_out = $PAS_file;
				 	$PAS_annotate_out =~ s/PASoutput/PAS.annotateOut/g;
				 	 my $profile = "my_profile.yaml";
    		     print "the file is $PAS_annotate_out \n";
		    		my $cmd = "annotator -s $PAS_annotate_in -b hg19 -p $profile -o $PAS_annotate_out";
		    		print($cmd."\n");
		    		system($cmd);
                }
        }
                   
}

sub groupPASByGene
{
	my $dbh = DBconnect();
	 while (my ($lane, $cellHash) = each %dataset)
        {
             while (my ($cell, $cell_file) = each %$cellHash)
           	{
      				my $PAS_annotate_OUT=  $cell_file;
      				my $tissuePath = "analysis/tissueSpecific";
      				$PAS_annotate_OUT =~ s/analysis/$tissuePath/g;
                    $PAS_annotate_OUT =~ s/uniq/PAS.annotateOut/g;
                  
                    open IN, "<$PAS_annotate_OUT" or die "Can not open file: $PAS_annotate_OUT";
                    
                    my $PAS_Genes_file = $PAS_annotate_OUT;
                    $PAS_Genes_file =~ s/PAS.annotateOut/PAS.Genes/g;
                 
                    open OUT, ">$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                    my $count_noncoding = 0;
                    
                    my %PASs;
                    print ("$PAS_Genes_file\t");
                    my $line = <IN>;
				 	while ($line = <IN>)
				 	{
			          	chomp ($line);
			          	my @values=  split(/\s+/, $line);
			          	my $annotation = $values[4];
			          	my $gene_name;
			          	if($annotation  =~ /RNA/) #non coding
			          	{
			          		$count_noncoding ++;
			          	}
			          	else
			          	{
			          		$gene_name = $values[@values-2];
			          		#print ("$gene_name\t");
			          		if($gene_name =~ /ENST/) #look up the name in DB
			          		{
			          			my $ENST = $gene_name;
								$gene_name =~ s/\+//g;
								$gene_name =~ s/[[\d][\d]*[\]]//;
								$gene_name =~ s/\[//;
								$gene_name =~ s/\-//g;
								#print ("File:$gene_name\t");
			          			$gene_name = "\'$gene_name\'";
			          			my @result = DBselectGeneID($dbh, $gene_name);
			          			$gene_name = $result[0];
			          			#print ("DB:$gene_name\t");
			          			if(($gene_name eq "") ||($gene_name eq " ")) 
			          			{
			          				print ("not found in DB $ENST \t");
			          			}
			          		}
			          	
				          	my $trailingAnnotattion = $values[@values-1];
				          	my ($count, $mode, $modeCount, $peak) = split(/\,/, $trailingAnnotattion);
				          	if($peak =~ /NarrowPeak/) #select for NarrowPeak only
				          	{
				          		if(($gene_name ne "") &&($gene_name ne " ")) 
			          			{
					          		#Could be changed into $modeCount instead of $count
					          		$PASs{$gene_name}->{$mode}->{'PAS_count'}=$count;
					          		#$PASs{$gene_name}->{"totalCount"} += $count;
					          		$PASs{$gene_name}->{$mode}->{'PAS_start'} = $values[2];
					          		$PASs{$gene_name}->{$mode}->{'PAS_end'} = $values[3];
			          			}
				          	}     		
				         }
				 	}
			
			         my $noOfTotalGenes = keys( %PASs );
			         print ("total number of genes $noOfTotalGenes\t");
			         print OUT ("Gene,PAS,PAS_count,PAS_start,PAS_end\n");
			         while (my ($gene, $modeHash) = each %PASs)
        			{
        			#	print OUT ("$gene");
             			while (my ($mode, $hash) = each %$modeHash)
           				{
           					print OUT ("$gene,$mode,$PASs{$gene}->{$mode}->{'PAS_count'},$PASs{$gene}->{$mode}->{'PAS_start'},$PASs{$gene}->{$mode}->{'PAS_end'}\n");
           					#if($mode ne "totalCount")
           					#{
           						#my $vote = $PAS_count / $PASs{$gene}->{"totalCount"};
           						#my $observed_usage = $vote / $noOfTotalGenes;
							#	if(($mode ne "PAS_start")&&($mode ne "PAS_end"))
							#	print OUT ("$gene,$mode,$PAS_count,$total,\n");
           					#}
           				}
           			#	print OUT ("\n");
        			}
        			
        			close OUT;
           	#print some debuggin info
           	print ("No of non coding RNA $count_noncoding\n");    
           	
           	}
        }
	closeConnection ($dbh);          	
}


sub getPASMedians
{
		my $PAS_median_file ="../../analysis/tissueSpecific/PAS_median_3CellTypes";#PreFinal/LiverKidney/PAS_median";
	#my $median_usage_file ="median_usage";
	open OUT, ">$PAS_median_file" or die "Can not open file: $PAS_median_file";
	  my %PASs;
	 my $max=0;
	 my $max_variance =0;
	  my $m;
	  my $g;
		while (my ($lane, $cellHash) = each %dataset)
        {
             while (my ($cell, $cell_file) = each %$cellHash)
           	{
      				my $PAS_Genes_file=  $cell_file;
      				my $tissuePath = "analysis/tissueSpecific";
      				$PAS_Genes_file =~ s/analysis/$tissuePath/g;
                    $PAS_Genes_file =~ s/uniq/PAS.Genes/g;
                 
                    open IN, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                  
                    my $line = <IN>;
				 	while ($line = <IN>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 		
				 		if(exists $PASs{$gene})
				 		{
				 			my %genePASs = %{$PASs{$gene}};
				 			my $added= 'false';
				 	ADDED:{	for my $i ( keys %genePASs ) 
				 			{
				 				for (my $j=0;$j<@{$genePASs{$i}{'start'}}; $j++)
				 				{
				 					if(((($PAS_start >= $genePASs{$i}{'start'}[$j])&&($PAS_start <= $genePASs{$i}{'end'}[$j]))||(($PAS_end >= $genePASs{$i}{'start'}[$j])&&($PAS_end <= $genePASs{$i}{'end'}[$j])))||
				 					((($genePASs{$i}{'start'}[$j] >=$PAS_start)&&($genePASs{$i}{'start'}[$j] <=$PAS_end))||(($genePASs{$i}{'end'}[$j] >=$PAS_start)&&($genePASs{$i}{'end'}[$j] <= $PAS_end))))
				 					
				 					{
				 						
				 						push(@{$PASs{$gene}{$i}{'start'}}, $PAS_start);
				 						push(@{$PASs{$gene}{$i}{'end'}},$PAS_end);
				 						push(@{$PASs{$gene}{$i}{'mode'}}, $mode);
				 						$added='true';
				 						last ADDED;
				 						
				 					}
				 				}
				 			}
				 	}
				 			if($added eq 'false')
				 			{
				 				my $c = keys %genePASs;
				 				$c++;
				 				my $p = "PAS$c";
				 				push(@{$PASs{$gene}{$p}{'start'}}, $PAS_start);
				 				push(@{$PASs{$gene}{$p}{'end'}},$PAS_end);
				 				push(@{$PASs{$gene}{$p}{'mode'}}, $mode);
				 			}

				 		}
				 		else
				 		{
				 			#$PASs{$gene}{'count'}=0;
				 			push(@{$PASs{$gene}{'PAS1'}{'start'}}, $PAS_start);
				 			push(@{$PASs{$gene}{'PAS1'}{'end'}},$PAS_end);
				 			push(@{$PASs{$gene}{'PAS1'}{'mode'}}, $mode);
				 			
				 		}
				 	}
           	}
        }
					my $m2;
         	my $g2; 		
	while (my ($gene, $PASHash) = each %PASs)
	{
		my %genePASs = %{$PASs{$gene}};
		print OUT ("$gene");
		for my $i ( keys %genePASs )
		{
			my @start = @{$PASs{$gene}{$i}{'start'}};
			my @end = @{$PASs{$gene}{$i}{'end'}};
			my @mode = @{$PASs{$gene}{$i}{'mode'}};
			
			my @idxes =  sort { $start[$a] <=> $start[$b] } 0..$#end;
			my @start_sorted  = @start [ @idxes ];
			my @end_sorted = @end[ @idxes ];
			my @mode_sorted = @mode[@idxes];
			
			my $median;
         	my $mid=0;
         
         	my $total = @start;
         	#if($total %2 ==0)
         	{
         			$mid = int($total /2);
         		$median = $mode_sorted[$mid];	
       	 	}
       	 	#else
       	 	#{
       	 	#	$mid = int($total /2);
       	 	#	$median = int(($mode_sorted[$mid] + $mode_sorted[$mid - 1]) / 2);
        	# }
        	 print OUT (",$median");
        	 my $range = $end_sorted[@end_sorted-1] - $start_sorted[0];
        	 
        	 #calculate variance
        	  my $sum1=0;
        	  my $size = @mode_sorted;
        	 for(my $i=0; $i<@mode_sorted; $i++) 
        	 {
 				$sum1= $sum1+ $mode_sorted[$i];
        	 }
        	 my $mean = $sum1/$size;
        	 
        	 my $sum=0;
        	 for(my $i=0; $i<@mode_sorted; $i++) 
        	 {
 				$sum=$sum + (($mode_sorted[$i]-$mean)**2);
        	 }
        	 my $variance = $sum/$size;
        	 my $std_dev = sqrt($variance);
        	 
        	
        	 if($std_dev > $max_variance)
        	 {
        	 	$max_variance = $std_dev;
        	 	$g2 = $gene;
        	 	$m2 = $mode[0];
        	 }
        	 
        	 if($range>$max)
        	 {
        	 	$max = $range;
        	 	$m=$median;
        	 	$g = $gene;
        	 }
		}
		print OUT ("\n");
	}	 
print ("the max range is $g,$m,$max\n");
print("The max variance is $g2,$m2,$max_variance\n");
}

sub AdjustedCount_wrong
{
	
		while (my ($lane, $cellHash) = each %dataset)
        {
             while (my ($cell, $cell_file) = each %$cellHash)
           	{my %PASs =();
           		my $t=0;
           		my $total=0;
      				my $PAS_Genes_file=  $cell_file;
      				my $tissuePath = "analysis/tissueSpecific";
      				$PAS_Genes_file =~ s/analysis/$tissuePath/g;
                   $PAS_Genes_file =~ s/uniq/PAS.Genes/g;
                 # my  $PAS_Genes_file = "../../analysis/tissueSpecific/shannonin.txt";
                    open IN, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                  	
                  	my $PAS_Genes_Adjusted = "$PAS_Genes_file.Adjusted_2";
                  	open OUT, ">$PAS_Genes_Adjusted" or die "Can not open file: $PAS_Genes_Adjusted";
                  
                  	my $total_count =0;
                    my $line = <IN>;
				 	while ($line = <IN>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 		$PASs{$gene}+=$PAS_count;
						$total_count += $PAS_count;
						$total++;
				 	}
				 	close IN;
				 	
				 	open IN2, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                    $line = <IN2>;
				 	while ($line = <IN2>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 		
				 		my $adjusted_count = (($PAS_count +1)/($total+$total_count)) ;#* (10**8);
				 		
				 		print OUT ("$gene,$mode,$PAS_count,$adjusted_count,$PAS_start,$PAS_end\n");
				 		$t += $adjusted_count;
				 	}
         			close OUT;
         			close IN2;
         	
         print ("$t\n");
           	}
           	
        }           	

}

sub FilterOnlyGenesInThree
{
	my %PASs =();
	while (my ($lane, $cellHash) = each %dataset)
    {
         while (my ($cell, $cell_file) = each %$cellHash)
         {
       		my $PAS_Genes_file=  $cell_file;
      		my $tissuePath = "analysis/tissueSpecific";
      		$PAS_Genes_file =~ s/analysis/$tissuePath/g;
            $PAS_Genes_file =~ s/uniq/PAS.Genes/g;
                    #$PAS_Genes_file = "../../analysis/tissueSpecific/shannonin.txt";
            open IN, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
            
            my $line = <IN>;
			while ($line = <IN>)
			{
				chomp ($line);
				my ($gene,$mode,$PAS_count,$PAS_start,$PAS_end ) = split(/\,/, $line);
				$PASs{$cell}{$gene}{$mode}= $line;
			}
		}
    }
    #my $PASs_2 = dclone(\%PASs);
    while (my ($lane, $cellHash) = each %dataset)
    {
         while (my ($cell, $cell_file) = each %$cellHash)
         {
         
           		#my $t=0;
      		my $PAS_Genes_file_out=  $cell_file;
      		my $tissuePath = "analysis/tissueSpecific";
      		$PAS_Genes_file_out =~ s/analysis/$tissuePath/g;
            $PAS_Genes_file_out =~ s/uniq/PAS.Genes.ExpressedInThree/g;
                    #$PAS_Genes_file = "../../analysis/tissueSpecific/shannonin.txt";
            open OUT, ">$PAS_Genes_file_out" or die "Can not open file: $PAS_Genes_file_out";
            print OUT ("Gene,PAS,PAS_count,PAS_start,PAS_end\n");
            
            
         		while (my ($gene, $modehash) = each %{$PASs{$cell}})
         		{
         			if((exists $PASs{'Liver'}{$gene} )&&(exists $PASs{'Brain'}{$gene})&&(exists $PASs{'Kidney'}{$gene} )) #expressed in both tissues
        			{
	         			while (my ($mode,$line) = each %{$modehash})
	         			{
	         				print OUT ("$line\n");
	         			}
        			}
         		}
    		
		}
    }
    close OUT;
    #read median
    my $PAS_median_file ="../../analysis/tissueSpecific/PAS_median_3CellTypes";
		open IN, "<$PAS_median_file" or die "Can not open file: $PAS_median_file";
        my %PAS_median;
        my $numberOfPas =0;
		while (my $line = <IN>)
		{
			chomp ($line);
			my @values = split(/\,/, $line);;
			for(my $i=1;$i<@values;$i++)
			{
				$PAS_median{$values[0]}{$values[$i]}=1;
				
			}
			
		}
		close IN;
		
		my $PAS_median_out = "../../analysis/tissueSpecific/PAS_median_3CellTypes.ExpressedInThree";
		open OUT, ">$PAS_median_out" or die "Can not open file: $PAS_median_out";
        while (my ($key,$hash) = each %PAS_median)
        {
        	if((exists $PASs{'Liver'}{$key} )&&(exists $PASs{'Brain'}{$key})&&(exists $PASs{'Kidney'}{$key} )) #expressed in both tissues
        	{
        		print OUT ("$key");
        		while (my ($mode,$nothing) = each %$hash)
        		{
        			print OUT (",$mode");
        		}
        		print OUT ("\n");
        	}
        }
    
	close IN;
	close OUT;
}
sub AdjustedCount_review
{
	#read medians
		my $PAS_median_file ="../../analysis/tissueSpecific/PAS_median_3CellTypes.ExpressedInThree";
		open IN, "<$PAS_median_file" or die "Can not open file: $PAS_median_file";
        my %PAS_median;
        my $numberOfPas =0;
		while (my $line = <IN>)
		{
			chomp ($line);
			my @values = split(/\,/, $line);;
			for(my $i=1;$i<@values;$i++)
			{
				$PAS_median{$values[0]}{$values[$i]}=1;
				
			}
			
		}
		close IN;
		
		my %PASs =();
		my %PAS_modes =();
		my %PASs_allCells =();
		my %relative_exp=();
		my %PAS_all=();
		while (my ($lane, $cellHash) = each %dataset)
        {
             while (my ($cell, $cell_file) = each %$cellHash)
           	{
           			
      				my $PAS_Genes_file=  $cell_file;
      				my $tissuePath = "analysis/tissueSpecific";
      				$PAS_Genes_file =~ s/analysis/$tissuePath/g;
                    $PAS_Genes_file =~ s/uniq/PAS.Genes.ExpressedInThree/g;
                    open IN, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                  	my $line = <IN>;
				 	while ( $line = <IN>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 		#if in median +-10, refer to it as median, else, put it as it is
				 		
					 	Exist:
					 	{
					 		for(my $k=0;$k<=10;$k++)
					 		{
						 		if(exists $PAS_median{$gene}{$mode+$k})
						 		{
						 			$PASs{$cell}{$gene}{$mode+$k} =$PAS_count;
						 			$PASs_allCells{$gene}{$mode+$k}+=$PAS_count;
						 			$PAS_all{$gene}{$mode+$k}=1;
						 			if ( grep { $_ eq $mode+$k } @{$PAS_modes{$gene}} ) {
						 			}
						 			else
						 			{push(@{$PAS_modes{$gene}},$mode+$k);}
									last Exist;
						 		}
					 		}
					 		for(my $k=10;$k>0;$k--)
					 		{
						 		if(exists $PAS_median{$gene}{$mode-$k})
						 		{
						 			$PASs{$cell}{$gene}{$mode-$k} =$PAS_count;
						 			$PASs_allCells{$gene}{$mode-$k}+=$PAS_count;
						 			$PAS_all{$gene}{$mode-$k}=1;
						 			if ( grep { $_ eq $mode-$k } @{$PAS_modes{$gene}} ) {
						 			}
						 			else
									{	push(@{$PAS_modes{$gene}},$mode-$k);}
						 			
									last Exist;
						 		}
					 		}
					 	
				 			$PASs{$cell}{$gene}{$mode} =$PAS_count;
				 			$PASs_allCells{$gene}{$mode}+=$PAS_count;
				 			$PAS_all{$gene}{$mode}=1;
							if ( grep { $_ eq $mode } @{$PAS_modes{$gene}} ) {
						 	}
							else
				 			{push(@{$PAS_modes{$gene}},$mode);}
							
			
					 	}
				 	}
				 	close IN;
				 	
           	}
        }
        
		my $PAS_median2 = dclone(\%PAS_all);
		my $PAS_median3 = dclone(\%PAS_all);
	
		while (my ($key, $value) = each %$PAS_median3)
		{
			while (my ($v, $nothing) = each %$value)
			{
				$numberOfPas ++;
			}
			
		}
		
		while (my ($lane, $cellHash) = each %dataset)
        {
             while (my ($cell, $cell_file) = each %$cellHash)
           	{my %PASs =();
           		#my $t=0;
      				my $PAS_Genes_file=  $cell_file;
      				my $tissuePath = "analysis/tissueSpecific";
      				$PAS_Genes_file =~ s/analysis/$tissuePath/g;
                    $PAS_Genes_file =~ s/uniq/PAS.Genes.ExpressedInThree/g;
                     open IN, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                  	
                  	#Adjust gene level
                  	my $total_count =0;
                    my $line = <IN>;
				 	while ($line = <IN>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 		$PASs{$gene}{"PAScount"}+=$PAS_count;
				 		
						#$total_count += $PAS_count;
				 	}
				 	close IN;
				 	
				 	
				 	my $PAS_Genes_Adjusted_genelevel = "$PAS_Genes_file.Adjusted_genelevel";
                  	open OUT, ">$PAS_Genes_Adjusted_genelevel" or die "Can not open file: $PAS_Genes_Adjusted_genelevel";
				 	print OUT ("Gene,Mode,PAScount,PAscountAdjustedGeneLevel,PASstart,PASend\n");
				 	
				 	open IN2, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                    $line = <IN2>;
                    my %PAS_covered = ();
				 	while ($line = <IN2>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 		my $numberOfPasInGene = keys %{$PAS_all{$gene}};
				 		#my $adjusted_count = (($PAS_count)/($PASs{$gene}{"PAScount"}));
				 		#my $adjusted_count = (($PAS_count+1)/($PASs{$gene}{"PAScount"}+$numberOfPasInGene));
				 		my $adjusted_count = ($PAS_count+(1/$numberOfPasInGene))/($PASs{$gene}{"PAScount"}+1);
				 		
				 		print OUT ("$gene,$mode,$PAS_count,$adjusted_count,$PAS_start,$PAS_end\n");
				 		$PAS_covered{$gene}{$mode} =1;
				 	}
			
				 	
				 	
				 	while (my ($genex, $medianHash) = each %$PAS_median2)
	 				{
	 	 				while (my ($pasx,$nothing) = each %$medianHash)
	 	 				{
	 	 					my $ex = "f";
				 			#for(my $k=0;$k<=10;$k++)
					 		{
						 		if(exists $PAS_covered{$genex}{$pasx})
						 		{
				 					$ex = "t";
						 		}
					 		}
					 		#for(my $k=10;$k>0;$k--)
					 		{
						 		if(exists $PAS_covered{$genex}{$pasx})
						 		{
				 					$ex = "t";
						 		}
					 		}
					 		if($ex eq "f") 			
           					{
           						my $numberOfPasInGene = keys %{$PAS_all{$genex}};
           						
						 		my $adjusted_count = (1/$numberOfPasInGene)/($PASs{$genex}{"PAScount"}+1);
						 	#	print ("will adjust\n");
						 		print OUT ("$genex,$pasx,0,$adjusted_count,0,0\n");
						 		
           					}
           				}
        			}
         			close OUT;
         			close IN2;
         			
         			#Adjust library level
         			
         			open IN, "<$PAS_Genes_Adjusted_genelevel" or die "Can not open file: $PAS_Genes_Adjusted_genelevel";
         			$line = <IN>;
         			my %PASs_adjusted=();
         			
         			my $summation_count=0;
         			while ($line = <IN>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$PAS_count_adjusted_genelevel,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 		$summation_count +=$PAS_count_adjusted_genelevel;
				 		
				 	
				 	}
         			close IN;
         			open IN2, "<$PAS_Genes_Adjusted_genelevel" or die "Can not open file: $PAS_Genes_Adjusted_genelevel";
         			$line = <IN2>;
         			
         			my $PAS_Genes_Adjusted = "$PAS_Genes_file.Adjusted_Final";
                  	open OUT, ">$PAS_Genes_Adjusted" or die "Can not open file: $PAS_Genes_Adjusted";
				 	print OUT ("Gene,Mode,PAScount,PAscountAdjustedGeneLevel,PASAdjusted,PASstart,PASend\n");
				 	#	print ("$PAS_Genes_Adjusted\n");
         			while ($line = <IN2>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$PAS_count_adjusted_genelevel,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 	
         				
         				#print ("$summation_count\n");
				 		#my $adjusted_count = (($PAS_count_adjusted_genelevel+1)/($summation_count + $numberOfPas)); #has no effect
				 		my $adjusted_count = (($PAS_count_adjusted_genelevel)/($summation_count ));
				 		
				 		print OUT ("$gene,$mode,$PAS_count,$PAS_count_adjusted_genelevel,$adjusted_count,$PAS_start,$PAS_end\n");
				 	}
				 	close IN2;
				 	close OUT;
           	}
           	
        }           	

}
sub AdjustedCount_running
{
	
		while (my ($lane, $cellHash) = each %dataset)
        {
             while (my ($cell, $cell_file) = each %$cellHash)
           	{my %PASs =();
           		my $t=0;
      				my $PAS_Genes_file=  $cell_file;
      				my $tissuePath = "analysis/tissueSpecific";
      				$PAS_Genes_file =~ s/analysis/$tissuePath/g;
                    $PAS_Genes_file =~ s/uniq/PAS.Genes/g;
                    #$PAS_Genes_file = "../../analysis/tissueSpecific/shannonin.txt";
                    open IN, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                  	
                  	my $PAS_Genes_Adjusted = "$PAS_Genes_file.Adjusted";
                  	#open OUT, ">$PAS_Genes_Adjusted" or die "Can not open file: $PAS_Genes_Adjusted";
                  
                  	my $total_count =0;
                    my $line = <IN>;
				 	while ($line = <IN>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 		$PASs{$gene}+=$PAS_count;
						$total_count += $PAS_count;
				 	}
				 	close IN;
				 	
				 	open IN2, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                    $line = <IN2>;
				 	while ($line = <IN2>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 		
				 		my $adjusted_count = ($PAS_count/($PASs{$gene}+$total_count)) * (10**8);
				 		
				 		#print OUT ("$gene,$mode,$PAS_count,$adjusted_count,$PAS_start,$PAS_end\n");
				 		$t += $adjusted_count;
				 	}
         			close OUT;
         			close IN2;
         	
         print ("$cell\t$t\n");
           	}
           	
        }           	

}
sub NormalizeCount
{
	
		while (my ($lane, $cellHash) = each %dataset)
        {
             while (my ($cell, $cell_file) = each %$cellHash)
           	{my %PASs =();
      				my $PAS_Genes_file=  $cell_file;
      				my $tissuePath = "analysis/tissueSpecific";
      				$PAS_Genes_file =~ s/analysis/$tissuePath/g;
                    $PAS_Genes_file =~ s/uniq/PAS.Genes/g;
                    open IN, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                  	
                  	my $PAS_Genes_Adjusted = "$PAS_Genes_file.Normalized";
                  	open OUT, ">$PAS_Genes_Adjusted" or die "Can not open file: $PAS_Genes_Adjusted";
                  
                  	my $total_count =0;
                    my $line = <IN>;
				 	while ($line = <IN>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 		$PASs{$gene}+=$PAS_count;
						$total_count += $PAS_count;
				 	}
				 	close IN;
				 	
				 	open IN2, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                    $line = <IN2>;
				 	while ($line = <IN2>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 		
				 		my $adjusted_count = ($PAS_count/($PASs{$gene}+$total_count)) * (10**8);
				 		if($adjusted_count <1)
				 		{
				 			print ("$cell, $adjusted_count\t"); 
				 		}
				 		print OUT ("$gene,$mode,$PAS_count,$adjusted_count,$PAS_start,$PAS_end\n");
				 	}
         			close OUT;
         			close IN2;
         
         
           	}
        }           	

}

