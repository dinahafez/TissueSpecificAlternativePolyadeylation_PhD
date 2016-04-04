#! /usr/bin/perl -w
use strict;
#use Switch;
use DBI;
use DBHandler;
 
my $UTRtail = 0;
my $samples = 10;

my @train=();
my @model =();

my %model_h;
my %train_h;
#getPositivePos();
#getPositiveExamples();

#less strict
#getNegativePos();

#more strict
#getNegativePosAfterStopCodon();

#getNegativeExamples();

#getSignificantSequences();
#writeFilesinExpectedFormat("Liver", 336, "neg");
#writeFilesinExpectedFormat("Brain", 233, "pos");
#getNegativePosAfterStopCodon();
#writeFilesinExpectedFormat("Brain", 233, "pos", "/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/SVM/Brain_positive.fa");
#writeFilesinExpectedFormat("Brain", 2330, "neg", "/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/SVM/Brain_negative.fa");

#writeFilesinExpectedFormat("Liver", 336, "pos", "/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/SVM/Liver_positive.fa");
#writeFilesinExpectedFormat("Liver", 3360, "neg", "/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/SVM/Liver_negative.fa");

getSignificantSequences_3CellTypes("../../analysis/tissueSpecific/Liver_Specific_3CellTypes_PAS_Original","../../analysis/tissueSpecific/SignificantPositionsIn_3CellTypes_Liver.fa","../../analysis/tissueSpecific/SignificantPositionsOut_3CellTypes_Liver.fa");
getSignificantSequences_3CellTypes("../../analysis/tissueSpecific/Brain_Specific_3CellTypes_PAS_Original","../../analysis/tissueSpecific/SignificantPositionsIn_3CellTypes_Brain.fa","../../analysis/tissueSpecific/SignificantPositionsOut_3CellTypes_Brain.fa");
getSignificantSequences_3CellTypes("../../analysis/tissueSpecific/Kidney_Specific_3CellTypes_PAS_Original","../../analysis/tissueSpecific/SignificantPositionsIn_3CellTypes_Kidney.fa","../../analysis/tissueSpecific/SignificantPositionsOut_3CellTypes_Kidney.fa");

sub getSignificantSequences_3CellTypes
{
	my ($significantFile,$faSeq_in, $faSeq_out) = @_;
	my $dbh = DBconnect();
	
	open IN, "<$significantFile" or die "Can not open file :$significantFile";
	
	#my $faSeq_Liver = "../../analysis/tissueSpecific/SignificantPositionsIn_Liver.fa";
	#my $seqOut_Liver = "../../analysis/tissueSpecific/SignificantPositionsOut_Liver.fa";
	#my $faSeq_Brain = "../../analysis/tissueSpecific/SignificantPositionsIn_Brain.fa";
	#my $seqOut_Brain = "../../analysis/tissueSpecific/SignificantPositionsOut_Brain.fa";
	
	open OUT, ">$faSeq_in" or die "Can not open file: $faSeq_in";
	#open OUT_Brain, ">$faSeq_Brain" or die "Can not open file: $faSeq_Brain";
	
	my $line = <IN>;
	while($line = <IN>)
	{
		chomp($line);
		my ($gene,$median,$pas,$dif) = split(/\,/, $line);
		my $sqlGene = "'$gene'";
		my @result = DBselectChrStrand($dbh, $sqlGene);
		my $chr = $result[0];
		my $strand = $result[1];
	
		
		my $start=0;
		my $end=0;
		if(($chr eq "") ||($chr eq " ")) 
		{
			 print ("not found in DB $gene \t");
		}
		else
		{
			if($strand eq "1")
			{	
				$start = $pas  - 101;
				$end = $pas + 100;
			}
			else
			{	
				$start = $pas  - 101;
				$end = $pas + 100;
			}	
		
		}
		 print OUT ("chr$chr:$start-$end\n");
	
	}
	
	my $cmd = "/nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/twoBitToFa /nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/hg19.2bit -seqList=$faSeq_in $faSeq_out";
	print ("$cmd\n");
	system ($cmd);
	
	
closeConnection ($dbh);          	
}

sub getSignificantSequences
{
	my $dbh = DBconnect();
	
	my $significantFile = "../../analysis/tissueSpecific/significant_final.csv";
	open IN, "<$significantFile" or die "Can not open file :$significantFile";
	
	my $faSeq_Liver = "../../analysis/tissueSpecific/SignificantPositionsIn_Liver.fa";
	my $seqOut_Liver = "../../analysis/tissueSpecific/SignificantPositionsOut_Liver.fa";
	my $faSeq_Brain = "../../analysis/tissueSpecific/SignificantPositionsIn_Brain.fa";
	my $seqOut_Brain = "../../analysis/tissueSpecific/SignificantPositionsOut_Brain.fa";
	
	open OUT_Liver, ">$faSeq_Liver" or die "Can not open file: $faSeq_Liver";
	open OUT_Brain, ">$faSeq_Brain" or die "Can not open file: $faSeq_Brain";
	
	my $line = <IN>;
	while($line = <IN>)
	{
		chomp($line);
		my ($gene,$pas,$dif,$liver,$brain) = split(/\,/, $line);
		my $sqlGene = "'$gene'";
		my @result = DBselectChrStrand($dbh, $sqlGene);
		my $chr = $result[0];
		my $strand = $result[1];
	
		
		my $start=0;
		my $end=0;
		if(($chr eq "") ||($chr eq " ")) 
		{
			 print ("not found in DB $gene \t");
		}
		else
		{
			#MIssing: get the original PAS not the grouped one???
			if($strand eq "1")
			{	
				$start = $pas  - 101;
				$end = $pas + 100;
			}
			else
			{	
				$start = $pas  - 101;
				$end = $pas + 100;
			}	
		
		
			#determine whether it is brain or liver specific
			if($liver > $brain)
			{
				print OUT_Liver ("chr$chr:$start-$end\n");
			}		
			elsif($brain > $liver)
			{
				print OUT_Brain ("chr$chr:$start-$end\n");
			}
			else
			{
				print ("no diff$gene, $pas\n");
			}
		}
	
	}
	
	my $cmd = "/nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/twoBitToFa /nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/hg19.2bit -seqList=$faSeq_Liver $seqOut_Liver";
	print ("$cmd\n");
	system ($cmd);
	
	my $cmd = "/nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/twoBitToFa /nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/hg19.2bit -seqList=$faSeq_Brain $seqOut_Brain";
	print ("$cmd\n");
	system ($cmd);
	
closeConnection ($dbh);          	
}

sub writeFilesinExpectedFormat
{
	my ($cell, $totalNoOfExample, $positiveOrNegative, $inFile) = @_;
#	my $in = "/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/SVM/SignificantPositionsOut_".$cell.".fa";
#	open IN, "<$in" or die "Can not open file :$in";
	
	my $out = "/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/SVM/".$cell."_".$positiveOrNegative.".fa";
	open OUT, ">$out" or die "Can not open file: $out";
	
	open IN, "<$inFile" or die "Can not open file :$inFile";
	
	my $Model = "/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/SVM/".$cell."_".$positiveOrNegative."_model.fa";
	open MODEL, ">$Model" or die "Can not open file: $Model";
	
	my $TrainTest = "/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/SVM/".$cell."_".$positiveOrNegative."_train_test.fa";
	open TRAIN, ">$TrainTest" or die "Can not open file: $TrainTest";
	
	getPartitionedSet($totalNoOfExample);
	
	my $labl;
	if($positiveOrNegative eq "pos")
	{
		$labl = "1";
	}
	else
	{
		$labl = "-1";
	}
	
	my $line;
	my $i=-1;
	while($line = <IN>)
	{
		chomp ($line);
		if($line =~ /^>/)
		{
			print OUT ("$line label=$labl\n");
			$i++;
		}
		else
		{
			print OUT ("$line\n");
		}
		
		#partition for model selection
		#if (grep /^$i/, @model)
		if(exists $model_h{$i})
		{
			if($line =~ /^>/)
			{
				print MODEL ("$line label=$labl\n");
			
			}
			else
			{
				print MODEL ("$line\n");
			}
		}	
		else
		{	
			if($line =~ /^>/)
			{
				print TRAIN ("$line label=$labl\n");
			
			}
			else
			{
				print TRAIN ("$line\n");
			}
		}
		
		
	}
	
	close IN;
	close OUT;
	close TRAIN;
	close MODEL;
}

sub getPartitionedSet
{
	my ($total) = @_;
	srand(time() ^ $$ ^ unpack "%32L*", `ps axww | gzip`);#set the random number seed for the rand operator
	
	
	my $model_total = int(0.2 * $total);
	my $train_total =  $total - $model_total;
	
	print ("$model_total\t$train_total\n");
	
	#my @ipartition = sample( $total); # random sampling
    
 	my @ipartition = sort {-1+int rand 3} (0..$total-1);
	
	%model_h=();
	%train_h = ();
    for(my $i=0; $i<$model_total; $i++)
    {
   		$model_h{$ipartition[$i]} = 1;
    }
    for(my $i=$model_total; $i<$total; $i++)
    {
  		$train_h{$ipartition[$i]} = 1;
    }
    
 
}

sub getNegativePosAfterStopCodon
 {
 	
 	my $dbh = DBconnect();
	
	my $significantFile = "../../analysis/tissueSpecific/significant_final.csv";
	open IN, "<$significantFile" or die "Can not open file :$significantFile";
	
	my $faSeq_Liver = "../../analysis/tissueSpecific/Liver_negative_in.fa";
	my $seqOut_Liver = "../../analysis/tissueSpecific/Liver_negative_out.fa";
	
	my $faSeq_Brain = "../../analysis/tissueSpecific/Brain_negative_in.fa";
	my $seqOut_Brain = "../../analysis/tissueSpecific/Brain_negative_out.fa";
	
	open OUT_Liver, ">$faSeq_Liver" or die "Can not open file: $faSeq_Liver";
	open OUT_Brain, ">$faSeq_Brain" or die "Can not open file: $faSeq_Brain";
	my $stop = 0;
 	my $empty =0;
 	my $rangeBad = 0;
	my $line = <IN>;
	
	while($line = <IN>)
	{
		chomp($line);
		
		my ($gene,$polya,$dif,$liver,$brain) = split(/\,/, $line);
		my $sqlGene = "'$gene'";
		my @result = DBselectChrStrand($dbh, $sqlGene);
		my $chr = $result[0];
		my $strand = $result[1];
	
		
		my $start=0;
		my $end=0;
		if(($chr eq "") ||($chr eq " ")) 
		{
			 print ("not found in DB $gene \t");
		}
		else
		{
			if($strand eq "1")
			{	
				$strand = '+';
			}
			else
			{	
				$strand = '-';
			}	
			
		}
	
 		
 		#select stop Codon position
 		my $correctStop = "true";
 	
 		my @rows_selected = DBselectUpdatedGeneStopCodonFromGene_max($dbh,$sqlGene,$strand);
 		my $stopCodon = $rows_selected[0];
 		if($stopCodon)
 		{
 			$stopCodon = $rows_selected[0];
 			
	 		#sanity check
	 		if($strand eq "+")
			{	
	 			if($stopCodon >=$polya)
	 			{
	 				@rows_selected = DBselectUpdatedGeneStopCodonFromGene_avg($dbh,$sqlGene,$strand);
 					$stopCodon = int($rows_selected[0]);
			 		if($stopCodon >=$polya)
			 		{
			 			@rows_selected = DBselectUpdatedGeneStopCodonFromGene_min($dbh,$sqlGene,$strand);
		 				$stopCodon = $rows_selected[0];
					 	if($stopCodon)
					 	{
					 		$stopCodon = $rows_selected[0];
					 			
						 	if($stopCodon >=$polya)
						 	{
						 		$correctStop = "false";
						 		print ("stop codon after polya $gene $stopCodon, $polya\n");
						 	}
								
						}
				 		
				 	}
						
				}
	 		}
		
			else
			{	
	 			if($stopCodon <=$polya)
	 			{
	 				@rows_selected = DBselectUpdatedGeneStopCodonFromGene_avg($dbh,$sqlGene,$strand);
 					$stopCodon = int($rows_selected[0]);
			 		if($stopCodon <=$polya)
			 		{
			 			@rows_selected = DBselectUpdatedGeneStopCodonFromGene_min($dbh,$sqlGene,$strand);
	 					$stopCodon = $rows_selected[0];
				 		if($stopCodon)
				 		{
				 			$stopCodon = $rows_selected[0];
				 			
					 		if($stopCodon <=$polya)
					 		{
					 			$correctStop = "false";
					 			print ("stop codon after polya $gene\n");
					 		}
							
						}
					}
	 			}
			}
		}
		else
		{
			print ("$gene does not exist\n");
		}
		
		srand(time() ^ $$ ^ unpack "%32L*", `ps axww | gzip`);#set the random number seed for the rand operator

		my $range;
		if($strand eq "+")
		{
			if($correctStop eq "true")		
			{
			 	$range =   $polya - 100 - $stopCodon +1;
			 	if($range < 100)
			 	{
			 		$rangeBad++;
		 			#In this case, we claim that the stop codon falls in polya -500
		 			$stopCodon = $polya - 500;
		 			$range =   $polya - 100 - $stopCodon +1;
				}
			}
			else
			{
				$stopCodon = $polya - 500;
				$range =   $polya - 100 - $stopCodon +1;
			}
 			#randomly pick 10 points
 			
		 	my $random_pos =0;
		 	for (my $i=0;$i<$samples;$i++)
		 	{ 
		 		$random_pos = int(rand($range)) + $stopCodon;
		 		my $upstream = $random_pos-101;
	 			my $downstream = $random_pos + 100;
	 			#determine whether it is brain or liver specific
				if($liver > $brain)
				{
					print OUT_Liver ("chr$chr:$upstream-$downstream\n");
				}		
				elsif($brain > $liver)
				{
					print OUT_Brain ("chr$chr:$upstream-$downstream\n");
				}
			
	 			
		 	}
		 }
		 else
		 {
		 	if($correctStop eq "true")		
			{
				$range =   $stopCodon - ($polya + 100) +1;
		 		if($range < 100)
		 		{
		 			#print ("range <100: stop = $stopCodon, poly: $polya \n");
		 			$rangeBad++;
		 			$stopCodon = $polya + 500;
		 			$range =   $stopCodon - ($polya + 100) +1;
		 		}
			}
			else
			{
				$stopCodon = $polya + 500;
		 		$range =   $stopCodon - ($polya + 100) +1;
			}
 			#randomly pick 10 points
		 	my $random_pos =0;
		 	for (my $i=0;$i<$samples;$i++)
		 	{ 
		 		$random_pos = int(rand($range)) + $polya + 100;
		 		my $upstream = $random_pos-101;
	 			my $downstream = $random_pos + 100;
	 			if($liver > $brain)
				{
					print OUT_Liver ("chr$chr:$upstream-$downstream\n");
				}		
				elsif($brain > $liver)
				{
					print OUT_Brain ("chr$chr:$upstream-$downstream\n");
				}
				
	 			
		 	}
		 }
	}
	closeConnection($dbh);
 	
 	my $cmd = "/nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/twoBitToFa /nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/hg19.2bit -seqList=$faSeq_Liver $seqOut_Liver";
	print ("$cmd\n");
	system ($cmd);

 	$cmd = "/nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/twoBitToFa /nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/hg19.2bit -seqList=$faSeq_Brain $seqOut_Brain";
	print ("$cmd\n");
	system ($cmd);
	
	printLabel($seqOut_Liver,"Liver");
	printLabel($seqOut_Brain,"Brain");
	
	print("number of genes that have stop > polya $stop\n");
	print ("number of gnees that do not have stop $rangeBad\n");
 }
 

sub printLabel
{
	my ($inFile, $cell)=@_; 
	open IN, "<$inFile" or die "Can not open file :$inFile";
	
	my $out = "../../analysis/tissueSpecific/".$cell."_negative.fa";	
	open OUT, ">$out" or die "Can not open file: $out";
	
	while( my $line = <IN>)
	{
		chomp ($line);
		if($line =~ /^>/)
		{
			print OUT ("$line label=-1\n");
		}
		else
		{
			print OUT ("$line\n");
		}
	}
		
}

#OLD CODE

sub getPositivePos

{
	my $constitutive = "../../analysis/all_Cells.cons";
	my $faSeq = "../../analysis/PolyaPositivePositionsIn.fa";
	my $seqOut = "../../analysis/PolyaPositivePositionsOut.fa";
	
	open IN, "<$constitutive" or die "Can not open file :$constitutive";
	open OUT, ">$faSeq" or die "Can not open file: $faSeq";
	
	while (my $line = <IN>){
		chomp ($line);
				
		my @values = split(',', $line);  	#cut on space or tab
		my $name = $values[0];	
		my $chr = $values[1];
		my $strand = $values[2];
		my $pos_5 = $values[3];
		my $pos_3 = $values[4];
		my $mode_count = $values[5];
		my $count = $values[6];
		
		my $start = 0;
		my $end = 0;
		if($strand eq "+")
		{	
			$start = $pos_3  - 101;
			$end = $pos_3 + 100;
		}
		else
		{	
			$start = $pos_5  - 101;
			$end = $pos_5 + 100;
		}		
		print OUT ("chr$chr:$start-$end\n");	
	}	
	my $cmd = "/nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/twoBitToFa /nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/hg19.2bit -seqList=$faSeq $seqOut";
	print ("$cmd\n");
	system ($cmd);
	
}

sub sample
{
	my ($n) = @_;
 my ( @nums, $iters, $i, $k, );

      #  seed the random number generator
    srand;  


      #  initialize with the identity permutation
    @nums = 1 .. $n;  

      #  create the permutation
    $iters = 12 * $n**3 * log($n) + 1;
    for ( $i = 1; $i <= $iters; $i++)
        {
           if (rand(1) <= .5)   #  Flip a coin, and if heads swap
                                # a random adjacent pair of elements.
           {
               $k = int( rand($n-1) );
               ( $nums [$k], $nums [$k + 1] ) = ( $nums [$k + 1], $nums [$k] );
           }
        }

    return @nums;
}


sub getPositiveExamples
 {
 	my $constitutive = "../../analysis/all_Cells.cons";
 	my $seqOut = "../../analysis/PolyaPositivePositionsOut.fa";
 	my $positiveExamples = "../../analysis/PolyaPositiveExamplesTrain2.fa";
	my $positiveExamples2 = "../../analysis/PolyaPositiveExamplesTest2.fa";
	my $positiveExamples3 = "../../analysis/PolyaPositiveExamplesModel2.fa";
 	
 	open IN_1, "<$constitutive" or die "Can not open file :$constitutive";
 	open IN_2, "<$seqOut" or die "Can not open file :$seqOut";
	open OUT, ">$positiveExamples" or die "Can not open file: $positiveExamples";
	open OUT_2, ">$positiveExamples2" or die "Can not open file: $positiveExamples2";
	open OUT_3, ">$positiveExamples3" or die "Can not open file: $positiveExamples3";
	
	@train=();
	@model =();
	%train_h=();
	%model_h=();
	my $totalNoOfExample = 6972;
	getPartitionedSet($totalNoOfExample);
	
	my $line_2 = <IN_2>;
	$line_2 =~/^>(\S+):(\d+)-(\d+)/ ||die $line_2;
	my $i=-1;
	while (my $line = <IN_1>){
		$i++;
		chomp ($line);
		my @values = split(',', $line);  	#cut on space or tab
		my $name = $values[0];	
		my $chr = $values[1];
		my $strand = $values[2];
		my $pos_5 = $values[3];
		my $pos_3 = $values[4];
		my $mode_count = $values[5];
		my $count = $values[6];
		
		#read the sequence
		my $seq="";
		while (($line_2 = <IN_2>) && ($line_2 !~/^>(\S+):(\d+)-(\d+)/ ))
		{
			$seq = $seq . $line_2;
		}
		
		my $polya= 0;
		if($strand eq "+")
		{	
			$polya = $pos_3 ;
		}
		else
		{	
			$polya = $pos_5 ;
		}	
		if (grep /^$i/, @model)
		{
			print OUT_3 (">chr$chr,$strand,$polya,$mode_count,$count label=1 $name\n");	
			print OUT_3 ($seq);	
		}	
		
		else
		{	
			print OUT (">chr$chr,$strand,$polya,$mode_count,$count label=1 $name\n");	
			print OUT ($seq);
		}
	}	
	close IN_1;
	close IN_2;
	close OUT;
	close OUT_3;
	close OUT_2;
 }

sub getNegativePosAfterStopCodon_old
 {
 	my $dbh = DBconnect();
 	my $constitutive = "../../analysis/all_Cells.cons";
 	my $faSeq = "../../analysis/PolyaNegativePositionsIn.fa";
	my $seqOut = "../../analysis/PolyaNegativePositionsOut.fa";
	
 	open IN_1, "<$constitutive" or die "Can not open file :$constitutive";	
 	open OUT, ">$faSeq" or die "Can not open file: $faSeq";
 	
 	my $stop = 0;
 	my $empty =0;
 	my $rangeBad = 0;
 	while (my $line = <IN_1>){
		chomp ($line);
		my @values = split(',', $line);  	#cut on space or tab
		my $name = "\'$values[0]\'";	
		my $chr = $values[1];
		my $strand = $values[2];
		my $pos_5 = $values[3];
		my $pos_3 = $values[4];
		my $mode_count = $values[5];
		my $count = $values[6];	
 	
 		#select start and end of the gene from the DB
 		#should use the polya site directly
 		my $polya;
 		if($strand eq "+")
		{	
			$polya = $pos_3 
		}
		else
		{	
			$polya = $pos_5;
		}
		
 		
 		#select stop Codon position
 		my @rows_selected = DBselectGeneStopCodon($dbh,$name,$strand);
 		my $stopCodon = $rows_selected[0];
 		my $correctStop = "true";
 		
		if(!$stopCodon)
 		{
 			@rows_selected = DBselectGeneStopCodonFromGene($dbh,$name,$strand);
 			$stopCodon = $rows_selected[0];	
 			if(!$stopCodon)
 			{
 				$correctStop = "false";
 			}
 		}
 		if($correctStop eq "true")		
 		{
 			$stopCodon = $rows_selected[0];
 			
	 		#sanity check
	 		if($strand eq "+")
			{	
	 			if($stopCodon >=$polya)
	 			{
	 				$correctStop = "false";
	 			}
			}
			else
			{
				if($stopCodon <=$polya)
	 			{
	 				$correctStop = "false";
	 			}
			}
		}
		
		srand(time() ^ $$ ^ unpack "%32L*", `ps axww | gzip`);#set the random number seed for the rand operator

		my $range;
		if($strand eq "+")
		{
			if($correctStop eq "true")		
			{
		 		$range =   $polya - 100 - $stopCodon +1;
		 		if($range < 100)
		 		{
		 			$rangeBad++;
		 			#In this case, we claim that the stop codon falls in polya -500
		 			$stopCodon = $polya - 500;
		 			$range =   $polya - 100 - $stopCodon +1;
		 		}
			}
			else
			{
				$stopCodon = $polya - 500;
		 		$range =   $polya - 100 - $stopCodon +1;
			}
 			#randomly pick 10 points
		 	my $random_pos =0;
		 	for (my $i=0;$i<$samples;$i++)
		 	{ 
		 		$random_pos = int(rand($range)) + $stopCodon;
		 		my $upstream = $random_pos-101;
	 			my $downstream = $random_pos + 100;
	 			print OUT ("chr$chr:$upstream-$downstream\n");
		 	}
		 }
		 else
		 {
		 	if($correctStop eq "true")		
			{
				$range =   $stopCodon - ($polya + 100) +1;
		 		if($range < 100)
		 		{
		 			#print ("range <100: stop = $stopCodon, poly: $polya \n");
		 			$rangeBad++;
		 			$stopCodon = $polya + 500;
		 			$range =   $stopCodon - ($polya + 100) +1;
		 		}
			}
			else
			{
				$stopCodon = $polya + 500;
		 		$range =   $stopCodon - ($polya + 100) +1;
			}
 			#randomly pick 10 points
		 	my $random_pos =0;
		 	for (my $i=0;$i<$samples;$i++)
		 	{ 
		 		$random_pos = int(rand($range)) + $polya + 100;
		 		my $upstream = $random_pos-101;
	 			my $downstream = $random_pos + 100;
	 			print OUT ("chr$chr:$upstream-$downstream\n");
		 	}
		 }
 		
 	}
 	
 	my $cmd = "/nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/twoBitToFa /nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/hg19.2bit -seqList=$faSeq $seqOut";
	print ("$cmd\n");
	system ($cmd);
	closeConnection($dbh);
	print("number of genes that have stop > polya $stop\n");
	print ("number of gnees that do not have stop $rangeBad\n");
 }
 

 sub getNegativePos
 {
 	my $dbh = DBconnect();
 	my $constitutive = "../../analysis/all_Cells.cons";
 	my $faSeq = "../../analysis/PolyaNegativePositionsIn.fa";
	my $seqOut = "../../analysis/PolyaNegativePositionsOut.fa";
	
 	open IN_1, "<$constitutive" or die "Can not open file :$constitutive";	
 	open OUT, ">$faSeq" or die "Can not open file: $faSeq";
 	while (my $line = <IN_1>){
		chomp ($line);
		my @values = split(',', $line);  	#cut on space or tab
		my $name = $values[0];	
		my $chr = $values[1];
		my $strand = $values[2];
		my $pos_5 = $values[3];
		my $pos_3 = $values[4];
		my $mode_count = $values[5];
		my $count = $values[6];	
 	
 		my $start =0;
 		my $end = 0;
 		#select start and end of the gene from the DB
 		$name = "\'$name\'";
 		my @rows_selected = DBselectGeneStartEnd($dbh, 'geneinfo',$name);
		foreach (@rows_selected)
 		{
 			if($_ eq "")
 			{
 				print ("Errorrr $name\n");
 			}
 			else
 			{
 				$start = $rows_selected[0];
 				$end = $rows_selected[1];
 				
 				#remove the last 400 nt
 				if($strand eq "+")
				{	
 					$end = $end -$UTRtail;
				}
				else
				{
					$start = $start - $UTRtail;
				}
 			}
 		}
 		
 		srand(time() ^ $$ ^ unpack "%32L*", `ps axww | gzip`);#set the random number seed for the rand operator
 		my $range =  $end - $start+1;
 		#randomly pick 10 points
 		for (my $i=0;$i<$samples;$i++)
 		{
 			my $random_pos = 0;
 			if($strand eq "+")
 			{	 
 				$random_pos = int(rand($range)) + $start;
 			}
 			else
 			{
 				$random_pos = int(rand($range)) + $end;
 			}
 			
 			my $upstream = $random_pos-101;
 			my $downstream = $random_pos + 100;
 			print OUT ("chr$chr:$upstream-$downstream\n");
 		}
 	}
 	
 	my $cmd = "/nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/twoBitToFa /nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/hg19.2bit -seqList=$faSeq $seqOut";
	print ("$cmd\n");
	system ($cmd);
	closeConnection($dbh);
 }
 
 sub getNegativeExamples
 {
 
 	my $constitutive = "../../analysis/all_Cells.cons";
 	my $seqOut = "../../analysis/PolyaNegativePositionsOut.fa";
 	my $negativeExamples = "../../analysis/PolyaNegativeExamplesTrain2.fa";
	my $negativeExamples2 = "../../analysis/PolyaNegativeExamplesTest2.fa";
	my $negativeExamples3 = "../../analysis/PolyaNegativeExamplesModel2.fa";
 	
 	open IN_1, "<$constitutive" or die "Can not open file :$constitutive";
 	open IN_2, "<$seqOut" or die "Can not open file :$seqOut";
	open OUT, ">$negativeExamples" or die "Can not open file: $negativeExamples";
	open OUT_2, ">$negativeExamples2" or die "Can not open file: $negativeExamples2";
	open OUT_3, ">$negativeExamples3" or die "Can not open file: $negativeExamples3";
	
	%train_h=();

	%model_h =();
	my $totalNoOfExample = (69720 * $samples);
	getPartitionedSet($totalNoOfExample);
	
	#my %model_h = map { $_->[0] => 1 } @model;
	#my %test_h = map { $_->[0] => 1 } @test;

	my $line_2 = <IN_2>;
	$line_2 =~/^>(\S+):(\d+)-(\d+)/ ||die $line_2;
	my $i=-1;
	while (my $line = <IN_1>){
	
		chomp ($line);
		my @values = split(',', $line);  	
		my $name = $values[0];	
		my $chr = $values[1];
		my $strand = $values[2];
		my $pos_5 = $values[3];
		my $pos_3 = $values[4];
		my $mode_count = $values[5];
		my $count = $values[6];
		
		for(my $j=0; $j<$samples; $j++)
		{
			$i++;
			chomp ($line_2);
			my @values2 = split('-', $line_2);  	
			my $randomPos = $values2[1]-100;	
			
			#read the sequence
			my $seq="";
			while (($line_2 = <IN_2>) && ($line_2 !~/^>(\S+):(\d+)-(\d+)/ ))
			{
				$seq = $seq . $line_2;
			}
			
			if( exists $model_h{$i})
			{
				print OUT_3 (">chr$chr,$strand,$randomPos label=-1 $name\n");
				print OUT_3 ($seq);	
			}	
		
			else
			{	
				print OUT (">chr$chr,$strand,$randomPos label=-1 $name\n");
				print OUT ($seq);
			}
			
			
			
		}
	}	
	close IN_1;
	close IN_2;
	close OUT;
	close OUT_2;
	close OUT_3;
 }
 
