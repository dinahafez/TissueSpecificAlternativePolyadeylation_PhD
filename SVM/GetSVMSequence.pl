#! /usr/local/bin/perl -w
use strict;
# use Switch;
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
#writeFilesinExpectedFormat2("Liver", 336, "neg");
%model_h = ();
%train_h = ();
#writeFilesinExpectedFormat2("Brain", 233, "pos");
#getNegativePosAfterStopCodon();
#writeFilesinExpectedFormat("Brain", 233, "pos", "/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/SVM_rc/SignificantPositionsOut_Brain2_rc.fa ");
#writeFilesinExpectedFormat("Brain", 2330, "neg", "/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/SVM_rc/Brain_negative_out_rc.fa");

#writeFilesinExpectedFormat("Liver", 336, "pos", "/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/SVM_rc/SignificantPositionsOut_Liver2_rc.fa");
#writeFilesinExpectedFormat("Liver", 3360, "neg", "/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/SVM_rc/Liver_negative_out_rc.fa");

#Trying to retrieve the original pas not the median and compare
#getSignificantSequencesForOriginalPAS("/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/Liver_Specific_PAS_Original", "../../analysis/tissueSpecific/SignificantPositionsIn_Liver2.fa", "../../analysis/tissueSpecific/SignificantPositionsOut_Liver2.fa");
#getSignificantSequencesForOriginalPAS("/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/Brain_Liver_common_PAS", "../../analysis/tissueSpecific/Brain_Liver_common_PositionsIn.fa", "../../analysis/tissueSpecific/Brain_Liver_common_out.fa");
#getSignificantSequencesForOriginalPAS("/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/Brain_Liver_constitutive_mapped_to_Brain_PAS", "../../analysis/tissueSpecific/Brain_Liver_costitutive_Brain_PositionsIn.fa", "../../analysis/tissueSpecific/Brain_Liver_constitutive_mapped_to_Brain_out.fa");
#getSignificantSequencesForOriginalPAS("/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/Brain_Liver_constitutive_mapped_to_Liver_PAS", "../../analysis/tissueSpecific/Brain_Liver_costitutive_Liver_PositionsIn.fa", "../../analysis/tissueSpecific/Brain_Liver_constitutive_mapped_to_Liver_out.fa");
#getSignificantSequencesForOriginalPAS_rc("/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/Brain_Liver_constitutive_mapped_to_Brain_PAS", "../../analysis/tissueSpecific/Brain_Liver_costitutive_Brain_PositionsIn.fa", "../../analysis/tissueSpecific/Brain_Liver_constitutive_mapped_to_Brain_out.fa","../../analysis/tissueSpecific/Brain_Liver_constitutive_mapped_to_Brain_out.strand" );
#getSignificantSequencesForOriginalPAS_rc("/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/Brain_Liver_constitutive_mapped_to_Liver_PAS", "../../analysis/tissueSpecific/Brain_Liver_costitutive_Liver_PositionsIn.fa", "../../analysis/tissueSpecific/Brain_Liver_constitutive_mapped_to_Liver_out.fa","../../analysis/tissueSpecific/Brain_Liver_constitutive_mapped_to_Liver_out.strand" );

#getSignificantSequencesForOriginalPAS_rc("/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/Liver_Specific_PAS_Original", "../../analysis/tissueSpecific/SignificantPositionsIn_Liver2.fa", "../../analysis/tissueSpecific/SignificantPositionsOut_Liver2.fa", "../../analysis/tissueSpecific/SignificantPositionsOut_Liver2.strand");
#getSignificantSequencesForOriginalPAS_rc("/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/PAS_Constitutive_mapped_closest", "/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/PAS_Constitutive_IN2.fa", "/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/PAS_Constitutive_out2.fa", "/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/PAS_Constitutive2.strand");


#getNegativePosAfterStopCodon_common_cons ("/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/Brain_Liver_common_PAS", "../../analysis/tissueSpecific/Brain_Liver_common_PositionsIn.fa","../../analysis/tissueSpecific/Brain_Liver_common_negative.fa", "Brain_Liver_common");
#getNegativePosAfterStopCodon_common_cons ("../../analysis/tissueSpecific/Brain_Liver_constitutive_mapped_to_Brain_PAS_allinDB", "../../analysis/tissueSpecific/Brain_Liver_cons_Brain_PositionsIn_allinDB.fa","../../analysis/tissueSpecific/Brain_Liver_constitutive_mapped_to_Brain_negative_allinDB.fa", "Brain_Liver_constitutive_mapped_to_Brain") ;
#getNegativePosAfterStopCodon_common_cons ("../../analysis/tissueSpecific/Brain_Liver_constitutive_mapped_to_Liver_PAS_allinDB", "../../analysis/tissueSpecific/Brain_Liver_cons_Liver_PositionsIn_allinDB.fa","../../analysis/tissueSpecific/Brain_Liver_constitutive_mapped_to_Liver_negative_allinDB.fa", "Brain_Liver_constitutive_mapped_to_Liver");

#reverseComplementNegativeGenes();

#writeFilesinExpectedFormat("Brain_Liver_common",1264,"pos","/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/SVM//Brain_Liver_common_out.fa");
#writeFilesinExpectedFormat("Brain_Liver_common",12640,"neg","/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/SVM//Brain_Liver_common_negative.fa",);
#writeFilesinExpectedFormat("Brain_Liver_cons_Brain",2945,"pos","../../analysis/tissueSpecific/Brain_Liver_constitutive_mapped_to_Brain_out.fa");
#writeFilesinExpectedFormat("Brain_Liver_cons_Brain",29450,"neg","../../analysis/tissueSpecific/Brain_Liver_constitutive_mapped_to_Brain_negative.fa",);
#writeFilesinExpectedFormat("Brain_Liver_cons_Liver",2945,"pos","../../analysis/tissueSpecific/Brain_Liver_constitutive_mapped_to_Liver_out.fa");
#writeFilesinExpectedFormat("Brain_Liver_cons_Liver",29450,"neg","../../analysis/tissueSpecific/Brain_Liver_constitutive_mapped_to_Liver_negative.fa",);

#writeFilesinExpectedFormat("Brain_Liver_cons_Brain",2945,"pos","/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/SVM_rc/Brain_Liver_constitutive_mapped_to_Brain_out_rc.fa");
#writeFilesinExpectedFormat("Brain_Liver_cons_Brain",29420,"neg","/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/SVM_rc/Brain_Liver_constitutive_mapped_to_Brain_out_rc.fa");

#getSignificantSequences_3CellTypes("../../analysis/tissueSpecific/Liver_Specific_3CellTypes_PAS_Original","../../analysis/tissueSpecific/SignificantPositionsIn_3CellTypes_Liver.fa","../../analysis/tissueSpecific/SignificantPositionsOut_3CellTypes_Liver.fa","../../analysis/tissueSpecific/SignificantPositionsOut_3CellTypes_Liver_rc.fa");
#getSignificantSequences_3CellTypes("../../analysis/tissueSpecific/Brain_Specific_3CellTypes_PAS_Original","../../analysis/tissueSpecific/SignificantPositionsIn_3CellTypes_Brain.fa","../../analysis/tissueSpecific/SignificantPositionsOut_3CellTypes_Brain.fa","../../analysis/tissueSpecific/SignificantPositionsOut_3CellTypes_Brain_rc.fa");
#getSignificantSequences_3CellTypes("../../analysis/tissueSpecific/Kidney_Specific_3CellTypes_PAS_Original","../../analysis/tissueSpecific/SignificantPositionsIn_3CellTypes_Kidney.fa","../../analysis/tissueSpecific/SignificantPositionsOut_3CellTypes_Kidney.fa","../../analysis/tissueSpecific/SignificantPositionsOut_3CellTypes_Kidney_rc.fa");

#getSignificantSequences_3CellTypes("../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Brain_Liver_Specific_3CellTypes_PAS_Original","../../analysis/tissueSpecific/HighlyExpressedInTwoCells/SignificantPositionsIn_3CellTypes_Brain_Liver.fa","../../analysis/tissueSpecific/HighlyExpressedInTwoCells/SignificantPositionsOut_3CellTypes_Brain_Liver.fa","../../analysis/tissueSpecific/HighlyExpressedInTwoCells/SignificantPositionsOut_3CellTypes_Brain_Liver_rc.fa");
#getSignificantSequences_3CellTypes("../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Brain_Kidney_Specific_3CellTypes_PAS_Original","../../analysis/tissueSpecific/HighlyExpressedInTwoCells/SignificantPositionsIn_3CellTypes_Brain_Kidney.fa","../../analysis/tissueSpecific/HighlyExpressedInTwoCells/SignificantPositionsOut_3CellTypes_Brain_Kidney.fa","../../analysis/tissueSpecific/HighlyExpressedInTwoCells/SignificantPositionsOut_3CellTypes_Brain_Kidney_rc.fa");
#getSignificantSequences_3CellTypes("../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Kidney_Liver_Specific_3CellTypes_PAS_Original","../../analysis/tissueSpecific/HighlyExpressedInTwoCells/SignificantPositionsIn_3CellTypes_Kidney_Liver.fa","../../analysis/tissueSpecific/HighlyExpressedInTwoCells/SignificantPositionsOut_3CellTypes_Kidney_Liver.fa","../../analysis/tissueSpecific/HighlyExpressedInTwoCells/SignificantPositionsOut_3CellTypes_Kidney_Liver_rc.fa");

#print ("KIDNEY\n");
#getSignificantSequences_3CellTypes("../../analysis/tissueSpecific/3CellTypes/constitutive_mapped_to_kidney_PAS","../../analysis/tissueSpecific/3CellTypes/constitutive_Kidney_In.fa","../../analysis/tissueSpecific/3CellTypes/constitutive_Kidney_out","../../analysis/tissueSpecific/3CellTypes/constitutive_Kidney_rc.fa");
#print ("BRAIN\n");
#getSignificantSequences_3CellTypes("../../analysis/tissueSpecific/3CellTypes/constitutive_mapped_to_Brain_PAS","../../analysis/tissueSpecific/3CellTypes/constitutive_Brain_In.fa","../../analysis/tissueSpecific/3CellTypes/constitutive_Brain_out","../../analysis/tissueSpecific/3CellTypes/constitutive_Brain_rc.fa");
#print ("LIVER\n");
#getSignificantSequences_3CellTypes("../../analysis/tissueSpecific/3CellTypes/constitutive_mapped_to_Liver_PAS","../../analysis/tissueSpecific/3CellTypes/constitutive_Liver_In.fa","../../analysis/tissueSpecific/3CellTypes/constitutive_Liver_out","../../analysis/tissueSpecific/3CellTypes/constitutive_Liver_rc.fa");

getNegativePosAfterStopCodon_common_cons2 ("/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/3CellTypes/constitutive_mapped_to_Brain_PAS", "../../analysis/tissueSpecific/3CellTypes/constitutive_mapped_to_Brain_PAS_PositionsIn.fa","../../analysis/tissueSpecific/3CellTypes/constitutive_mapped_to_Brain_negative.fa", "3Cells_Constitutive");
writeFilesinExpectedFormat("constitutive", 21340, "neg", "../../analysis/tissueSpecific/3CellTypes/constitutive_mapped_to_Brain_negative_rc.fa");


#writeFilesinExpectedFormat("Brain", 71, "pos", "../../../analysis/tissueSpecific/HighlyExpressedInOneCell/SignificantPositionsOut_3CellTypes_Brain_rc.fa");
#writeFilesinExpectedFormat("Liver", 69, "pos", "../../../analysis/tissueSpecific/HighlyExpressedInOneCell/SignificantPositionsOut_3CellTypes_Liver_rc.fa");
#writeFilesinExpectedFormat("Kidney", 52, "pos", "../../../analysis/tissueSpecific/HighlyExpressedInOneCell/SignificantPositionsOut_3CellTypes_Kidney_rc.fa");

#writeFilesinExpectedFormat("Brain_Liver", 89, "pos", "../../../analysis/tissueSpecific/HighlyExpressedInTwoCells/SignificantPositionsOut_3CellTypes_Brain_Liver_rc.fa");
#writeFilesinExpectedFormat("Brain_Kidney", 33, "pos", "../../../analysis/tissueSpecific/HighlyExpressedInTwoCells/SignificantPositionsOut_3CellTypes_Brain_Kidney_rc.fa");
#writeFilesinExpectedFormat("Kidney_Liver", 120, "pos", "../../../analysis/tissueSpecific/HighlyExpressedInTwoCells/SignificantPositionsOut_3CellTypes_Kidney_Liver_rc.fa");

#writeFilesinExpectedFormat("constitutive_Kidney", 2134, "neg", "../../analysis/tissueSpecific/3CellTypes/constitutive_Kidney_rc.fa");
#writeFilesinExpectedFormat("constitutive_Brain", 2134, "neg", "../../analysis/tissueSpecific/3CellTypes/constitutive_Brain_rc.fa");
#writeFilesinExpectedFormat("constitutive_Liver", 2134, "neg", "../../analysis/tissueSpecific/3CellTypes/constitutive_Liver_rc.fa");

#NEEL
#compareTranscriptLength("../../analysis/tissueSpecific/145_Liver_hg19.PAS.Genes","../../analysis/tissueSpecific/Brain_PA_hg19.PAS.Genes","../../analysis/tissueSpecific/145_Kidney_hg19.PAS.Genes","../../analysis/tissueSpecific/Kidney_Specific_3CellTypes_PAS_Original");  


	
sub getSignificantSequences_3CellTypes
{
	my ($significantFile,$faSeq_in, $faSeq_out,$faSeq_out_rc) = @_;
	my $dbh = DBconnect();
	
	open IN, "<$significantFile" or die "Can not open file :$significantFile";
	
	open OUT, ">$faSeq_in" or die "Can not open file: $faSeq_in";
	
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
			 print ("$gene not found in DB $gene \t");
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
	
	reverseComplementNegativeGenes_3CellTypes($significantFile, $faSeq_out, $faSeq_out_rc);
	
closeConnection ($dbh);          	
}

sub reverseComplementNegativeGenes_3CellTypes
{
	my ($significantFile,$seqOut, $seqOut_rc) = @_;
	my $dbh = DBconnect();
	
	open IN, "<$significantFile" or die "Can not open file :$significantFile";
		
	open IN_Liver, "<$seqOut" or die "Can not open file: $seqOut";
	
	open OUT_Liver, ">$seqOut_rc" or die "Can not open file: $seqOut_rc";
	
	my $seq = "";
	my $L_label = "";
	my $B_label="";
	my $stop = "false";
	my $L_l;
	my $B_l;		
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
			$seq = "";
			
			#read the sequence from either brain or liver
			if($L_label eq "")
			{
				$L_label  = <IN_Liver>; 
				chomp($L_label);
			}
			else
			{
				$L_label = $L_l;
			}
			
			LOOP:{
				while ($L_l = <IN_Liver>)
				{
					chomp($L_l);
					if($L_l =~ m/^\>/)
					{
						last LOOP;						
					}
					else
					{
					
						$seq = $seq . $L_l; 
					}
					
				}
			}
		}
		
		#check strand, reverse complement if -ve
		if($strand eq "1")
		{	
			print OUT_Liver ("$L_label\n$seq\n");
			
		}
		else
		{	
			# reverse the DNA sequence
	        my $revcomp = reverse($seq);
	
			# complement the reversed DNA sequence
	        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
			
			print OUT_Liver ("$L_label\n$revcomp\n");
			
		}
			
	}
	
	
	
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
			
			#if(($liver > $brain)&&($pas{$gene}{$pas}{'liver'} <$pas{$gene}{$pas}{'brain'}))
			#{
		#		print  ("CHECK $gene,$pas,$liver,$brain,$pas{$gene}{$pas}{'liver'},$pas{$gene}{$pas}{'brain'}\n");
		#	}		
		#	elsif(($brain > $liver)&&($pas{$gene}{$pas}{'liver'} >$pas{$gene}{$pas}{'brain'}))
		#	{
		#		print  ("CHECK $gene,$pas,$liver,$brain,$pas{$gene}{$pas}{'liver'},$pas{$gene}{$pas}{'brain'}\n");
		#	}
			
			
			
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

sub reverseComplementNegativeGenes
{
	my $dbh = DBconnect();
	my $significantFile = "../../analysis/tissueSpecific/significant_final.csv";
	open IN, "<$significantFile" or die "Can not open file :$significantFile";
	
	my $seqOut_Liver = "../../analysis/tissueSpecific/SignificantPositionsOut_Liver.fa";
	my $seqOut_Brain = "../../analysis/tissueSpecific/SignificantPositionsOut_Brain.fa";
	
	my $seqOut_Liver_rc = "../../analysis/tissueSpecific/SignificantPositionsOut_Liver_rc.fa";
	my $seqOut_Brain_rc = "../../analysis/tissueSpecific/SignificantPositionsOut_Brain_rc.fa";
	
	open IN_Liver, "<$seqOut_Liver" or die "Can not open file: $seqOut_Liver";
	open IN_Brain, "<$seqOut_Brain" or die "Can not open file: $seqOut_Brain";
	open OUT_Liver, ">$seqOut_Liver_rc" or die "Can not open file: $seqOut_Liver_rc";
	open OUT_Brain, ">$seqOut_Brain_rc" or die "Can not open file: $seqOut_Brain_rc";
	
	my $seq = "";
	my $L_label = "";
	my $B_label="";
	my $stop = "false";
	my $L_l;
	my $B_l;		
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
			$seq = "";
			
			#read the sequence from either brain or liver
			if($liver > $brain) #liver
			{
				if($L_label eq "")
				{
					$L_label  = <IN_Liver>; 
					chomp($L_label);
				}
				else
				{
					$L_label = $L_l;
				}
				
				LOOP:{
					while ($L_l = <IN_Liver>)
					{
						chomp($L_l);
						if($L_l =~ m/^\>/)
						{
							last LOOP;						
						}
						else
						{
						
							$seq = $seq . $L_l; 
						}
						
					}
				}
			}		
			elsif($brain > $liver) #brain
			{
				
				if($B_label eq "")
				{
					$B_label  = <IN_Brain>; 
					chomp($B_label);
				}
				else
				{
					$B_label = $B_l;
				}
				
				LOOP2:{
					while ($B_l = <IN_Brain>)
					{
						chomp($B_l);
						if($B_l =~ m/^\>/)
						{
							last LOOP2;
						}
						else
						{
							$seq = $seq . $B_l; 
						}
						
					}
				}
			}
		
			#check strand, reverse complement if -ve
			if($strand eq "1")
			{	
				if($liver > $brain)
				{
					print OUT_Liver ("$L_label\n$seq\n");
				}		
				elsif($brain > $liver)
				{
					print OUT_Brain ("$B_label\n$seq\n");
				}
			}
			else
			{	
				# reverse the DNA sequence
		        my $revcomp = reverse($seq);
		
				# complement the reversed DNA sequence
		        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
				
				
				if($liver > $brain)
				{
					print OUT_Liver ("$L_label\n$revcomp\n");
				}		
				elsif($brain > $liver)
				{
					print OUT_Brain ("$B_label\n$revcomp\n");
				}
			}
			
		}
	
	}
	
closeConnection ($dbh);          	
}

sub getSignificantSequencesForOriginalPAS_rc
{
	my ($infile, $faSeq, $seqOut, $strandFile) = @_;
	
	my %genes = ();
	my $dbh = DBconnect();
	
	open IN, "<$infile" or die "Can not open file :$infile";
	
	open OUT, ">$faSeq" or die "Can not open file: $faSeq";
	open STRAND_OUT, ">$strandFile" or die "Can not open file: $strandFile";
	
	my $line = <IN>;
	while($line = <IN>)
	{
		chomp($line);
		my ($gene,$pas_median,$pas) = split(/\,/, $line);
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
		elsif($chr eq "MT")
		{
			print ("chr MT $gene\n");
		}
		elsif($chr eq "M")
		{
			print ("chr M $gene\n");
		}
		else
		{
			if(exists $genes{$gene})
			{
				print("this gene already exist $gene\n");
			}else
			{
				$genes{$gene} = 1;
			}
			if($strand eq "1")
			{	
				$start = $pas  - 101;
				$end = $pas + 100;
				print STRAND_OUT ("1\n");
			}
			else
			{	
				$start = $pas  - 101;
				$end = $pas + 100;
				print STRAND_OUT ("-1\n");
			}	
		
			print OUT ("chr$chr:$start-$end\n");			
		}
	
	}
	
	close (OUT);
	print ("done mapping\n");
	my $cmd = "/nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/twoBitToFa /nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/hg19.2bit -seqList=$faSeq $seqOut";
	print ("$cmd\n");
	system ($cmd);
	
	close (STRAND_OUT);
	#close (OUT)
	
	
	closeConnection ($dbh);   

	#reverse complement negative genes
	my $seq_file = $seqOut;
	$seq_file =~ s/.fa/_rc.fa/g;
	
	open IN, "<$seqOut" or die "Can not open file :$seqOut";
	open STRAND, "<$strandFile" or die "Can not open file: $strandFile";
	open OUT, ">$seq_file" or die "Can not open file :$seq_file";
    
    my $label = <IN>;
    chomp($label);
    my $seq;
    my $strand;
    my $Line;
	while ($strand = <STRAND>)
	{
		LOOP:{
			$seq = "";
			while ($Line = <IN>)
			{
				chomp($Line);
				if($Line =~ m/^\>/)
				{
					last LOOP;						
				}
				else
				{
					$seq = $seq . $Line; 
				}
			}
		}
		chomp ($strand);
		if($strand eq "1")
		{
			print OUT ("$label\n$seq\n");
		}
		else
		{
			# reverse the DNA sequence
		    my $revcomp = reverse($seq);
		
			# complement the reversed DNA sequence
		    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
		    
		    print OUT ("$label\n$revcomp\n");
		}
		
		
		$label = $Line;
	}
	close IN;
	close STRAND;
	close OUT;
}

sub getSignificantSequencesForOriginalPAS
{
	my ($infile, $faSeq,$seqOut) = @_;
	
	my $dbh = DBconnect();
	
	open IN, "<$infile" or die "Can not open file :$infile";
	
	open OUT, ">$faSeq" or die "Can not open file: $faSeq";
	
	my $line = <IN>;
	while($line = <IN>)
	{
		chomp($line);
		my ($gene,$pas_median,$pas) = split(/\,/, $line);
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
		
			print OUT ("chr$chr:$start-$end\n");			
		}
	
	}
	
	my $cmd = "/nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/twoBitToFa /nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/hg19.2bit -seqList=$faSeq $seqOut";
	print ("$cmd\n");
	system ($cmd);
	
#	my $cmd = "/nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/twoBitToFa /nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/hg19.2bit -seqList=$faSeq_Brain $seqOut_Brain";
#	print ("$cmd\n");
#	system ($cmd);
	
closeConnection ($dbh);          	
}




sub writeFilesinExpectedFormat
{
	my ($cell, $totalNoOfExample, $positiveOrNegative, $inFile) = @_;
	
	my $out = "../../analysis/tissueSpecific/".$cell."_".$positiveOrNegative.".fa";
	open OUT, ">$out" or die "Can not open file: $out";
	
	open IN, "<$inFile" or die "Can not open file :$inFile";
	
	my $Model = "../../analysis/tissueSpecific/".$cell."_".$positiveOrNegative."_model.fa";
	open MODEL, ">$Model" or die "Can not open file: $Model";
	
	my $TrainTest = "../../analysis/tissueSpecific/".$cell."_".$positiveOrNegative."_train_test.fa";
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

sub writeFilesinExpectedFormat2
{
	my ($cell, $totalNoOfExample, $positiveOrNegative) = @_;
		
	my $in = "../../analysis/tissueSpecific/SignificantPositionsOut_".$cell."2.fa";
	open IN, "<$in" or die "Can not open file :$in";

	
	my $Model = "/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/SVM/".$cell."_model.fa";
	open MODEL, ">$Model" or die "Can not open file: $Model";
	
	my $TrainTest = "/nfs/central/home/dmh31/RIP/PA_all_2/analysis/tissueSpecific/SVM/".$cell."_train_test.fa";
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
			
			$i++;
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
	
	#open OUT_Liver, ">$faSeq_Liver" or die "Can not open file: $faSeq_Liver";
	#open OUT_Brain, ">$faSeq_Brain" or die "Can not open file: $faSeq_Brain";
	
	my $Brain_strand = "../../analysis/tissueSpecific/Brain_negative_strand.fa";
	my $Liver_strand = "../../analysis/tissueSpecific/Liver_negative_strand.fa";
	
	open STRAND_Liver, ">$Liver_strand" or die "Can not open file: $Liver_strand";
	open STRAND_Brain, ">$Brain_strand" or die "Can not open file: $Brain_strand";
	
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
					#print OUT_Liver ("chr$chr:$upstream-$downstream\n");
				}		
				elsif($brain > $liver)
				{
					#print OUT_Brain ("chr$chr:$upstream-$downstream\n");
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
				#	print OUT_Liver ("chr$chr:$upstream-$downstream\n");
				}		
				elsif($brain > $liver)
				{
				#	print OUT_Brain ("chr$chr:$upstream-$downstream\n");
				}
				
	 			
		 	}
		 }
		 if($liver > $brain)
		{
			print STRAND_Liver ("$strand\n");
		}		
		elsif($brain > $liver)
		{
			print STRAND_Brain ("$strand\n");
		}
	}
	closeConnection($dbh);
 	
 	my $cmd = "/nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/twoBitToFa /nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/hg19.2bit -seqList=$faSeq_Liver $seqOut_Liver";
	print ("$cmd\n");
	#system ($cmd);

 	$cmd = "/nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/twoBitToFa /nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/hg19.2bit -seqList=$faSeq_Brain $seqOut_Brain";
	print ("$cmd\n");
	#system ($cmd);
	
	#printLabel($seqOut_Liver,"Liver");
	#printLabel($seqOut_Brain,"Brain");
	
	#print("number of genes that have stop > polya $stop\n");
	#print ("number of gnees that do not have stop $rangeBad\n");
	
	
	#reverse complement negative genes
	my $seq_file_Liver = $seqOut_Liver;
	$seq_file_Liver =~ s/.fa/_rc.fa/g;
	
	my $seq_file_Brain = $seqOut_Brain;
	$seq_file_Brain =~ s/.fa/_rc.fa/g;
	
	open IN, "<$seqOut_Liver" or die "Can not open file :$seqOut_Liver";
	open STRAND, "<$Liver_strand" or die "Can not open file: $Liver_strand";
	open OUT, ">$seq_file_Liver" or die "Can not open file :$seq_file_Liver";
    
    my $label = <IN>;
    chomp($label);
    my $seq;
    my $strand;
    my $Line;
	while ($strand = <STRAND>)
	{
		for (my $i=0;$i<$samples;$i++)
		{
			LOOP:{
				$seq = "";
				while ($Line = <IN>)
				{
					chomp($Line);
					if($Line =~ m/^\>/)
					{
						last LOOP;						
					}
					else
					{
						$seq = $seq . $Line; 
					}
				}
			}
			chomp ($strand);
			if($strand eq "+")
			{
				print OUT ("$label\n$seq\n");
			}
			else
			{
				# reverse the DNA sequence
			    my $revcomp = reverse($seq);
			
				# complement the reversed DNA sequence
			    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
			    
			    print OUT ("$label\n$revcomp\n");
			}	
			$label = $Line;
		}
	}
	for (my $i=0;$i<$samples;$i++)
		{
			LOOP:{
				$seq = "";
				while ($Line = <IN>)
				{
					chomp($Line);
					if($Line =~ m/^\>/)
					{
						last LOOP;						
					}
					else
					{
						$seq = $seq . $Line; 
					}
				}
			}
			chomp ($strand);
			if($strand eq "+")
			{
				print OUT ("$label\n$seq\n");
			}
			else
			{
				# reverse the DNA sequence
			    my $revcomp = reverse($seq);
			
				# complement the reversed DNA sequence
			    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
			    
			    print OUT ("$label\n$revcomp\n");
			}	
			$label = $Line;
		}
	close IN;
	close STRAND;
	close OUT;
	
	#BRAIN
	open IN, "<$seqOut_Brain" or die "Can not open file :$seqOut_Brain";
	open STRAND, "<$Brain_strand" or die "Can not open file: $Brain_strand";
	open OUT, ">$seq_file_Brain" or die "Can not open file :$seq_file_Brain";
    
    my $label = <IN>;
    chomp($label);
    my $seq;
    my $strand;
    my $Line;
	while ($strand = <STRAND>)
	{
		print ($strand);
		for (my $i=0;$i<$samples;$i++)
		{
			LOOP:{
				$seq = "";
				while ($Line = <IN>)
				{
					chomp($Line);
					if($Line =~ m/^\>/)
					{
						last LOOP;						
					}
					else
					{
						$seq = $seq . $Line; 
					}
				}
			}
			chomp ($strand);
			if($strand eq "+")
			{
				print OUT ("$label\n$seq\n");
			}
			else
			{
				# reverse the DNA sequence
			    my $revcomp = reverse($seq);
			
				# complement the reversed DNA sequence
			    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
			    
			    print OUT ("$label\n$revcomp\n");
			}
			
			
			$label = $Line;
		}
	}
	
	close IN;
	close STRAND;
	close OUT;
 }
 
sub getNegativePosAfterStopCodon_common_cons2
 {
 	my ($inFile, $faSeq_in,$faSeq_out, $name) = @_;
 	my $dbh = DBconnect();
	my $counter = 0;
	
	srand(time() ^ $$ ^ unpack "%32L*", `ps axww | gzip`);#set the random number seed for the rand operator
	
	
	my $strand_file = $faSeq_in;
	$strand_file =~ s/.fa/_strand.fa/g;
	open IN, "<$inFile" or die "Can not open file :$inFile";
	open STRAND, ">$strand_file" or die "Can not open file: $strand_file";
	open OUT, ">$faSeq_in" or die "Can not open file: $faSeq_in";
	my $count = 0;
	my $stop = 0;
 	my $empty =0;
 	my $rangeBad = 0;
	my $line = <IN>;

	my $f  = "f";
	my $stopPolya = 0;

	while($line = <IN>)
	{
		my $flag = "t";
		chomp($line);
		
		my ($gene,$median,$polya,$diff) = split(/\,/, $line);
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
			
	 		#select stop Codon position
	 		my $correctStop = "true";
	 		my $range = 0;
	 		
	 		my @rows_selected = DBselectUpdatedGeneStopCodonFromGene_max($dbh,$sqlGene,$strand);
	 		my $stopCodon = $rows_selected[0];
	 		if($stopCodon)
	 		{
	 			$stopCodon = $rows_selected[0];
	 			
		 		if($strand eq "+")
				{
					if($stopCodon >=$polya)
		 			{
		 				$stopCodon = $polya - 500;
						$range =   $polya - 100 - $stopCodon +1;	
					}
					else
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
		 		}
				else  #negative strand
				{	
			 		if($stopCodon <=$polya)		
			 		{
						$stopCodon = $polya + 500;
			 			$range =   $stopCodon - ($polya + 100) +1;
					}
					else
					{
						$range =   $stopCodon - ($polya + 100) +1;
				 		if($range < 100)
				 		{
				 			$rangeBad++;
				 			$stopCodon = $polya + 500;
				 			$range =   $stopCodon - ($polya + 100) +1;
				 		}
					}
				}	
			}
			else
			{
				#no stop codon selected
				if($strand eq "+")
				{
					$stopCodon = $polya - 500;
					$range =   $polya - 100 - $stopCodon +1;	
				}
				else
				{
					$stopCodon = $polya + 500;
		 			$range =   $stopCodon - ($polya + 100) +1;
				}
			}
			
 			#randomly pick 10 points
 			
		 	my $random_pos =0;
		 	for (my $i=0;$i<$samples;$i++)
		 	{ 
		 		$random_pos = int(rand($range)) + $stopCodon;
		 		my $upstream = $random_pos-101;
	 			my $downstream = $random_pos + 100;
	 			#determine whether it is brain or liver specific
				print OUT ("chr$chr:$upstream-$downstream\n");	 			
		 	}
		 	
		 	print STRAND ("$strand\n");
		 	
		 	if($range == 0)
		 	{
		 		print ("somethign is missing\n");
		 	}
		} #end else
	}#end while 
	
	print ("count is $count\n");
	closeConnection($dbh);
 	
 	my $cmd = "/nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/twoBitToFa /nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/hg19.2bit -seqList=$faSeq_in $faSeq_out";
	print ("$cmd\n");
	system ($cmd);

	#printLabel($faSeq_out,$name);
	print("number of genes that have stop > polya $stop\n");
	print ("number of gnees that do not have stop $rangeBad\n");
	
	close STRAND;
	close IN;
	close OUT;
	
	#reverse complement negative genes
	
	my $faSeq_rc = $faSeq_out;
	$faSeq_rc =~ s/.fa/_rc.fa/g;
	
	open STRAND, "<$strand_file" or die "Can not open file: $strand_file";
	open IN, "<$faSeq_out" or die "Can not open file: $faSeq_out";
	open OUT, ">$faSeq_rc" or die "Can not open file: $faSeq_rc";

	my $label = <IN>;
    chomp($label);
    my $seq;
    my $strand;
    my $Line;
	while ($strand = <STRAND>)
	{
		for (my $i=0;$i<$samples;$i++)
		{ 
			LOOP:{
				$seq = "";
				while ($Line = <IN>)
				{
					chomp($Line);
					if($Line =~ m/^\>/)
					{
						last LOOP;						
					}
					else
					{
						$seq = $seq . $Line; 
					}
				}
			}
			chomp ($strand);
			if($strand eq "+")
			{
				print OUT ("$label\n$seq\n");
			}
			else
			{
				# reverse the DNA sequence
			    my $revcomp = reverse($seq);
			
				# complement the reversed DNA sequence
			    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
			    
			    print OUT ("$label\n$revcomp\n");
			}
			
			
			$label = $Line;
		}
	}
	
	close IN;
	close STRAND;
	close OUT;
	print ("number of genes that have no stop codon = $counter\n");
 }
sub getNegativePosAfterStopCodon_common_cons
 {
 	my ($inFile, $faSeq_in,$faSeq_out, $name) = @_;
 	my $dbh = DBconnect();
	
	my $strand_file = $faSeq_in;
	$strand_file =~ s/.fa/_strand.fa/g;
	open IN, "<$inFile" or die "Can not open file :$inFile";
	open STRAND, ">$strand_file" or die "Can not open file: $strand_file";
	open OUT, ">$faSeq_in" or die "Can not open file: $faSeq_in";
	my $count = 0;
	my $stop = 0;
 	my $empty =0;
 	my $rangeBad = 0;
	my $line = <IN>;

	my $f  = "f";
	#if($f eq "t"){	
	while($line = <IN>)
	{
		chomp($line);
		
		my ($gene,$median,$polya,$diff) = split(/\,/, $line);
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
		 $count++;
		 print STRAND ("$strand\n");
	#}
	}
	print ("count is $count\n");
	closeConnection($dbh);
 	
 	my $cmd = "/nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/twoBitToFa /nfs/labs/ohlerlab/sata/data/dina/PA_all_2/genome/hg19.2bit -seqList=$faSeq_in $faSeq_out";
	print ("$cmd\n");
	system ($cmd);

	#printLabel($faSeq_out,$name);
	print("number of genes that have stop > polya $stop\n");
	print ("number of gnees that do not have stop $rangeBad\n");
	
	close STRAND;
	close IN;
	close OUT;
	
	#reverse complement negative genes
	
	my $faSeq_rc = $faSeq_out;
	$faSeq_rc =~ s/.fa/_rc.fa/g;
	
	open STRAND, "<$strand_file" or die "Can not open file: $strand_file";
	open IN, "<$faSeq_out" or die "Can not open file: $faSeq_out";
	open OUT, ">$faSeq_rc" or die "Can not open file: $faSeq_rc";

	my $label = <IN>;
    chomp($label);
    my $seq;
    my $strand;
    my $Line;
	while ($strand = <STRAND>)
	{
		for (my $i=0;$i<$samples;$i++)
		{ 
			LOOP:{
				$seq = "";
				while ($Line = <IN>)
				{
					chomp($Line);
					if($Line =~ m/^\>/)
					{
						last LOOP;						
					}
					else
					{
						$seq = $seq . $Line; 
					}
				}
			}
			chomp ($strand);
			if($strand eq "+")
			{
				print OUT ("$label\n$seq\n");
			}
			else
			{
				# reverse the DNA sequence
			    my $revcomp = reverse($seq);
			
				# complement the reversed DNA sequence
			    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
			    
			    print OUT ("$label\n$revcomp\n");
			}
			
			
			$label = $Line;
		}
	}
	
	close IN;
	close STRAND;
	close OUT;
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
 
 
 ############

 
 sub compareTranscriptLength   #Neel
{
	my($liverF,$brainF,$kidneyF,$specificFile) = @_;
	my $dbh = DBconnect();
	my %PAS_Liver = ();
	my %PAS_Brain = ();
	my %PAS_Kidney = ();
	my $max_dif = 0;
	
	open IN, "<$liverF" or die "Can not open file: $liverF";
	my $line = <IN>;
	while ($line = <IN>)
	{
		chomp ($line);
		my ($gene,$mode,$PAS_count,$Total,$PAS_start,$PAS_end) = split(/\,/, $line);
		$PAS_Liver{$gene}{$mode}=$PAS_count;
	}
	close IN;
	
	open IN, "<$brainF" or die "Can not open file: $brainF";
	$line = <IN>;
	while ($line = <IN>)
	{
		chomp ($line);
		my ($gene,$mode,$PAS_count,$Total,$PAS_start,$PAS_end) = split(/\,/, $line);
		$PAS_Brain{$gene}{$mode}=$PAS_count;
	}
	close IN;
	
	open IN, "<$kidneyF" or die "Can not open file: $kidneyF";
	my $line = <IN>;
	while ($line = <IN>)
	{
		chomp ($line);
		my ($gene,$mode,$PAS_count,$Total,$PAS_start,$PAS_end) = split(/\,/, $line);
		$PAS_Kidney{$gene}{$mode}=$PAS_count;
	}
	close IN;
	
	open SPEC, "<$specificFile" or die "Can not open file: $specificFile";
	$line = <SPEC>;
	my $shortCount =0;
	my $longCount = 0;
	my $midCount = 0;
	
	my $total=0;
	while ($line = <SPEC>)
	{
		chomp ($line);
	
		my ($gene, $PAS_median,$PAS_original, $dif) = split(/\,/, $line);
		
		my $sqlGene = "'$gene'";
		my @result = DBselectChrStrand($dbh, $sqlGene);
		my $chr = $result[0];
		my $strand = $result[1];
		my $set = "0";
		
		if($PAS_original eq "196769432")
				{
					print ("here");
				}
			if((isLongest($strand,$PAS_original,$PAS_Brain{$gene}) eq "true")&&(isLongest($strand,$PAS_original,$PAS_Liver{$gene}) eq "true"))
			{
				$longCount++;
				$set = "1";
			}
			if((isShortest($strand,$PAS_original,$PAS_Brain{$gene}) eq "true")&&(isShortest($strand,$PAS_original,$PAS_Liver{$gene}) eq "true"))
			{
				$shortCount ++;
				
				
					print ("$gene,$PAS_original,$strand\t");
				
			}
			#else
			#{
			#	$midCount++;
			#}
			
			$total ++;
		
	}
	
	print ("Total number is $total\n");
	print ("number of shortest specific transcripts:$shortCount\n");
	print ("number of Longest specific transcripts:$longCount\n");
	print ("number of middle specific transcripts:$midCount\n");
	
	closeConnection($dbh);
	
}

sub isShortest
{
	my ($strand,$PAS,$PAS_other) = @_;

	my $short = "true";
		if($strand eq "1")
		{
			while (my ($other, $count) = each (%$PAS_other))
			{
				if($other < $PAS )
				{
					$short = "false";
					
				}
			}
		}
		else
		{
			while (my ($other, $count) = each (%$PAS_other))
			{
				
				if($other > $PAS )
				{
					$short = "false";
				
				}
			}
		}
	
	return $short;
	
}


sub isLongest
{
	my ($strand,$PAS,$PAS_other) = @_;
	my $long = "true";
	
	
	if($strand eq "1")
	{
		while (my ($other, $count) = each (%$PAS_other))
		{
			if(($other > $PAS ))
			{
				$long = "false";
				
			}
		}
	}
	else
	{
		while (my ($other, $count) = each (%$PAS_other))
		{
			if($PAS eq "196769432")
				{
					print ("Long $other\t");
				}
			if(($other < $PAS ))
			{
				$long = "false";
				
			}
		}
	}
	return $long;
	

}
 
