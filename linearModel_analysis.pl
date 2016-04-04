use strict;
use Switch;
use List::Util qw(sum);
#use Statistics::R ;
use List::Util qw(shuffle);
use List::Util qw(sum);
use Math::Random qw(random_multinomial);
use Math::Random qw(random_permutation);
use List::MoreUtils qw(firstidx);
 use List::MoreUtils qw(indexes);
 use Statistics::Zscore qw();
 #use DBHandler;
 use POSIX;
 use Storable qw(dclone);


	
#This script performs the downstream analysis of the linear model
#It first runs the old school analysis methodologies like ChiSquare and Shannon Entropy
#Then it generates the files required for running the mixed model.
#the linear model itself is written in a separate  R script
#After the output of linear model is generated, A followup analysis is performed on the output
#Including pvale calculation 

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



#Analysis of the generated PAS counts
#Signinificance test using ChiSquare and Shannon Entropy (old school)
ChiSquareTest();
calculateShanonEntropy_review();

#My linear model
MixedModel();
#MixedModel_allCellTypes ();


#Analysing the output
RemoveSeq2();
compareWithMac();
prepareBedFileIGV_3CellTypes();
countA();
getSigPvalues();
calculateShanonEntropyForSignificantLME();
pvalue_bed();
prepareFisherExact_Daniel();
parseFisherFile();

computeDifferences_original_3cellTypes();
computeDifferences_allCellTypes();
getSignificantDifferences_allCellTypes();

computeDifferences_permuted();
computeDifferences_original();
retrieveOriginalPASnotMedian ();
retrieveOriginalPASnotMedianFromConstitutive2();
retrieveOriginalPASnotMedianFromCommon();
subStringSequences("../../analysis/tissueSpecific/SVM_rc/Liver_pos.fa", "../../analysis/tissueSpecific/SVM_rc/Liver_pos_rc_upstream.fa", "../../analysis/tissueSpecific/SVM_rc/Liver_pos_rc_downstream.fa");
subStringSequences("../../analysis/tissueSpecific/SVM_rc/Brain_Liver_constitutive_mapped_to_Brain_out_rc.fa", "../../analysis/tissueSpecific/SVM_rc/Brain_Liver_constitutive_mapped_to_Brain_rc_upstream.fa", "../../analysis/tissueSpecific/SVM_rc/Brain_Liver_constitutive_mapped_to_Brain_rc_downstream.fa");
subStringSequences2("../../analysis/tissueSpecific/SVM_3CellTypes/Constitutive/constitutive_Brain_neg.fa", "../../analysis/tissueSpecific/SeqLogo/constitutive_Brain_40.fa");

#Downstream analysis for figures in the paper
prepareBedFileIGV();
 analyze_residual_all_cells_types();

calculateMedianUsage();
calculateDistanceToMedianUsage();

computeZscore_original_3cellTypes();
computeDifferences_permuted_3cellTypes();

#Higly expressed in one cell type -- significance test
computeMaxDifFromMedian_original_3cellTypes();
computeMaxDifFromMedian_permuted_3cellTypes();
prepareBedFileIGV_3CellTypes();


#Highly expressed in two cell types
computeHighlyExpressedinTwoCells_original();
computeHighlyExpressedinTwoCells_permuted();

getSpecificExclusive2();

retrieveOriginalPASnotMedian_3cellTypes_highestInOneCell();
retrieveOriginalPASnotMedian_3cellTypes_highestInTwoCells();

compareTranscriptLength("../../analysis/tissueSpecific/PAS_median_3CellTypes", "../../analysis/tissueSpecific/HighlyExpressedInOneCell/Brain_Specific_3CellTypes_PAS_Original");



sub ChiSquareTest
{
	#read medians
		my $PAS_median_file ="../../analysis/tissueSpecific/PAS_median";
		open IN, "<$PAS_median_file" or die "Can not open file: $PAS_median_file";
        my %PAS_median;
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
		while (my ($lane, $cellHash) = each %dataset)
        {
             while (my ($cell, $cell_file) = each %$cellHash)
           	{
           			
      				my $PAS_Genes_file=  $cell_file;
      				my $tissuePath = "analysis/tissueSpecific";
      				$PAS_Genes_file =~ s/analysis/$tissuePath/g;
                    $PAS_Genes_file =~ s/uniq/PAS.Genes.Adjusted/g;
                    open IN, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                  	
				 	while (my $line = <IN>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$PAS_count_normalized,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 		#if in median +-10, refer to it as median, else, put it as it is
				 			
					 	Exist:
					 	{
					 		for(my $k=0;$k<=10;$k++)
					 		{
						 		if(exists $PAS_median{$gene}{$mode+$k})
						 		{
						 			$PASs{$gene}{$cell}{$mode+$k} =$PAS_count_normalized;
						 			
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
						 			$PASs{$gene}{$cell}{$mode-$k} =$PAS_count_normalized;
						 			if ( grep { $_ eq $mode-$k } @{$PAS_modes{$gene}} ) {
						 			}
						 			else
									{	push(@{$PAS_modes{$gene}},$mode-$k);}
						 			
									last Exist;
						 		}
					 		}
					 	
				 			$PASs{$gene}{$cell}{$mode} =$PAS_count_normalized;
							if ( grep { $_ eq $mode } @{$PAS_modes{$gene}} ) {
						 	}
							else
				 			{push(@{$PAS_modes{$gene}},$mode);}
							
			
					 	}
				 	}
				 	close IN;
           	}
        }
        
        #output Chi test files
        my $totalGenes =0;
        my $genes_alternative =0;
        my $PAS_Genes_Chi = "../../analysis/tissueSpecific/PAS_Chi";
       	open OUT_F, ">$PAS_Genes_Chi" or die "Can not open file: $PAS_Genes_Chi";
       	
       	my $R = Statistics::R->new();
		$R->start_sharedR ;
        
        my $count_sig = 0;
        my $count_total = 0;
        
        my $All_Genes_PAS = "../../analysis/tissueSpecific/ALL_Genes_PAS";
       	open OUT_PAS, ">$All_Genes_PAS" or die "Can not open file: $All_Genes_PAS";
       	
        while (my ($gene, $modeHash) = each %PAS_modes)
        {
        	my $chi = "chi.csv";
       		open OUT, ">$chi" or die "Can not open file: $chi";
       	
        	$totalGenes ++;
        	my @modes = sort(@{$PAS_modes{$gene}});
        	if(@modes >1) #only print out genes that have more than 1 polyA site
        	{
        		print OUT_PAS ("$gene\n");
        		$genes_alternative ++;
	            for(my $i=0; $i< @modes; $i++)
	            {
	            	if($i==@modes-1)
	            	{
	            		print OUT ("$modes[$i]");
	            		print OUT_PAS ("$modes[$i]");
	            	}
	            	else
	            	{
	            		print OUT ("$modes[$i],");
	            		print OUT_PAS ("$modes[$i],");
	            	}
	            }
	           	print OUT ("\n");
	           	print OUT_PAS ("\n");
	           	
	           	
	           	while (my ($cell, $modeHash) = each %{$PASs{$gene}})
	           	{
	           		print OUT ("$cell,");
	           		print OUT_PAS ("$cell,");
	           		for(my $i=0; $i< @modes; $i++)
	            	{
	            		if($i==@modes-1)
	            		{
		            		if(exists $PASs{$gene}{$cell}{$modes[$i]})
							{
								print OUT ("$PASs{$gene}{$cell}{$modes[$i]}");
								print OUT_PAS ("$PASs{$gene}{$cell}{$modes[$i]}");
							}            		
							else
							{
								print OUT ("0");
								print OUT_PAS ("0");
							}
	            		}
	            		else
	            		{
		            		if(exists $PASs{$gene}{$cell}{$modes[$i]})
							{
								print OUT ("$PASs{$gene}{$cell}{$modes[$i]},");
								print OUT_PAS ("$PASs{$gene}{$cell}{$modes[$i]},");
							}            		
							else
							{
								print OUT ("0,");
								print OUT_PAS ("0,");
					
							}
	            		}
	            		
	           		}
	           		print OUT ("\n");
	           		print OUT_PAS ("\n");
	           	}
        	
        		close OUT;
        		
        		
        		#perform Chi square test
        		$R -> run('d <- read.csv ("chi.csv")');
        		$R -> run('options(warn=-1)');
				$R->run('chisq.test(d)');
				my $result = $R-> read();
				#$result2 =~ s/\$p.value//g;
				#$result2 =~ s/\n//g; 
				#$result2 =~ s/\[1]//g;
				#$result2 =~ s/Warning message:In chisq.test(d) : Chi-squared approximation may be incorrect//g;
				#print OUT_F ("$gene\t$result2\n");
				chomp ($result);
				my @values = split(/\s+/, $result);
				#print  ("$values[@values-1]\n");
				my $pvalue = $values[@values-1];
				print OUT_F ("$gene\t$pvalue\n");
				if ($pvalue <= 0.01)
				{
					$count_sig++;
				}
				$count_total++;
        	}
        }
	close OUT_PAS;
	print ("number of significant $count_sig / $count_total\n");
	
	my $perc =  $genes_alternative/ $totalGenes;
	print ("total number of genes alternative  $genes_alternative /$totalGenes = $perc\n");
	exit;
}

sub calculatePvalues
{
	#read original
		my $PAS_original ="Original1_header.csv";
		open NAME, "<$PAS_original" or die "Can not open file: $PAS_original";
		
		my $PAS_original_beta ="Original1_beta.csv";
		open BETA, "<$PAS_original_beta" or die "Can not open file: $PAS_original_beta";
        
        my %PAS_betas;
        
        my $name_line = <NAME>;
        my $beta_line = <BETA>;
        
        my @name_values = ();
		my @beta_values = ();
		
        my %PAS_betas;
		while ($name_line = <NAME>)
		{
			$beta_line = <BETA>;
			
			chomp ($name_line);
			chomp($beta_line);
			
			$name_line =~ s/"//g;
			$beta_line =~ s/"//g;
			
			@name_values = split(/\,/, $name_line);
			@beta_values = split(/\,/, $beta_line);
			
			for(my $i=0;$i<@name_values;$i++)
			{
				$PAS_betas{$name_values[$i]}{'beta'} = $beta_values[$i];
			}
			
		}
		close NAME;
		close BETA;
		
		my $start = 1;
		my $end = 10;
		my $flag = "true";
		while ($flag eq "true")
		{
			
		for(my $i=$start; $i<=$end; $i++)
		{
			@name_values = ();
			@beta_values = ();
			
			my $PAS_permute ="beta_header$i.csv";
			open NAME, "<$PAS_permute" or die  "Can not open file: $PAS_permute";;
			
			my $PAS_permute_beta ="beta$i.csv";
			open BETA, "<$PAS_permute_beta" or die "Can not open file: $PAS_permute_beta";;
	        
	        my $name_line = <NAME>;
	        my $beta_line = <BETA>;
	        
	       while ($name_line = <NAME>)
			{
				$beta_line = <BETA>;
				
				chomp ($name_line);
				chomp($beta_line);
				
				$name_line =~ s/"//g;
				$beta_line =~ s/"//g;
				
				@name_values = split(/\,/, $name_line);
				@beta_values = split(/\,/, $beta_line);
				
				for(my $i=0;$i<@name_values;$i++)
				{
					if(abs($PAS_betas{$name_values[$i]}{'beta'}) <= abs($beta_values[$i]))
					{
						$PAS_betas{$name_values[$i]}{'pvalue_count'}++;
					}
					
				}
				
			}
			close NAME;
			close BETA;	
		}
		if($start == 1)
		{
			$start = 51;
			$end = 60;
		}
		elsif($start == 51)
		{
			$start = 601;
			$end = 612;
			
		}
		elsif($start == 601)
		{
			$start = 671;
			$end = 677;
			
		}
		elsif($start == 671)
		{
			$start = 741;
			$end = 752;
			
		}
		elsif($start == 741)
		{
			$start = 801;
			$end = 811;
			
		}
		elsif($start == 801)
		{
			$start = 871;
			$end = 881;
			
		}
		elsif($start == 871)
		{
			$start = 126;
			$end = 600;
			
		}
		else
		{
			$flag = "false";
		}
			print("start is $start\n");
		}

		my $outfile_sig = "pvalues_significant";
		open OUTSIG, ">$outfile_sig" or die "Can not open file: $outfile_sig";
		
		my $outfile_nonsig = "pvalues_not_significant";
		open OUTNONSIG, ">$outfile_nonsig" or die "Can not open file: $outfile_nonsig";
		
		my $noOfPermutation = 10-1+1 + 60-51+1+ 612-601+1+ 677-671+1+ 752-741+1 + 811-801+1 + 881-871+1 + 600 - 126+1; #$end - $start +1;
		print ("Number of permutation is $noOfPermutation\n");
		
		my $specific05=0;
		my $specific01 = 0;
		my $total = 0;
		while (my ($pas, $value_hash) = each %PAS_betas)
		{
		   		my $pvalue = $PAS_betas{$pas}{'pvalue_count'} / $noOfPermutation;
		  		
		  		if(($pvalue <0.05)|| ($pvalue ==0))
		  		{
		  			
		  			$specific05++;
		  		}
		  		else
		  		{
		  			print OUTNONSIG ("$pas,$pvalue\n");
		  		}
		  		if(($pvalue <0.01)|| ($pvalue ==0))
		  		{print OUTSIG ("$pas,$pvalue\n");
		  			$specific01++;
		  		}
		  		$total++;
		}
		close OUT;
		
		 my $perc05 = $specific05 / $total;
		 my $perc01 = $specific01 / $total;
		
		print ("pvalue < 0.05 $specific05 / $total = $perc05 of number of permutation\n");
		print ("pvalue < 0.01 $specific01 / $total = $perc01 of number of permutation\n");
		
}
sub getSigPvalues
{
	#read original
		my $PAS_original ="../../analysis/tissueSpecific/pvalue_corrected.csv";
		open NAME, "<$PAS_original" or die "Can not open file: $PAS_original";
		
		
		my $outfile_sig = "../../analysis/tissueSpecific/pvalues_significant";
		open OUTSIG, ">$outfile_sig" or die "Can not open file: $outfile_sig";
		
		my $outfile_nonsig = "../../analysis/tissueSpecific/pvalues_not_significant";
		open OUTNONSIG, ">$outfile_nonsig" or die "Can not open file: $outfile_nonsig";
		
	  my %PAS_betas;
        
        my $name_line = <NAME>;
       
        my @name_values = ();
		my @beta_values = ();
		my $fdr_spec =0;
		my $specific05 = 0;
		my $total = 0;
		
        my %PAS_betas;
		while ($name_line = <NAME>)
		{
			chomp ($name_line);
			
			$name_line =~ s/"//g;
			
			my ($gene, $pas,$pvalue,$fdr) =  split(/\,/, $name_line);
			
			if($pvalue <=0.01)
		  	{
		  		print OUTSIG ("$gene,$pas,$pvalue\n");	
		  		$specific05++;
		  	}
		  	else
	  		{
	  			print OUTNONSIG ("$gene,$pas,$pvalue\n");
	  		}
		  	if($fdr <= 0.02)
		  	{
		  		$fdr_spec++;
		  		
		  	}
	  		
	  		$total++;
		}
		close NAME;
		close BETA;
		
		
		
		
		 my $perc05 = $specific05 / $total *100;
		 my $percFDR = $fdr_spec / $total *100;
		
		print ("pvalue < 0.05 $specific05 / $total = $perc05 \n");
		print ("pvalue < 0.01 $fdr_spec / $total = $fdr_spec \n");
		
}

sub getSigPvalues_permute
{
	#read original
		my $PAS_original ="../../analysis/tissueSpecific/PAS_model.csv";
		open NAME, "<$PAS_original" or die "Can not open file: $PAS_original";
		
		my $PAS_original_beta ="../../analysis/tissueSpecific/PAS_model_palues.csv";
		open BETA, "<$PAS_original_beta" or die "Can not open file: $PAS_original_beta";
        
        my %PAS_betas;
        
        my $name_line = <NAME>;
        my $beta_line = <BETA>;
        
        my @name_values = ();
		my @beta_values = ();
		
        my %PAS_betas;
		while ($name_line = <NAME>)
		{
			$beta_line = <BETA>;
			
			chomp ($name_line);
			chomp($beta_line);
			
			$name_line =~ s/"//g;
			$beta_line =~ s/"//g;
			
			@name_values = split(/\,/, $name_line);
			@beta_values = split(/\,/, $beta_line);
			
			for(my $i=0;$i<@name_values;$i++)
			{
				$PAS_betas{$name_values[$i]}{'beta'} = $beta_values[$i];
			}
			
		}
		close NAME;
		close BETA;
		
		
		my $outfile_sig = "../../analysis/tissueSpecific/pvalues_significant";
		open OUTSIG, ">$outfile_sig" or die "Can not open file: $outfile_sig";
		
		my $outfile_nonsig = "../../analysis/tissueSpecific/pvalues_not_significant";
		open OUTNONSIG, ">$outfile_nonsig" or die "Can not open file: $outfile_nonsig";
		
		my $specific05=0;
		my $specific01 = 0;
		my $total = 0;
		while (my ($pas, $value_hash) = each %PAS_betas)
		{
		   		my $pvalue = $PAS_betas{$pas}{'beta'} ;
		   		
		  		if((abs($pvalue) <0.05)|| ($pvalue ==0))
		  		{
		  			
		  			$specific05++;
		  		}
		  		else
		  		{
		  			print OUTNONSIG ("$pas,$pvalue\n");
		  		}
		  		if((abs($pvalue) < 0.01)|| ($pvalue ==0))
		  		{
		  			print OUTSIG ("$pas,$pvalue\n");
		  			$specific01++;
		  		}
		  		$total++;
		}
		close OUT;
		
		 my $perc05 = $specific05 / $total *100;
		 my $perc01 = $specific01 / $total *100;
		
		print ("pvalue < 0.05 $specific05 / $total = $perc05 \n");
		print ("pvalue < 0.01 $specific01 / $total = $perc01 \n");
		
}


sub calculateShanonEntropy_review
{
	#read medians
		my $PAS_median_file ="../../analysis/tissueSpecific/PAS_median_3CellTypes.ExpressedInThree";
		open IN, "<$PAS_median_file" or die "Can not open file: $PAS_median_file";
        my %PAS_median;
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
		while (my ($lane, $cellHash) = each %dataset)
        {
             while (my ($cell, $cell_file) = each %$cellHash)
           	{
           			
      				my $PAS_Genes_file=  $cell_file;
      				my $tissuePath = "analysis/tissueSpecific";
      				$PAS_Genes_file =~ s/analysis/$tissuePath/g;
                    $PAS_Genes_file =~ s/uniq/PAS.Genes.ExpressedInThree.Adjusted_Final/g;
                    
                    open IN, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                  	my $line = <IN>;
				 	while ( $line = <IN>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$PAS_count_adjusted_genelevel,$PAS_count_normalized,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 		#if in median +-10, refer to it as median, else, put it as it is
				 		
					 	Exist:
					 	{
					 		for(my $k=0;$k<=10;$k++)
					 		{
						 		if(exists $PAS_median{$gene}{$mode+$k})
						 		{
						 			$PASs{$cell}{$gene}{$mode+$k} =$PAS_count_normalized;
						 			$PASs_allCells{$gene}{$mode+$k}+=$PAS_count_normalized;
						 			
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
						 			$PASs{$cell}{$gene}{$mode-$k} =$PAS_count_normalized;
						 			$PASs_allCells{$gene}{$mode-$k}+=$PAS_count_normalized;
						 			if ( grep { $_ eq $mode-$k } @{$PAS_modes{$gene}} ) {
						 			}
						 			else
									{	push(@{$PAS_modes{$gene}},$mode-$k);}
						 			
									last Exist;
						 		}
					 		}
					 	
				 			$PASs{$cell}{$gene}{$mode} =$PAS_count_normalized;
				 			$PASs_allCells{$gene}{$mode}+=$PAS_count_normalized;
							if ( grep { $_ eq $mode } @{$PAS_modes{$gene}} ) {
						 	}
							else
				 			{push(@{$PAS_modes{$gene}},$mode);}
							
			
					 	}
				 	}
				 	close IN;
				 	
           	}
        }
        
        #calcluate shanon entropy
        #calculate relative expression
        while (my ($cell, $geneHash) = each %PASs)
        {
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{
        		my @modes = @{$PAS_modes{$gene}};
        	
	        	#if(@modes >1) #only print out genes that have more than 1 polyA site
	        	{
	        	#	if((exists $PASs{'Liver'}{$gene} )&&(exists $PASs{'Brain'}{$gene})&&(exists $PASs{'Kidney'}{$gene} )) #expressed in both tissues
        		{
	        		while (my ($PAS,$normalized_count) = each %$PAShash )
	        		{
	        		#	print("$gene,$PAS\n");
	        			$relative_exp{$cell}{$gene}{$PAS} = $normalized_count / $PASs_allCells{$gene}{$PAS};
	        		}
	        	}
	        	}
        	
        	}
        }
        
        #calculate entropy
        my %entropy_per_cell = ();
        while (my ($cell, $geneHash) = each %relative_exp)
        {
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{
        		while (my ($PAS,$rel_exp) = each %$PAShash )
        		{
        			$entropy_per_cell{$cell}{$gene}{$PAS} =$rel_exp * (log($rel_exp)/log(2));
        		}
        	}
        }
        
        my %entropy = ();
        while (my ($gene, $array) = each %PAS_modes)
        {
        	my @modes = @{$PAS_modes{$gene}};
        	
	        
        	#if(@modes >1) #only print out genes that have more than 1 polyA site
        	{
        	
	        	for(my $i=0; $i < @{$array}; $i++)
	        	{
	        		my @PASarray = @{$array};
	        		my $PAS = $PASarray[$i];
	        		#for all cells
	        		while (my ($lane, $cellHash) = each %dataset)
			        {
			             while (my ($cell, $cell_file) = each %$cellHash)
			           	{
			           		
			           		if (exists $entropy_per_cell{$cell}{$gene}{$PAS})
		        			{
	        					$entropy{$gene}{$PAS} += $entropy_per_cell{$cell}{$gene}{$PAS};
		        			}
		        			else
			        		{
			        			print ("$cell,$gene,$PAS\n");
			        			$entropy{$gene}{$PAS} += 0.001 * (log(0.001)/log(2));
			        		}
		        			
			        	}
			        }
        		}
        
        	}
        }
            
        #calculate categorical tissue specificity
        my %tissue_specificity = ();
         while (my ($cell, $geneHash) = each %PASs)
        {
        	
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{
        	#	if((exists $PASs{'Liver'}{$gene} )&&(exists $PASs{'Brain'}{$gene})&&(exists $PASs{'Kidney'}{$gene} )) #expressed in both tissues
        		{
        		my @modes = @{$PAS_modes{$gene}};
        		
	        
	       # 	if(@modes >1) #only print out genes that have more than 1 polyA site
	        	{
	        	
	        		while (my ($PAS,$normalized_count) = each %$PAShash )
	        		{
	        			$tissue_specificity{$cell}{$gene}{$PAS} = (-1 * $entropy{$gene}{$PAS}) - (log($relative_exp{$cell}{$gene}{$PAS})/log(2));
	        		}
	        	}
        		}
        	
        	}
        }
        
        #print to files
        
        my $entropy_file = "../../analysis/tissueSpecific/Entropy_Hvalue";
		open OUT, ">$entropy_file" or die "Can not open file: $entropy_file";	
        print OUT ("PAS,entropy\n");
        
        while (my ($gene, $PAShash) = each %entropy)
        {
        		while (my ($PAS,$spec) = each %$PAShash )
        		{
        			my $ent = -1 * $entropy{$gene}{$PAS} ;
        			print OUT ("$PAS,$ent\n");
        		
        	}
        }
        	close OUT;

        
        while (my ($cell, $geneHash) = each %tissue_specificity)
        {
        	my $relative_exp_file = "../../analysis/tissueSpecific/$cell.Entropy_Final2";
		 	open OUT, ">$relative_exp_file" or die "Can not open file: $relative_exp_file";	
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{
        		while (my ($PAS,$spec) = each %$PAShash )
        		{#added here trial
        			#if((exists $tissue_specificity{'Liver'}{$gene} )&&(exists $tissue_specificity{'Brain'}{$gene})&&(exists $tissue_specificity{'Kidney'}{$gene} )) #expressed in both tissues
        			{
        			print OUT ("$gene,$PAS,$tissue_specificity{$cell}{$gene}{$PAS}\n");
        			}
        		}
        	}
        	close OUT;
        }
        
        my $relative_exp_file = "../../analysis/tissueSpecific/Entropy_3cells_Final";
		open OUT, ">$relative_exp_file" or die "Can not open file: $relative_exp_file";	
        print OUT ("PAS,entropy\n");
       
       
        my $geneFile = "../../analysis/tissueSpecific/GeneChrMap.txt";
		open IN, "<$geneFile" or die "Can not open file: $geneFile";
	    
	    my %geneChrMap=();
	    my $line = <IN>;   	
		while ( $line = <IN>)
		{
			chomp ($line);
			my ($geneID, $transcriptID, $chr) = split(/\s+/, $line);
			$geneChrMap{$geneID} = $chr;
			
		}
		close IN;		 
	
				
		my $EntropyBed = "../../analysis/tissueSpecific/Entropy_3Cells_Final.bed";
		open OUT_2, ">$EntropyBed" or die "Can not open file: $EntropyBed";	
        
        while (my ($cell, $geneHash) = each %tissue_specificity)
        {
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{
        		while (my ($PAS,$spec) = each %$PAShash )
        		{
        			my @modes = @{$PAS_modes{$gene}};
        		
	        
	        	if(@modes >1) #only print out genes that have more than 1 polyA site
	      { 
	       
        		#	if((exists $tissue_specificity{'Liver'}{$gene} )&&(exists $tissue_specificity{'Brain'}{$gene})&&(exists $tissue_specificity{'Kidney'}{$gene} )) #expressed in both tissues
        		{
        			print OUT ("$PAS,$tissue_specificity{$cell}{$gene}{$PAS}\n");
        		 
        			#if($tissue_specificity{$cell}{$gene}{$PAS} <=0.5)
        			{
	        			if(exists $geneChrMap{$gene})
						{
							my $start = $PAS -1;
							my $num = eval sprintf('%.2f', $tissue_specificity{$cell}{$gene}{$PAS});
							print OUT_2 ("chr$geneChrMap{$gene}\t$start\t$PAS\t$num\n");
							
						}
						else {print ("$gene\t");}
        			}
        		}
        	
        			
        		}}
        	}
        	
        }
        
       
        	 	
}


sub prepareFisherExact_Daniel
{
	#read medians
		my $PAS_median_file ="../../analysis/tissueSpecific/PAS_median_3CellTypes.ExpressedInThree";
		open IN, "<$PAS_median_file" or die "Can not open file: $PAS_median_file";
        my %PAS_median;
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
		
		#get totals
		my $kidney_total=0;
		my $liver_total=0;
		my $brain_total=0;
		
		while (my ($lane, $cellHash) = each %dataset)
        {
             while (my ($cell, $cell_file) = each %$cellHash)
           	{
           			
      				my $PAS_Genes_file=  $cell_file;
      				my $tissuePath = "analysis/tissueSpecific";
      				$PAS_Genes_file =~ s/analysis/$tissuePath/g;
                    $PAS_Genes_file =~ s/uniq/PAS.Genes.ExpressedInThree.Adjusted_Final/g;
                    open IN, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                  	
				 	while (my $line = <IN>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$other ) = split(/\,/, $line);
				 		if($cell eq "Kidney")
				 		{
				 			$kidney_total+=$PAS_count;
				 		}
				 		elsif($cell eq "Liver")
				 		{
				 			$liver_total+=$PAS_count;
				 		}
				 		elsif($cell eq "Brain")
				 		{
				 			$brain_total += $PAS_count;
				 		}
				 		else
				 		{
				 			print ("error\n");
				 		}
				 	
				 	}
           	}
        }
        
        
        
		while (my ($lane, $cellHash) = each %dataset)
        {
             while (my ($cell, $cell_file) = each %$cellHash)
           	{
           			
      				my $PAS_Genes_file=  $cell_file;
      				my $tissuePath = "analysis/tissueSpecific";
      				$PAS_Genes_file =~ s/analysis/$tissuePath/g;
                    $PAS_Genes_file =~ s/uniq/PAS.Genes.ExpressedInThree.Adjusted_Final/g;
                    open IN, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                  	
				 	while (my $line = <IN>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$adjustedGene, $adjusted_lib,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 		#if in median +-10, refer to it as median, else, put it as it is
				 		
					 	Exist:
					 	{
					 		for(my $k=0;$k<=10;$k++)
					 		{
						 		if(exists $PAS_median{$gene}{$mode+$k})
						 		{
						 			$PASs{$cell}{$gene}{$mode+$k} =$PAS_count;
						 			$PASs_allCells{$gene}{$mode+$k}+=$PAS_count;
						 			
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
						 			if ( grep { $_ eq $mode-$k } @{$PAS_modes{$gene}} ) {
						 			}
						 			else
									{	push(@{$PAS_modes{$gene}},$mode-$k);}
						 			
									last Exist;
						 		}
					 		}
					 	
				 			$PASs{$cell}{$gene}{$mode} =$PAS_count;
				 			$PASs_allCells{$gene}{$mode}+=$PAS_count;
							if ( grep { $_ eq $mode } @{$PAS_modes{$gene}} ) {
						 	}
							else
				 			{push(@{$PAS_modes{$gene}},$mode);}
							
			
					 	}
				 	}
				 	close IN;
				 	
           	}
        }
        
        #calcluate Fisher
       	
        
         my $fisher_in = "../../analysis/tissueSpecific/Fisher_in.csv";
		open OUT, ">$fisher_in" or die "Can not open file: $fisher_in";	
        
        print OUT ("cell,gene,pas1,pas2,count1,count2,count1_normalized,count2_normalized,count1other,count2other\n");
        while (my ($cell, $geneHash) = each %PASs)
        {
        	my $other_to_normalize ;
	       	if($cell eq "Kidney")
			{
				$other_to_normalize = int(($liver_total+$brain_total)/$kidney_total);
			}
			 elsif($cell eq "Liver")
			{
				$other_to_normalize = int(($kidney_total+$brain_total)/$liver_total);
			}
			 elsif($cell eq "Brain")
			 {
				$other_to_normalize = int(($liver_total+$kidney_total)/$brain_total);
			}
			
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{	
				my $pas1="0";
				my $pas2="0";
				my $count1="0";
				my $count2="0";
				my $count1other="0";
				my $count2other="0";
				my $count1_normalized = 0;
				my $count2_normalized = 0;
				
				#for  my $key (sort { $a <=> $b } keys(%$PAShash) )
				foreach (sort { ($$PAShash{$b} <=> $$PAShash{$a}) || ($a <=> $b) } keys %$PAShash) 
				{
					
					if($pas1 eq "0")
					{
						$pas1 =$_;
						$count1 = $$PAShash{$_}; #*$other_to_normalize
						$count1_normalized = $count1 * $other_to_normalize;
						$count1other = $PASs_allCells{$gene}{$pas1} - $count1;
					}
					elsif($pas2 eq "0")
					{
						$pas2 = $_;
						$count2 = $$PAShash{$_};
						$count2other = $PASs_allCells{$gene}{$pas2} - $count2;
						$count2_normalized = $count2 * $other_to_normalize;
						
						print OUT ("$cell,$gene,$pas1,$pas2,$count1,$count2,$count1_normalized,$count2_normalized,$count1other,$count2other\n");
						$pas1=$pas2;
						$pas2="0";
						$count1=$count2;
						$count2="0";
						$count1other=$count2other;
						$count2other="0";
						$count1_normalized = $count2_normalized;
					}
				}
			
        	}
        }
        close OUT;
        
        
}
sub parseFisherFile
{
	my $fisher_file ="../../analysis/tissueSpecific/file.txt";
	
     my $geneFile = "../../analysis/tissueSpecific/GeneChrMap.txt";
		open IN, "<$geneFile" or die "Can not open file: $geneFile";
	    
	    my %geneChrMap=();
	    my $line = <IN>;   	
		while ( $line = <IN>)
		{
			chomp ($line);
			my ($geneID, $transcriptID, $chr) = split(/\s+/, $line);
			$geneChrMap{$geneID} = $chr;
			
		}
		close IN;	
		my $FisherBed = "../../analysis/tissueSpecific/Fisher.bed";
		open OUT, ">$FisherBed" or die "Can not open file: $FisherBed";	
		
		my $FisherBedsig = "../../analysis/tissueSpecific/Fisher_sig.bed";
		open OUT2, ">$FisherBedsig" or die "Can not open file: $FisherBedsig";	
		
        open IN, "<$fisher_file" or die "Can not open file: $fisher_file";
    
		
		while (my $line = <IN>)
		{
			chomp ($line);
			my ($num,$cell,$gene,$PAS1,$PAS2,$count1,$count2,$normalized_count1,$normalized_count2,$count1other,$count2other,$p)  = split(/\s+/, $line);
			
		
				$gene =~s/"//g;
				#print $gene;
				if(exists $geneChrMap{$gene})
				{
					my $start;
					my $end;
					if($normalized_count1 >= $normalized_count2)
					{
						$start = $PAS1 -1;
						$end = $PAS1;
						
					}
					else
					{
						$start = $PAS2-1;
						$end = $PAS2;
						print ("should not be\n");
					}
					
					
					print OUT ("chr$geneChrMap{$gene}\t$start\t$end\t$p\n");
					
				
				if($p < 0.05)
				{
					print OUT2 ("chr$geneChrMap{$gene}\t$start\t$end\t$p\n");
				}
			
        			
				}
			
		}
		close OUT;
		close OUT2;
} 


sub calculateShanonEntropy
{
	#read medians
		my $PAS_median_file ="../../analysis/tissueSpecific/PAS_median_3CellTypes";
		open IN, "<$PAS_median_file" or die "Can not open file: $PAS_median_file";
        my %PAS_median;
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
		while (my ($lane, $cellHash) = each %dataset)
        {
             while (my ($cell, $cell_file) = each %$cellHash)
           	{
           			
      				my $PAS_Genes_file=  $cell_file;
      				my $tissuePath = "analysis/tissueSpecific";
      				$PAS_Genes_file =~ s/analysis/$tissuePath/g;
                    $PAS_Genes_file =~ s/uniq/PAS.Genes.Adjusted_2/g;
                    open IN, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                  	
				 	while (my $line = <IN>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$PAS_count_normalized,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 		#if in median +-10, refer to it as median, else, put it as it is
				 		
					 	Exist:
					 	{
					 		for(my $k=0;$k<=10;$k++)
					 		{
						 		if(exists $PAS_median{$gene}{$mode+$k})
						 		{
						 			$PASs{$cell}{$gene}{$mode+$k} =$PAS_count_normalized;
						 			$PASs_allCells{$gene}{$mode+$k}+=$PAS_count_normalized;
						 			
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
						 			$PASs{$cell}{$gene}{$mode-$k} =$PAS_count_normalized;
						 			$PASs_allCells{$gene}{$mode-$k}+=$PAS_count_normalized;
						 			if ( grep { $_ eq $mode-$k } @{$PAS_modes{$gene}} ) {
						 			}
						 			else
									{	push(@{$PAS_modes{$gene}},$mode-$k);}
						 			
									last Exist;
						 		}
					 		}
					 	
				 			$PASs{$cell}{$gene}{$mode} =$PAS_count_normalized;
				 			$PASs_allCells{$gene}{$mode}+=$PAS_count_normalized;
							if ( grep { $_ eq $mode } @{$PAS_modes{$gene}} ) {
						 	}
							else
				 			{push(@{$PAS_modes{$gene}},$mode);}
							
			
					 	}
				 	}
				 	close IN;
				 	
           	}
        }
        
        #calcluate shanon entropy
        #calculate relative expression
        while (my ($cell, $geneHash) = each %PASs)
        {
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{
        		my @modes = @{$PAS_modes{$gene}};
        	
	        	if(@modes >1) #only print out genes that have more than 1 polyA site
	        	{
	        		if((exists $PASs{'Liver'}{$gene} )&&(exists $PASs{'Brain'}{$gene})&&(exists $PASs{'Kidney'}{$gene} )) #expressed in both tissues
        		{
	        		while (my ($PAS,$normalized_count) = each %$PAShash )
	        		{
	        			$relative_exp{$cell}{$gene}{$PAS} = $normalized_count / $PASs_allCells{$gene}{$PAS};
	        		}
	        	}
	        	}
        	
        	}
        }
        
        #calculate entropy
        my %entropy_per_cell = ();
        while (my ($cell, $geneHash) = each %relative_exp)
        {
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{
        		while (my ($PAS,$rel_exp) = each %$PAShash )
        		{
        			$entropy_per_cell{$cell}{$gene}{$PAS} =$rel_exp * (log($rel_exp)/log(2));
        		}
        	}
        }
        
        my %entropy = ();
        while (my ($gene, $array) = each %PAS_modes)
        {
        	my @modes = @{$PAS_modes{$gene}};
        	
	        
        	if(@modes >1) #only print out genes that have more than 1 polyA site
        	{
        	
	        	for(my $i=0; $i < @{$array}; $i++)
	        	{
	        		my @PASarray = @{$array};
	        		my $PAS = $PASarray[$i];
	        		#for all cells
	        		while (my ($lane, $cellHash) = each %dataset)
			        {
			             while (my ($cell, $cell_file) = each %$cellHash)
			           	{
			           		
			           		if (exists $entropy_per_cell{$cell}{$gene}{$PAS})
		        			{
	        					$entropy{$gene}{$PAS} += $entropy_per_cell{$cell}{$gene}{$PAS};
		        			}
		        			else
			        		{
			        			$entropy{$gene}{$PAS} += 0.001 * (log(0.001)/log(2));
			        		}
		        			
			        	}
			        }
        		}
        
        	}
        }
            
        #calculate categorical tissue specificity
        my %tissue_specificity = ();
         while (my ($cell, $geneHash) = each %PASs)
        {
        	
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{
        		if((exists $PASs{'Liver'}{$gene} )&&(exists $PASs{'Brain'}{$gene})&&(exists $PASs{'Kidney'}{$gene} )) #expressed in both tissues
        		{
        		my @modes = @{$PAS_modes{$gene}};
        		
	        
	        	if(@modes >1) #only print out genes that have more than 1 polyA site
	        	{
	        	
	        		while (my ($PAS,$normalized_count) = each %$PAShash )
	        		{
	        			$tissue_specificity{$cell}{$gene}{$PAS} = (-1 * $entropy{$gene}{$PAS}) - (log($relative_exp{$cell}{$gene}{$PAS})/log(2));
	        		}
	        	}
        		}
        	
        	}
        }
        
        #print to files
        while (my ($cell, $geneHash) = each %tissue_specificity)
        {
        	my $relative_exp_file = "../../analysis/tissueSpecific/$cell.Entropy_2";
		 	open OUT, ">$relative_exp_file" or die "Can not open file: $relative_exp_file";	
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{
        		while (my ($PAS,$spec) = each %$PAShash )
        		{
        			print OUT ("$gene,$PAS,$tissue_specificity{$cell}{$gene}{$PAS}\n");
        		}
        	}
        	close OUT;
        }
        
        my $relative_exp_file = "../../analysis/tissueSpecific/Entropy_3cells_2";
		open OUT, ">$relative_exp_file" or die "Can not open file: $relative_exp_file";	
        print OUT ("PAS,entropy\n");
       
       
        my $geneFile = "../../analysis/tissueSpecific/GeneChrMap.txt";
		open IN, "<$geneFile" or die "Can not open file: $geneFile";
	    
	    my %geneChrMap=();
	    my $line = <IN>;   	
		while ( $line = <IN>)
		{
			chomp ($line);
			my ($geneID, $transcriptID, $chr) = split(/\s+/, $line);
			$geneChrMap{$geneID} = $chr;
			
		}
		close IN;		 
	
				
		my $EntropyBed = "../../analysis/tissueSpecific/Entropy_3Cells_2.bed";
		open OUT_2, ">$EntropyBed" or die "Can not open file: $EntropyBed";	
        
        while (my ($cell, $geneHash) = each %tissue_specificity)
        {
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{
        		while (my ($PAS,$spec) = each %$PAShash )
        		{
        			if((exists $tissue_specificity{'Liver'}{$gene} )&&(exists $tissue_specificity{'Brain'}{$gene})&&(exists $tissue_specificity{'Kidney'}{$gene} )) #expressed in both tissues
        		{
        			print OUT ("$PAS,$tissue_specificity{$cell}{$gene}{$PAS}\n");
        		 
        			#if($tissue_specificity{$cell}{$gene}{$PAS} <=3.2)
        			{
	        			if(exists $geneChrMap{$gene})
						{
							my $start = $PAS -1;
							my $num = eval sprintf('%.2f', $tissue_specificity{$cell}{$gene}{$PAS});
							print OUT_2 ("chr$geneChrMap{$gene}\t$start\t$PAS\t$num\n");
							
						}
						else {print ("$gene\t");}
        			}
        		}
        			
        		}
        	}
        	
        }
        
        my $EntropyTry = "../../analysis/tissueSpecific/Entropy_3Cells_try";
		open OUT_3, ">$EntropyTry" or die "Can not open file: $EntropyTry";	
        
       	while (my ($gene,$PAShash) = each %entropy)
        	{
        		while (my ($PAS,$spec) = each %$PAShash )
        		{
        			if((exists $tissue_specificity{'Liver'}{$gene} )&&(exists $tissue_specificity{'Brain'}{$gene})&&(exists $tissue_specificity{'Kidney'}{$gene} )) #expressed in both tissues
        		{
        			print OUT_3 ("$PAS,$entropy{$gene}{$PAS}\n");
        		 
        			#if($tissue_specificity{$cell}{$gene}{$PAS} <=3.2)
        			{
	        			if(exists $geneChrMap{$gene})
						{
							my $start = $PAS -1;
						#	my $num = eval sprintf('%.2f', $tissue_specificity{$cell}{$gene}{$PAS});
						#	print OUT_2 ("chr$geneChrMap{$gene}\t$start\t$PAS\t$num\n");
						}
						else {print ("$gene\t");}
        			}
        	
        			
        		}
        	}
        	
        }
        
        close OUT;
        close OUT_2;
        	 	
}

sub calculateShanonEntropy_old
{
	#read medians
		my $PAS_median_file ="../../analysis/tissueSpecific/PAS_median_3CellTypes";
		open IN, "<$PAS_median_file" or die "Can not open file: $PAS_median_file";
        my %PAS_median;
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
		while (my ($lane, $cellHash) = each %dataset)
        {
             while (my ($cell, $cell_file) = each %$cellHash)
           	{
           			
      				my $PAS_Genes_file=  $cell_file;
      				my $tissuePath = "analysis/tissueSpecific";
      				$PAS_Genes_file =~ s/analysis/$tissuePath/g;
                    $PAS_Genes_file =~ s/uniq/PAS.Genes.Adjusted_2/g;
                    open IN, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                  	
				 	while (my $line = <IN>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$PAS_count_normalized,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 		#if in median +-10, refer to it as median, else, put it as it is
				 		
					 	Exist:
					 	{
					 		for(my $k=0;$k<=10;$k++)
					 		{
						 		if(exists $PAS_median{$gene}{$mode+$k})
						 		{
						 			$PASs{$cell}{$gene}{$mode+$k} =$PAS_count_normalized;
						 			$PASs_allCells{$gene}{$mode+$k}+=$PAS_count_normalized;
						 			
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
						 			$PASs{$cell}{$gene}{$mode-$k} =$PAS_count_normalized;
						 			$PASs_allCells{$gene}{$mode-$k}+=$PAS_count_normalized;
						 			if ( grep { $_ eq $mode-$k } @{$PAS_modes{$gene}} ) {
						 			}
						 			else
									{	push(@{$PAS_modes{$gene}},$mode-$k);}
						 			
									last Exist;
						 		}
					 		}
					 	
				 			$PASs{$cell}{$gene}{$mode} =$PAS_count_normalized;
				 			$PASs_allCells{$gene}{$mode}+=$PAS_count_normalized;
							if ( grep { $_ eq $mode } @{$PAS_modes{$gene}} ) {
						 	}
							else
				 			{push(@{$PAS_modes{$gene}},$mode);}
							
			
					 	}
				 	}
				 	close IN;
				 	
           	}
        }
        
        #calcluate shanon entropy
        #calculate relative expression
        while (my ($cell, $geneHash) = each %PASs)
        {
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{
        		my @modes = @{$PAS_modes{$gene}};
        	
	        	if(@modes >1) #only print out genes that have more than 1 polyA site
	        	{
	        		if((exists $PASs{'Liver'}{$gene} )&&(exists $PASs{'Brain'}{$gene})&&(exists $PASs{'Kidney'}{$gene} )) #expressed in both tissues
        		{
	        		while (my ($PAS,$normalized_count) = each %$PAShash )
	        		{
	        			$relative_exp{$cell}{$gene}{$PAS} = $normalized_count / $PASs_allCells{$gene}{$PAS};
	        		}
	        	}
	        	}
        	
        	}
        }
        
        #calculate entropy
        my %entropy_per_cell = ();
        while (my ($cell, $geneHash) = each %relative_exp)
        {
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{
        		while (my ($PAS,$rel_exp) = each %$PAShash )
        		{
        			$entropy_per_cell{$cell}{$gene}{$PAS} =$rel_exp * (log($rel_exp)/log(2));
        		}
        	}
        }
        
        my %entropy = ();
        while (my ($gene, $array) = each %PAS_modes)
        {
        	my @modes = @{$PAS_modes{$gene}};
        	
	        
        	if(@modes >1) #only print out genes that have more than 1 polyA site
        	{
        	
	        	for(my $i=0; $i < @{$array}; $i++)
	        	{
	        		my @PASarray = @{$array};
	        		my $PAS = $PASarray[$i];
	        		#for all cells
	        		while (my ($lane, $cellHash) = each %dataset)
			        {
			             while (my ($cell, $cell_file) = each %$cellHash)
			           	{
			           		
			           		if (exists $entropy_per_cell{$cell}{$gene}{$PAS})
		        			{
	        					$entropy{$gene}{$PAS} += $entropy_per_cell{$cell}{$gene}{$PAS};
		        			}
		        			else
			        		{
			        			$entropy{$gene}{$PAS} += 0.001 * (log(0.001)/log(2));
			        		}
		        			
			        	}
			        }
        		}
        
        	}
        }
            
        #calculate categorical tissue specificity
        my %tissue_specificity = ();
         while (my ($cell, $geneHash) = each %PASs)
        {
        	
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{
        		if((exists $PASs{'Liver'}{$gene} )&&(exists $PASs{'Brain'}{$gene})&&(exists $PASs{'Kidney'}{$gene} )) #expressed in both tissues
        		{
        		my @modes = @{$PAS_modes{$gene}};
        		
	        
	        	if(@modes >1) #only print out genes that have more than 1 polyA site
	        	{
	        	
	        		while (my ($PAS,$normalized_count) = each %$PAShash )
	        		{
	        			$tissue_specificity{$cell}{$gene}{$PAS} = (-1 * $entropy{$gene}{$PAS}) - (log($relative_exp{$cell}{$gene}{$PAS})/log(2));
	        		}
	        	}
        		}
        	
        	}
        }
        
        #print to files
        while (my ($cell, $geneHash) = each %tissue_specificity)
        {
        	my $relative_exp_file = "../../analysis/tissueSpecific/$cell.Entropy_2";
		 	open OUT, ">$relative_exp_file" or die "Can not open file: $relative_exp_file";	
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{
        		while (my ($PAS,$spec) = each %$PAShash )
        		{
        			print OUT ("$gene,$PAS,$tissue_specificity{$cell}{$gene}{$PAS}\n");
        		}
        	}
        	close OUT;
        }
        
        my $relative_exp_file = "../../analysis/tissueSpecific/Entropy_3cells_2";
		open OUT, ">$relative_exp_file" or die "Can not open file: $relative_exp_file";	
        print OUT ("PAS,entropy\n");
       
       
        my $geneFile = "../../analysis/tissueSpecific/GeneChrMap.txt";
		open IN, "<$geneFile" or die "Can not open file: $geneFile";
	    
	    my %geneChrMap=();
	    my $line = <IN>;   	
		while ( $line = <IN>)
		{
			chomp ($line);
			my ($geneID, $transcriptID, $chr) = split(/\s+/, $line);
			$geneChrMap{$geneID} = $chr;
			
		}
		close IN;		 
	
				
		my $EntropyBed = "../../analysis/tissueSpecific/Entropy_3Cells_2.bed";
		open OUT_2, ">$EntropyBed" or die "Can not open file: $EntropyBed";	
        
        while (my ($cell, $geneHash) = each %tissue_specificity)
        {
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{
        		while (my ($PAS,$spec) = each %$PAShash )
        		{
        			if((exists $tissue_specificity{'Liver'}{$gene} )&&(exists $tissue_specificity{'Brain'}{$gene})&&(exists $tissue_specificity{'Kidney'}{$gene} )) #expressed in both tissues
        		{
        			print OUT ("$PAS,$tissue_specificity{$cell}{$gene}{$PAS}\n");
        		 
        			#if($tissue_specificity{$cell}{$gene}{$PAS} <=3.2)
        			{
	        			if(exists $geneChrMap{$gene})
						{
							my $start = $PAS -1;
							my $num = eval sprintf('%.2f', $tissue_specificity{$cell}{$gene}{$PAS});
							print OUT_2 ("chr$geneChrMap{$gene}\t$start\t$PAS\t$num\n");
							
						}
						else {print ("$gene\t");}
        			}
        		}
        			
        		}
        	}
        	
        }
        
        my $EntropyTry = "../../analysis/tissueSpecific/Entropy_3Cells_try";
		open OUT_3, ">$EntropyTry" or die "Can not open file: $EntropyTry";	
        
       	while (my ($gene,$PAShash) = each %entropy)
        	{
        		while (my ($PAS,$spec) = each %$PAShash )
        		{
        			if((exists $tissue_specificity{'Liver'}{$gene} )&&(exists $tissue_specificity{'Brain'}{$gene})&&(exists $tissue_specificity{'Kidney'}{$gene} )) #expressed in both tissues
        		{
        			print OUT_3 ("$PAS,$entropy{$gene}{$PAS}\n");
        		 
        			#if($tissue_specificity{$cell}{$gene}{$PAS} <=3.2)
        			{
	        			if(exists $geneChrMap{$gene})
						{
							my $start = $PAS -1;
						#	my $num = eval sprintf('%.2f', $tissue_specificity{$cell}{$gene}{$PAS});
						#	print OUT_2 ("chr$geneChrMap{$gene}\t$start\t$PAS\t$num\n");
						}
						else {print ("$gene\t");}
        			}
        	
        			
        		}
        	}
        	
        }
        
        close OUT;
        close OUT_2;
        	 	
}

sub calculateShanonEntropyForSignificantLME
{
	#read medians
		my $PAS_median_file ="../../analysis/tissueSpecific/PAS_median_3CellTypes";
		open IN, "<$PAS_median_file" or die "Can not open file: $PAS_median_file";
        my %PAS_median;
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
		while (my ($lane, $cellHash) = each %dataset)
        {
             while (my ($cell, $cell_file) = each %$cellHash)
           	{
           			
      				my $PAS_Genes_file=  $cell_file;
      				my $tissuePath = "analysis/tissueSpecific";
      				$PAS_Genes_file =~ s/analysis/$tissuePath/g;
                    $PAS_Genes_file =~ s/uniq/PAS.Genes.Adjusted_Final/g;
                    open IN, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                  	
				 	while (my $line = <IN>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$PAS_genelevel,$PAS_count_normalized,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 		#if in median +-10, refer to it as median, else, put it as it is
				 		
					 	Exist:
					 	{
					 		for(my $k=0;$k<=10;$k++)
					 		{
						 		if(exists $PAS_median{$gene}{$mode+$k})
						 		{
						 			$PASs{$cell}{$gene}{$mode+$k} =$PAS_count_normalized;
						 			$PASs_allCells{$gene}{$mode+$k}+=$PAS_count_normalized;
						 			
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
						 			$PASs{$cell}{$gene}{$mode-$k} =$PAS_count_normalized;
						 			$PASs_allCells{$gene}{$mode-$k}+=$PAS_count_normalized;
						 			if ( grep { $_ eq $mode-$k } @{$PAS_modes{$gene}} ) {
						 			}
						 			else
									{	push(@{$PAS_modes{$gene}},$mode-$k);}
						 			
									last Exist;
						 		}
					 		}
					 	
				 			$PASs{$cell}{$gene}{$mode} =$PAS_count_normalized;
				 			$PASs_allCells{$gene}{$mode}+=$PAS_count_normalized;
							if ( grep { $_ eq $mode } @{$PAS_modes{$gene}} ) {
						 	}
							else
				 			{push(@{$PAS_modes{$gene}},$mode);}
							
			
					 	}
				 	}
				 	close IN;
				 	
           	}
        }
        
        #calcluate shanon entropy
        #calculate relative expression
        while (my ($cell, $geneHash) = each %PASs)
        {
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{
        		my @modes = @{$PAS_modes{$gene}};
        		if((exists $PASs{'Liver'}{$gene} )&&(exists $PASs{'Brain'}{$gene})&&(exists $PASs{'Kidney'}{$gene} )) #expressed in both tissues
        		{
	        
	        	if(@modes >1) #only print out genes that have more than 1 polyA site
	        	{
	        		while (my ($PAS,$normalized_count) = each %$PAShash )
	        		{
	        			$relative_exp{$cell}{$gene}{$PAS} = $normalized_count / $PASs_allCells{$gene}{$PAS};
	        		}
	        	}
        		}
	        	
        	}
        }
        
        #calculate entropy
        my %entropy_per_cell = ();
        while (my ($cell, $geneHash) = each %relative_exp)
        {
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{
        		while (my ($PAS,$rel_exp) = each %$PAShash )
        		{
        			$entropy_per_cell{$cell}{$gene}{$PAS} =$rel_exp * (log($rel_exp)/log(2));
        		}
        	}
        }
        
        my %entropy = ();
        while (my ($gene, $array) = each %PAS_modes)
        {
        	my @modes = @{$PAS_modes{$gene}};
        	
	        
        	if(@modes >1) #only print out genes that have more than 1 polyA site
        	{
        	
	        	for(my $i=0; $i < @{$array}; $i++)
	        	{
	        		my @PASarray = @{$array};
	        		my $PAS = $PASarray[$i];
	        		#for all cells
	        		while (my ($lane, $cellHash) = each %dataset)
			        {
			             while (my ($cell, $cell_file) = each %$cellHash)
			           	{
			           		
			           		if (exists $entropy_per_cell{$cell}{$gene}{$PAS})
		        			{
	        					$entropy{$gene}{$PAS} += $entropy_per_cell{$cell}{$gene}{$PAS};
		        			}
		        			else
			        		{
			        			$entropy{$gene}{$PAS} += 0.001 * (log(0.001)/log(2));
			        		}
		        			
			        	}
			        }
	        	}	
        	}
        }
            
        #calculate categorical tissue specificity
        my %tissue_specificity = ();
         while (my ($cell, $geneHash) = each %PASs)
        {
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{
        		my @modes = @{$PAS_modes{$gene}};
        		if(@modes >1) #only print out genes that have more than 1 polyA site
	        	{
	        	if((exists $PASs{'Liver'}{$gene} )&&(exists $PASs{'Brain'}{$gene})&&(exists $PASs{'Kidney'}{$gene} )) #expressed in both tissues
        		{
	        		while (my ($PAS,$normalized_count) = each %$PAShash )
	        		{
	        			$tissue_specificity{$cell}{$gene}{$PAS} = (-1 * $entropy{$gene}{$PAS}) - (log($relative_exp{$cell}{$gene}{$PAS})/log(2));
	        		}
        		}
	        	}
        	}
        }
        
        #read the significant PASs
        my $pvaluesSig = "../../analysis/tissueSpecific/HighlyExpressedInOneCell/Brain_Liver_Kidney_significant.csv";
        open IN, "<$pvaluesSig" or die "Can not open file: $pvaluesSig";	
        my %significant_pvalues=();
        while (my $line = <IN>)
		{
			chomp ($line);
			my ($gene_p,$pas_p,$cell,$dif,$l,$k,$b ) = split(/\,/, $line);
			$significant_pvalues{$gene_p}{$pas_p} = $dif;
			
			#chomp ($line);
			#my @values = split(/\,/, $line);
			#my ($cell_p,$gene_p,$pas_p) = split(/\//, $values[0]);
			#$significant_pvalues{$cell_p}{$gene_p}{$pas_p} = $values[1];
		}
			
        
        
        #print to files
#        while (my ($cell, $geneHash) = each %tissue_specificity)
#        {
#        	my $relative_exp_file = "../../analysis/tissueSpecific/$cell.Entropy";
#		 	open OUT, ">$relative_exp_file" or die "Can not open file: $relative_exp_file";	
#        	while (my ($gene,$PAShash) = each %$geneHash)
#        	{
#        		while (my ($PAS,$spec) = each %$PAShash )
#        		{
#        			print OUT ("$gene,$PAS,$tissue_specificity{$cell}{$gene}{$PAS}\n");
#        		}
#        	}
#        	close OUT;
#        }

#TEST
my $geneFile = "../../analysis/tissueSpecific/GeneChrMap.txt";
		open IN, "<$geneFile" or die "Can not open file: $geneFile";
	    
	    my %geneChrMap=();
	    my $line = <IN>;   	
		while ( $line = <IN>)
		{
			chomp ($line);
			my ($geneID, $transcriptID, $chr) = split(/\s+/, $line);
			$geneChrMap{$geneID} = $chr;
			
		}
		close IN;	

#####
        
        my $sig_file = "../../analysis/tissueSpecific/Entropy_sig_3Cells_Final";
        my $nonsig_file = "../../analysis/tissueSpecific/Entropy_nonsig_3Cells_Final";
		
		open OUTSIG, ">$sig_file" or die "Can not open file: $sig_file";	
		open OUTNONSIG, ">$nonsig_file" or die "Can not open file: $nonsig_file";
			
        print OUTSIG ("PAS,entropy\n");
        print OUTNONSIG ("PAS,entropy\n");
        my $count=0;
        my $count_sig =0;
        my $count_shannon_sig = 0;
        while (my ($cell, $geneHash) = each %tissue_specificity)
        {
       		
       		my $geneHash = $tissue_specificity{$cell};
       		
        	while (my ($gene,$PAShash) = each %$geneHash)
        	{
        		while (my ($PAS,$spec) = each %$PAShash )
        		{
        			
        			if(exists $significant_pvalues{$gene}{$PAS})
        			{
        			
        				#if(exists $PASs{$cell}{$gene}{$PAS})
        				print OUTSIG ("$PAS,$tissue_specificity{$cell}{$gene}{$PAS}\n");
        				#if((exists $tissue_specificity{'Liver'}{$gene}{$PAS} )&&(exists $tissue_specificity{'Brain'}{$gene}{$PAS})&&(exists $tissue_specificity{'Kidney'}{$gene}{$PAS} )) #expressed in both tissues
        				{
        				if($tissue_specificity{$cell}{$gene}{$PAS} >1.5)
        				{
        					$count_shannon_sig ++;
        					#print ("$cell\t$gene\t$geneChrMap{$gene}\t$PAS\t$tissue_specificity{$cell}{$gene}{$PAS}\t$significant_pvalues{$gene}{$PAS}\n");
        				}
        				$count_sig++;
        				}
        			}
        			
        			else
        			{print OUTNONSIG ("$cell,$gene,$PAS,$tissue_specificity{$cell}{$gene}{$PAS}\n");
        				if($tissue_specificity{$cell}{$gene}{$PAS} <1)
        				{if((exists $tissue_specificity{'Liver'}{$gene}{$PAS} )&&(exists $tissue_specificity{'Brain'}{$gene}{$PAS})&&(exists $tissue_specificity{'Kidney'}{$gene}{$PAS} )) #expressed in both tissues
        				
        				{
        					#print ("$cell\t$gene\t$geneChrMap{$gene}\t$PAS,$tissue_specificity{$cell}{$gene}{$PAS}\n");
        				}
        				}
        			}
        		}
        			$count++;
        	
        	}
        	
        }
        	
        
        close OUT;
        
        print ("CountSig = $count_shannon_sig/$count_sig/$count\n");
       
       
       
				 	 	
}

sub MixedModel_allCellTypes  #wrong function 
{
	#read medians
		my $PAS_median_file ="../../analysis/tissueSpecific/PAS_median_all_cells_noCellLine";
		open IN, "<$PAS_median_file" or die "Can not open file: $PAS_median_file";
        my %PAS_median;
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
	#	my %PAS_modes =();
		while (my ($lane, $cellHash) = each %dataset)
        {
             while (my ($cell, $cell_file) = each %$cellHash)
           	{
           			
      				my $PAS_Genes_file=  $cell_file;
      				my $tissuePath = "analysis/tissueSpecific";
      				$PAS_Genes_file =~ s/analysis/$tissuePath/g;
                    $PAS_Genes_file =~ s/uniq/PAS.Genes.Adjusted/g;
                    open IN, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                  	
				 	while (my $line = <IN>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$PAS_count_normalize,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 		#if in median +-10, refer to it as median, else, put it as it is
				 			
					 	Exist:
					 	{
					 		for(my $k=0;$k<=10;$k++)
					 		{
						 		if(exists $PAS_median{$gene}{$mode+$k})
						 		{
						 			$PASs{$gene}{$cell}{$mode+$k} =$PAS_count;
						 			
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
						 			$PASs{$gene}{$cell}{$mode-$k} =$PAS_count;
						 			if ( grep { $_ eq $mode-$k } @{$PAS_modes{$gene}} ) {
						 			}
						 			else
									{	push(@{$PAS_modes{$gene}},$mode-$k);}
						 			
									last Exist;
						 		}
					 		}
					 	
				 			$PASs{$gene}{$cell}{$mode} =$PAS_count;
							if ( grep { $_ eq $mode } @{$PAS_modes{$gene}} ) {
						 	}
							else
				 			{push(@{$PAS_modes{$gene}},$mode);}
							
			
					 	}
				 	}
				 	close IN;
           	}
        }
        
        my $in2cell =0;
        my $in1cell = 0;
        my %geneIn1cell= ();
        my $brain=0;
        my $liver =0;
        
        #Add pseudocounts 0/1 for each PAS that does not exist in the data
        #If the gene is not expressed at all in one of the tissues, do not include *********************
        my %PAS_toprint=();
        my %PAS_constitutive = ();
        my $count = keys %PASs;
        
        my $gene_cons=0;
        my $gene_alternative=0;
        
        while (my ($gene, $cellHash) = each %PASs)
        {
        	#if(keys %{$PASs{$gene}} ==2) #expressed in both tissues
        	{
        		my @modes = sort(@{$PAS_modes{$gene}});
        		if(@modes >1) #only print out genes that have more than 1 polyA site
        		{
        			################################
        			while (my ($lane, $cellHash_dataset) = each %dataset)
        			{
             			while (my ($cell_dataset, $cell_file_dataset) = each %$cellHash_dataset)
           				{
	   		           		for(my $i=0; $i< @modes; $i++)
			            	{
			              		if(exists $PASs{$gene}{$cell_dataset}{$modes[$i]})
								{
									$PAS_toprint{$gene}{$cell_dataset}{$modes[$i]}=$PASs{$gene}{$cell_dataset}{$modes[$i]};
								}            		
								else
								{
									$PAS_toprint{$gene}{$cell_dataset}{$modes[$i]} = 0;
									
								}
			           		}
           				}
           				
        			}
		           
        		}
        	}
        }        	
        		
        #output mixed model files - original
        my $All_Genes_PAS = "../../analysis/tissueSpecific/PAS_all_cells_noCellLine_try";
       	#my $All_Genes_PAS = "../../analysis/tissueSpecific/gene";
       	printMixedModelFile2(\%PAS_modes,\%PAS_toprint,$All_Genes_PAS);
       	       	
       	#my $constitutive_PAS = "../../analysis/tissueSpecific/PAS_Brain_Liver_constitutive";
       	#printMixedModelFile2(\%PAS_modes,\%PAS_constitutive,$constitutive_PAS);
       	     
       	#permute
       	#permute_head(\%PAS_toprint);
}

sub MixedModel
{
	#read medians
		my $PAS_median_file ="../../analysis/3cells/PAS_median_Brain";#PreFinal/LiverKidney/PAS_median";
		
		#my $PAS_median_file ="../../analysis/tissueSpecific/PAS_median_3CellTypes";#PreFinal/LiverKidney/PAS_median";
		#my $PAS_median_file ="../../analysis/tissueSpecific/PAS_median_all_cells_noCellLine";
		open IN, "<$PAS_median_file" or die "Can not open file: $PAS_median_file";
        my %PAS_median;
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
	#	my %PAS_modes =();
		while (my ($lane, $cellHash) = each %dataset)
        {
             while (my ($cell, $cell_file) = each %$cellHash)
           	{
           			
      				my $PAS_Genes_file=  $cell_file;
      				my $tissuePath = "analysis/tissueSpecific";
      				$PAS_Genes_file =~ s/analysis/$tissuePath/g;
                    $PAS_Genes_file =~ s/uniq/PAS.Genes.Adjusted/g;
                    open IN, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                  	
				 	while (my $line = <IN>)
				 	{
				 		chomp ($line);
				 		my ($gene,$mode,$PAS_count,$PAS_count_normalize,$PAS_start,$PAS_end ) = split(/\,/, $line);
				 		#if in median +-10, refer to it as median, else, put it as it is
				 			
					 	Exist:
					 	{
					 		for(my $k=0;$k<=10;$k++)
					 		{
					 			
						 		if(exists $PAS_median{$gene}{$mode+$k})
						 		{
						 			$PASs{$gene}{$cell}{$mode+$k} =$PAS_count;
						 			
						 			if ( grep { $_ eq $mode+$k } @{$PAS_modes{$gene}} ) {
						 			}
						 			else
						 			{push(@{$PAS_modes{$gene}},$mode+$k);}
									last Exist;
						 		}
					 		}
					 		#print ("\n minus: ");
					 		for(my $k=10;$k>0;$k--)
					 		{
					 			
						 		if(exists $PAS_median{$gene}{$mode-$k})
						 		{
						 			$PASs{$gene}{$cell}{$mode-$k} =$PAS_count;
						 			if ( grep { $_ eq $mode-$k } @{$PAS_modes{$gene}} ) {
						 			}
						 			else
									{	push(@{$PAS_modes{$gene}},$mode-$k);}
						 			
									last Exist;
						 		}
					 		}
					 	
				 			$PASs{$gene}{$cell}{$mode} =$PAS_count;
							if ( grep { $_ eq $mode } @{$PAS_modes{$gene}} ) {
						 	}
							else
				 			{push(@{$PAS_modes{$gene}},$mode);}
							
			
					 	}
				 	}
				 	close IN;
           	}
        }
        
        my $in2cell =0;
        my $in1cell = 0;
        my %geneIn1cell= ();
        my $brain=0;
        my $liver =0;
        
        #Add pseudocounts 0/1 for each PAS that does not exist in the data
        #If the gene is not expressed at all in one of the tissues, do not include *********************
        my %PAS_toprint=();
        my %PAS_constitutive = ();
        my $count = keys %PASs;
        
        my $gene_cons=0;
        my $gene_alternative=0;
        
        my $cons_rest=0;
     
     	my %PAS_try=();
     	
     	my %gene_count_total=();#
        my %gene_count_cons = ();
        my %gene_count_cons_in_3cell=();
        my %gene_count_alt=();
        my %gene_count_alt_in_3cell=();#
        my %gene_count_all_in_3cell=();#
        
        while (my ($gene, $cellHash) = each %PASs)
        {
        	my @modes = sort(@{$PAS_modes{$gene}});
        	
        	if(!exists $gene_count_total{$gene})
			{
				$gene_count_total{$gene}= 1;
			}
			
        	if(keys %{$PASs{$gene}} ==3) #expressed in all tissues
        	{
        		if(!exists $gene_count_all_in_3cell{$gene})
				{
					$gene_count_all_in_3cell{$gene}= 1;
				}	
        		if(@modes >1) #only print out genes that have more than 1 polyA site
        		{
        			if(!exists $gene_count_alt_in_3cell{$gene})
					{
						$gene_count_alt_in_3cell{$gene}= 1;
					}
		
        			 while (my ($cell, $pas_hash) = each %$cellHash)
		           	{
		           		for(my $i=0; $i< @modes; $i++)
		            	{
		              		if(exists $PASs{$gene}{$cell}{$modes[$i]})
							{
								$PAS_toprint{$gene}{$cell}{$modes[$i]}=$PASs{$gene}{$cell}{$modes[$i]};
							}            		
							else
							{
								$PAS_toprint{$gene}{$cell}{$modes[$i]} = 0;
								
							}
		           		}
		           	}
		           	$gene_alternative++;
        		}
        	}
        	if(@modes >1)
        	{
        		if(!exists $gene_count_alt{$gene})
				{
					$gene_count_alt{$gene}= 1;
				}
        	}
        	
        	#constitutive sites doe not have to be expressed in the three cell types, they only have to have one mode
        	if(@modes == 1)
        		{
        			if(!exists $gene_count_cons{$gene})
					{
						$gene_count_cons{$gene}= 1;
					}
					if(keys %{$PASs{$gene}} ==3)
					{
						if(!exists $gene_count_cons_in_3cell{$gene})
						{
							$gene_count_cons_in_3cell{$gene}= 1;
						}
					}
        			#this should be the set of PASs that are constitutives
        			 while (my ($cell, $pas_hash) = each %$cellHash)
		           	{
		           		
	              		if(exists $PASs{$gene}{$cell}{$modes[0]})
						{
							if(!(exists $PAS_constitutive{$gene}{$modes[0]}))
							{
								$PAS_constitutive{$gene}{$modes[0]}=$PASs{$gene}{$cell}{$modes[0]};
								$gene_cons++;
								
								
							}
						}            		
		           	}
        			
        		}
        		
        		if(@modes == 1)
        		{
        			#this should be the set of PASs that are constitutives and expresses d in the 3 cells
        			if(keys %{$PASs{$gene}} ==3) #expressed in all tissues
        			{
        				$cons_rest++;
        				$PAS_try{$gene}{$modes[0]}=$PASs{$gene}{'Brain'}{$modes[0]};
        			}
        			
        		}
        		
        		
        	
        
        	
        	
        }   
       
       print ("number of constitutive everywhere $gene_cons\n");
       print ("number of restricted cons $cons_rest\n");
       

        
        my $count_total = keys %gene_count_total;
        my $count_total_in3 = keys %gene_count_all_in_3cell;
        my $count_cons = keys %gene_count_cons;
        my $count_cons_3cell = keys %gene_count_cons_in_3cell;
        my $count_alt= keys %gene_count_alt;
        my $count_alt_3cell = keys %gene_count_alt_in_3cell;
        
        
        print ("number of gene count total = $count_total\n");
       	print ("number of gene count total expressed in 3 cell = $count_total_in3 \n");
       	print ("number of gene constitutive $count_cons \n");
       	print ("number of gene constitutive in 3 cell = $count_cons_3cell \n");
       	print ("number of gene count alternative = $count_alt \n");
       	print ("number of gene count alternative in 3 cell = $count_alt_3cell \n");
   
       
       
        #output mixed model files - original
        my $All_Genes_PAS = "../../analysis/tissueSpecific/PAS_try";#PreFinal/LiverKidney/PAS_Alternative";
       	#my $All_Genes_PAS = "../../analysis/tissueSpecific/gene";
       	printMixedModelFile2(\%PAS_modes,\%PAS_toprint,$All_Genes_PAS);
       	       	
       	my $constitutive_PAS = "../../analysis/tissueSpecific/PAS_cons_try";
       	#printMixedModelFile2_constitutive(\%PAS_modes,\%PAS_try,$constitutive_PAS);
       	     
       	#permute
       	#If two cell line (Brain and liver)
       	#permute_head_2cellTypes(\%PAS_toprint);
       	#if more than 2, all cell types
       	#If the dataset is too large, we need to sample the PASs?
       	#permute_head_allCellTypes(\%PAS_toprint);
}
sub permute_head_allCellTypes
{
	my ($PASs) = @_;	
	
       	#permute
       	my $noOfPermutation = 1000;
       	#put the cell types in anarray
       	my @tissues;
       	while (my ($cell, $array_type_Hash) = each %cellTypes_hash)
       	{
       		push(@tissues, $cell);
       	}
       	
       	
       	srand (time ^ $$ ^ unpack "%L*", `ps axww | gzip`);
       	my $newPASs;
       	for (my $i=1; $i<= $noOfPermutation; $i++)
       	{
       		#Permute tissues (generate permuted assignment)
       		my $permutedTissues = permuteTissues(\@tissues);
       		
       		#select PAS samples to use
       		my $PAS_to_permute = $PASs;   #samplePas($PASs);
       		
       		$newPASs={};	
       		my $file = "../../analysis/tissueSpecific/permutedFiles/PAS.$i";
       		
       		while (my ($gene, $cellHash) = each %$PAS_to_permute)
    		{
    			for(my $i=0; $i<@tissues;$i++)
    			{
    				    			
	    			permute ($tissues[$i],$permutedTissues->[$i], $PAS_to_permute, $newPASs, $gene);	
    			}		
    		}     		
	       	#print
		
	       	printMixedModelFile2(\%PAS_modes,$newPASs,$file);
	       	
       		
       
      }
}
sub permuteTissues()
{
	my ($tissues) = @_;
       	
	my @permuted_tissues = random_permutation(@$tissues);
	#print (@$tissues);
	#print ("\n");
	#print (@permuted_tissues);
	return \@permuted_tissues;
}

sub samplePas()
{
	#$PAS_toprint{$gene}{$cell}{$modes[$i]} = 0;
	
	my ($PASs) = @_;
	
	srand(time() ^ $$ ^ unpack "%32L*", `ps axww | gzip`);#set the random number seed for the rand operator
	
	my $total = keys %$PASs;
	
	my $permute_total = int(0.3 * $total);
	my $other_total =  $total - $permute_total;
	
#	print ("$total\t$permute_total\t$other_total\n");
 
 	my @ipartition = sort {-1+int rand 3} (0..$total-1);
#	print ("partition @ipartition \n");
	
	my %permute_h=();
	
	for(my $i=0; $i<$permute_total; $i++)
    {
   		$permute_h{$ipartition[$i]} = 1;
    }
	

	my $count =0;
	my $PAS_sample =  ();
    for my $key (keys %$PASs)
    {
    	#if ( grep { $_ eq $count} @ipartition)
    	if (exists $permute_h{$count} )
    	{
    		$PAS_sample->{$key} = $PASs->{$key};
    	}
    	$count++;
    }
    
    $total = keys %$PAS_sample;
 #   print ("total sample $total\n");
    return $PAS_sample;
}


sub permute_head_2cellTypes
{
	my ($PASs) = @_;	

		          
       	#permute
       	my $noOfPermutation = 1000;
       	
       	srand (time ^ $$ ^ unpack "%L*", `ps axww | gzip`);
       	my $newPASs;
       	for (my $i=1; $i<= $noOfPermutation; $i++)
       	{
       		$newPASs={};	
       		my $file = "../../analysis/tissueSpecific/PreFinal/LiverKidney/permutedFiles/PAS.$i";
       		
       		#randomly select a number, if it is <0.5, each tissue will be assigned polya from its distributopn, else from the other tissue's destribution
	    	while (my ($gene, $cellHash) = each %$PASs)
    		{
	    		my $random_number = rand();
	    
	    		if($random_number < 0.5)
	    		{
	    			permute("Kidney", "Kidney", $PASs, $newPASs, $gene);	
	    			permute("Liver", "Liver", $PASs, $newPASs,$gene);		
	    		}
	    		else
	    		{
	    			permute("Kidney", "Liver", $PASs, $newPASs,$gene);	
	    			permute("Liver", "Kidney", $PASs, $newPASs,$gene);		
	    		}
						
    		}
       		
       	#print
	
       	printMixedModelFile2(\%PAS_modes,$newPASs,$file);
       		
       
      }
}
      
sub permute
{
	my ($tissue1,$tissue2, $PAS_hash, $newPAS_hash, $gene)  = @_;
	
	#get parameters of the Liver
	 #draw a random dist from the liver, based on the count of the brain
	 #exchange
	 
	 #random assign
	 my @tissue1_count = values %{$PAS_hash ->{$gene}{$tissue1}};
	 my $tissue1_total = sum(@tissue1_count);
	 
	 #get parameters of the multinomial for tissue2, sort according to the sorting in modes
	 my @tissue2_parameters=();
	 
	 my @modes = sort(@{$PAS_modes{$gene}});
	
	 for (my $m =0; $m <@modes; $m++)
	 {
	 	push( @tissue2_parameters, $PAS_hash->{$gene}{$tissue2}{$modes[$m]});
	 }
	 
	 my $tissue2_total = sum(@tissue2_parameters);
	#print ("Tissue2 total $tissue2_total\n");	
		
		
	
	 for(my $f =0; $f<@tissue2_parameters; $f++)
	 {
	 	$tissue2_parameters[$f] = $tissue2_parameters[$f] / $tissue2_total;
	 	
	 }
	#print ("$tissue1\t$tissue2\t$gene\t$tissue1_total\t@tissue2_parameters \n");
	
	my $max = @tissue2_parameters;
	my $lastIsZero = 0;
	while($tissue2_parameters[$max-1] ==0)
	{
		$lastIsZero ++;
		$max --;
	}
	if($lastIsZero == 0)
	{	 my @permuted_PAS =  random_multinomial($tissue1_total,@tissue2_parameters);
		 #print ("$gene,$tissue1,$tissue2,$tissue1_total,@permuted_PAS\n");
		
		for (my $m =0; $m <@modes; $m++)
		{
		 	$newPAS_hash->{$gene}{$tissue1}{$modes[$m]}=$permuted_PAS[$m];
		 	
		}
	}
	else
	{
	#	print ("in else $max  $lastIsZero\n");
		my @new_tissue2_parameters = splice(@tissue2_parameters, 0,$max);
		if(@new_tissue2_parameters == 1)
		{
			for (my $m =0; $m <@modes; $m++)
			{
				if($m < @new_tissue2_parameters)
			 	{
			 		$newPAS_hash->{$gene}{$tissue1}{$modes[$m]}=$tissue1_total;
			 	}
			 	else
			 	{$newPAS_hash->{$gene}{$tissue1}{$modes[$m]}=0;}
			}
			
		}
		else
		{
			my @permuted_PAS =  random_multinomial($tissue1_total,@new_tissue2_parameters);
			 #print ("$gene,$tissue1,$tissue2,$tissue1_total,@permuted_PAS\n");
		#	print ("new param @new_tissue2_parameters \n");
			for (my $m =0; $m <@modes; $m++)
			{
				if($m < @new_tissue2_parameters)
			 	{$newPAS_hash->{$gene}{$tissue1}{$modes[$m]}=$permuted_PAS[$m];
			 	}
			 	else
			 	{$newPAS_hash->{$gene}{$tissue1}{$modes[$m]}=0;}
			}
		}
	}

}    
    
sub printMixedModelFile2_constitutive()
{
	my ($PAS_modes_hash,$PASs_hash,$filename) = @_;
	my %PAS_modes = %$PAS_modes_hash;
	my %PASs = %$PASs_hash;
	
	open OUT_PAS, ">$filename" or die "Can not open file: $filename";
      
    print OUT_PAS ("gene,pas\n");
       	
       	my $count =0;
        while (my ($gene, $modeHash) = each %PASs)
        {
        	my @modes = sort(@{$PAS_modes{$gene}});
        	print OUT_PAS ("$gene,$modes[0]\n");
        	
        }
	close OUT_PAS;
}
sub printMixedModelFile2()
{
	my ($PAS_modes_hash,$PASs_hash,$filename) = @_;
	my %PAS_modes = %$PAS_modes_hash;
	my %PASs = %$PASs_hash;
	
	open OUT_PAS, ">$filename" or die "Can not open file: $filename";
      
    print OUT_PAS ("tissue,gene,pas,count\n");
       	
       	my $count =0;
        while (my ($gene, $modeHash) = each %PASs)
        {
        	my @modes = sort(@{$PAS_modes{$gene}});
        	#if(@modes >1) #only print out genes that have more than 1 polyA site
        	{
        			
        		for(my $i=0; $i< @modes; $i++)
	            {
		           	while (my ($cell, $modeHash) = each %{$PASs{$gene}})
		           	{
		           		
		            		print OUT_PAS ("$cell,$gene,$modes[$i],");
		              		if(exists $PASs{$gene}{$cell}{$modes[$i]})
							{
								if($PASs{$gene}{$cell}{$modes[$i]} == 0)
								{
									print OUT_PAS ("1\n");
									#print ("$gene is zero\n");
								}
								else
								{
									my $count_pesudo = $PASs{$gene}{$cell}{$modes[$i]}+1;
									
									print OUT_PAS ("$count_pesudo\n");
								}
							}  
							    		
							
		           	}
		         }
        	
        	}
        }
	close OUT_PAS;
}

sub printMixedModelFile()
{
	my ($PAS_modes_hash,$PASs_hash,$filename) = @_;
	my %PAS_modes = %$PAS_modes_hash;
	my %PASs = %$PASs_hash;
	
	open OUT_PAS, ">$filename" or die "Can not open file: $filename";
      
    print OUT_PAS ("tissue,gene,pas,count\n");
       	
        while (my ($gene, $modeHash) = each %PAS_modes)
        {
        	my @modes = sort(@{$PAS_modes{$gene}});
        	if(@modes >1) #only print out genes that have more than 1 polyA site
        	{
	           	while (my ($cell, $modeHash) = each %{$PASs{$gene}})
	           	{
	           		for(my $i=0; $i< @modes; $i++)
	            	{
	            		print OUT_PAS ("$cell,$gene,$modes[$i],");
	              		if(exists $PASs{$gene}{$cell}{$modes[$i]})
						{
							print OUT_PAS ("$PASs{$gene}{$cell}{$modes[$i]}\n");
						}            		
						else
						{
							print OUT_PAS ("0\n");			
						}
	           		}
	           	}
        	}
        }
	close OUT_PAS;
}

sub computeDifferences_allCellTypes
{
	
	my $infile_path = "resultFiles/PAS";
	
	my $permutation_start = 1;
	my $permutation_end = 1000;
	my %differences=(); 	
	for (my $i=$permutation_start; $i<$permutation_end;$i++)
	{
		my $differenceFile = $infile_path . ".$i.residual";
		open DIFF, "<$differenceFile" or die "Can not open file: $differenceFile";
		my $line1 = <DIFF>;  
		while ( $line1 = <DIFF>)
		{
			chomp ($line1);
			my $dif = $line1;
			#my $line2 = <DIFF>;
			#chomp($line2);
			#my $dif = abs($line1-$line2);
			$differences{$dif} ++;
		}

		close DIFF;
	}
	
	
	
	my $diffOut = "permuted_dif_allCellTypes.csv";
	open OUT, ">$diffOut" or die "Can not open file: $diffOut";
	print OUT ("Diff\n");	
	
	my $sig = 0;
	my $total=0;
	
	while (my ($dif, $count) = each %differences)
    {
    	for(my $c = 0; $c <$count; $c++)
    	{
    		print OUT ("$dif\n");
    	}	
    }
       	
    close OUT;
   # my $perc = $sig / $total;
    
   # print ("sig pas = $sig/$total = $perc\n")
}

sub computeDifferences_permuted
{
	
	my $infile_path = "resultFiles/PAS";
	
	my $permutation_start = 1;
	my $permutation_end = 1000;
	my %differences=(); 	
	for (my $i=$permutation_start; $i<$permutation_end;$i++)
	{
		my $differenceFile = $infile_path . ".$i.residual";
		open DIFF, "<$differenceFile" or die "Can not open file: $differenceFile";
		my $line1 = <DIFF>;  
		while ( $line1 = <DIFF>)
		{
			chomp ($line1);
			my $line2 = <DIFF>;
			chomp($line2);
			my $dif = abs($line1-$line2);
			$differences{$dif} ++;
		}

		close DIFF;
	}
	
	
	
	my $diffOut = "permuted_dif.csv";
	open OUT, ">$diffOut" or die "Can not open file: $diffOut";
	print OUT ("Diff\n");	
	
	my $sig = 0;
	my $total=0;
	
	while (my ($dif, $count) = each %differences)
    {
    	for(my $c = 0; $c <$count; $c++)
    	{
    		print OUT ("$dif\n");
    	}	
    }
       	
    close OUT;
   # my $perc = $sig / $total;
    
   # print ("sig pas = $sig/$total = $perc\n")
}

sub getSignificantDifferences_allCellTypes
{
	#For downregulated < -2.125, UPregulated > 2.75
	my $differenceFile = "../../analysis/tissueSpecific/PAS_expressed_genes_only_out";
	open DIFF, "<$differenceFile" or die "Can not open file: $differenceFile";
	
	
	my $PasFile = "../../analysis/tissueSpecific/PAS_expressed_genes_only";
	open PAS, "<$PasFile" or die "Can not open file: $PasFile";

	my $pasLine = <PAS>;
	my $line = <DIFF>;  
	my %differences=(); 
	my %geneStats = ();	
	
	while ( $line = <DIFF>)
	{
		chomp ($line);
		$pasLine = <PAS>;
		my ($cell, $gene, $polya, $count) = split(/\,/, $pasLine);
		
		$differences{$cell}{$gene}{$polya}{"dif"} = $line;
		$geneStats{$gene}{"total"} += $line;
		$geneStats{$gene}{"count"}++;
		$geneStats{$gene}{$cell}{$polya} = $line;
	}
	close PAS;
	close DIFF; 
	
	my $geneFile = "../../analysis/tissueSpecific/GeneChrMap.txt";
	open IN, "<$geneFile" or die "Can not open file: $geneFile";
    my %geneChrMap=();
    my $line = <IN>;   	
	while ( $line = <IN>)
	{
		chomp ($line);
		my ($geneID, $transcriptID, $chr) = split(/\s+/, $line);
		$geneChrMap{$geneID} = $chr;
		
	}
	close IN;
	
	
	while (my ($cell, $geneHash) = each %differences)
	{
		my $file = "../../analysis/tissueSpecific/".$cell."_PAS_significant.csv";
		open OUT_sig, ">$file" or die "Can not open file: $file";
		
		my $file_bed = "../../analysis/tissueSpecific/".$cell."_PAS_significant.bed";
		open OUT_bed, ">$file_bed" or die "Can not open file: $file_bed";
	
		while (my ($gene, $pasHash) = each %$geneHash)
	    {
	    	while (my ($pas, $difHash) = each %$pasHash)
	       	{ 
       			if($pasHash->{$pas}{'dif'} >2.75)
       			{
       				if(exists $geneChrMap{$gene})
					{
						my $start = $pas -1;
						print OUT_bed ("chr$geneChrMap{$gene}\t$start\t$pas\t$pasHash->{$pas}{'dif'}\n");
						
					}
					
	       			print OUT_sig ("$gene,$pas,$pasHash->{$pas}{'dif'}\n");
       			}
	       		
	       	}
	    }
	    close OUT_sig;
	    close OUT_bed;
	}
	
	calculateZscore(\%geneStats);
	
	
}


sub getZscore_original_3cellTypes
{
	my $differenceFile = "../../analysis/tissueSpecific/PAS_Brain_Liver_Kidney_out";
	open DIFF, "<$differenceFile" or die "Can not open file: $differenceFile";
	
	
	my $PasFile = "../../analysis/tissueSpecific/PAS_Brain_Liver_Kidney";
	open PAS, "<$PasFile" or die "Can not open file: $PasFile";

	my $pasLine = <PAS>;
	my $line = <DIFF>;  
	my %differences=(); 
	my %geneStats = ();	
	
	
	
	while ( $line = <DIFF>)
	{
		chomp ($line);
		$pasLine = <PAS>;
		my ($cell, $gene, $polya, $count) = split(/\,/, $pasLine);
		
		$differences{$cell}{$gene}{$polya}{"dif"} = $line;
		$geneStats{$gene}{$polya}{"total"} += $line;
		$geneStats{$gene}{$polya}{"count"}++;
		$geneStats{$gene}{$cell}{$polya} = $line;
	}
	close PAS;
	close DIFF; 
	
	my $geneFile = "../../analysis/tissueSpecific/GeneChrMap.txt";
	open IN, "<$geneFile" or die "Can not open file: $geneFile";
    my %geneChrMap=();
    my $line = <IN>;   	
	while ( $line = <IN>)
	{
		chomp ($line);
		my ($geneID, $transcriptID, $chr) = split(/\s+/, $line);
		$geneChrMap{$geneID} = $chr;
		
	}
	close IN;
	
	
	while (my ($cell, $geneHash) = each %differences)
	{
		my $file = "../../analysis/tissueSpecific/".$cell."_PAS_significant.csv";
		open OUT_sig, ">$file" or die "Can not open file: $file";
		
		my $file_bed = "../../analysis/tissueSpecific/".$cell."_PAS_significant.bed";
		open OUT_bed, ">$file_bed" or die "Can not open file: $file_bed";
	
		while (my ($gene, $pasHash) = each %$geneHash)
	    {
	    	while (my ($pas, $difHash) = each %$pasHash)
	       	{ 
       			if($pasHash->{$pas}{'dif'} >2.75)
       			{
       				if(exists $geneChrMap{$gene})
					{
						my $start = $pas -1;
						print OUT_bed ("chr$geneChrMap{$gene}\t$start\t$pas\t$pasHash->{$pas}{'dif'}\n");
						
					}
					
	       			print OUT_sig ("$gene,$pas,$pasHash->{$pas}{'dif'}\n");
       			}
	       		
	       	}
	    }
	    close OUT_sig;
	    close OUT_bed;
	}
	
	calculateZscore(\%geneStats);
	
	
}


sub calculateZscore
{
	#mu = 	 (sum  x) / n
	#z = (x-mu)/sigma
	#to get standard deviation of each gene:
	#subtract the mean from every number to get the list of deviations. 
	#Create a list of these numbers. It's OK to get negative numbers here. 
	#Next, square the resulting list of numbers
	#Add up all of the resulting squares to get their total sum. 
	#Divide your result by one less than the number of items in the list.
	
	my ($geneStats) = @_;
	
	while (my ($gene, $meanHash) = each %$geneStats)
	{
		$geneStats->{$gene}{"mean"} = $geneStats->{$gene}{"total"} / $geneStats->{$gene}{"count"};
		
	}
	my %deviation;
	while (my ($gene, $cellHash) = each %$geneStats)
	{
		#print ("GENE $gene\t");
		while (my ($cell, $pasHash) = each %$cellHash)
		{
			
			if(($cell ne "mean")&&($cell ne "total")&&($cell ne "count")&&($cell ne "sumDeviation"))
			{
				while (my ($pas, $residHash) = each %$pasHash)
	       		{
	       			$deviation{$gene}{$cell}{$pas} = ($pasHash->{$pas} - $geneStats->{$gene}{"mean"}) ** 2; 
	       			$geneStats->{$gene}{"sumDeviation"} += $deviation{$gene}{$cell}{$pas} ;
	       			
			
	       		}			
			}
		}
	}
	while (my ($gene, $meanHash) = each %$geneStats)
	{
		$geneStats->{$gene}{"variance"} = $geneStats->{$gene}{"sumDeviation"} / ($geneStats->{$gene}{"count"}-1);
		$geneStats->{$gene}{"standardDeviation"} = sqrt($geneStats->{$gene}{"variance"}); 
	}
	
	my $m =  $geneStats->{"ENSG00000221869"}{"mean"};
	print ("MEAN: $m\n");
	my $var = $geneStats->{"ENSG00000221869"}{"standardDeviation"};
	print ("variance $var\n");
	print ("Total: $geneStats->{\"ENSG00000221869\"}{\"total\"}\n");
	print ("COUNT: $geneStats->{\"ENSG00000221869\"}{\"count\"}\n");
}

sub by_number {
    if ($a < $b){ -1 } elsif ($a > $b) { 1 } else { 0 }
}

sub getMaxDifFromMedian
{
	my ($handle,$original,$maxCell,@values) = @_;
	my $median;
	my $mid = int @values/2;

	my @cellTypes = ("Liver", "Kidney", "Brain");
	my @permutation = sort { $values[$a] <=> $values[$b] } (0..$#values);
	
	my @sorted_values = sort by_number @values;
	my $medianCell;
	@cellTypes = @cellTypes[@permutation];
	
	
	if (@values % 2) {
    $median = $sorted_values[ $mid ];
	$medianCell = $cellTypes[$mid];
	} 
#	else {
 #   $median = ($sorted_values[$mid-1] + $sorted_values[$mid])/2;
  #  $medianCell = $cellTypes[1];
	#} 
	
	my $maxDif = abs($sorted_values[2] - $median);
	$$maxCell = $cellTypes[2];
	

	if($original eq "1")
	{
	print $handle ("$$maxCell,$medianCell,$maxDif\n");
	}
	else
	{
		print $handle ("$maxDif\n");
		
	}
	return $maxDif;
}

sub getMinCell
{
	my ($handle,$original,$minCell,@values) = @_;
	my $median;
	my $mid = int @values/2;

	my @cellTypes = ("Liver", "Kidney", "Brain");
	my @permutation = sort { $values[$a] <=> $values[$b] } (0..$#values);
	
	my @sorted_values = sort by_number @values;
	@cellTypes = @cellTypes[@permutation];
	
	#check signs for different conditions
	my $mean = ($sorted_values[2]+$sorted_values[1])/2;
	my $dif = abs($mean - $sorted_values[0]);
	
	$$minCell = $cellTypes[0];
	

	if($original eq "1")
	{
		print $handle ("$$minCell,$dif\n");
	}
	else
	{
		print $handle ("$dif\n");
		
	}
	return $dif;
}


sub computeDifferences_original
{
	my $differenceFile = "../../analysis/tissueSpecific/original_residuals.csv";#original_residuals.csv";
	open DIFF, "<$differenceFile" or die "Can not open file: $differenceFile";
	
	my $sigFile = "../../analysis/tissueSpecific/significant_final.csv";#original_residuals.csv";
	#open OUT_sig, ">$sigFile" or die "Can not open file: $sigFile";
	
	my $commonFile = "../../analysis/tissueSpecific/nonsignificant_final.csv";#original_residuals.csv";
	open OUT_common, ">$commonFile" or die "Can not open file: $commonFile";
	
	my $PasFile = "../../analysis/tissueSpecific/PAS_Brain_Liver";
	open PAS, "<$PasFile" or die "Can not open file: $PasFile";
	#line1 liver, line2 brain
	my $pasLine1 = <PAS>;
	my $pasLine2;
	my $line1 = <DIFF>;  
	my %differences=(); 	
	while ( $line1 = <DIFF>)
	{
		chomp ($line1);
		my $line2 = <DIFF>;
		chomp($line2);
		my $dif = abs($line1-$line2);
		$pasLine1 = <PAS>;
		my ($cell1, $gene1, $polya1, $count1) = split(/\,/, $pasLine1);
		$pasLine2 = <PAS>;
		my ($cell2, $gene2, $polya2, $count2) = split(/\,/, $pasLine2);
		
		$differences{$gene1}{$polya1}{"dif"} = $dif;
		$differences{$gene1}{$polya1}{"Brain"} = $line2;
		$differences{$gene1}{$polya1}{"Liver"} = $line1;		
	}
	close PAS;
	close DIFF;
	
	my $diffOut = "../../analysis/tissueSpecific/original_dif.csv";
	open OUT, ">$diffOut" or die "Can not open file: $diffOut";
	print OUT ("Gene,PAS,Diff\n");	
	#print OUT_sig ("Gene,PAS,Diff,Liver,Brain\n");	
	print OUT_common ("Gene,PAS,Diff,Liver,Brain\n");	
	my $sig = 0;
	my $total=0;
	
	my %significant_genes=();
	while (my ($gene, $pasHash) = each %differences)
    {
    	while (my ($pas, $difHash) = each %$pasHash)
       	{ 
       		print OUT ("$gene,$pas,$pasHash->{$pas}{'dif'}\n");
       		if($pasHash->{$pas}{'dif'} >2.475)
       		{
       			#print OUT_sig ("$gene,$pas,$pasHash->{$pas}{'dif'},$pasHash->{$pas}{'Liver'},$pasHash->{$pas}{'Brain'},\n");
       			
       			
       			if(! exists $significant_genes{$gene})
       			{$significant_genes{$gene} = 1;
       				$sig ++;
       			}
       		}
       		else#if($pasHash->{$pas}{'dif'} <0.4)
       		{
       			print OUT_common ("$gene,$pas,$pasHash->{$pas}{'dif'},$pasHash->{$pas}{'Liver'},$pasHash->{$pas}{'Brain'},\n");
       			
       		}
       		$total++;
       	}
       	
    }
    close OUT;
    my $perc = $sig / $total;
    
    print ("sig pas = $sig/$total = $perc\n")
}

sub computeMaxDifFromMedian_original_3cellTypes
#HighlyExpressedInOneCell
{
	my $differenceFile = "../../analysis/tissueSpecific/PAS_Brain_Liver_Kidney_out";
	open DIFF, "<$differenceFile" or die "Can not open file: $differenceFile";
		
	my $PasFile = "../../analysis/tissueSpecific/PAS_Brain_Liver_Kidney";
	open PAS, "<$PasFile" or die "Can not open file: $PasFile";
	
	my $sigFile = "../../analysis/tissueSpecific/HighlyExpressedInOneCell/Brain_Liver_Kidney_significant.csv";
	open SIG, ">$sigFile" or die "Can not open file: $sigFile";
	print SIG ("Gene,PAS,maxCell,maxDif,Liver,Kidney,Brain\n");
	
	my $all = "../../analysis/tissueSpecific/HighlyExpressedInOneCell/Heart_Liver_Kidney.csv";
	open ALL, ">$all" or die "Can not open file: $all";
	print ALL ("Gene,PAS,maxCell,maxDif,Liver,Kidney,Brain\n");
	
	#line1 liver, line2 kidney, line3 brain
	my $pasLine = <PAS>;
	my $pasLine2;
	my $pasLine3;
	my $line1 = <DIFF>; 	
	
	my %residual = ();
	my $noOfCells = 3;
	my @zscores=();
	
	#print dif in a file
	my $diffOut = "../../analysis/tissueSpecific/HighlyExpressedInOneCell/Brain_Liver_Kidney_Diff_computed_original.csv";
	open my $OUT, ">$diffOut" or die "Can not open file: $diffOut";
	print $OUT ("maxCell,MedianCell,Diff\n");	
	
	while ( $line1 = <DIFF>)
	{
		my @resid=();
		chomp ($line1);
		push @resid,$line1;
		
		for(my $i=0; $i < $noOfCells-1; $i++)
		{
			$line1=<DIFF>;
			chomp($line1);
			push @resid,$line1;
			
		}
		my @cellTypes = ("Liver", "Kidney", "Brain");
		my $type = "1";
		my $maxCell;
		my $maxDif = getMaxDifFromMedian($OUT, $type,\$maxCell ,@resid);
		
		$pasLine = <PAS>;
		my ($cell, $gene, $polya, $count) = split(/\,/, $pasLine);
		$pasLine=<PAS>;
		$pasLine=<PAS>;
		if($maxDif > 2.503)
		{
			print SIG ("$gene,$polya,$maxCell,$maxDif,$resid[0],$resid[1],$resid[2]\n");
		}
		print ALL ("$gene,$polya,$maxCell,$maxDif,$resid[0],$resid[1],$resid[2]\n");
	}
	close $OUT;
	close DIFF;
	close PAS;
	close SIG;
}

sub computeMaxDifFromMedian_permuted_3cellTypes
{
	
	my $infile_path = "../../analysis/tissueSpecific/resultFiles/PAS";
	
	my $permutation_start = 1;
	my $permutation_end = 1000;
	my %differences=(); 
	
	my $diffOut = "../../analysis/tissueSpecific/HighlyExpressedInOneCell/Brain_Liver_Kidney_Diff_computed_permuted.csv";
	open my $OUT, ">$diffOut" or die "Can not open file: $diffOut";
	print $OUT ("Diff\n");	
	my $noOfCells = 3;
	
	for (my $i=$permutation_start; $i<$permutation_end;$i++)
	{
		my $differenceFile = $infile_path . ".$i.residual";
		open DIFF, "<$differenceFile" or die "Can not open file: $differenceFile";
		my $line1 = <DIFF>; 
		
		while ( $line1 = <DIFF>)
		{
			my @resid=();
			chomp ($line1);
			push @resid,$line1;
			
			for(my $i=0; $i < $noOfCells-1; $i++)
			{
				$line1=<DIFF>;
				chomp($line1);
				push @resid,$line1;
				
			}
			my $type = "0";
			my $maxCell;
			my $maxDif = getMaxDifFromMedian($OUT,$type,$maxCell,@resid);
		
		
		}

		close DIFF;
	}
	close $OUT;
	
	
	

}

sub computeHighlyExpressedinTwoCells_original
{
	my $differenceFile = "../../analysis/tissueSpecific/PAS_Brain_Liver_Kidney_out";
	open DIFF, "<$differenceFile" or die "Can not open file: $differenceFile";
		
	my $PasFile = "../../analysis/tissueSpecific/PAS_Brain_Liver_Kidney";
	open PAS, "<$PasFile" or die "Can not open file: $PasFile";
	
	my $sigFile = "../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Brain_Liver_Kidney_significant.csv";
	open SIG, ">$sigFile" or die "Can not open file: $sigFile";
	print SIG ("Gene,PAS,maxCell,maxDif,Liver,Kidney,Brain\n");
	
	#line1 liver, line2 kidney, line3 brain
	my $pasLine = <PAS>;
	my $pasLine2;
	my $pasLine3;
	my $line1 = <DIFF>; 	
	
	my %residual = ();
	my $noOfCells = 3;
	
	#print dif in a file
	my $diffOut = "../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Brain_Liver_Kidney_Diff_computed_original.csv";
	open my $OUT, ">$diffOut" or die "Can not open file: $diffOut";
	print $OUT ("minCell,Diff\n");	
	
	my $liver_kidney_BG  = "../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Kidney_Liver_BG_original.csv";
	my $brain_kidney_BG  = "../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Brain_Kidney_BG_original.csv";
	my $brain_liver_BG  = "../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Brain_Liver_BG_original.csv";
	
	my $liver_kidney_BG_Sig  = "../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Kidney_Liver_BG_significant.csv";
	my $brain_kidney_BG_Sig  = "../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Brain_Kidney_BG_significant.csv";
	my $brain_liver_BG_Sig  = "../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Brain_Liver_BG_significant.csv";
	
	
	open my $LK, ">$liver_kidney_BG" or die "Can not open file: $liver_kidney_BG";
	print $LK ("Diff\n");
	
	open my $BK, ">$brain_kidney_BG" or die "Can not open file: $brain_kidney_BG";
	print $BK ("Diff\n");
	
	open my $BL, ">$brain_liver_BG" or die "Can not open file: $brain_liver_BG";
	print $BL ("Diff\n");
	
	open my $LKS, ">$liver_kidney_BG_Sig" or die "Can not open file: $liver_kidney_BG_Sig";
	print $LKS ("Gene,PAS,maxCell,maxDif,Liver,Kidney,Brain\n");
	
	open my $BKS, ">$brain_kidney_BG_Sig" or die "Can not open file: $brain_kidney_BG_Sig";
	print $BKS ("Gene,PAS,maxCell,maxDif,Liver,Kidney,Brain\n");
	
	open my $BLS, ">$brain_liver_BG_Sig" or die "Can not open file: $brain_liver_BG_Sig";
	print $BLS ("Gene,PAS,maxCell,maxDif,Liver,Kidney,Brain\n");
	
	while ( $line1 = <DIFF>)
	{
		my @resid=();
		chomp ($line1);
		push @resid,$line1;
		
		for(my $i=0; $i < $noOfCells-1; $i++)
		{
			$line1=<DIFF>;
			chomp($line1);
			push @resid,$line1;
			
		}
		my @cellTypes = ("Liver", "Kidney", "Brain");
		my $type = "1";
		my $minCell;
		my $dif = getMinCell($OUT, $type,\$minCell ,@resid);
		
		$pasLine = <PAS>;
		my ($cell, $gene, $polya, $count) = split(/\,/, $pasLine);
		$pasLine=<PAS>;
		$pasLine=<PAS>;
		if($dif > 2.8)
		{
			print SIG ("$gene,$polya,$minCell,$dif,$resid[0],$resid[1],$resid[2]\n");
		}
		#BrainLiver: 2.807
		#BrainKidney:2.846
		#LiverKidney:2.7625
		
		if($minCell eq "Brain")
		{
			print $LK ("$dif\n");
			if($dif > 2.7625)
			{
				print $LKS ("$gene,$polya,$minCell,$dif,$resid[0],$resid[1],$resid[2]\n");
			}
		}
		elsif($minCell eq "Liver")
		{
			print $BK ("$dif\n");
			if($dif > 2.846)
			{
				print $BKS("$gene,$polya,$minCell,$dif,$resid[0],$resid[1],$resid[2]\n");
			}
		}
		elsif($minCell eq "Kidney")
		{
			print $BL ("$dif\n");
			if($dif > 2.807)
			{
				print $BLS ("$gene,$polya,$minCell,$dif,$resid[0],$resid[1],$resid[2]\n");
			}
		}
		else
		{
			print ("error\n");
		}
	}
	close $OUT;
	close DIFF;
	close PAS;
	close SIG;
}

sub computeHighlyExpressedinTwoCells_permuted
{
	
	my $infile_path = "../../analysis/tissueSpecific/resultFiles/PAS";
	
	my $permutation_start = 1;
	my $permutation_end = 1000;
	my %differences=(); 
	
	my $diffOut = "../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Brain_Liver_Kidney_Diff_computed_permuted.csv";
	open my $OUT, ">$diffOut" or die "Can not open file: $diffOut";
	print $OUT ("Diff\n");	
	my $noOfCells = 3;
	
	
	my $liver_kidney_BG  = "../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Kidney_Liver_BG_permuted.csv";
	my $brain_kidney_BG  = "../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Brain_Kidney_BG_permuted.csv";
	my $brain_liver_BG  = "../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Brain_Liver_BG_permuted.csv";
	
	open my $LK, ">$liver_kidney_BG" or die "Can not open file: $liver_kidney_BG";
	print $LK ("Diff\n");
	
	open my $BK, ">$brain_kidney_BG" or die "Can not open file: $brain_kidney_BG";
	print $BK ("Diff\n");
	
	open my $BL, ">$brain_liver_BG" or die "Can not open file: $brain_liver_BG";
	print $BL ("Diff\n");
	
			
	for (my $i=$permutation_start; $i<$permutation_end;$i++)
	{
		my $differenceFile = $infile_path . ".$i.residual";
		open DIFF, "<$differenceFile" ;#or die "Can not open file: $differenceFile";
		my $line1 = <DIFF>; 
		
		while ( $line1 = <DIFF>)
		{
			my @resid=();
			chomp ($line1);
			push @resid,$line1;
			
			for(my $i=0; $i < $noOfCells-1; $i++)
			{
				$line1=<DIFF>;
				chomp($line1);
				push @resid,$line1;
				
			}
			my $type = "0";
			my $minCell;
			my $dif = getMinCell($OUT,$type,\$minCell,@resid);
		
			
			
			if($minCell eq "Brain")
			{
				print $LK ("$dif\n");
			}
			elsif($minCell eq "Liver")
			{
				print $BK ("$dif\n");
			}
			elsif($minCell eq "Kidney")
			{
				print $BL ("$dif\n");
			}
			else
			{
				print ("error\n");
			}
		
		}

		close DIFF;
	}
	close $OUT;
	
	
	

}

sub getSpecificExclusive
{
	my $specific_2cells = "../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Brain_Liver_Kidney_significant.csv";
	open CELL2, "<$specific_2cells" or die "Can not open file: $specific_2cells";
		
	my $specific_1cell = "../../analysis/tissueSpecific/HighlyExpressedInOneCell/Brain_Liver_Kidney_significant.csv";
	open CELL1, "<$specific_1cell" or die "Can not open file: $specific_1cell";
		
	my $cell2 = <CELL2>;
	my $cell1 = <CELL1>;
	
	my %Kidney_Liver = ();
	my %Brain_Liver = ();
	my %Brain_Kidney = ();
	 
	while ($cell2 = <CELL2>)
	{
		my ($gene,$PAS,$minCell,$maxDif,$Liver,$Kidney,$Brain) = split(/\,/, $cell2);
		if($minCell eq "Brain")
		{
			$Kidney_Liver{$gene}{$PAS} = 1;
		}
		elsif($minCell eq "Liver")
		{
			$Brain_Kidney{$gene}{$PAS} = 1;
		}
		elsif($minCell eq "Kidney")
		{
			$Brain_Liver{$gene}{$PAS} = 1;
		}
		else
		{print ("ERROR");}	
	}
	close CELL2;
	
	my $Brain_Liver_count=0;
	my $Brain_Kidney_count=0;
	my $Brain_only_count=0;
	
	my $Liver_Brain_count=0;
	my $Liver_Kidney_count=0;
	my $Liver_only_count=0;
	
	my $Kidney_Liver_count=0;
	my $Kidney_Brain_count=0;
	my $Kidney_only_count=0;
	
	my %Brain_Liver_counted=();
	
	while ($cell1 = <CELL1>)
	{
		my ($gene,$PAS,$maxCell,$maxDif,$Liver,$Kidney,$Brain) = split(/\,/, $cell1);
		
		if($maxCell eq "Brain")
		{
			if (exists ($Brain_Liver{$gene}) && exists ($Brain_Liver{$gene}{$PAS}))
			{
				$Brain_Liver_count ++;
				print ("Brain_Liver:$gene,$PAS\n");
			}
			elsif (exists($Brain_Kidney{$gene}) && exists ($Brain_Kidney{$gene}{$PAS}))
			{
				$Brain_Kidney_count++;
				print ("Brain_Kidney:$gene,$PAS\n");
			}
			else
			{
				$Brain_only_count ++;
			}
		}
		elsif($maxCell eq "Liver")
		{
			if (exists ($Brain_Liver{$gene}) && exists ($Brain_Liver{$gene}{$PAS}))
			{
					$Liver_Brain_count ++;
					print ("Liver_Brain:$gene,$PAS\n");
			}
			elsif (exists($Kidney_Liver{$gene}) && exists ($Kidney_Liver{$gene}{$PAS}))
			{
				$Liver_Kidney_count++;
				print ("Liver_Kidney:$gene,$PAS\n");
			}
			else
			{
				$Liver_only_count ++;
			}
		}
		elsif($maxCell eq "Kidney")
		{
			if (exists ($Brain_Kidney{$gene}) && exists ($Brain_Kidney{$gene}{$PAS}))
			{
				$Kidney_Brain_count ++;
				print ("Kidney_Brain:$gene,$PAS\n");
			}
			elsif (exists($Kidney_Liver{$gene}) && exists ($Kidney_Liver{$gene}{$PAS}))
			{
				$Kidney_Liver_count++;
				print ("Kidney_Liver:$gene,$PAS\n");
			}
			else
			{
				$Kidney_only_count ++;
			}
		}
		
	}
	
	print ("Numner of Brain_Liver: $Brain_Liver_count\n");
	print ("Number of Brain_Kidney: $Brain_Kidney_count\n");
	print ("Number of Brain only specific: $Brain_only_count\n");
	
	print ("Numner of Brain_Liver: $Liver_Brain_count\n");
	print ("Number of Liver_Kidney: $Liver_Kidney_count\n");
	print ("Number of Liver only specific: $Liver_only_count\n");
	
	print ("Numner of Kidney_Brain: $Kidney_Brain_count\n");
	print ("Number of Kidney_Liver: $Kidney_Liver_count\n");
	print ("Number of Kidney only specific: $Kidney_only_count\n");
	
}

sub getSpecificExclusive2
{
	my $specific_2cells = "../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Brain_Liver_Kidney_significant.csv";
	open CELL2, "<$specific_2cells" or die "Can not open file: $specific_2cells";
		
	my $specific_1cell = "../../analysis/tissueSpecific/HighlyExpressedInOneCell/Brain_Liver_Kidney_significant.csv";
	open CELL1, "<$specific_1cell" or die "Can not open file: $specific_1cell";
		
	my $cell2 = <CELL2>;
	my $cell1 = <CELL1>;
	
	my %Kidney_Liver = ();
	my %Brain_Liver = ();
	my %Brain_Kidney = ();
	 
	my %Brain_only=();
	my %Liver_only=();
	my %Kidney_only = ();
	 
	while ($cell1 = <CELL1>)
	{
		my ($gene,$PAS,$maxCell,$maxDif,$Liver,$Kidney,$Brain) = split(/\,/, $cell1);
		
		if($maxCell eq "Brain")
		{
			$Brain_only{$gene}{$PAS} = 1;
		}
		elsif($maxCell eq "Liver")
		{
			$Liver_only{$gene}{$PAS} = 1;
		}
		elsif($maxCell eq "Kidney")
		{
			$Kidney_only{$gene}{$PAS} = 1;
		}
	} 
	my $Liver_count=0;
	my $Kidney_count=0;
	
	while ($cell2 = <CELL2>)
	{
		my ($gene,$PAS,$minCell,$maxDif,$Liver,$Kidney,$Brain) = split(/\,/, $cell2);
		if($minCell eq "Brain")
		{
			if(exists($Liver_only{$gene}{$PAS}))
			{
				#$Liver_count++;
				print("$gene,$PAS\n");
			}
			elsif(exists($Kidney_only{$gene}{$PAS}))
			{
				#$Kidney_count++;
				print ("$gene,$PAS\n");
			}
			else
			{
				$Kidney_Liver{$gene}{$PAS}=1;
			}
		}
		elsif($minCell eq "Liver")
		{
			if(exists($Kidney_only{$gene}{$PAS}))
			{
				#$Kidney_count++;
				print("$gene,$PAS\n");
			}
			elsif(exists($Brain_only{$gene}{$PAS}))
			{
				#$Liver_count++;
				print ("$gene,$PAS\n");
			}
			else
			{
				$Brain_Kidney{$gene}{$PAS}=1;
			}
		}
		if($minCell eq "Kidney")
		{
			if(exists($Liver_only{$gene}{$PAS}))
			{
				$Liver_count++;
				print("$gene,$PAS\n");
			}
			elsif(exists($Brain_only{$gene}{$PAS}))
			{
				$Kidney_count++;
				print ("$gene,$PAS\n");
			}
			else
			{
				$Brain_Liver{$gene}{$PAS}=1;
			}
		}
		
	}
	close CELL2;
	
	print ("Not in Brain: Liver $Liver_count\n");
	print ("Not in Brain: Kidney $Kidney_count\n");
	my $Brain_Liver_count=0;
	my $Brain_Kidney_count=0;
	my $Brain_only_count=0;
	
	my $Liver_Brain_count=0;
	my $Liver_Kidney_count=0;
	my $Liver_only_count=0;
	
	my $Kidney_Liver_count=0;
	my $Kidney_Brain_count=0;
	my $Kidney_only_count=0;
	
	my %Brain_Liver_counted=();
	
	
	
	print ("Numner of Brain_Liver: $Brain_Liver_count\n");
	print ("Number of Brain_Kidney: $Brain_Kidney_count\n");
	print ("Number of Brain only specific: $Brain_only_count\n");
	
	print ("Numner of Brain_Liver: $Liver_Brain_count\n");
	print ("Number of Liver_Kidney: $Liver_Kidney_count\n");
	print ("Number of Liver only specific: $Liver_only_count\n");
	
	print ("Numner of Kidney_Brain: $Kidney_Brain_count\n");
	print ("Number of Kidney_Liver: $Kidney_Liver_count\n");
	print ("Number of Kidney only specific: $Kidney_only_count\n");
	
}
sub computeZscore_original_3cellTypes
{
	my $differenceFile = "../../analysis/tissueSpecific/PAS_Brain_Liver_Kidney_out";#original_residuals.csv";
	open DIFF, "<$differenceFile" or die "Can not open file: $differenceFile";
		
	my $PasFile = "../../analysis/tissueSpecific/PAS_Brain_Liver_Kidney";
	open PAS, "<$PasFile" or die "Can not open file: $PasFile";
	
	#line1 liver, line2 kidney, line3 brain
	my $pasLine = <PAS>;
	my $pasLine2;
	my $pasLine3;
	my $line1 = <DIFF>; 	
	
	my %residual = ();
	my $noOfCells = 3;
	my @zscores=();
	
	#print zscroes in a file
	my $diffOut = "../../analysis/tissueSpecific/Brain_Liver_Kidney_original_zscore.csv";
	open OUT, ">$diffOut" or die "Can not open file: $diffOut";
	print OUT ("Gene,PAS,Diff\n");	
	
	while ( $line1 = <DIFF>)
	{
		my @resid=();
		chomp ($line1);
		push @resid,$line1;
		
		for(my $i=0; $i < $noOfCells-1; $i++)
		{
			$line1=<DIFF>;
			chomp($line1);
			push @resid,$line1;
			
		}
		
		my $z = Statistics::Zscore->new;
		my $resid_zscore = $z->standardize( \@resid );
		
		
	
		for (my $i=0; $i<$noOfCells; $i++)
		{
			$pasLine = <PAS>;
			my ($cell, $gene, $polya, $count) = split(/\,/, $pasLine);
			print OUT ("$cell,$gene,$polya,$$resid_zscore[$i]\n");
		}
	}
	close OUT;
	close DIFF;
	#close PAS;
	
}

sub computeZscore_permuted_3cellTypes
{
	
	my $infile_path = "../../analysis/tissueSpecific/resultFiles/PAS";
	
	my $permutation_start = 1;
	my $permutation_end = 1000;
	my %differences=(); 
	
	my $diffOut = "../../analysis/tissueSpecific/Brain_Liver_Kidney_permuted_zscore.csv";
	open OUT, ">$diffOut" or die "Can not open file: $diffOut";
	print OUT ("Diff\n");	
	my $noOfCells = 3;
	
	for (my $i=$permutation_start; $i<$permutation_end;$i++)
	{
		my $differenceFile = $infile_path . ".$i.residual";
		open DIFF, "<$differenceFile" ;#or die "Can not open file: $differenceFile";
		my $line1 = <DIFF>; 
		
		while ( $line1 = <DIFF>)
		{
			my @resid=();
			chomp ($line1);
			push @resid,$line1;
			
			for(my $i=0; $i < $noOfCells-1; $i++)
			{
				$line1=<DIFF>;
				chomp($line1);
				push @resid,$line1;
				
			}
			
			my $z = Statistics::Zscore->new;
			my $resid_zscore = $z->standardize( \@resid );
		
			for (my $i=0; $i<$noOfCells; $i++)
			{
				print OUT ("$$resid_zscore[$i]\n");
			}
		}

		close DIFF;
	}
	
	

}


sub getOriginalPAS
{
	my($InHandle,$OutHandle, $PAS_cell) = @_;
	
	my %PASs = ();
	my $max_dif = 0;
	my $line = <$InHandle>;
	while ($line = <$InHandle>)
	{
		chomp ($line);
		my ($gene, $PAS,$count, $start, $end) = split(/\,/, $line);
		$PASs{$gene}{$PAS}=$count;
	}
	
	
	#find closest pas

	my %closest_pases = ();
	
	 while (my ($gene, $medianHash) = each %$PAS_cell)
	 {
	 	 while (my ($pas_median,$nothing) = each %$medianHash)
	 	 {
	 	 	 while (my ($pas_original, $cellHash) = each %{$PASs{$gene}})
	 	 	 {
	 	 	 	my $dif = abs($pas_original - $pas_median);
	 	 	 	if(exists $closest_pases{$gene}{$pas_median})
	 	 	 	{
	 	 	 		if ($dif < $closest_pases{$gene}{$pas_median}{"DIF"})
	 	 	 		{
	 	 	 			$closest_pases{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	 	 			$closest_pases{$gene}{$pas_median}{"DIF"} = $dif;
	 	 	 		}
	 	 	 	}
	 	 	 	else
	 	 	 	{
	 				$closest_pases{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	  			$closest_pases{$gene}{$pas_median}{"DIF"} = $dif;
	 	 		}
	 	 	 	
	 	 	 }
	 	 }
	 }
	
	print $OutHandle ("Gene,median,original,differece\n");
	while (my ($gene, $medianHash) = each %closest_pases)
	 {
	 	 while (my ($pas_median,$nothing) = each %$medianHash)
	 	 {
	 	 	print $OutHandle ("$gene,$pas_median,$closest_pases{$gene}{$pas_median}{'PAS'},$closest_pases{$gene}{$pas_median}{'DIF'}\n");
	 	 	if($closest_pases{$gene}{$pas_median}{'DIF'} > $max_dif)
	 	 	{
	 	 		$max_dif = $closest_pases{$gene}{$pas_median}{'DIF'};
	 	 	}
	 	 }
	 }
	 
	 print ("Max dif is $max_dif\n");
	 
}

sub getOriginalPAS_closestOne
{
	my($InHandle1, $InHandle2,$OutHandle, $PAS_cell) = @_;
	
	my %PASs1 = ();
	my $max_dif = 0;
	my $line = <$InHandle1>;
	while ($line = <$InHandle1>)
	{
		chomp ($line);
		my ($gene, $PAS,$count, $start, $end) = split(/\,/, $line);
		$PASs1{$gene}{$PAS}=$count;
	}
	
	
	my %PASs2 = ();
	$line = <$InHandle2>;
	while ($line = <$InHandle2>)
	{
		chomp ($line);
		my ($gene, $PAS,$count, $start, $end) = split(/\,/, $line);
		$PASs2{$gene}{$PAS}=$count;
	}
	
	#find closest pas

	my %closest_pases1 = ();
	my %closest_pases2 = ();
	
	 while (my ($gene, $medianHash) = each %$PAS_cell)
	 {
	 	 while (my ($pas_median,$nothing) = each %$medianHash)
	 	 {
	 	 	 while (my ($pas_original, $cellHash) = each %{$PASs1{$gene}})
	 	 	 {
	 	 	 	my $dif = abs($pas_original - $pas_median);
	 	 	 	if(exists $closest_pases1{$gene}{$pas_median})
	 	 	 	{
	 	 	 		if ($dif < $closest_pases1{$gene}{$pas_median}{"DIF"})
	 	 	 		{
	 	 	 			$closest_pases1{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	 	 			$closest_pases1{$gene}{$pas_median}{"DIF"} = $dif;
	 	 	 		}
	 	 	 	}
	 	 	 	else
	 	 	 	{
	 				$closest_pases1{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	  			$closest_pases1{$gene}{$pas_median}{"DIF"} = $dif;
	 	 		}
	 	 	 	
	 	 	 }
	 	 	 while (my ($pas_original, $cellHash) = each %{$PASs2{$gene}})
	 	 	 {
	 	 	 	my $dif = abs($pas_original - $pas_median);
	 	 	 	if(exists $closest_pases2{$gene}{$pas_median})
	 	 	 	{
	 	 	 		if ($dif < $closest_pases2{$gene}{$pas_median}{"DIF"})
	 	 	 		{
	 	 	 			$closest_pases2{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	 	 			$closest_pases2{$gene}{$pas_median}{"DIF"} = $dif;
	 	 	 		}
	 	 	 	}
	 	 	 	else
	 	 	 	{
	 				$closest_pases2{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	  			$closest_pases2{$gene}{$pas_median}{"DIF"} = $dif;
	 	 		}
	 	 	 	
	 	 	 }
	 	 }
	 }
	
	print $OutHandle ("Gene,median,original,differece\n");
	while (my ($gene, $medianHash) = each %closest_pases1)
	 {
	 	 while (my ($pas_median,$nothing) = each %$medianHash)
	 	 {
	 	 	if($closest_pases1{$gene}{$pas_median}{'DIF'} <=$closest_pases2{$gene}{$pas_median}{'DIF'})
	 	 	{
	 	 		print $OutHandle ("$gene,$pas_median,$closest_pases1{$gene}{$pas_median}{'PAS'},$closest_pases1{$gene}{$pas_median}{'DIF'}\n");
	 	 	}
	 	 	else
	 	 	{
	 	 		print $OutHandle ("$gene,$pas_median,$closest_pases2{$gene}{$pas_median}{'PAS'},$closest_pases2{$gene}{$pas_median}{'DIF'}\n");
	 	 	}
	 	 	
	 	 }
	 }
	 
}

sub retrieveOriginalPASnotMedian_3cellTypes_highestInOneCell
{
	
	my $file = "../../analysis/tissueSpecific/HighlyExpressedInOneCell/Brain_Liver_Kidney_significant.csv";
	open IN, "<$file" or die "Can not open file: $file";
    my %PAS_Liver = ();
    my %PAS_Brain =();
    my %PAS_Kidney=();
    my $line = <IN>;
        
    my $max_dif=0;
	
	while ( $line = <IN>)
	{
		chomp ($line);
		my ($gene, $PAS_median,$maxCell,$dif,$liver,$kidney,$brain) = split(/\,/, $line);
		if($maxCell eq "Liver")
		{
			$PAS_Liver{$gene}{$PAS_median}=1;
		}
		elsif($maxCell eq "Brain")
		{
			$PAS_Brain{$gene}{$PAS_median}=1;
		}
		elsif($maxCell eq "Kidney")
		{
			$PAS_Kidney{$gene}{$PAS_median}=1;
		}
		else
		{print ("none\n");}
	}
	
	close IN;
	
	my $brain_out  = "../../analysis/tissueSpecific/HighlyExpressedInOneCell/Brain_Specific_3cellTypes_PAS_Original";
	my $liver_out  = "../../analysis/tissueSpecific/HighlyExpressedInOneCell/Liver_Specific_3CellTypes_PAS_Original";
	my $kidney_out  = "../../analysis/tissueSpecific/HighlyExpressedInOneCell/Kidney_Specific_3CellTypes_PAS_Original";
	
	
	my $brain_in  = "../../analysis/tissueSpecific/Brain_PA_hg19.PAS.Genes";
	my $liver_in  = "../../analysis/tissueSpecific/145_Liver_hg19.PAS.Genes";
	my $kidney_in  = "../../analysis/tissueSpecific/145_Kidney_hg19.PAS.Genes";
		
	#LIVER
	open my $InHandle, "<$liver_in" or die "Can not open file: $liver_in";
	open my $OutHandle, ">$liver_out" or die "Can not open file: $liver_out";
	getOriginalPAS($InHandle, $OutHandle,\%PAS_Liver);
	close $InHandle;
	close $OutHandle;
	
	 #BRAIN
	 open  $InHandle, "<$brain_in" or die "Can not open file: $brain_in";
	open  $OutHandle, ">$brain_out" or die "Can not open file: $brain_out";
	getOriginalPAS($InHandle, $OutHandle,\%PAS_Brain);
	close $InHandle;
	close $OutHandle;
	
	#Kidney
	open  $InHandle, "<$kidney_in" or die "Can not open file: $kidney_in";
	open  $OutHandle, ">$kidney_out" or die "Can not open file: $kidney_out";
	getOriginalPAS($InHandle, $OutHandle,\%PAS_Kidney);
	close $InHandle;
	close $OutHandle;
	
	
}
sub retrieveOriginalPASnotMedian_3cellTypes_highestInTwoCells
{#get the original PAS that is closer to the median
	
	my $file = "../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Brain_Liver_Kidney_significant.csv";
	open IN, "<$file" or die "Can not open file: $file";
    my %PAS_Kidney_Liver = ();
    my %PAS_Brain_Liver =();
    my %PAS_Brain_Kidney=();
    my $line = <IN>;
        
    my $max_dif=0;
	
	while ( $line = <IN>)
	{
		chomp ($line);
		my ($gene, $PAS_median,$minCell,$dif,$liver,$kidney,$brain) = split(/\,/, $line);
		if($minCell eq "Liver")
		{
			$PAS_Brain_Kidney{$gene}{$PAS_median}=1;
		}
		elsif($minCell eq "Brain")
		{
			$PAS_Kidney_Liver{$gene}{$PAS_median}=1;
		}
		elsif($minCell eq "Kidney")
		{
			$PAS_Brain_Liver{$gene}{$PAS_median}=1;
		}
		else
		{print ("none\n");}
	}
	
	close IN;
	
	my $brain_liver_out  = "../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Brain_Liver_Specific_3CellTypes_PAS_Original";
	my $kidney_liver_out  = "../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Kidney_Liver_Specific_3CellTypes_PAS_Original";
	my $brain_kidney_out  = "../../analysis/tissueSpecific/HighlyExpressedInTwoCells/Brain_Kidney_Specific_3CellTypes_PAS_Original";
	
	
	my $brain_in  = "../../analysis/tissueSpecific/Brain_PA_hg19.PAS.Genes";
	my $liver_in  = "../../analysis/tissueSpecific/145_Liver_hg19.PAS.Genes";
	my $kidney_in  = "../../analysis/tissueSpecific/145_Kidney_hg19.PAS.Genes";
		
	#Brain_LIVER
	open my $InHandle1, "<$liver_in" or die "Can not open file: $liver_in";
	open my $InHandle2, "<$brain_in" or die "Can not open file: $brain_in";
	open my $OutHandle, ">$brain_liver_out" or die "Can not open file: $brain_liver_out";
	getOriginalPAS_closestOne($InHandle1,$InHandle2, $OutHandle,\%PAS_Brain_Liver);
	close $InHandle1;
	close $InHandle2;
	close $OutHandle;
	
	 #BRAIN_kidney
	 open  $InHandle1, "<$brain_in" or die "Can not open file: $brain_in";
	 open  $InHandle2, "<$kidney_in" or die "Can not open file: $kidney_in";
	open  $OutHandle, ">$brain_kidney_out" or die "Can not open file: $brain_kidney_out";
	getOriginalPAS_closestOne($InHandle1,$InHandle2, $OutHandle,\%PAS_Brain_Kidney);
	close $InHandle1;
	close $InHandle2;
	close $OutHandle;
	
	#Kidney_liver
	open  $InHandle1, "<$kidney_in" or die "Can not open file: $kidney_in";
	open  $InHandle2, "<$liver_in" or die "Can not open file: $liver_in";
	open  $OutHandle, ">$kidney_liver_out" or die "Can not open file: $kidney_liver_out";
	getOriginalPAS_closestOne($InHandle1,$InHandle2, $OutHandle,\%PAS_Kidney_Liver);
	close $InHandle1;
	close $InHandle2;
	close $OutHandle;
	
	
}

sub retrieveOriginalPASnotMedian_BrainLiver
{
	
	my $file = "../../analysis/tissueSpecific/significant_final.csv";
	open IN, "<$file" or die "Can not open file: $file";
    my %PAS_Liver = ();
    my %PAS_Brain =();
        my $line = <IN>;
        
        my $max_dif=0;
	while ( $line = <IN>)
	{
		chomp ($line);
		my ($gene, $PAS_median,$dif,$liver_dif,$brain_dif) = split(/\,/, $line);
		if($liver_dif > $brain_dif)
		{
			$PAS_Liver{$gene}{$PAS_median}=1;
		}
		elsif($liver_dif < $brain_dif)
		{
			$PAS_Brain{$gene}{$PAS_median}=1;
		}
		else
		{print ("none\n");}
	}
	
	close IN;
	
	my $brain_out  = "../../analysis/tissueSpecific/Brain_Specific_PAS_Original";
	my $liver_out  = "../../analysis/tissueSpecific/Liver_Specific_PAS_Original";
	
	my $brain_in  = "../../analysis/tissueSpecific/Brain_PA_hg19.PAS.Genes";
	my $liver_in  = "../../analysis/tissueSpecific/145_Liver_hg19.PAS.Genes";
		#LIVER
		
	my %PASs = ();
	open IN, "<$liver_in" or die "Can not open file: $liver_in";
	my $line = <IN>;
	while ($line = <IN>)
	{
		chomp ($line);
		my ($gene, $PAS,$count, $start, $end) = split(/\,/, $line);
		$PASs{$gene}{$PAS}=$count;
	}
	close IN;
	
	#find closest pas

	my %closest_pases = ();
	
	 while (my ($gene, $medianHash) = each %PAS_Liver)
	 {
	 	 while (my ($pas_median,$nothing) = each %$medianHash)
	 	 {
	 	 	 while (my ($pas_original, $cellHash) = each %{$PASs{$gene}})
	 	 	 {
	 	 	 	my $dif = abs($pas_original - $pas_median);
	 	 	 	if(exists $closest_pases{$gene}{$pas_median})
	 	 	 	{
	 	 	 		if ($dif < $closest_pases{$gene}{$pas_median}{"DIF"})
	 	 	 		{
	 	 	 			$closest_pases{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	 	 			$closest_pases{$gene}{$pas_median}{"DIF"} = $dif;
	 	 	 		}
	 	 	 	}
	 	 	 	else
	 	 	 	{
	 				$closest_pases{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	  			$closest_pases{$gene}{$pas_median}{"DIF"} = $dif;
	 	 		}
	 	 	 	
	 	 	 }
	 	 }
	 }
	open OUT, ">$liver_out" or die "Can not open file: $liver_out";
	print OUT ("Gene,median,original,differece\n");
	while (my ($gene, $medianHash) = each %closest_pases)
	 {
	 	 while (my ($pas_median,$nothing) = each %$medianHash)
	 	 {
	 	 	print OUT ("$gene,$pas_median,$closest_pases{$gene}{$pas_median}{'PAS'},$closest_pases{$gene}{$pas_median}{'DIF'}\n");
	 	 	if($closest_pases{$gene}{$pas_median}{'DIF'} > $max_dif)
	 	 	{
	 	 		$max_dif = $closest_pases{$gene}{$pas_median}{'DIF'};
	 	 	}
	 	 }
	 }
	 close OUT;
	 
	 #BRAIN
	%PASs = ();
	open IN, "<$brain_in" or die "Can not open file: $brain_in";
	my $line = <IN>;
	while ($line = <IN>)
	{
		chomp ($line);
		my ($gene, $PAS,$count, $start, $end) = split(/\,/, $line);
		$PASs{$gene}{$PAS}=$count;
	}
	close IN;
	 %closest_pases = ();
	
	 while (my ($gene, $medianHash) = each %PAS_Brain)
	 {
	 	 while (my ($pas_median,$nothing) = each %$medianHash)
	 	 {
	 	 	 while (my ($pas_original, $cellHash) = each %{$PASs{$gene}})
	 	 	 {
	 	 	 	my $dif = abs($pas_original - $pas_median);
	 	 	 	if(exists $closest_pases{$gene}{$pas_median})
	 	 	 	{
	 	 	 		if ($dif < $closest_pases{$gene}{$pas_median}{"DIF"})
	 	 	 		{
	 	 	 			$closest_pases{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	 	 			$closest_pases{$gene}{$pas_median}{"DIF"} = $dif;
	 	 	 		}
	 	 	 	}
	 	 	 	else
	 	 	 	{
	 				$closest_pases{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	  			$closest_pases{$gene}{$pas_median}{"DIF"} = $dif;
	 	 		}
	 	 	 	
	 	 	 }
	 	 }
	 }
	open OUT, ">$brain_out" or die "Can not open file: $brain_out";
	print OUT ("Gene,median,original,differece\n");
	while (my ($gene, $medianHash) = each %closest_pases)
	 {
	 	 while (my ($pas_median,$nothing) = each %$medianHash)
	 	 {
	 	 	print OUT ("$gene,$pas_median,$closest_pases{$gene}{$pas_median}{'PAS'},$closest_pases{$gene}{$pas_median}{'DIF'}\n");
	 	 	if($closest_pases{$gene}{$pas_median}{'DIF'} > $max_dif)
	 	 	{
	 	 		$max_dif = $closest_pases{$gene}{$pas_median}{'DIF'};
	 	 	}
	 	 }
	 }
	 close OUT;
	 print ("max dif is $max_dif\n");
}

sub retrieveOriginalPASnotMedianFromConstitutive2
#const may not be expressed in all tissues
{
	
	my %PAS_median = ();
	my $brain_out;
	my $liver_out;
	my $kidney_out;
	  
	my $file = "../../analysis/tissueSpecific/Final/PAS_Constitutive";#SVM_3CellTypes/PAS_Brain_Liver_Kidney_constitutive";
		open IN, "<$file" or die "Can not open file: $file";
	    %PAS_median = ();
	 
	   	my $line = <IN>;   
		while ( $line = <IN>)
		{
			chomp ($line);
			my ($gene,$PAS_median) = split(/\,/, $line);
			$PAS_median{$gene}{$PAS_median}=1;
		}
		
		close IN;
		
	my $cons_out  = "../../analysis/tissueSpecific/Final/PAS_Constitutive_mapped_closest";
		
	my $brain_in  = "../../analysis/tissueSpecific/Brain_PA_hg19.PAS.Genes";
	my $liver_in  = "../../analysis/tissueSpecific/145_Liver_hg19.PAS.Genes";
	my $kidney_in  = "../../analysis/tissueSpecific/145_Kidney_hg19.PAS.Genes";
	
	my $max_dif = 0;
		
	my %PAS_liver = ();
	my %PAS_kidney = ();
	my %PAS_brain = ();
	
	#READ PASs
	open IN, "<$liver_in" or die "Can not open file: $liver_in";
	my $line = <IN>;
	while ($line = <IN>)
	{
		chomp ($line);
		my ($gene, $PAS,$count, $start, $end) = split(/\,/, $line);
		$PAS_liver{$gene}{$PAS}=$count;
	}
	close IN;
	
	open IN, "<$brain_in" or die "Can not open file: $brain_in";
	 $line = <IN>;
	while ($line = <IN>)
	{
		chomp ($line);
		my ($gene, $PAS,$count, $start, $end) = split(/\,/, $line);
		$PAS_brain{$gene}{$PAS}=$count;
	}
	close IN;
	
	open IN, "<$kidney_in" or die "Can not open file: $kidney_in";
	 $line = <IN>;
	while ($line = <IN>)
	{
		chomp ($line);
		my ($gene, $PAS,$count, $start, $end) = split(/\,/, $line);
		$PAS_kidney{$gene}{$PAS}=$count;
	}
	close IN;
	
	
	
	#find closest pas

	my %closest_pases = ();
	
	 while (my ($gene, $medianHash) = each %PAS_median)
	 {
	 	 while (my ($pas_median,$nothing) = each %$medianHash)
	 	 {
	 	 	 while (my ($pas_original, $cellHash) = each %{$PAS_brain{$gene}})
	 	 	 {
	 	 	 	my $dif = abs($pas_original - $pas_median);
	 	 	 	if(exists $closest_pases{$gene}{$pas_median})
	 	 	 	{
	 	 	 		if ($dif < $closest_pases{$gene}{$pas_median}{"DIF"})
	 	 	 		{
	 	 	 			$closest_pases{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	 	 			$closest_pases{$gene}{$pas_median}{"DIF"} = $dif;
	 	 	 		}
	 	 	 	}
	 	 	 	else
	 	 	 	{
	 				$closest_pases{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	  			$closest_pases{$gene}{$pas_median}{"DIF"} = $dif;
	 	 		}
	 	 	 	
	 	 	 }
	 	 	 ####
	 	 	 while (my ($pas_original, $cellHash) = each %{$PAS_liver{$gene}})
	 	 	 {
	 	 	 	my $dif = abs($pas_original - $pas_median);
	 	 	 	if(exists $closest_pases{$gene}{$pas_median})
	 	 	 	{
	 	 	 		if ($dif < $closest_pases{$gene}{$pas_median}{"DIF"})
	 	 	 		{
	 	 	 			$closest_pases{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	 	 			$closest_pases{$gene}{$pas_median}{"DIF"} = $dif;
	 	 	 		}
	 	 	 	}
	 	 	 	else
	 	 	 	{
	 				$closest_pases{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	  			$closest_pases{$gene}{$pas_median}{"DIF"} = $dif;
	 	 		}
	 	 	 	
	 	 	 }
	 	 	 ###
	 	 	 while (my ($pas_original, $cellHash) = each %{$PAS_kidney{$gene}})
	 	 	 {
	 	 	 	my $dif = abs($pas_original - $pas_median);
	 	 	 	if(exists $closest_pases{$gene}{$pas_median})
	 	 	 	{
	 	 	 		if ($dif < $closest_pases{$gene}{$pas_median}{"DIF"})
	 	 	 		{
	 	 	 			$closest_pases{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	 	 			$closest_pases{$gene}{$pas_median}{"DIF"} = $dif;
	 	 	 		}
	 	 	 	}
	 	 	 	else
	 	 	 	{
	 				$closest_pases{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	  			$closest_pases{$gene}{$pas_median}{"DIF"} = $dif;
	 	 		}
	 	 	 	
	 	 	 }
	 	 }
	 }
	open OUT, ">$cons_out" or die "Can not open file: $cons_out";
	print OUT ("Gene,median,original,differece\n");
	while (my ($gene, $medianHash) = each %closest_pases)
	 {
	 	 while (my ($pas_median,$nothing) = each %$medianHash)
	 	 {
	 	 	print OUT ("$gene,$pas_median,$closest_pases{$gene}{$pas_median}{'PAS'},$closest_pases{$gene}{$pas_median}{'DIF'}\n");
	 	 	if($closest_pases{$gene}{$pas_median}{'DIF'} > $max_dif)
	 	 	{
	 	 		$max_dif = $closest_pases{$gene}{$pas_median}{'DIF'};
	 	 	}
	 	 }
	 }
	 close OUT;
	 
		 print ("max dif is $max_dif\n");
	
}
sub retrieveOriginalPASnotMedianFromConstitutive
{
	
	my %PAS_median = ();
	my $brain_out;
	my $liver_out;
	my $kidney_out;
	  
	my $file = "../../analysis/tissueSpecific/Final/PAS_Constitutive";#SVM_3CellTypes/PAS_Brain_Liver_Kidney_constitutive";
		open IN, "<$file" or die "Can not open file: $file";
	    %PAS_median = ();
	 
	   	my $line = <IN>;   
		while ( $line = <IN>)
		{
			chomp ($line);
			my ($gene,$PAS_median) = split(/\,/, $line);
			$PAS_median{$gene}{$PAS_median}=1;
		}
		
		close IN;
		
		$brain_out  = "../../analysis/tissueSpecific/Final/Brain_Liver_constitutive_mapped_to_Brain_PAS";
		$liver_out  = "../../analysis/tissueSpecific/Final/Brain_Liver_constitutive_mapped_to_Liver_PAS";
		$kidney_out  = "../../analysis/tissueSpecific/Final/Brain_Liver_constitutive_mapped_to_kidney_PAS";
	
	my $brain_in  = "../../analysis/tissueSpecific/Brain_PA_hg19.PAS.Genes";
	my $liver_in  = "../../analysis/tissueSpecific/145_Liver_hg19.PAS.Genes";
	my $kidney_in  = "../../analysis/tissueSpecific/145_Kidney_hg19.PAS.Genes";
	#LIVER
	my $max_dif = 0;
	
	
		
	my %PASs = ();
	open IN, "<$liver_in" or die "Can not open file: $liver_in";
	my $line = <IN>;
	while ($line = <IN>)
	{
		chomp ($line);
		my ($gene, $PAS,$count, $start, $end) = split(/\,/, $line);
		$PASs{$gene}{$PAS}=$count;
	}
	close IN;
	
	#find closest pas

	my %closest_pases = ();
	
	 while (my ($gene, $medianHash) = each %PAS_median)
	 {
	 	 while (my ($pas_median,$nothing) = each %$medianHash)
	 	 {
	 	 	 while (my ($pas_original, $cellHash) = each %{$PASs{$gene}})
	 	 	 {
	 	 	 	my $dif = abs($pas_original - $pas_median);
	 	 	 	if(exists $closest_pases{$gene}{$pas_median})
	 	 	 	{
	 	 	 		if ($dif < $closest_pases{$gene}{$pas_median}{"DIF"})
	 	 	 		{
	 	 	 			$closest_pases{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	 	 			$closest_pases{$gene}{$pas_median}{"DIF"} = $dif;
	 	 	 		}
	 	 	 	}
	 	 	 	else
	 	 	 	{
	 				$closest_pases{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	  			$closest_pases{$gene}{$pas_median}{"DIF"} = $dif;
	 	 		}
	 	 	 	
	 	 	 }
	 	 }
	 }
	open OUT, ">$liver_out" or die "Can not open file: $liver_out";
	print OUT ("Gene,median,original,differece\n");
	while (my ($gene, $medianHash) = each %closest_pases)
	 {
	 	 while (my ($pas_median,$nothing) = each %$medianHash)
	 	 {
	 	 	print OUT ("$gene,$pas_median,$closest_pases{$gene}{$pas_median}{'PAS'},$closest_pases{$gene}{$pas_median}{'DIF'}\n");
	 	 	if($closest_pases{$gene}{$pas_median}{'DIF'} > $max_dif)
	 	 	{
	 	 		$max_dif = $closest_pases{$gene}{$pas_median}{'DIF'};
	 	 	}
	 	 }
	 }
	 close OUT;
	 
	 #BRAIN
	%PASs = ();
	open IN, "<$brain_in" or die "Can not open file: $brain_in";
	my $line = <IN>;
	while ($line = <IN>)
	{
		chomp ($line);
		my ($gene, $PAS,$count, $start, $end) = split(/\,/, $line);
		$PASs{$gene}{$PAS}=$count;
	}
	close IN;
	 %closest_pases = ();
	
	 while (my ($gene, $medianHash) = each %PAS_median)
	 {
	 	 while (my ($pas_median,$nothing) = each %$medianHash)
	 	 {
	 	 	 while (my ($pas_original, $cellHash) = each %{$PASs{$gene}})
	 	 	 {
	 	 	 	my $dif = abs($pas_original - $pas_median);
	 	 	 	if(exists $closest_pases{$gene}{$pas_median})
	 	 	 	{
	 	 	 		if ($dif < $closest_pases{$gene}{$pas_median}{"DIF"})
	 	 	 		{
	 	 	 			$closest_pases{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	 	 			$closest_pases{$gene}{$pas_median}{"DIF"} = $dif;
	 	 	 		}
	 	 	 	}
	 	 	 	else
	 	 	 	{
	 				$closest_pases{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	  			$closest_pases{$gene}{$pas_median}{"DIF"} = $dif;
	 	 		}
	 	 	 	
	 	 	 }
	 	 }
	 }
	open OUT, ">$brain_out" or die "Can not open file: $brain_out";
	print OUT ("Gene,median,original,differece\n");
	while (my ($gene, $medianHash) = each %closest_pases)
	 {
	 	 while (my ($pas_median,$nothing) = each %$medianHash)
	 	 {
	 	 	print OUT ("$gene,$pas_median,$closest_pases{$gene}{$pas_median}{'PAS'},$closest_pases{$gene}{$pas_median}{'DIF'}\n");
	 	 	if($closest_pases{$gene}{$pas_median}{'DIF'} > $max_dif)
	 	 	{
	 	 		$max_dif = $closest_pases{$gene}{$pas_median}{'DIF'};
	 	 	}
	 	 }
	 }
	 close OUT;
	 
	  #Kidney
	%PASs = ();
	open IN, "<$kidney_in" or die "Can not open file: $kidney_in";
	my $line = <IN>;
	while ($line = <IN>)
	{
		chomp ($line);
		my ($gene, $PAS,$count, $start, $end) = split(/\,/, $line);
		$PASs{$gene}{$PAS}=$count;
	}
	close IN;
	 %closest_pases = ();
	
	 while (my ($gene, $medianHash) = each %PAS_median)
	 {
	 	 while (my ($pas_median,$nothing) = each %$medianHash)
	 	 {
	 	 	 while (my ($pas_original, $cellHash) = each %{$PASs{$gene}})
	 	 	 {
	 	 	 	my $dif = abs($pas_original - $pas_median);
	 	 	 	if(exists $closest_pases{$gene}{$pas_median})
	 	 	 	{
	 	 	 		if ($dif < $closest_pases{$gene}{$pas_median}{"DIF"})
	 	 	 		{
	 	 	 			$closest_pases{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	 	 			$closest_pases{$gene}{$pas_median}{"DIF"} = $dif;
	 	 	 		}
	 	 	 	}
	 	 	 	else
	 	 	 	{
	 				$closest_pases{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	  			$closest_pases{$gene}{$pas_median}{"DIF"} = $dif;
	 	 		}
	 	 	 	
	 	 	 }
	 	 }
	 }
	open OUT, ">$kidney_out" or die "Can not open file: $kidney_out";
	print OUT ("Gene,median,original,differece\n");
	while (my ($gene, $medianHash) = each %closest_pases)
	 {
	 	 while (my ($pas_median,$nothing) = each %$medianHash)
	 	 {
	 	 	print OUT ("$gene,$pas_median,$closest_pases{$gene}{$pas_median}{'PAS'},$closest_pases{$gene}{$pas_median}{'DIF'}\n");
	 	 	if($closest_pases{$gene}{$pas_median}{'DIF'} > $max_dif)
	 	 	{
	 	 		$max_dif = $closest_pases{$gene}{$pas_median}{'DIF'};
	 	 	}
	 	 }
	 }
	 close OUT;
	 print ("max dif is $max_dif\n");
	
}

sub retrieveOriginalPASnotMedianFromCommon
{
	
	my %PAS_median = ();
	my $brain_out;
	my $liver_out;
	  
	my $file = "../../analysis/tissueSpecific/common_final.csv";
		open IN, "<$file" or die "Can not open file: $file";
	    %PAS_median = ();
	 
	   	my $line = <IN>;   
		while ( $line = <IN>)
		{
			chomp ($line);
			my ($gene,$PAS_median,$Diff,$Liver,$Brain) = split(/\,/, $line);
			$PAS_median{$gene}{$PAS_median}=1;
		}
		
		close IN;
		
	my	$out  = "../../analysis/tissueSpecific/Brain_Liver_common_PAS";
	
	my $brain_in  = "../../analysis/tissueSpecific/Brain_PA_hg19.PAS.Genes";
	my $liver_in  = "../../analysis/tissueSpecific/145_Liver_hg19.PAS.Genes";
	#LIVER
	my $max_dif = 0;
		
	my %PASs = ();
	open IN, "<$liver_in" or die "Can not open file: $liver_in";
	$line = <IN>;
	while ($line = <IN>)
	{
		chomp ($line);
		my ($gene, $PAS,$count, $start, $end) = split(/\,/, $line);
		$PASs{$gene}{$PAS}=$count;
	}
	close IN;
	
	#find closest pas

	my %closest_pases_liver = ();
	
	 while (my ($gene, $medianHash) = each %PAS_median)
	 {
	 	 while (my ($pas_median,$nothing) = each %$medianHash)
	 	 {
	 	 	 while (my ($pas_original, $cellHash) = each %{$PASs{$gene}})
	 	 	 {
	 	 	 	my $dif = abs($pas_original - $pas_median);
	 	 	 	if(exists $closest_pases_liver{$gene}{$pas_median})
	 	 	 	{
	 	 	 		if ($dif < $closest_pases_liver{$gene}{$pas_median}{"DIF"})
	 	 	 		{
	 	 	 			$closest_pases_liver{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	 	 			$closest_pases_liver{$gene}{$pas_median}{"DIF"} = $dif;
	 	 	 		}
	 	 	 	}
	 	 	 	else
	 	 	 	{
	 				$closest_pases_liver{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	  			$closest_pases_liver{$gene}{$pas_median}{"DIF"} = $dif;
	 	 		}
	 	 	 	
	 	 	 }
	 	 }
	 }
	
	 
	 #BRAIN
	%PASs = ();
	open IN, "<$brain_in" or die "Can not open file: $brain_in";
	my $line = <IN>;
	while ($line = <IN>)
	{
		chomp ($line);
		my ($gene, $PAS,$count, $start, $end) = split(/\,/, $line);
		$PASs{$gene}{$PAS}=$count;
	}
	close IN;
	my %closest_pases_brain = ();
	
	 while (my ($gene, $medianHash) = each %PAS_median)
	 {
	 	 while (my ($pas_median,$nothing) = each %$medianHash)
	 	 {
	 	 	 while (my ($pas_original, $cellHash) = each %{$PASs{$gene}})
	 	 	 {
	 	 	 	my $dif = abs($pas_original - $pas_median);
	 	 	 	if(exists $closest_pases_brain{$gene}{$pas_median})
	 	 	 	{
	 	 	 		if ($dif < $closest_pases_brain{$gene}{$pas_median}{"DIF"})
	 	 	 		{
	 	 	 			$closest_pases_brain{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	 	 			$closest_pases_brain{$gene}{$pas_median}{"DIF"} = $dif;
	 	 	 		}
	 	 	 	}
	 	 	 	else
	 	 	 	{
	 				$closest_pases_brain{$gene}{$pas_median}{"PAS"} = $pas_original;
	 	  			$closest_pases_brain{$gene}{$pas_median}{"DIF"} = $dif;
	 	 		}
	 	 	 	
	 	 	 }
	 	 }
	 }
	 
	open OUT, ">$out" or die "Can not open file: $out";

	print OUT ("Gene,median,original,differece\n");
	#print the closest from either brain or liver
	while (my ($gene, $medianHash) = each %closest_pases_brain)
	 {
	 	 while (my ($pas_median,$nothing) = each %$medianHash)
	 	 {
	 	 	if($closest_pases_brain{$gene}{$pas_median}{'DIF'} < $closest_pases_liver{$gene}{$pas_median}{'DIF'} )
	 	 	{print OUT ("$gene,$pas_median,$closest_pases_brain{$gene}{$pas_median}{'PAS'},$closest_pases_brain{$gene}{$pas_median}{'DIF'}\n");
	 	 		if($closest_pases_brain{$gene}{$pas_median}{'DIF'} > $max_dif)
		 	 	{
		 	 		$max_dif = $closest_pases_brain{$gene}{$pas_median}{'DIF'};
		 	 	}
	 	 	}
	 	 	else
	 	 	{
	 	 		{print OUT ("$gene,$pas_median,$closest_pases_liver{$gene}{$pas_median}{'PAS'},$closest_pases_liver{$gene}{$pas_median}{'DIF'}\n");}
	 	 	
	 	 		if($closest_pases_liver{$gene}{$pas_median}{'DIF'} > $max_dif)
	 	 		{
	 	 			$max_dif = $closest_pases_liver{$gene}{$pas_median}{'DIF'};
	 	 		}
	 	 	}
	 	 }
	 }
	 close OUT;
	 print ("max dif is $max_dif\n");
	
}
sub pvalue_bed
{
	my $inputFile = "../../analysis/tissueSpecific/pvalue_final_5.csv";
	my $bedfile_sig = "../../analysis/tissueSpecific/pvalue_5_sig.bed";
	
	my $geneFile = "../../analysis/tissueSpecific/GeneChrMap.txt";
	open IN, "<$geneFile" or die "Can not open file: $geneFile";
    
    my %geneChrMap=();
    my $line = <IN>;   	
	while ( $line = <IN>)
	{
		chomp ($line);
		my ($geneID, $transcriptID, $chr) = split(/\s+/, $line);
		$geneChrMap{$geneID} = $chr;
		
	}
	close IN;
	
	open OUT, ">$bedfile_sig" or die "Can not open file: $bedfile_sig";
	open IN, "<$inputFile" or die "Can not open file: $inputFile";
	
	my $line ; 
	my $sig=0;
	my $total=0;	
	while ( $line = <IN>)
	{
		chomp ($line);
		my ($gene,$pas,$value) = split(/\,/, $line); 	
		$gene =~ s/"//g;
		$pas =~ s/"//g;
		if(exists $geneChrMap{$gene})
		{
			if($value < 0.0001)
			{
				my $start = $pas -1;
				print OUT ("chr$geneChrMap{$gene}\t$start\t$pas\t$value\n");
				$sig++;
			}
			$total++;
		}
		else {print ("$gene\t");}
	}
	close IN;
	close OUT;
	my $perc = $sig/$total*100;
	
	print ("significant pvalues = $sig / $total = $perc\n");
}

sub prepareBedFileIGV_3CellTypes
{
	my $inputFile = "../../analysis/tissueSpecific/HighlyExpressedInOneCell/Brain_Liver_Kidney_Significant.csv";
	my $bedfile_sig = "../../analysis/tissueSpecific/Brain_Liver_Kidney_significant.bed";
	
	#my $inputFile = "../../analysis/tissueSpecific/pvalues_not_significant";
	#my $bedfile_nonsig = "../../analysis/tissueSpecific/non_significant_pas.bed";
	
	my $geneFile = "../../analysis/tissueSpecific/GeneChrMap.txt";
	open IN, "<$geneFile" or die "Can not open file: $geneFile";
    
    my %geneChrMap=();
    my $line = <IN>;   	
	while ( $line = <IN>)
	{
		chomp ($line);
		my ($geneID, $transcriptID, $chr) = split(/\s+/, $line);
		$geneChrMap{$geneID} = $chr;
		
	}
	close IN;
	open OUT, ">$bedfile_sig" or die "Can not open file: $bedfile_sig";
	open IN, "<$inputFile" or die "Can not open file: $inputFile";
	
	my $line ; 	
	while ( $line = <IN>)
	{
		chomp ($line);
		my ($gene,$pas,$cell,$dif,$Liver,$kidney,$Brain) = split(/\,/, $line); 	
	
		if(exists $geneChrMap{$gene})
		{
			my $start = $pas -1;
			my $num = eval sprintf('%.2f', $dif);
			print OUT ("chr$geneChrMap{$gene}\t$start\t$pas\t$num\n");
			if($gene eq 'ENSG00000062282')
			{print  ("chr$geneChrMap{$gene}\t$start\t$pas\t$num\n");
				
			}
		}
		
	}
	close IN;
	close OUT;


	while (my ($lane, $cellHash) = each %dataset)
        {
             while (my ($cell, $cell_file) = each %$cellHash)
           	{
           		my $total = 0;		
      				my $PAS_Genes_file=  $cell_file;
      				my $tissuePath = "analysis/tissueSpecific";
      				$PAS_Genes_file =~ s/analysis/$tissuePath/g;
                    $PAS_Genes_file =~ s/uniq/PAS.Genes/g;
                    open IN, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                  	
                  	my $output = $PAS_Genes_file.".bed";
                  	open OUT, ">$output" or die "Can not open file: $output";
				 	while (my $line = <IN>)
				 	{
				 		chomp ($line);
				 		my ($gene,$PAS,$PAS_count,$PAS_start,$PAS_end ) = split(/\,/, $line);
						if(exists $geneChrMap{$gene})
						{
							my $start = $PAS -1;
							
							print OUT ("chr$geneChrMap{$gene}\t$start\t$PAS\t$PAS_count\n");
							$total =$total + $PAS_count;
							
						}
						
				 	}
				 	close OUT;
				 	close IN;
				 	#print ("$PAS_Genes_file $total\n");
           	}
        }
	
}
sub prepareBedFileIGV_2Cell
{
	my $inputFile = "../../analysis/tissueSpecific/Brain_Liver_dif.csv";
	my $bedfile_sig = "../../analysis/tissueSpecific/Brain_Liver_dif.bed";
	
	#my $inputFile = "../../analysis/tissueSpecific/pvalues_not_significant";
	#my $bedfile_nonsig = "../../analysis/tissueSpecific/non_significant_pas.bed";
	
	my $geneFile = "../../analysis/tissueSpecific/GeneChrMap.txt";
	open IN, "<$geneFile" or die "Can not open file: $geneFile";
    
    my %geneChrMap=();
    my $line = <IN>;   	
	while ( $line = <IN>)
	{
		chomp ($line);
		my ($geneID, $transcriptID, $chr) = split(/\s+/, $line);
		$geneChrMap{$geneID} = $chr;
		
	}
	close IN;
	open OUT, ">$bedfile_sig" or die "Can not open file: $bedfile_sig";
	open IN, "<$inputFile" or die "Can not open file: $inputFile";
	
	my $line ; 	
	while ( $line = <IN>)
	{
		chomp ($line);
		my ($gene,$pas,$dif,$Liver,$Brain) = split(/\,/, $line); 	
	
		if(exists $geneChrMap{$gene})
		{
			my $start = $pas -1;
			my $num = eval sprintf('%.2f', $dif);
			print OUT ("chr$geneChrMap{$gene}\t$start\t$pas\t$num\n");
		}
		else {print ("$gene\t");}
	}
	close IN;
	close OUT;

	
	while (my ($lane, $cellHash) = each %dataset)
        {
             while (my ($cell, $cell_file) = each %$cellHash)
           	{
           		my $total = 0;		
      				my $PAS_Genes_file=  $cell_file;
      				my $tissuePath = "analysis/tissueSpecific";
      				$PAS_Genes_file =~ s/analysis/$tissuePath/g;
                    $PAS_Genes_file =~ s/uniq/PAS.Genes/g;
                    open IN, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
                  	
                  	my $output = $PAS_Genes_file.".bed";
                  	open OUT, ">$output" or die "Can not open file: $output";
				 	while (my $line = <IN>)
				 	{
				 		chomp ($line);
				 		my ($gene,$PAS,$PAS_count,$PAS_start,$PAS_end ) = split(/\,/, $line);
						if(exists $geneChrMap{$gene})
						{
							my $start = $PAS -1;
							print OUT ("chr$geneChrMap{$gene}\t$start\t$PAS\t$PAS_count\n");
							$total =$total + $PAS_count;
						}
						else {print ("$gene\t");}
				 	}
				 	close OUT;
				 	close IN;
				 	print ("$PAS_Genes_file $total\n");
           	}
        }
	
}
#ALL cell types
sub analyze_residual_all_cells_types
{
	my $inputFile = "../../analysis/tissueSpecific/PAS_expressed_genes_only_out";
	my $PAS_File ="../../analysis/tissueSpecific/PAS_expressed_genes_only";
	
	my $colon = "../../analysis/tissueSpecific/colon_residual.bed";
	my $kidney = "../../analysis/tissueSpecific/kidney145_residual.bed";
	my $brest = "../../analysis/tissueSpecific/brest_residual.bed";
	my $skeletal = "../../analysis/tissueSpecific/skeletal_residual.bed";
	my $heart = "../../analysis/tissueSpecific/heart_residual.bed";
	my $lung = "../../analysis/tissueSpecific/lung_residual.bed";
	my $testis = "../../analysis/tissueSpecific/testis_residual.bed";
	my $pancreas = "../../analysis/tissueSpecific/pancreas_residual.bed";
	my $liver = "../../analysis/tissueSpecific/liver_residual.bed";
	my $fetal = "../../analysis/tissueSpecific/fetal_residual.bed";
	my $brain = "../../analysis/tissueSpecific/brain_residual.bed";
	my $spleen = "../../analysis/tissueSpecific/spleen_residual.bed";
	my $prostate = "../../analysis/tissueSpecific/prostate_residual.bed";
		
	my $geneFile = "../../analysis/tissueSpecific/GeneChrMap.txt";
	open IN, "<$geneFile" or die "Can not open file: $geneFile";
    
    my %geneChrMap=();
    my $line = <IN>;   	
	while ( $line = <IN>)
	{
		chomp ($line);
		my ($geneID, $transcriptID, $chr) = split(/\s+/, $line);
		$geneChrMap{$geneID} = $chr;
		
	}
	close IN;
	

	open PAS, "<$PAS_File" or die "Can not open file: $PAS_File";
	open resid, "<$inputFile" or die "Can not open file: $inputFile";
	
	open COLON, ">$colon" or die "Can not open file: $colon";
	open KIDNEY, ">$kidney" or die "Can not open file: $kidney";
	open BREST, ">$brest" or die "Can not open file: $brest";
	open SKELETAL, ">$skeletal" or die "Can not open file: $skeletal";
	open HEART, ">$heart" or die "Can not open file: $heart";
	open LUNG, ">$lung" or die "Can not open file: $lung";
	open TESTIS, ">$testis" or die "Can not open file: $testis";
	open PANCREAS, ">$pancreas" or die "Can not open file: $pancreas";
	open LIVER, ">$liver" or die "Can not open file: $liver";
	open FETAL, ">$fetal" or die "Can not open file: $fetal";
	open BRAIN, ">$brain" or die "Can not open file: $brain";
	open SPLEEN, ">$spleen" or die "Can not open file: $spleen";
	open PROSTATE, ">$prostate" or die "Can not open file: $prostate";
	
	my @colon_resid;
	my @kidney_resid;
	my @brest_resid;
	my @skeletal_resid;
	my @heart_resid;
	my @lung_resid;
	my @testis_resid;
	my @pancreas_resid;
	my @liver_resid;
	my @fetal_resid;
	my @brain_resid;
	my @spleen_resid;
	my @prostate_resid;
	
	
	my @gene_pas;

	my $pas_line = <PAS>;
	my $resid_line = <resid>;
	 	
	while ( $pas_line = <PAS>)
	{
		chomp ($pas_line);
		$resid_line = <resid>;
		chomp($resid_line);
		my ($tissue,$gene,$pas,$count) = split(/\,/, $pas_line); 	
	
		
		if(exists $geneChrMap{$gene})
		{
			my $start = $pas -1;
			my $num = eval sprintf('%.2f', $resid_line);
			switch ($tissue)
			{
				case "Colon" {print COLON ("chr$geneChrMap{$gene}\t$start\t$pas\t$num\n");
					my $apas = $gene."_".$pas;
					push (@gene_pas, $apas);
					push (@colon_resid, $num);
				}
				case "Kidney" {print KIDNEY ("chr$geneChrMap{$gene}\t$start\t$pas\t$num\n");
					push (@kidney_resid, $num);
				}
				case "Brest" {print BREST ("chr$geneChrMap{$gene}\t$start\t$pas\t$num\n");
					push (@brest_resid, $num);
				}
				case "Skeletal_Muscle" {print SKELETAL ("chr$geneChrMap{$gene}\t$start\t$pas\t$num\n");
					push (@skeletal_resid, $num);
				}
				case "Heart" {print HEART ("chr$geneChrMap{$gene}\t$start\t$pas\t$num\n");
					push (@heart_resid, $num);
				}
				case "Lung" {print LUNG ("chr$geneChrMap{$gene}\t$start\t$pas\t$num\n");
					push (@lung_resid, $num);}
				case "Testis" {print TESTIS ("chr$geneChrMap{$gene}\t$start\t$pas\t$num\n");
					push (@testis_resid, $num);
				}
				case "Pancreas" {print PANCREAS ("chr$geneChrMap{$gene}\t$start\t$pas\t$num\n");
					push (@pancreas_resid, $num);
				}
				case "Liver" {print LIVER ("chr$geneChrMap{$gene}\t$start\t$pas\t$num\n");
					push (@liver_resid, $num);
				}
				case "Fetal" {print FETAL ("chr$geneChrMap{$gene}\t$start\t$pas\t$num\n");
					push (@fetal_resid, $num);
				}
				case "Brain" {print BRAIN ("chr$geneChrMap{$gene}\t$start\t$pas\t$num\n");
					push (@brain_resid, $num);
					if($gene eq "ENSG00000104442")
					{
						print ("chr$geneChrMap{$gene}\t$start\t$pas\n");
						print ("$pas\t$num\n");
					}
				}
				case "Spleen" {print SPLEEN ("chr$geneChrMap{$gene}\t$start\t$pas\t$num\n");
					push (@spleen_resid, $num);
				}
				case "Prostate"	 {print PROSTATE ("chr$geneChrMap{$gene}\t$start\t$pas\t$num\n");
					push (@prostate_resid, $num);}
				
			}
			
		}
		else {print ("$gene\t");}
	}
	#sort descendingly
	
	 ( $brain_genes_sorted_ref,$brain_resid_sorted_ref ) = sort_arrays_according_resid_score (\@brain_resid, \@gene_pas);
#	 ( $colon_genes_sorted_ref,$colon_resid_sorted_ref ) = sort_arrays_according_resid_score (\@colon_resid, \@gene_pas);
#	 ( $kidney_genes_sorted_ref,$kidney_resid_sorted_ref ) = sort_arrays_according_resid_score (\@kidney_resid, \@gene_pas);
#	 ( $brest_genes_sorted_ref,$brest_resid_sorted_ref ) = sort_arrays_according_resid_score (\@brest_resid, \@gene_pas);
#	( $skeletal_genes_sorted_ref,$skeletal_resid_sorted_ref ) = sort_arrays_according_resid_score (\@skeletal_resid, \@gene_pas);
#	 ( $heart_genes_sorted_ref,$heart_resid_sorted_ref ) = sort_arrays_according_resid_score (\@heart_resid, \@gene_pas);
#	 ( $lung_genes_sorted_ref,$lung_resid_sorted_ref ) = sort_arrays_according_resid_score (\@lung_resid, \@gene_pas);
#	 ( $testis_genes_sorted_ref,$testis_resid_sorted_ref ) = sort_arrays_according_resid_score (\@testis_resid, \@gene_pas);
#	 ( $pancreas_genes_sorted_ref,$pancreas_resid_sorted_ref ) = sort_arrays_according_resid_score (\@pancreas_resid, \@gene_pas);
#	 ( $liver_genes_sorted_ref,$liver_resid_sorted_ref ) = sort_arrays_according_resid_score (\@liver_resid, \@gene_pas);
#	 ( $fetal_genes_sorted_ref,$fetal_resid_sorted_ref ) = sort_arrays_according_resid_score (\@fetal_resid, \@gene_pas);
#	 ( $spleen_genes_sorted_ref,$spleen_resid_sorted_ref ) = sort_arrays_according_resid_score (\@spleen_resid, \@gene_pas);
#	 ( $prostate_genes_sorted_ref,$prostate_resid_sorted_ref ) = sort_arrays_according_resid_score (\@prostate_resid, \@gene_pas);
	
	#get UR genes and UR_other (top 200)
 	#SHOULD FIND ANOTHER WAY TO SET THIS THREASHOLD	
	getURandURother("brain", $brain_genes_sorted_ref ,$brain_resid_sorted_ref);
	#getURandURother("liver", $liver_genes_sorted_ref ,$liver_resid_sorted_ref);
	
}

sub getURandURother_permutation
{
	my ($cell,$genes_sorted_ref ,$resid_sorted_ref)	=@_;
	my @UR_resid = @{$resid_sorted_ref}[0..199];
	my @UR_genes = @{$genes_sorted_ref}[0..199];
		
	#This ensured that this set of genes was purely comprised of genes that were up-regulated in other cell types and had lower than mean expression in the brain.
	#We first made a set comprising all UR genes from all other cell types
	my ($UR_other_genes_ref, $UR_other_resid_ref) = getUROtherGeneSet($cell);
	
	my $count = @{$UR_other_genes_ref};
	print ("count in UR other is $count\n");
	
	my @todelete; 
	#We then removed the brain UR genes from this global UR other gene set
	for(my $i=0; $i<@{$UR_other_genes_ref}; )
	{
		if ( grep { $_ eq $UR_other_genes_ref->[$i] } @UR_genes)
		{
			my @indexes = indexes { $_ eq $UR_other_genes_ref->[$i] } @{$UR_other_genes_ref};
			for(my $j = 0; $j< @indexes; $j++)
			{
				splice(@{$UR_other_genes_ref}, $indexes[$j]-$j, 1);
				splice(@{$UR_other_resid_ref}, $indexes[$j]-$j, 1);
			}
			$i = $i-@indexes;
		}
		else
		{
			$i++;	
		}
	
		
		#we further removed genes from this set that had expression z-score >= 0 in the brain
		my $index = firstidx { $_ eq $UR_other_genes_ref->[$i] } @{$genes_sorted_ref};
		
		if($resid_sorted_ref->[$index] >= 0)
		{
			my @indexes = indexes { $_ eq $UR_other_genes_ref->[$i] } @{$UR_other_genes_ref};
			for(my $j = 0; $j< @indexes; $j++)
			{
				splice(@{$UR_other_genes_ref}, $indexes[$j]-$j, 1);
				splice(@{$UR_other_resid_ref}, $indexes[$j]-$j, 1);
			}
			$i=$i-@indexes;
		}
		else
		{
			$i++;
		}
	}
	
	my $count = @{$UR_other_genes_ref};
	print ("count in UR other is $count\n");
	
	writeBedFile (\@UR_genes,\@UR_resid, "../../analysis/tissueSpecific/".$cell."_UR.bed");
	writeBedFile ($UR_other_genes_ref,$UR_other_resid_ref, "../../analysis/tissueSpecific/".$cell."_UR_other.bed");
	
}    
  
sub getURandURother
{
	my ($cell,$genes_sorted_ref ,$resid_sorted_ref)	=@_;
	my @UR_resid = @{$resid_sorted_ref}[0..199];
	my @UR_genes = @{$genes_sorted_ref}[0..199];
		
	#This ensured that this set of genes was purely comprised of genes that were up-regulated in other cell types and had lower than mean expression in the brain.
	#We first made a set comprising all UR genes from all other cell types
	my ($UR_other_genes_ref, $UR_other_resid_ref) = getUROtherGeneSet($cell);
	
	my $count = @{$UR_other_genes_ref};
	print ("count in UR other is $count\n");
	
	my @todelete; 
	#We then removed the brain UR genes from this global UR other gene set
	for(my $i=0; $i<@{$UR_other_genes_ref}; )
	{
		if ( grep { $_ eq $UR_other_genes_ref->[$i] } @UR_genes)
		{
			my @indexes = indexes { $_ eq $UR_other_genes_ref->[$i] } @{$UR_other_genes_ref};
			for(my $j = 0; $j< @indexes; $j++)
			{
				splice(@{$UR_other_genes_ref}, $indexes[$j]-$j, 1);
				splice(@{$UR_other_resid_ref}, $indexes[$j]-$j, 1);
			}
			$i = $i-@indexes;
		}
		else
		{
			$i++;	
		}
	
		
		#we further removed genes from this set that had expression z-score >= 0 in the brain
		my $index = firstidx { $_ eq $UR_other_genes_ref->[$i] } @{$genes_sorted_ref};
		
		if($resid_sorted_ref->[$index] >= 0)
		{
			my @indexes = indexes { $_ eq $UR_other_genes_ref->[$i] } @{$UR_other_genes_ref};
			for(my $j = 0; $j< @indexes; $j++)
			{
				splice(@{$UR_other_genes_ref}, $indexes[$j]-$j, 1);
				splice(@{$UR_other_resid_ref}, $indexes[$j]-$j, 1);
			}
			$i=$i-@indexes;
		}
		else
		{
			$i++;
		}
	}
	
	my $count = @{$UR_other_genes_ref};
	print ("count in UR other is $count\n");
	
	writeBedFile (\@UR_genes,\@UR_resid, "../../analysis/tissueSpecific/".$cell."_UR.bed");
	writeBedFile ($UR_other_genes_ref,$UR_other_resid_ref, "../../analysis/tissueSpecific/".$cell."_UR_other.bed");
	
}  
sub writeBedFile
{
	my ($genes, $resid, $file) = @_;
	open OUT, ">$file" or die "Can not open file: $file";
		
	my $geneFile = "../../analysis/tissueSpecific/GeneChrMap.txt";
	open IN, "<$geneFile" or die "Can not open file: $geneFile";
    
    my %geneChrMap=();
    my $line = <IN>;   	
	while ( $line = <IN>)
	{
		chomp ($line);
		my ($geneID, $transcriptID, $chr) = split(/\s+/, $line);
		$geneChrMap{$geneID} = $chr;
		
	}
	close IN;
	
	for(my $i=0; $i<@{$genes}; $i++)
	{
		my ($gene, $pas) = split(/_/, $genes->[$i]);
		if(exists $geneChrMap{$gene})
		{
			my $start = $pas -1;
			print OUT ("chr$geneChrMap{$gene}\t$start\t$pas\t$resid->[$i]\t$genes->[$i]\n");
		}
	}
		

}


sub getUROtherGeneSet
{
	my ($currentCell) = @_;
	my @UR_other_resid =();
	my @UR_other_genes = ();
	
	
	while (my ($cell, $array_type_Hash) = each %cellTypes_hash)
	{
		if($cell ne $currentCell)
		{
			while (my ($array_type, $array) = each %$array_type_Hash)
			{
				#print ("arrat $$array->[0] \n");
				
				if($array_type eq "resid")
				{
					push(@UR_other_resid ,@{$$array}[0..199]);
				}
				elsif ($array_type eq "genes")
				{
					push(@UR_other_genes ,@{$$array}[0..199]);
				}
				
			}
		}
	}
		
	return (\@UR_other_genes, \@UR_other_resid);
	
}
sub sort_arrays_according_resid_score
{
	my ($resid_ref, $gene_pas_ref) = @_;
	
	my @residual = @{$resid_ref};
	my @gene_pas = @{$gene_pas_ref};
	
	my @permutation = sort { $residual[$b] <=> $residual[$a] } (0..$#residual);
	#print ("permutation $permutation[0], $permutation[1]\n");
	my @genes_sorted_array = @gene_pas[@permutation];
	my @resid_sorted = @residual[@permutation];
	
	my @indexes = indexes { $_ eq "ENSG00000104442_66515528" } @genes_sorted_array;
	print ("index of this gene is @indexes\n");
	print ("$genes_sorted_array[0],$resid_sorted[0]\n$genes_sorted_array[70],$resid_sorted[70]\n");
	return (\@genes_sorted_array,\@resid_sorted);
}
##################################################################
#old functions incorrect
sub calculateMedianUsage
{
	my $median_usage_file ="../../analysis/tissueSpecific/median";
	#my $median_usage_file ="median_usage";
	open OUT, ">$median_usage_file" or die "Can not open file: $median_usage_file";
	  my %PASs;
	 
	  
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
				 		my @values  = split(/\,/, $line);
				 		my $gene = $values[0];
				 		for(my $i=1; $i<@values; $i+=4)
				 		{
				 			my $mode = $values[$i];
				 			if(exists $PASs{$gene}{$mode})
				 			{
				 				$mode = $mode;
				 			}
				 			elsif(exists $PASs{$gene}{$mode+1})
				 			{
				 				$mode = $mode +1;
				 			}
				 			elsif(exists $PASs{$gene}{$mode+2})
				 			{
				 				$mode = $mode +2;
				 			}
				 			elsif(exists $PASs{$gene}{$mode+3})
				 			{
				 				$mode = $mode +3;
				 			}
				 			elsif(exists $PASs{$gene}{$mode+4})
				 			{
				 				$mode = $mode +4;
				 			}
				 			elsif(exists $PASs{$gene}{$mode+5})
				 			{
				 				$mode = $mode +5;
				 			}
				 			elsif(exists $PASs{$gene}{$mode-1})
				 			{
				 				$mode = $mode -1;
				 			}
				 			elsif(exists $PASs{$gene}{$mode-2})
				 			{
				 				$mode = $mode -2;
				 			}
				 			elsif(exists $PASs{$gene}{$mode-3})
				 			{
				 				$mode = $mode -3;
				 			}
				 			elsif(exists $PASs{$gene}{$mode-4})
				 			{
				 				$mode = $mode -4;
				 			}
				 			elsif(exists $PASs{$gene}{$mode-5})
				 			{
				 				$mode = $mode -5;
				 			}
				 		
				 			
				 			push(@{$PASs{$gene}{$mode}},$values[$i+3]); 
				 		}	
				 	}
				 	close IN;
           	}
        }
        
        print OUT ("Gene,mode,totalCountAllTissues,median_usage\n");
        my $median_usage;
         while (my ($gene, $modeHash) = each %PASs)
         {
        	#print OUT ("$gene");
           	foreach my $mode (sort keys %{$PASs{$gene}}) {    	
         		my @array = sort(@{$PASs{$gene}{$mode}});
         		my $total = @array;
         		my $median;
         		my $mid=0;
         		if($total %2 ==0)
         	 	{
         	 		$mid = int($total /2);
         			$median = $array[$mid];	
         	 	}
         	 	else
         	 	{
         	 		$mid = int($total /2);
					$median = ($array[$mid] + $array[$mid - 1]) / 2;
        	 	}
         	 	print OUT ("$gene,$mode,$total,$median\n");
				$median_usage->{$gene}->{$mode}= $median;
         		
          	}
         #	print OUT ("\n");
         }
         
     calculateDistanceToMedianUsage($median_usage);
         
}	

sub calculateDistanceToMedianUsage
{
	my ($median_usage) = @_;
	#my %median_usage = %$href;

	while (my ($lane, $cellHash) = each %dataset)
    {
     	while (my ($cell, $cell_file) = each %$cellHash)
	  	{
	  		my $PAS_Genes_file=  $cell_file;
      		my $tissuePath = "analysis/tissueSpecific";
      		$PAS_Genes_file =~ s/analysis/$tissuePath/g;
            $PAS_Genes_file =~ s/uniq/PAS.Genes/g;
                 
            open IN, "<$PAS_Genes_file" or die "Can not open file: $PAS_Genes_file";
              
            my $PAS_distance = $PAS_Genes_file;
            $PAS_distance =~ s/PAS.Genes/PAS.Genes.Distance/g;
            open OUT, ">$PAS_distance" or die "Can not open file: $PAS_distance";      
            print OUT ("Gene,mode,PAScount,PASvote,PAS_observed_usage,median,distanceToMedian\n");
           
       
            my $line = <IN>;
			while ($line = <IN>)
			{
				chomp ($line);
				 my @values  = split(/\,/, $line);
				 my $gene = $values[0];
				 my $med_use =0;
				 #for(my $i=1; $i<@values; $i+=4)
				 #{
				 		my $mode = $values[1];
				 		
				 		if(exists $median_usage->{$gene}->{$mode})
				 		{
				 			$med_use = $median_usage->{$gene}->{$mode};
				 			
				 		}
				 		elsif(exists $median_usage->{$gene}->{$mode+1})
				 		{
				 			$med_use = $median_usage->{$gene}->{$mode+1};
				 		}
				 		elsif(exists $median_usage->{$gene}->{$mode+2})
				 		{
				 			$med_use = $median_usage->{$gene}->{$mode+2};
				 		}
				 		elsif(exists $median_usage->{$gene}->{$mode+3})
				 		{
				 			$med_use = $median_usage->{$gene}->{$mode+3};
				 		}
				 		elsif(exists $median_usage->{$gene}->{$mode+4})
				 		{
				 			$med_use = $median_usage->{$gene}->{$mode+4};
				 		}
				 		elsif(exists $median_usage->{$gene}->{$mode+5})
				 		{
				 			$med_use = $median_usage->{$gene}->{$mode+5};
				 		}
				 		elsif(exists $median_usage->{$gene}->{$mode-1})
				 		{
				 			$med_use = $median_usage->{$gene}->{$mode-1};
				 		}
				 		elsif(exists $median_usage->{$gene}->{$mode-2})
				 		{
				 			$med_use = $median_usage->{$gene}->{$mode-2};
				 		}
				 		elsif(exists $median_usage->{$gene}->{$mode-3})
				 		{
				 			$med_use = $median_usage->{$gene}->{$mode-3};
				 		}
				 		elsif(exists $median_usage->{$gene}->{$mode-4})
				 		{
				 			$med_use = $median_usage->{$gene}->{$mode-4};
				 		}
				 		elsif(exists $median_usage->{$gene}->{$mode-5})
				 		{
				 			$med_use = $median_usage->{$gene}->{$mode-5};
				 		}
				 		
				 		#print OUT ("Gene,mode,PAScount,PASvote,PAS_observed_usage\n");	
				 	my $observed_usage = $values[4];
				 	if($med_use ==0)
				 	{
				 		print ("$PAS_Genes_file: med_use=0, $gene,$values[4]\n");
				 	}
				 	else
				 	{
				 		my $distance = ($observed_usage - $med_use) / $med_use;
				 		print OUT ("$gene,$values[1],$values[2],$values[3],$values[4],$med_use,$distance\n");
				 	}
				# }
			}
	  	}
    }
}
				     		
sub calculateMedianUsage2
{
	my $median_usage_file ="../../analysis/tissueSpecific/median_usage";
	open OUT, ">$median_usage_file" or die "Can not open file: $median_usage_file";
	  my %PASs;
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
				 		my @values  = split(/\,/, $line);
				 		my $gene = $values[0];
				 		for(my $i=1; $i<@values; $i+=4)
				 		{
				 			if (exists $PASs{$gene}{$values[$i]})
				 			{
				 			}
				 			else
				 			{
				 				$PASs{$gene}{$values[$i]};
				 			}
				 			push(@{$PASs{$gene}{$values[$i]}},$values[$i+3]); 
				 		}	
				 	}
				 	close IN;
           	}
        }
        
        print OUT ("Gene,mode,totalCountAllTissues,median_usage\n");
        my %median_usage;
         while (my ($gene, $modeHash) = each %PASs)
         {
        	print OUT ("$gene");
         	#while (my ($mode, $usage_array_ref) = each %$modeHash)
         	foreach my $mode (sort keys %{$PASs{$gene}}) {
    
         	
        # 	for my $mode  (sort { @{$PASs{$gene}{$a}} <=> @{$PASs{$gene}{$b}} }keys %{$PASs{$gene}} )
         #	{
         		#my @array = sort(@$usage_array_ref);
        		my @array = sort(@{$PASs{$gene}{$mode}});
         		my $total = @array;
         		my $median;
         		my $mid=0;
         		if($total %2 ==0)
         	 	{
         	 		$mid = int($total /2);
         			$median = $array[$mid];	
         	 	}
         	 	else
         	 	{
         	 		$mid = int($total /2);
					$median = ($array[$mid] + $array[$mid - 1]) / 2;
        	 	}
         	 	print OUT (",$mode,$total,$median");
				$median_usage{$gene}{$mode}= {$median};
         		
          	}
         	print OUT ("\n");
         }
}	         		
         		
sub subStringSequences
{
	my ($in, $up , $down) = @_;
	
	open IN, "<$in" or die "Can not open file: $in";
	open OUT_up, ">$up" or die "Can not open file: $up";
	open OUT_down, ">$down" or die "Can not open file: $down";
	my $line = <IN>;
	print OUT_up ("$line");
	print OUT_down ("$line");
	my $seq="";
	while ($line = <IN>)
	{
		chomp $line;
		if($line =~ /^>/)
		{
			my $upstream = substr($seq,0,105);
			my $downstream = substr($seq, 95,105);
			
			print OUT_up ("$upstream\n");
			print OUT_down ("$downstream\n");
			print OUT_up ("$line\n");
			print OUT_down ("$line\n");
			$seq= "";
			
		}
		else
		{
			$seq = $seq.$line;
		}
		
	}
	
	my $upstream = substr($seq,0,105);
			my $downstream = substr($seq, 95,105);
			
			print OUT_up ("$upstream\n");
			print OUT_down ("$downstream\n");
			
	close OUT_up;
	close OUT_down;
	close IN;
	 
	
}         		
       
sub subStringSequences2
{
	my ($in, $out) = @_;
	
	open IN, "<$in" or die "Can not open file: $in";
	open OUT, ">$out" or die "Can not open file: $out";
	
	my $line = <IN>;
	chomp($line);
	print OUT ("$line \n");
	my $seq="";
	while ($line = <IN>)
	{
		chomp $line;
		if($line =~ /^>/)
		{
			my $upstream = substr($seq,60,81);
		
			if($upstream !~ m/(A|C|G|T|c|t|a|g)*/)
            {print ("Error $upstream \n");}
            
            $upstream = uc($upstream);
			print OUT ("$upstream\n");
			print OUT ("$line \n");
			$seq= "";
			
		}
		else
		{
			$seq = $seq.$line;
		}
		
	}
	
	my $upstream = substr($seq,60,81);
	$upstream = uc($upstream);
    print OUT ("$upstream\n");
    #print OUT ("$upstream\n");
			
	close OUT;
	close IN;
	 
	
}        
################################
#These two functions generate 90% tss

sub filter2Uniq5Read_2
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
			}
			close IN;
			
			while (my ($chr, $strand_hash) = each %reads) 
			{	while (my ($strand, $pos_hash) = each %$strand_hash)
				{
					while(my ($pos,$hash) = each %$pos_hash)
					{
						if($reads{$chr}->{$strand}->{$pos}->{"count5"} >=2)
						{
							for(my $i=0; $i<$reads{$chr}->{$strand}->{$pos}->{"noOfReads"}; $i++)
							{
								if($strand eq "+")
								{
									print OUT ("$chr+\t$pos\t$pos\t3_utr\t100\t$strand\n");
								}
								else
								{
								print OUT ("$chr-\t$pos\t$pos\t3_utr\t100\t$strand\n");
								}
							}
						}
					}
				}
			}
			close OUT;
			my $sorted = $filteredUniq5."sorted";
			my $cmd = "sort -t$\'\t\' -k1r,1 -k 2n,2 $filteredUniq5  > $sorted";
			print ($cmd);
    		system($cmd);
    		
			
			my $cellType = $cell_file;
			my $path = '../../analysis/';
			$cellType =~ s/$path//g;
		
			my $cmd = "mkdir ../../analysis/tissueSpecific/$cellType.fseq4";
			system($cmd);
		
			my $fseq = "../../analysis/tissueSpecific/$cellType.fseq4";
		
			$cmd = "/home/dc97/biotools/fseq/bin/fseq -f 0 -l 33 -of wig -d ../../analysis/tissueSpecific -o $fseq $sorted";
	    	print($cmd."\n");
	    	system($cmd);
	   
	    	#my $cutoff = sample_wig($cellType);
	    	
		}
	}
}
sub PAS
{
	
	while (my ($lane, $cellHash) = each %dataset)
	{
		while (my ($cell, $cell_file) = each %$cellHash)
		{
			my %reads=();
			my $filtered_InternalPriming = $cell_file;
			my $tissuePath = "../analysis/tissueSpecific";
			$filtered_InternalPriming =~ s/analysis/$tissuePath/g;
			$filtered_InternalPriming =~ s/uniq/filtered.InternalPriming/g;
			open IN, "<$filtered_InternalPriming" or die "Can not open file :$filtered_InternalPriming";
			
			my $PASin = $filtered_InternalPriming;
			$PASin =~ s/filtered.InternalPriming/PASin/;
			open OUT, ">$PASin" or die "Can not open file: $PASin";
			my $line;
			while ( $line = <IN>)
			{
					chomp ($line);
					$line =~ s/^\s+//; #trim leading space
				
					my ($uniq_count,$chr, $strand,$pos_5, $pos_3) = split(/\s+/, $line);
			 		$reads{$chr}->{$strand}->{$pos_3}->{"count5"} ++;
			 		$reads{$chr}->{$strand}->{$pos_3}->{"noOfReads"} += $uniq_count;
			}
			close IN;
			
			while (my ($chr, $strand_hash) = each %reads) 
			{	while (my ($strand, $pos_hash) = each %$strand_hash)
				{
					while(my ($pos,$hash) = each %$pos_hash)
					{
						if($reads{$chr}->{$strand}->{$pos}->{"count5"} >=2)
						{
								print OUT ("$chr,$strand,$pos,$reads{$chr}->{$strand}->{$pos}->{'noOfReads'}\n");
						}
					}
				}
			}
			close OUT;
			my $sorted = $PASin."sorted";
			my $cmd = "sort -t, -k1r,1 -k 2r,2 -k3n,3 $PASin  > $sorted";
			print ($cmd);
    		system($cmd);
    		
			my $cellType = $cell_file;
			my $path = '../../analysis/';
			$cellType =~ s/$path//g;
			
	                        
	       	my $PAS_out =  $PASin;
	        $PAS_out =~ s/PASin/PASout/g;
	        my $cmd = "java GenerateTSSs -s 0 -w ../../../analysis/tissueSpecific/$cellType.fseq4/ -c $sorted -o $PAS_out";
	        print($cmd."\n");
	        system($cmd);
			
	    	
		}
	}
}
#################################################



sub RemoveSeq2
{
	
	my $file = "../../analysis/tissueSpecific/HighlyExpressedInOneCell/Brain_pos.fa";
	open INseq, "<$file" or die "Can not open file :$file";
	
	my $out = "Brain.fa";
	open OUT, ">$out" or die "Can not open file :$out";
	my $count=0;
	 
	while (my $mark = <INseq> ){
		
		my $seq =<INseq>;
		my $sub = substr($seq,100,25);
		if($sub =~ m/AAAAAA/i)
		{
			print ("$mark$seq\n");
			$count++;
		}
		else
		{
			print OUT ("$mark$seq");
		}
	}
	print ("number with priming $count\n");
}

sub compareWithMac
{
	
	my $file = "Brain_Fwd.bed";#./../analysis/tissueSpecific/HighlyExpressedInOneCell/Brain_pos.fa";
	open INseq, "<$file" or die "Can not open file :$file";
	my %brain=();
	
	while (my $line = <INseq> )
	{
		chomp ($line);
		my ($chr,$start,$end,$score) = split(/\s+/,$line);
		#$brain{$chr}{$start} = 1;
		
		$brain{$chr}{$end} = 1;
	}	
	
	close INseq;
	
	my $file = "Brain_Rev.bed";#./../analysis/tissueSpecific/HighlyExpressedInOneCell/Brain_pos.fa";
	open INseq, "<$file" or die "Can not open file :$file";
	
	
	while (my $line = <INseq> )
	{
		chomp ($line);
		my ($chr,$start,$end,$score) = split(/\s+/,$line);
		#$brain{$chr}{$start} = 1;
		$brain{$chr}{$end} = 1;
	}	
	close INseq;
	
	
	my $file = "../../analysis/tissueSpecific/HighlyExpressedInOneCell/Brain_pos.fa";
	open INseq, "<$file" or die "Can not open file :$file";
	my $count =0;
	my $total=0;
	
	while (my $mark = <INseq> )
	{
		chomp ($mark);
		my ($chr,$nu) = split(/:/,$mark);
		$chr =~ s/>//;
		my ($start,$end) = split(/-/,$nu);
		my $poly = $start + 101;
		#print ("$chr\t$poly\n");
		for (my $i=0; $i<=5; $i++)
		{
			my $p = $poly +$i;
			if(exists $brain{$chr}{$p})
			{
				$count++;
				
			}
		}
		for (my $i=5; $i>0; $i--)
		{	my $p = $poly -$i;
			if(exists $brain{$chr}{$p})
			{
				$count++;
				
			}
		}
		
		$total++;
		my $seq = <INseq>;
	}
		print ("found $count/$total\n");
		

}

sub countA
{
	
	my $file = "../../analysis/tissueSpecific/HighlyExpressedInOneCell/Liver_pos.fa";
	open INseq, "<$file" or die "Can not open file :$file";
	my $count =0;
	my $total=0;
	my $count12=0;
	my $count10=0;
	my $count8=0;
	my $count6= 0 ;
	
	my $mark = <INseq> ;
	$mark = <INseq>;
	
	#while (my $mark = <INseq> )
	{
		chomp ($mark);
		my $string = uc(substr($mark,100,25));
		print ("$string\n");
		for (my $i=0; $i<=17; $i++)
		{	my $s = substr($string,$i,8);
			print ("$s\n");
			if($s =~ m/^A+[G|C|T|A]A*[G|C|T|A]A+$/i)
			{
				$count++;
				$i=25;
			}
		#print ("$string\n");
		}
		
		
		
		 $mark = <INseq> ;
		$total++;
	}

print ("Count 8 with mismatch: $count\t$total\n");	

}
#kidney 34/52 (+/-2)
#37(+/-5)
#liver 40/69, 45
#Brain 14/71
