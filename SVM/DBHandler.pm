#! usr/bin/perl -w
use strict;
use DBI;

#This script is mainly to access sql database on the server
#created on Feb 23 2011
#by Dina

$| = 1;    #If set to nonzero, forces a flush after every write or print

sub DBconnect
{
	my $dsn = 'dbi:mysql:database=d_h;host=mysql-sandbox.igsp.duke.edu';
	my $username = 'd_hafez';
	my $password = 'Dodo3';
	
	#Establish connection to database
	my $dbh = DBI->connect($dsn, $username, $password, { RaiseError=>1, AutoCommit=>0, ShowErrorStatement=>1 }) or die "Can't connect to database: $DBI::errstr\n";

return $dbh;

}

sub closeConnection
{
	my ($dbh) = @_;
	
	#disconnect
	$dbh->disconnect;
}

sub DBselect 
{
	my ($dbh,$table, $ORF) = @_;
	
	my $sthSelect = $dbh->prepare('select * from '.$table.' where Systematic_Name = ?');
	
	#run the query
	$sthSelect->execute(($ORF));
	
	#iterate over the set of selected rows
	#while (my @record = $sthSelect->fetchrow_array())
	#{
	#  print "systematic Name= " . $record[0] . "\t Chromosome = " . $record[1] . "\t start = " . $record[2] . "\t end = " . $record[3] . "\n";
	#}
	
	return $sthSelect->fetchrow_array()
	
}

sub DBselectIN 
{
	my ($dbh,$table, $ORF) = @_;
	
	my $sthSelect = $dbh->prepare('select DISTINCT wiki from '.$table.' where Ensembl IN ('.$ORF.')');
	
	#run the query
	$sthSelect->execute();
	
	#iterate over the set of selected rows
	#while (my @record = $sthSelect->fetchrow_array())
	#{
	#  print "systematic Name= " . $record[0] . "\t Chromosome = " . $record[1] . "\t start = " . $record[2] . "\t end = " . $record[3] . "\n";
	#}
	
	return $sthSelect->fetchrow_array()
	
}

sub DBselectGeneID
{
	my ($dbh, $ORF) = @_;
	
	my $sthSelect = $dbh->prepare('select DISTINCT GeneID from GeneMap where TranscriptID IN ('.$ORF.')');
	
	#run the query
	$sthSelect->execute();
	
	#iterate over the set of selected rows
	#while (my @record = $sthSelect->fetchrow_array())
	#{
	#  print "systematic Name= " . $record[0] . "\t Chromosome = " . $record[1] . "\t start = " . $record[2] . "\t end = " . $record[3] . "\n";
	#}
	
	return $sthSelect->fetchrow_array()
	
}

sub DBselectChrStrand
{
	my ($dbh, $ORF) = @_;
	
	my $sthSelect = $dbh->prepare('select DISTINCT Chromosome,strand from EnsembleTranscript where GeneID IN ('.$ORF.')');
	
	#run the query
	$sthSelect->execute();
	
	#iterate over the set of selected rows
	#while (my @record = $sthSelect->fetchrow_array())
	#{
	#  print "systematic Name= " . $record[0] . "\t Chromosome = " . $record[1] . "\t start = " . $record[2] . "\t end = " . $record[3] . "\n";
	#}
	
	return $sthSelect->fetchrow_array();
	
}


sub DBselectIN_2 
{
	my ($dbh,$table, $ORF) = @_;
	
	my $sthSelect = $dbh->prepare('select DISTINCT Refseq from '.$table.' where Ensembl IN ('.$ORF.')');
	
	#run the query
	$sthSelect->execute();
	
	return $sthSelect->fetchrow_array()
	
}

sub DBselectGeneStartEnd 
{
	my ($dbh,$table, $ORF) = @_;
	
	my $sthSelect = $dbh->prepare('select DISTINCT genestart, geneend from '.$table.' where wiki = '.$ORF.' OR refseq ='.$ORF);
	
	#run the query
	$sthSelect->execute();
	
	
	#iterate over the set of selected rows
	#while (my @record = $sthSelect->fetchrow_array())
	#{
	#  print "systematic Name= " . $record[0] . "\t Chromosome = " . $record[1] . "\t start = " . $record[2] . "\t end = " . $record[3] . "\n";
	#}
	
	return $sthSelect->fetchrow_array()
	
}

sub DBselectGeneStopCodon 
{
	my ($dbh, $ORF, $strand) = @_;
	my $sthSelect ;
	if($strand eq '+')
	{
		 $sthSelect = $dbh->prepare('select min(end) from stopcodons where wiki=' .$ORF. ' OR refseq ='.$ORF);
	}
	else
	{
		 $sthSelect = $dbh->prepare('select max(start) from stopcodons where wiki=' .$ORF. ' OR refseq ='.$ORF);
	}

	
	#run the query
	$sthSelect->execute();
	
	
	return $sthSelect->fetchrow_array()
	
}

sub DBselectGeneStopCodonFromGene 
{
	my ($dbh, $ORF, $strand) = @_;
	my $sthSelect ;
	if($strand eq '+')
	{
		 $sthSelect = $dbh->prepare('select min(stopcodons.end) from  stopcodons left join geneinfo on stopcodons.refseq=  geneinfo.refseq where geneinfo.wiki=' .$ORF);
	}
	else
	{
		 $sthSelect = $dbh->prepare('select min(stopcodons.start) from  stopcodons left join geneinfo on stopcodons.refseq=  geneinfo.refseq where geneinfo.wiki=' .$ORF);
	}

	
	#run the query
	$sthSelect->execute();
	
	
	return $sthSelect->fetchrow_array()
	
}
 