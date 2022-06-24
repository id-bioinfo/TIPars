#####Copyright Owned by Marcus Shum, 01/01/2019######

#####Edited by Yongtao Ye, 20/06/2022######
##1) be parallel to run PastML
##2) change to python3 for add-node script(TREEMANUPULATION_AddInnodeNameToTreeByArgument.py)
##3) apply another script (splitEachColumn using C++) to generating statefiles for PastML


use threads;
use strict;
use Data::Dumper qw(Dumper);
use Time::HiRes qw(gettimeofday);

sub chomp2{ $_[0] =~ s/[\n\r]//gi; }

##Definition of variables##
my $line;			##Line reader for files, String
my $name;			##Pointer for fasta sequence name, String
my %sequence;		##Hash table for the fasta sequence of each taxa, HashTable
my $i;				##dummy variable for counting, Integer
my %names;			##Hash table for memorizing the name of the taxa based on the appearing sequence in the Fasta file input, HashTable
my $length;			##Record the length of the MSA, Integer
my $new;			##Variable for file read-write, FilePointer
my $cmd;			##To store the bash command to run, String
my $cpd;			##To perform $cmd, String
my $ancseq_file;	##Setting the name for the Ancestral Sequence File, String
my $anc;			##Variable for file read-write, FilePointer
my $h;				##Dummy variable for For Loop, Integer
my $ff4;			##Variable for grabbing PASTML result tab file, String
my @sep;			##Array for storage, Array
my @sep2;			##Array for storage, Array
my %seq;			##Hash Table for storing the rescontructed sequence from PastML result, HashTable
my $key;			##Dummy Variable for hash table looping, String
my $click;			##Dummy Variable for click checking
my $k;				##Dummy Variable for looping;
my $nucl;			##Variable to hold nucleotide character in creating the file for PASTML
###########################
my $numArgs = $#ARGV + 1;
if ($numArgs != 4) {
    print "\nUsage: ReconstrcutAncestralSeqByPastML.pl tree.nwk taxaMSA.fas outdir numberthreads\n";
    exit;
}

print "PLESAE CONFIRMED THAT YOU HAVE INSTALLED:\n1. PastML (https://pastml.pasteur.fr/)\n2. ETE3 (http://etetoolkit.org/download/)\n3. Python3 and Perl5 (or above)\n".
"4. The in-house script for adding internal node names (TREEMANUPULATION) \n".
"5. The in-house script for generating statefiles for PastML (splitEachColumn) \n".
"Both above in-house scripts should be under the SAME directory.\n\n";

my $tree = $ARGV[0];
chomp2 $tree;
open(TREE, "$tree");
my $fasta = $ARGV[1];
chomp2 $fasta;
open(FASTA, "$fasta");
my $outdir = $ARGV[2];
chomp2 $outdir;
my $threads = int($ARGV[3]);    ##number of threads new add on 20220122 by yyt

##get alignment length";
#$line = <FASTA>;
#$line = <FASTA>;
#chomp2 $line;
#$length = length($line);
$length=0;
while($line = <FASTA>){
	chomp2 $line;
	if ($line =~ />/ && $length != 0) {
		last;
	}
	if($line =~ />/)
	{
		$length=0;
	}
	else
	{
		$length=$length + length($line);
	}
}

print "Step1: GENERATING TABLE FILE FOR PASTML... \n";

my $start_time=gettimeofday;
$cmd = "./splitEachColumn ".$fasta." ".$outdir." ".$threads;
$cpd = `$cmd`;
my $end_time=gettimeofday;
my $used_time=$end_time-$start_time; 
print "Processing used ".$used_time." seconds\n";
print "FILE GENERATION FOR PASTML COMPLETED!!!\n\n";

print "Step2: ADDITION OF INNODE NAME TO TREE FILE...\n";
$start_time=gettimeofday;
$cmd = "python3 TREEMANUPULATION_AddInnodeNameToTreeByArgument.py ".$tree;
$cpd = `$cmd`;
$end_time=gettimeofday;
$used_time=$end_time-$start_time; 
print "Processing used ".$used_time." seconds\n";
print "Addition of Name Completed!!!\n\n";

print "Step3: START RUNNING PASTML......\n";
$start_time=gettimeofday;
my $interations = int($length/$threads);
my $iter = 0;
for($iter=0; $iter<$threads; $iter++){
	my $beginIdx = $iter*$interations;
	my $endIdx = ($iter+1)*$interations;
	if($iter == $threads-1){
		$endIdx = $length;
	}
	async {
		for($i=$beginIdx;$i<$endIdx;$i++){
			$cmd = "pastml --tree ".$tree."_InnodeNameAdded --data ".$outdir."/position_".$i.".txt --columns NUCL --prediction_method ACCTRAN --work_dir ".$outdir."/position_".$i."_generated.txt";
			$cpd = `$cmd`;
			##print "PastML for ".$i;
		}
	};
}
$_->join() for threads->list;

$end_time=gettimeofday;
$used_time=$end_time-$start_time; 
print "Processing used ".$used_time." seconds\n";
print "Finished Running PASTML!!!\n\n";

print "Step4: START TO COMBINE PASTML RESULT...\n";
$start_time=gettimeofday;

for($h=0;$h<$length;$h++){
	$cmd = 'sed -E "s/A\/T\/C\/G/N/g" '.$outdir.'/position_'.$h.'_generated.txt/combined_ancestral_states.tab > '.$outdir.'/position_'.$h.'_generated.txt/combined_ancestral_states.tab_renamed';
	$cpd = `$cmd`;
	#print $h."\n";
	my $ff4 = $outdir.'/position_'.$h.'_generated.txt/combined_ancestral_states.tab_renamed';chomp $ff4;
	open(FASTA, "$ff4");
	my %record;			##Hash table for storage of NUCL cahracter, HashTable
	$line = <FASTA>;
	while($line = <FASTA>){
		chomp2 $line;
		@sep = split /	/, $line;
		if(exists($record{$sep[0]})){
			if($sep[1] eq "A"){
				if($record{$sep[0]} eq "A"){
					$record{$sep[0]} = "A";
				}
				elsif($record{$sep[0]} eq "T"){
					$record{$sep[0]} = "W";
				}
				elsif($record{$sep[0]} eq "C"){
					$record{$sep[0]} = "M";
				}
				elsif($record{$sep[0]} eq "G"){
					$record{$sep[0]} = "R";
				}
				elsif($record{$sep[0]} eq "N"){
					$record{$sep[0]} = "A";
				}
				elsif($record{$sep[0]} eq "-"){
					$record{$sep[0]} = "A";
				}
				elsif($record{$sep[0]} eq "R"){
					$record{$sep[0]} = "R";
				}
				elsif($record{$sep[0]} eq "Y"){
					$record{$sep[0]} = "H";
				}
				elsif($record{$sep[0]} eq "K"){
					$record{$sep[0]} = "D";
				}
				elsif($record{$sep[0]} eq "M"){
					$record{$sep[0]} = "M";
				}
				elsif($record{$sep[0]} eq "S"){	
					$record{$sep[0]} = "V";
				}
				elsif($record{$sep[0]} eq "W"){
					$record{$sep[0]} = "W";
				}
				elsif($record{$sep[0]} eq "B"){
					$record{$sep[0]} = "N";
				}
				elsif($record{$sep[0]} eq "D"){
					$record{$sep[0]} = "D";
				}
				elsif($record{$sep[0]} eq "H"){
					$record{$sep[0]} = "H";
				}
				elsif($record{$sep[0]} eq "V"){
					$record{$sep[0]} = "V";
				}
				elsif($record{$sep[0]} eq "X"){
					$record{$sep[0]} = "A";
				}
			}
			elsif($sep[1] eq "T"){
				if($record{$sep[0]} eq "A"){
					$record{$sep[0]} = "W";
				}
				elsif($record{$sep[0]} eq "T"){
					$record{$sep[0]} = "T";
				}
				elsif($record{$sep[0]} eq "C"){
					$record{$sep[0]} = "Y";
				}
				elsif($record{$sep[0]} eq "G"){
					$record{$sep[0]} = "K";
				}
				elsif($record{$sep[0]} eq "N"){
					$record{$sep[0]} = "T";
				}
				elsif($record{$sep[0]} eq "-"){
					$record{$sep[0]} = "T";
				}
				elsif($record{$sep[0]} eq "R"){
					$record{$sep[0]} = "D";
				}
				elsif($record{$sep[0]} eq "Y"){
					$record{$sep[0]} = "Y";
				}
				elsif($record{$sep[0]} eq "K"){
					$record{$sep[0]} = "K";
				}
				elsif($record{$sep[0]} eq "M"){
					$record{$sep[0]} = "H";
				}
				elsif($record{$sep[0]} eq "S"){	
					$record{$sep[0]} = "B";
				}
				elsif($record{$sep[0]} eq "W"){
					$record{$sep[0]} = "W";
				}
				elsif($record{$sep[0]} eq "B"){
					$record{$sep[0]} = "B";
				}
				elsif($record{$sep[0]} eq "D"){
					$record{$sep[0]} = "D";
				}
				elsif($record{$sep[0]} eq "H"){
					$record{$sep[0]} = "H";
				}
				elsif($record{$sep[0]} eq "V"){
					$record{$sep[0]} = "N";
				}
				elsif($record{$sep[0]} eq "X"){
					$record{$sep[0]} = "T";
				}
			}
			elsif($sep[1] eq "C"){
				if($record{$sep[0]} eq "A"){
					$record{$sep[0]} = "M";
				}
				elsif($record{$sep[0]} eq "T"){
					$record{$sep[0]} = "Y";
				}
				elsif($record{$sep[0]} eq "C"){
					$record{$sep[0]} = "C";
				}
				elsif($record{$sep[0]} eq "G"){
					$record{$sep[0]} = "S";
				}
				elsif($record{$sep[0]} eq "N"){
					$record{$sep[0]} = "C";
				}
				elsif($record{$sep[0]} eq "-"){
					$record{$sep[0]} = "C";
				}
				elsif($record{$sep[0]} eq "R"){
					$record{$sep[0]} = "V";
				}
				elsif($record{$sep[0]} eq "Y"){
					$record{$sep[0]} = "Y";
				}
				elsif($record{$sep[0]} eq "K"){
					$record{$sep[0]} = "B";
				}
				elsif($record{$sep[0]} eq "M"){
					$record{$sep[0]} = "M";
				}
				elsif($record{$sep[0]} eq "S"){	
					$record{$sep[0]} = "S";
				}
				elsif($record{$sep[0]} eq "W"){
					$record{$sep[0]} = "H";
				}
				elsif($record{$sep[0]} eq "B"){
					$record{$sep[0]} = "B";
				}
				elsif($record{$sep[0]} eq "D"){
					$record{$sep[0]} = "N";
				}
				elsif($record{$sep[0]} eq "H"){
					$record{$sep[0]} = "H";
				}
				elsif($record{$sep[0]} eq "V"){
					$record{$sep[0]} = "V";
				}				
				elsif($record{$sep[0]} eq "X"){
					$record{$sep[0]} = "C";
				}
			}
			elsif($sep[1] eq "G"){
				if($record{$sep[0]} eq "A"){
					$record{$sep[0]} = "R";
				}
				elsif($record{$sep[0]} eq "T"){
					$record{$sep[0]} = "K";
				}
				elsif($record{$sep[0]} eq "C"){
					$record{$sep[0]} = "S";
				}
				elsif($record{$sep[0]} eq "G"){
					$record{$sep[0]} = "G";
				}
				elsif($record{$sep[0]} eq "N"){
					$record{$sep[0]} = "G";
				}
				elsif($record{$sep[0]} eq "-"){
					$record{$sep[0]} = "G";
				}
				elsif($record{$sep[0]} eq "R"){
					$record{$sep[0]} = "R";
				}
				elsif($record{$sep[0]} eq "Y"){
					$record{$sep[0]} = "B";
				}
				elsif($record{$sep[0]} eq "K"){
					$record{$sep[0]} = "K";
				}
				elsif($record{$sep[0]} eq "M"){
					$record{$sep[0]} = "V";
				}
				elsif($record{$sep[0]} eq "S"){	
					$record{$sep[0]} = "S";
				}
				elsif($record{$sep[0]} eq "W"){
					$record{$sep[0]} = "D";
				}
				elsif($record{$sep[0]} eq "B"){
					$record{$sep[0]} = "B";
				}
				elsif($record{$sep[0]} eq "D"){
					$record{$sep[0]} = "D";
				}
				elsif($record{$sep[0]} eq "H"){
					$record{$sep[0]} = "N";
				}
				elsif($record{$sep[0]} eq "V"){
					$record{$sep[0]} = "V";
				}				
				elsif($record{$sep[0]} eq "X"){
					$record{$sep[0]} = "G";
				}
			}
			elsif($sep[1] eq "N"){
				if($record{$sep[0]} eq "A"){
					$record{$sep[0]} = "A";
				}
				elsif($record{$sep[0]} eq "T"){
					$record{$sep[0]} = "T";
				}
				elsif($record{$sep[0]} eq "C"){
					$record{$sep[0]} = "C";
				}
				elsif($record{$sep[0]} eq "G"){
					$record{$sep[0]} = "G";
				}
				elsif($record{$sep[0]} eq "N"){
					$record{$sep[0]} = "N";
				}
				elsif($record{$sep[0]} eq "-"){
					$record{$sep[0]} = "N";
				}
				elsif($record{$sep[0]} eq "R"){
					$record{$sep[0]} = "R";
				}
				elsif($record{$sep[0]} eq "Y"){
					$record{$sep[0]} = "Y";
				}
				elsif($record{$sep[0]} eq "K"){
					$record{$sep[0]} = "K";
				}
				elsif($record{$sep[0]} eq "M"){
					$record{$sep[0]} = "M";
				}
				elsif($record{$sep[0]} eq "S"){	
					$record{$sep[0]} = "S";
				}
				elsif($record{$sep[0]} eq "W"){
					$record{$sep[0]} = "W";
				}
				elsif($record{$sep[0]} eq "B"){
					$record{$sep[0]} = "B";
				}
				elsif($record{$sep[0]} eq "D"){
					$record{$sep[0]} = "D";
				}
				elsif($record{$sep[0]} eq "H"){
					$record{$sep[0]} = "H";
				}
				elsif($record{$sep[0]} eq "V"){
					$record{$sep[0]} = "V";
				}
				elsif($record{$sep[0]} eq "X"){
					$record{$sep[0]} = "N";
				}
			}
			elsif($sep[1] eq "X"){
				if($record{$sep[0]} eq "A"){
					$record{$sep[0]} = "A";
				}
				elsif($record{$sep[0]} eq "T"){
					$record{$sep[0]} = "T";
				}
				elsif($record{$sep[0]} eq "C"){
					$record{$sep[0]} = "C";
				}
				elsif($record{$sep[0]} eq "G"){
					$record{$sep[0]} = "G";
				}
				elsif($record{$sep[0]} eq "N"){
					$record{$sep[0]} = "N";
				}
				elsif($record{$sep[0]} eq "-"){
					$record{$sep[0]} = "-";
				}
				elsif($record{$sep[0]} eq "R"){
					$record{$sep[0]} = "R";
				}
				elsif($record{$sep[0]} eq "Y"){
					$record{$sep[0]} = "Y";
				}
				elsif($record{$sep[0]} eq "K"){
					$record{$sep[0]} = "K";
				}
				elsif($record{$sep[0]} eq "M"){
					$record{$sep[0]} = "M";
				}
				elsif($record{$sep[0]} eq "S"){	
					$record{$sep[0]} = "S";
				}
				elsif($record{$sep[0]} eq "W"){
					$record{$sep[0]} = "W";
				}
				elsif($record{$sep[0]} eq "B"){
					$record{$sep[0]} = "B";
				}
				elsif($record{$sep[0]} eq "D"){
					$record{$sep[0]} = "D";
				}
				elsif($record{$sep[0]} eq "H"){
					$record{$sep[0]} = "H";
				}
				elsif($record{$sep[0]} eq "V"){
					$record{$sep[0]} = "V";
				}
				elsif($record{$sep[0]} eq "X"){
					$record{$sep[0]} = "X";
				}
			}
		}
		else{
			$record{$sep[0]} = $sep[1];
		}
	}
	
	for $key(keys %record){
		if(exists($seq{$key})){
			$seq{$key} = $seq{$key}.$record{$key};
		}
		else{
			$seq{$key} = $record{$key};
		}
	}
	close FASTA;
}

$end_time=gettimeofday;
$used_time=$end_time-$start_time; 
print "Processing used ".$used_time." seconds\n";
print "Printing Out the Ancestral Sequence...\n\n";

my $ancseq_file = $outdir.'/ancestral_sequence.fasta';
open(my $anc, '>', $ancseq_file) or die "Could not open '$ancseq_file' $!";

for $key (keys %seq){
	if(exists($sequence{$key})!=1){
		$seq{$key} =~ s/X/-/g;
		print $anc ">".$key."\n".$seq{$key}."\n";
	}
}

for($i=0;$i<$length;$i++){
	$cmd = "rm -rf ".$outdir."/position_".$i.".txt";
	$cpd = `$cmd`;
	$cmd = "rm -rf ".$outdir."/position_".$i."_generated.txt";
	$cpd = `$cmd`;
}
print "Finished All Process!!! Thank you for using!!!\n\n";




