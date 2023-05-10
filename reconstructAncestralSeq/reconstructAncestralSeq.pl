#####Copyright Owned by Marcus Shum, 01/01/2019######
use strict;
use Data::Dumper qw(Dumper);
sub chomp2{ $_[0] =~ s/[\n\r]//gi; }

##Definition of variables##
my $line;			##Line reader for files, String
my $name;			##Pointer for fasta sequence name, String
my %sequence;		##Hash table for the fasta sequence of each taxa, HashTable
my $i;				##dummy variable for counting, Integer
my %names;			##Hash table for memorizing the name of the taxa based on the appearing sequence in the Fasta file input, HashTable
my $length;			##Record the length of the MSA, Integer
my $y;				##Dummy variable for counting in file creation, Integer
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
my $m;				##Dummy Variable for looping;
my $nucl;			##Variable to hold nucleotide character in creating the file for PASTML
###########################
print $ARGV[0]."\n";
print $ARGV[1]."\n";
print $ARGV[2]."\n";
print $ARGV[3]."\n";
print "PLESAE CONFIRMED THAT YOU HAVE INSTALLED:\n1. PastML (https://pastml.pasteur.fr/)\n2. ETE3 (http://etetoolkit.org/download/)\n3. Python2 (NOT!!3!!!)\n4. The In-house python script of the TREEMANUPULATION is under the SAME directory\n\nBeforing using this script\n";
print "If you have confirmed that you have downloaded the above programs, please click \"Enter\" to proceed...\n";
# $click=<STDIN>;
# print "Please type in the name of/full path to the tree file(in newick format): \n";
my $tree = $ARGV[0];	##Defined tree file input
chomp2 $tree;
open(TREE, "$tree");
# print "Please type in the name of/full path to the sequence file(in single-lined fasta format): \n";
my $fasta = $ARGV[1];	##Defined fasta file input
chomp2 $fasta;
open(FASTA, "$fasta");

# get desired output file path
# print "Please type in the full file path(with name ends with tree) of the ouptut tree with innode: \n";
my $outtree = $ARGV[2];
# print "Please type in the full file lspath(with name ends with fasta) of the ouptut accestral sequence fasta: \n";
my $outanc = $ARGV[3];

print "Step1: READING IN FASTA FILE... \n\n";
$i=0;
while($line = <FASTA>){
	chomp2 $line;
	$name = substr $line, 1;
	$line = <FASTA>;
	chomp2 $line;
	$sequence{$name} = $line;
	$names{$i} = $name;
	$length = length($line);
	$i=$i+1;
}
print "FINISH READING IN FASTA FILE!!! \n\n";

print "Step2: NOW GENERATING TABLE FILE FOR PASTML... \n\n";
$y=1;
for($m =0; $m<$length;$m++){
	my $filename = "position_".$y.".txt";	##Output file creation
	open($new, '>', $filename) or die "Could not open '$filename' $!";
	
	### Generation of the input file for PASTML when ambiguous nucleotide is found ###
	print $new "TAXA	NUCL\n";
	for($k = 0; $k<$i;$k++){			
		$nucl = substr $sequence{$names{$k}}, $m, 1;
		$nucl = uc $nucl;
		if($nucl eq "-"){
			print $new $names{$k}."	X\n";
		}
		elsif($nucl eq "N"){
			print $new $names{$k}."	A\n";
			print $new $names{$k}."	T\n";
			print $new $names{$k}."	C\n";
			print $new $names{$k}."	G\n";
		}
		elsif($nucl eq "R"){
			print $new $names{$k}."	A\n";
			print $new $names{$k}."	G\n";
		}
		elsif($nucl eq "Y"){
			print $new $names{$k}."	T\n";
			print $new $names{$k}."	C\n";
		}
		elsif($nucl eq "K"){
			print $new $names{$k}."	G\n";
			print $new $names{$k}."	T\n";
		}
		elsif($nucl eq "M"){
			print $new $names{$k}."	A\n";
			print $new $names{$k}."	C\n";
		}
		elsif($nucl eq "S"){
			print $new $names{$k}."	G\n";
			print $new $names{$k}."	C\n";
		}
		elsif($nucl eq "W"){
			print $new $names{$k}."	A\n";
			print $new $names{$k}."	T\n";
		}
		elsif($nucl eq "B"){
			print $new $names{$k}."	G\n";
			print $new $names{$k}."	C\n";
			print $new $names{$k}."	T\n"
		}
		elsif($nucl eq "D"){
			print $new $names{$k}."	A\n";
			print $new $names{$k}."	G\n";
			print $new $names{$k}."	T\n"
		}
		elsif($nucl eq "H"){
			print $new $names{$k}."	A\n";
			print $new $names{$k}."	C\n";
			print $new $names{$k}."	T\n"
		}
		elsif($nucl eq "V"){
			print $new $names{$k}."	A\n";
			print $new $names{$k}."	C\n";
			print $new $names{$k}."	G\n"
		}
		else{
			print $new $names{$k}."	".$nucl."\n";
		}
	}
	print "position: ".$y."/".$length."\n";
	$y=$y+1;
	close $new;
}
print "FILE GENERATION FOR PASTML COMPLETED!!!\n\n";

print "Step3: ADDITION OF INNODE NAME TO TREE FILE...\n\n";
$cmd = "source activate ete3-py2 && python /code/TIPars/reconstructAncestralSeq/TREEMANUPULATION_AddInnodeNameToTreeByArgument.py ".$tree." && conda deactivate";
$cpd = `$cmd`;
print "Addition of Name Completed!!!\n\n";

print "Step4: START RUNNING PASTML......\n\n";

for($i=1;$i<$y;$i++){
	$cmd = "source activate pastml-py3-4 && pastml --tree ./".$tree."_InnodeNameAdded --data position_".$i.".txt --columns NUCL --prediction_method ACCTRAN --work_dir position_".$i."_generated.txt && conda deactivate";
	$cpd = `$cmd`;
}

print "Finished Running PASTML!!!\n\n";

print "Step5: START TO COMBINE PASTML RESULT...\n\n";

for($h=1;$h<$y;$h++){
	$cmd = 'sed -E "s/A\/T\/C\/G/N/g" position_'.$h.'_generated.txt/combined_ancestral_states.tab > position_'.$h.'_generated.txt/combined_ancestral_states.tab_renamed';
	$cpd = `$cmd`;
	print $h."\n";
	my $ff4 = 'position_'.$h.'_generated.txt/combined_ancestral_states.tab_renamed';chomp $ff4;
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

print "Printing Out the Ancestral Sequence...\n\n";

my $ancseq_file = 'ancestral_sequence.fasta';
open(my $anc, '>', $ancseq_file) or die "Could not open '$ancseq_file' $!";

for $key (keys %seq){
	if(exists($sequence{$key})!=1){
		$seq{$key} =~ s/X/-/g;
		print $anc ">".$key."\n".$seq{$key}."\n";
	}
}

for($i=0;$i<$y;$i++){
	$cmd = "rm -rf position_".$i.".txt";
	$cpd = `$cmd`;
	$cmd = "rm -rf position_".$i."_generated.txt";
	$cpd = `$cmd`;
}

# move output files to designated location
$cmd = "mv -f ./".$tree."_InnodeNameAdded ".$outtree;
$cpd = `$cmd`;
$cmd = "mv -f ./".$ancseq_file." ".$outanc;
$cpd = `$cmd`;


print "!!! PLEASE USE THE FILE \"\*TREE.FILE\*_InnodeNameAdded\"  AND \"ancestral_sequence.fasta\" FOR TIPars Step. !!!\n\n";
print "Finished All Process!!! Thank you for using!!!\n\n";
