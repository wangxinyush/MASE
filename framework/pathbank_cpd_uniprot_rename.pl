#！/usr/bin/perl
use warnings;
use strict;

my $dir = "E:/metabolism_prediction/Results/1_network";

#1. 获取代谢名到chEBI的映射信息, 1314个cpd，1089个到chebi
open CPD,"<","../pathbank_primary_pathway/pathbank_primary_metabolites.csv" or die $!;
my %cpd_chebi_hs = ();
while(<CPD>){
	s/[\r\n]//g;
	my @arr = split /"/;
	my @new_arr = ();
	foreach my $i(0..$#arr){
		if($i % 2 == 1){ #奇数 说明是被扩起来的
			push @new_arr, $arr[$i] if($arr[$i]);
		}
		else{
			my @temp = split /,/, $arr[$i], -1;
					
			shift @temp if $i != 0;
			pop @temp if $i != $#arr;
			
			foreach my $t(@temp){
				push @new_arr, $t; 
			}
		}
	}
	
	next if($new_arr[3] ne "Homo sapiens");
	#next if($new_arr[3] ne "Mus musculus");
	if($new_arr[8]){
		$cpd_chebi_hs{$new_arr[5]} = $new_arr[8];
	}
}
close CPD;

#1.1 获取chEBI对应的化合物名称信息
open CHEBI,"<","$dir/annotation/chEBI/compounds.tsv" or die $!;
<CHEBI>;
my %chebi_id_name_hs = ();
while(<CHEBI>){
	s/[\r\n]//g;
	my @arr = split /\t/;
	my ($id, $name) = ($arr[2], $arr[5]);
	$id =~ s/CHEBI://g;
	$chebi_id_name_hs{$id} = $name; #28879 => XXXX
}
close CHEBI;

#2. 获取hsa/mmu到UNIPROT的映射信息 
open PROTEIN,"<","../pathbank_primary_pathway/pathbank_primary_proteins.csv" or die $!;
my %protein_uniprot_hs = ();
while(<PROTEIN>){
	s/[\r\n]//g;
	
	my @arr = split /"/;
	my @new_arr = ();
	foreach my $i(0..$#arr){
		if($i % 2 == 1){ #奇数 说明是被扩起来的
			push @new_arr, $arr[$i] if($arr[$i]);
		}
		else{
			my @temp = split /,/, $arr[$i], -1;
					
			shift @temp if $i != 0;
			pop @temp if $i != $#arr;
			
			foreach my $t(@temp){
				push @new_arr, $t; 
			}
		}
	}
	
	next if($new_arr[3] ne "Homo sapiens");
	#next if($new_arr[3] ne "Mus musculus");
	if($new_arr[4]){
		$protein_uniprot_hs{$new_arr[5]} = $new_arr[4];
	}
}
close PROTEIN;

#2.1 获取UNIPROT对应的蛋白质名称信息
open UNI,"<","$dir/annotation/uniprot/uniprotkb_human_20230920.tsv" or die $!;
#open UNI,"<","$dir/annotation/uniprot/uniprotkb_mouse_20230920.tsv" or die $!;
<UNI>;
my %uniprot_id_name_hs = ();
while(<UNI>){
	s/[\r\n]//g;
	my @arr = split /\t/;
	my ($id, $entry_name, $protein_name, $gene_name) = @arr[0..3];
	$uniprot_id_name_hs{$id} = join "\t",($entry_name, $protein_name, $gene_name);
}
close UNI;

#3. Rename
open REACT,"<","Pathbank_HSA_clean_compound_enzyme.txt" or die $!;
open OUT,">","./rename/hsa_chebi_uniprot_reaction.txt" or die $!;

#open REACT,"<","Pathbank_MMU_clean_compound_enzyme.txt" or die $!;
#open OUT,">","./rename/mmu_chebi_uniprot_reaction.txt" or die $!;

print OUT "ChEBI_ID\tUniProt_ID\tType\tChEBI_Name\tUniProt_Name\tUniProt_Protein_Name\tGene_Name\t";

my @header_arr = qw(Pathway_id Reaction_id Reaction_name Compound Enzyme Type Compound_id Enzyme_id Compound_sboterm Enzyme_sboterm);
print OUT join "\t", @header_arr;

while(<REACT>){
	s/[\r\n]//g;
	my ($pathway_id, $reaction_id, $reaction_name, $cpd, $enzyme, $type, $cpd_id, $enzyme_id, $cpd_term, $enzyme_term) = split /\t/;
	my ($chebi_id, $uniprot_id, $chebi_name, $uniprot_str) = ("NA", "NA", "NA", "NA\tNA\tNA"); #uniprot_str:UniProt_Name\tUniProt_Protein_Name\tGene_Name
	
	if(exists $cpd_chebi_hs{$cpd}){
		$chebi_id = $cpd_chebi_hs{$cpd};
		
		if(exists $chebi_id_name_hs{$chebi_id}){ #has name
			$chebi_name = $chebi_id_name_hs{$chebi_id};
		}
	}
	
	if(exists $protein_uniprot_hs{$enzyme}){
		$uniprot_id = $protein_uniprot_hs{$enzyme};
		
		if(exists $uniprot_id_name_hs{$uniprot_id}){ #has name
			$uniprot_str = $uniprot_id_name_hs{$uniprot_id};
		}
	}
	
	print OUT join "\t", ($chebi_id, $uniprot_id, $type, $chebi_name, $uniprot_str);
	print OUT "\t$_\n";
}
close REACT;
close OUT;
