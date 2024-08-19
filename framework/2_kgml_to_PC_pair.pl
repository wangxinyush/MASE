#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;

#my $indir = "./hsa_kgml";
#open OUT,">","hsa_compound_gene_react.txt" or die $!;
#my $indir = "./mmu_kgml";
#open OUT,">","mmu_compound_gene_react.txt" or die $!;
my $indir = "./mcc_kgml";
open OUT,">","mcc_compound_gene_react.txt" or die $!;

my @kgml_files = glob "$indir/*.xml";
foreach my $kgml_file (sort @kgml_files){
	#my $kgml_file = "./test_hsa00270/hsa00270.xml";
	my ($ko_id,$path,$suffix) = fileparse($kgml_file, ".xml");

	my %reaction_gene_hs = (); #reaction id => gene1 gene2 ... (space splited)
	
	#获取该文件中所有的reaction以及其对应的gene
	open KGML,"<",$kgml_file or die $!;
	while(<KGML>){
		if(/<entry /){ #entry
			if(/type="gene"/){ #只取gene - reaction
				if(/id="(\d*)" name="(.*?)" type="gene" reaction="(.*?)"/){
					my ($gene_id, $react_id) = ($2, $3);
					
					if(exists $reaction_gene_hs{$react_id}){ 
						$reaction_gene_hs{$react_id} .= " $gene_id";
					}
					else{
						$reaction_gene_hs{$react_id} = $gene_id; #reaction id => gene1 gene2 ...
					}
				}
			}
		}
	}
	close KGML;

	#获取reaction相关的化合物
	my %reaction_K_hs = (); #rnXXX => 
	my %reaction_C_hs = ();
	open KGML,"<",$kgml_file or die $!;
	while(<KGML>){
		if(/<reaction/){ #reaction
			/name="(.*?)" type="(.*?)"/; #reaction id string
			my ($react_id, $react_type) = ($1, $2); #reaction id, type
			
			if(exists $reaction_gene_hs{$react_id}){ #有反应信息
				while(<KGML>){
					last if(/<\/reaction>/);

					#判断是product还是substrate
					if(/<substrate/){
						/name="(.*?)"/;
						my $substrate = $1;
						#R10992\tC00109\thsa:110 hsa:123\tsubstrate
						print OUT "$ko_id\t$react_id\t$substrate\t$reaction_gene_hs{$react_id}\tsubstrate\t$react_type\n";
					}
					elsif(/<product/){
						/name="(.*?)"/;
						my $product = $1;
						print OUT "$ko_id\t$react_id\t$product\t$reaction_gene_hs{$react_id}\tproduct\t$react_type\n";
					}
				}
			}
			else{
				warn "$react_id has no information in $kgml_file.\n";
			}
			
		}
	}
	
}



