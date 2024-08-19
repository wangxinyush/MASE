#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;

#INPUT & OUTPUT
#my $indir = "./01.homo_sapiens.3.1.sbml";
#open OUT,">","Reactome_HSA_compound_enzyme.txt" or die $!;

my $indir = "../00.rawData/mmu_sbml";
open OUT,">","./Reactome_MMU_compound_enzyme.txt" or die $!;

my @header_arr = qw(Pathway_id Reaction_id Reaction_name Compound Enzyme Type Compound_id Enzyme_id Compound_sboterm Enzyme_sboterm);
print OUT join "\t", @header_arr;
print OUT "\n";

#my @sbml_files = ("$indir/R-HSA-15869.sbml");
my @sbml_files = glob "$indir/*.sbml";
foreach my $sbml_file (sort @sbml_files){
	#my $sbml_file = "./test/R-HSA-5334118.sbml";
	my ($pathway_id,$path,$suffix) = fileparse($sbml_file, ".sbml");

	my %species_hs = (); #species id => species name
	my %species_sbo_hs = (); #species id => species sboTerm
	open SBML,"<",$sbml_file or die $!;
	while(<SBML>){
		if(/<species /){ #species
			my ($species_id, $species_name, $species_sboTerm) = (/id="(.*?)"/, /name="(.*?)"/, /sboTerm="(.*?)"/); #species id, name, sboTerm

			if(exists $species_hs{$species_id}){ 
				warn "$species_id, $species_name is not unique in $pathway_id!\n";
			}
			else{
				$species_hs{$species_id} = $species_name; #species id => species name
				$species_sbo_hs{$species_id} = $species_sboTerm; #species id => species sboTerm
			}
		}
	}
	close SBML;

	open SBML,"<",$sbml_file or die $!;
	while(<SBML>){
		if(/<reaction /){ #reaction start
			my ($reaction_id, $reaction_name) = (/id="(.*?)"/, /name="(.*?)"/);  #reaction id, name
			my @substrates_arr = ();
			my @products_arr = ();
			my @enzymes_arr = ();
			while(<SBML>){ #get substrates, products, enzymes for this reaction.
				last if(/<\/reaction/);
				if(/<listOfReactants/){ #substrates start
					while(<SBML>){
						last if (/<\/listOfReactants/);
						if(/<speciesReference/){ #substrates
							my ($species_id) = (/species="(.*?)"/);
							if(exists $species_hs{$species_id}){
								push @substrates_arr, $species_id;
							}
							else{
								warn "$species_id not existed in species!\n";
							}
						}
					}
				}
				if(/<listOfProducts/){ #products start
					while(<SBML>){
						last if (/<\/listOfProducts/);
						if(/<speciesReference/){ #substrates
							my ($species_id) = (/species="(.*?)"/);
							if(exists $species_hs{$species_id}){
								push @products_arr, $species_id;
							}
							else{
								warn "$species_id not existed in species!\n";
							}
						}
					}
				}
				if(/<listOfModifiers/){ #modifiers start
					while(<SBML>){
						last if (/<\/listOfModifiers/);
						if(/<modifierSpeciesReference/){ #substrates
							my ($species_id) = (/species="(.*?)"/);
							if(exists $species_hs{$species_id}){
								push @enzymes_arr, $species_id;
							}
							else{
								warn "$species_id not existed in species!\n";
							}
						}
					}
				}
			}
			
			foreach my $s (@substrates_arr){
				foreach my $e (@enzymes_arr){
					#my @header_arr = qw(Pathway_id Reaction_id Reaction_name Compound Enzyme Type Compound_id Enzyme_id Compound_sboterm Enzyme_sboterm);
					my ($s_sbo, $e_sbo) = ($species_sbo_hs{$s}, $species_sbo_hs{$e});
					$s_sbo = "NA" unless $s_sbo;
					$e_sbo = "NA" unless $e_sbo;
					print OUT "$pathway_id\t$reaction_id\t$reaction_name\t$species_hs{$s}\t$species_hs{$e}\tsubstrates\t$s\t$e\t$s_sbo\t$e_sbo\n";
				}
			}
			foreach my $s (@products_arr){
				foreach my $e (@enzymes_arr){
					my ($s_sbo, $e_sbo) = ($species_sbo_hs{$s}, $species_sbo_hs{$e});
					$s_sbo = "NA" unless $s_sbo;
					$e_sbo = "NA" unless $e_sbo;
					print OUT "$pathway_id\t$reaction_id\t$reaction_name\t$species_hs{$s}\t$species_hs{$e}\tproducts\t$s\t$e\t$s_sbo\t$e_sbo\n";
				}
			}
		}
	}
	close SBML;
}

close OUT;