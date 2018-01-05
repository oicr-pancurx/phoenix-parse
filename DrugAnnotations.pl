#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use List::MoreUtils qw/uniq any/;

my $outdir = shift;
my $sample = shift;

my $summaryFile = "$outdir/$sample.summary.csv";
my $variantsFile = "$outdir/$sample.variants.csv";
my $annotationHTML = "$outdir/$sample.onco.html";
my $trialsHTML = "$outdir/$sample.trials.html";
my $oncoKBFormattedFile = "/.mounts/labs/PCSI/users/azhang/OncoAnnotation/OncoKBFormatted.tsv";
my $GrainneFormattedFile = "/.mounts/labs/PCSI/users/azhang/OncoAnnotation/GrainneFormatted.tsv";
my $CivicFormattedFile = "/.mounts/labs/PCSI/users/azhang/OncoAnnotation/CivicFormatted.tsv";
my $ClinicalTrialsFormattedFile = "/.mounts/labs/PCSI/users/azhang/OncoAnnotation/ClinicalTrialsFormatted.tsv";
my $DrugBankFormattedFile = "/.mounts/labs/PCSI/users/azhang/OncoAnnotation/DrugBankFormatted.tsv";


my %vars;

my $refSeqFile = "/oicr/data/genomes/homo_sapiens_mc/refSeq/refSeq_genes.bed";
#my %genes = lookUpGenePositions($refSeqFile, \%vars);


#$finalTargetHash{$gene $var_id $match} where var_id = $vars->{$gene}{variants}{x}
#$finalTargetHash{$gene $var_id $match}{var_id, gene, variant_class, prot_variant, dbsnp, germline, match, ab_counts
#database_name is updated and temporary, all info from vars hash,
#oncokb_papers, civic_papers = ; separated flattened text, oncokb_link, civic _link = hyperlink
#oncokb_drugs, civic_drugs, grainne_drugs = ; separated flattened text, oncokb_database_version, civic_, grainne_
#oncokb_oncogenicity, oncokb_consequence
#drugs-oncokb{$drug}{$cancer}{oncokb_levels}, drugs-civic{$drug}{$cancer}{@civic_evidence, @civic_support, @civic_outcome, @civic_rating, @civic_summary},
#drugs-grainne{$drug}{$cancer}{trials_id, trials_drugclass, trials_status, trials_sponsor}
my %finalTargetHash;
my %ClinicalTrials;
	#$ClinicalTrials{$gene}{@variants, @aa_context, @nuc_context, @mutation_type, copy_number, @ab_count}
	#$ClinicalTrials{$gene}{trials}{$trial_id}{title, trial_status, trial_summary, trial_description, @cancers, @drugs, @sponsors, URL}
my %DrugBankMatches;
	#$DrugBankMatches{$gene}{@variants, @aa_context, @nuc_context, @mutation_type, copy_number, ab_counts}
	#$DrugBankMatches{$gene}{drugs}{$drug_id}{drugbank_status, drugbank_name, URL}

my %outputhash;

parseVariantsFile($variantsFile, \%vars);
addClinicalTrials($ClinicalTrialsFormattedFile, \%vars, 'ClinicalTrials.GOV');
addDrugBank($DrugBankFormattedFile, \%vars, 'DrugBank.CA');
printTrialsHTML($trialsHTML, \%ClinicalTrials, 1, 'Clinical Trial');
printTrialsHTML($trialsHTML, \%DrugBankMatches, 2, 'Drug Bank');

addOncoAnno($oncoKBFormattedFile, \%vars, 'oncokb');
addOncoAnno($GrainneFormattedFile,  \%vars, 'grainne');
addOncoAnno($CivicFormattedFile, \%vars, 'civic');

#print Dumper %finalTargetHash;
printOncoHTML($annotationHTML, \%vars);

sub addOncoAnno {
	
	#"Grainne Variants are matched to ANY uncommon variants with same gene name and some type of mutation_type unless otherwise specified";
	#"Currently cannot match OncoKB exon x mutations, splice, and a couple other variant types to Annovar";
	
	my $file = shift;
	my $vars = shift;
	my $database = shift;

	my %columns;
	
	my @allheaders;
	
	open (FILE, "<", $file) or die;
	while (my $l = <FILE>) {
		chomp $l;
		
		#get all headers and the columns they're in
		if ($l =~ /variant_class/) {
			@allheaders = split /\t/, $l;
			for (my $i = 0; $i < scalar @allheaders; $i++) {
				my $header = $allheaders[$i];
				$columns{$header} = $i;
			}
		}
	
		else {
			
			#for each line in the database, it can match to multiple variants in the sample if the match is ambiguous
		
					my @var = @{matchOncoDBtoVars($l, \%vars, \%columns)};
					if (scalar @var > 0) {
				
				for my $var (@var) {
					my $gene = $finalTargetHash{$var}{gene};
					my $varpair = $finalTargetHash{$var}{var_id};
					my $match;
					if ($finalTargetHash{$var}{match} =~ /ambiguous/) {
						$match = 0.5;
					}
					else {
						$match = 1;
					}
					my $header;
					
					my @info = split /\t/, $l;
					#these headers are already added when matching variants, should be used for variant identification, but do not need to reannotated
					my @filledheaders = qw/gene var_id match variant_class prot_variant dbsnp germline database_name ab_counts/;
					my ($drug, $cancer, @cancers);
							
					#for each line of database info, there is 1 drug, and ; spliced cancers, | spliced info within the same cancer context
					#add other annotations to $finalTargetHash{$var}{$drug}{$cancer} subhash
					for (my $i = 0; $i < scalar @allheaders; $i++) {
						$header = $allheaders[$i];
						if ($header =~ /drug(s)?$/i) {
							$drug = $info[$i];
							$drug =~ s/,/ /g;
						}
						if ($header =~ /cancer/i) {
							$cancer = $info[$i];
							@cancers = split /;/, $cancer;
						}
					}
					for (my $i = 0; $i < scalar @allheaders; $i++) {
						$header = $allheaders[$i];
						
						unless (grep /$header/, @filledheaders) {
							
							#for cancer-drug specific info, i.e. levels, evidence_statement etc., embed it into "drugs-database" hash of hashes
							#the headers listed below like database_name, drug and cancer are already included in the hash naming
							#oncokb_consequence and oncogenicity etc. are not specific to the variant, independent of drug/cancer
							unless ($header =~ /drug(s)?$|cancers|consequence|oncogenicity|papers|pubmed|database_version|database_name|link$/i) {
								my @cancer_info = split /;/, $info[$i];
								for (my $c = 0; $c < scalar @cancers; $c++) {
									$cancer = $cancers[$c];
									$cancer =~ s/,/ /g;
									my $cancerinfo = $cancer_info[$c];
									$cancerinfo =~ s/,/ /g;
									
									@{$finalTargetHash{"$gene\t$varpair\t$match"}{"drugs-$database"}{$drug}{$cancer}{$header}} = (split /\|/, $cancerinfo);
								}
							}
							
							#variant info from each database, oncogenicity, consequence, papers, link and drugs are flattened into text and added to finalTargetHash
							#and vars hash
							#NOTE: oncokb_papers in particular are / / separated because they're not associated with cancer, so they're not | or ; separated
							elsif ($header =~ /oncokb_papers/) {
								my @papers = (split / /, $info[$i]);
								$finalTargetHash{"$gene\t$varpair\t$match"}{$header} = (split / /, $info[$i]);
								for my $paper (@papers) {
									if (exists $vars->{$gene}{variants}{$varpair}{$header}) {
										unless (grep /^$paper$/, $vars->{$gene}{variants}{$varpair}{$header}) {
											$vars->{$gene}{variants}{$varpair}{$header} .= "; $paper";
										}
									}
									else {
										$vars->{$gene}{variants}{$varpair}{$header} = "$paper";
									}
								}
							}
							#this is the civic equivalent of papers, but these are cancer/drug associated, so they are split by ; for cancer, and | within the cancer context
							elsif ($header =~ /pubmed/) {
								my @papers = split /[;|]/, $info[$i];
								for my $paper (@papers) {
									if (exists $vars->{$gene}{variants}{$varpair}{$database."_papers"}) {
										unless (grep /^$paper$/, $vars->{$gene}{variants}{$varpair}{$database."_papers"}) {
											$vars->{$gene}{variants}{$varpair}{$database."_papers"} .= "; $paper";
										}
									}
									else {
										$vars->{$gene}{variants}{$varpair}{$database."_papers"} = "$paper";
									}
								}
							}
							elsif ($header =~ /(drugs|drug)$/i) {
								if (exists $vars->{$gene}{variants}{$varpair}{$database."_drugs"}) {
									unless (grep /^$drug$/i, $vars->{$gene}{variants}{$varpair}{$database."_drugs"}) {
										$vars->{$gene}{variants}{$varpair}{$database."_drugs"} .= "; $drug";
									}
								}
								else {
									$vars->{$gene}{variants}{$varpair}{$database."_drugs"} = "$drug";
								}
							}
							elsif ($header =~ /link|oncogenicity|consequence/) {
								$finalTargetHash{"$gene\t$varpair\t$match"}{$header} = $info[$i];
								$vars->{$gene}{variants}{$varpair}{$header} = $info[$i];
							}
							#database_version
							elsif ($header =~ /version/) {
								$finalTargetHash{"$gene\t$varpair\t$match"}{lc("$database\_$header")} = $info[$i];
								$vars->{$gene}{variants}{$varpair}{lc("$database\_$header")} = $info[$i];
							}
							#so far nothing is missed, but may eventually need this eventually for generic header info to be added from finalTargetHash to vars
							elsif (exists $finalTargetHash{"$gene\t$varpair\t$match"}{$header}) {
								$vars->{$gene}{variants}{$varpair}{$header} = $finalTargetHash{"$gene\t$varpair\t$match"}{$header};
								print "$header has not been added to vars yet and will be added now\n";
							}
					}
				}
			}
		}
	}		
	}
	close FILE;
	#print Dumper (\%finalTargetHash);
}

#ignore regional mutations and splice variants, cannot match
#only match to variants that are not common
#parses common annotations fields (germline, variant_class, prot_variants, or dbsnp fields) to find the matching variant in the sample
#if there is a match, common fields and all info already known about variant is added to finalTargetHash
#variant and its matched condition (confident, ambiguous) are returned

sub matchOncoDBtoVars {
	
	my $l = shift;
	my $vars = shift;
	#indicates which column info is stored in $l
	my $columns = shift;
	#genes affected in sample
	
	my @matchedvars;

	unless ($l =~ /variant_class/) {
				
		my @info = split /\t/, $l;

		my ($gene, @variant_class, @vars, $dbsnp, $germline, $database_name);
		$database_name = $info[$columns->{database_name}];
		
		my $match = 0;
		
			$gene = $info[$columns->{gene}];
		
			@variant_class = split /\|/, $info[$columns->{variant_class}];
			@vars = split /\|/, $info[$columns->{prot_variant}];
			$dbsnp = $info[$columns->{dbsnp}];
			$germline = $info[$columns->{germline}];

			#civic has variants which require multiple conditions to be met, e.g. specific fusion is then mutated and amplified
			#each variant is searched for, and if all these conditions are met, then it is a true match (very rough implementation)
			my @samevar;
			if (exists $vars->{$gene}) {
			
			for my $varpair (sort keys %{$vars->{$gene}{variants}}) {
				my $rarity;
	
				if (not exists $vars->{$gene}{variants}{$varpair}{rarity}) {
					$rarity = 'NA';
				}
				else {
					$rarity = $vars->{$gene}{variants}{$varpair}{rarity};
				}
				
				#first eliminate common variants, and make sure germline info is a match
				if ($rarity ne 'common' and
					 (($germline eq 'NA') or ($vars->{$gene}{variants}{$varpair}{mutation_class} =~/$germline/))) {
					
					#then for each of the variants
					for (my $i = 0; $i < scalar @vars; $i++) {
						my $var = $vars[$i];
						my $variant_class = $variant_class[$i];
						my $samevar = 0;
							
						#match by dbsnp first; if so, update variant_class and prot_variant fields accordingly
						if ($vars->{$gene}{variants}{$varpair}{dbsnp} ne 'NA' and ($vars->{$gene}{variants}{$varpair}{dbsnp} eq $var)
							 and ($vars->{$gene}{variants}{$varpair}{mutation_type} !~ /altered promoter/)) {
							
							if ($germline eq 'NA' or ($vars->{$gene}{variants}{$varpair}{mutation_class} =~ /$germline/)) {
								$samevar = 1;
								if ($germline eq 'NA') {
									if ($vars->{$gene}{variants}{$varpair}{mutation_class} =~ /germline/) {
										$germline = 'germline';
									}
									else {
										$germline = 'somatic';
									}
									$info[$columns->{germline}]  = $germline;
								}
								if ($vars->{$gene}{variants}{$varpair}{aa_context} =~ /p.(.*)$/) {
									my @aa_contexts = split /\|/, $1;
									$var = $aa_contexts[0];
									
								}
								else {
									$var = 'NA';
								}
								$info[$columns->{prot_variant}] = $var;
								$variant_class = $vars->{$gene}{variants}{$varpair}{mutation_type};
								$info[$columns->{variant_class}] = $variant_class;
							}
						}
							
						#otherwise, try to match by prot_variant for nonsyn, indels etc.
						elsif ($vars->{$gene}{variants}{$varpair}{aa_context} ne 'NA') {
							my @aa_contexts = split /\|/, $vars->{$gene}{variants}{$varpair}{aa_context};
												
							if ($variant_class =~ /nonsynonymous/) {
								#1.W802D                        W802D    ok test
								if ($var =~ /^[A-Z]([0-9]+)[A-Z]$/) {
									 if (grep /p.$var$/, @aa_contexts) {
											$samevar = 1;
									 }
								}
								#2.W802*                        W802[A-Z] or frameshift or indel, assume is nonsynonymous but can't tell
								elsif ($var =~ /^([A-Z])([0-9]+)\*$/) {
									 if (grep /p.$1$2[a-z_][a-z0-9_]*$/i, @aa_contexts) {
											$samevar = 0.5;
									 }
								}
							}
												
							elsif ($variant_class =~ /indel/) {
								#1.D26_L28del   26_28del ok no test for trunc
								if ($var =~ /^([A-Z])([0-9]+)_([A-Z])([0-9]+)del$/) {
									if (grep /p.$2_$4del$/, @aa_contexts) {
										$samevar = 1;
									}
								}
								
								#2.26_28del   26_28del ok no test for trunc
								elsif ($var =~ /^([0-9]+)_([0-9]+)del$/) {
									if (grep /p.$var$/, @aa_contexts) {
										$samevar = 1;
									}
								}
								
								#3. V559delinsNP                V559delinsNP
								elsif ($var =~ /^([A-Z])([0-9]+)delins([A-Z]+)$/) {
									 if (grep /p.$var$/, @aa_contexts) {
											$samevar = 1;
									 }
								}
								
								#4. L747_T751delinsP             L747_T751delinsP ok
								elsif ($var =~ /^([A-Z])([0-9]+)_([A-Z])([0-9]+)delins([A-Z]+)$/) {
									 if (grep /p.$var$/, @aa_contexts) {
											$samevar = 1;
									 }
								}
								
								#5. 559delinsTNP            [A-Z]559delinsTNP
								elsif ($var =~ /^([0-9]+)delins([A-Z]+)$/) {
									 if (grep /p.[A-Z]$var$/, @aa_contexts) {
											$samevar = 1;
									 }
								}
								
								#6. D71delinsD*            D71delinsD[A-Z]
								elsif ($var =~ /^([A-Z])([0-9]+)delins([A-Z]+)\*$/) {
									 if (grep /p.$1$2delins$3[A-Z]$/, @aa_contexts) {
											$samevar = 0.5;
									 }
								}
							}
							
							#NOTE: frameshift in oncoKB may be registered as frameshift or other by Annovar
							elsif ($variant_class =~ /frameshift/) {
								#1.E23KfsX20  					E23fs
								if ($var =~ /^([A-Z])([0-9]+)([A-Z]?)fs/) {
									if (grep /p.$1$2fs/, @aa_contexts) {
										$samevar = 1;
									}
								}
							}
							
							elsif ($variant_class =~ /duplication/) {
								#14. D171_L174dup                D171delinsD..LD or L174delinsD..LD..L ok test
								if ($var =~ /^([A-Z])([0-9]+)_([A-Z])([0-9]+)dup$/) {
									 my $length = $4 - $2 - 1;
									 if ($length < 0) {
											$length = 0;
									 }
									 if (grep /p.($1$2delins$1[A-Z]{$length}$3$1|$3$4delins($1[A-Z]{$length}$3)\1)/, @aa_contexts) {
											$samevar = 1;
									 }
								}
							}
							
							#variant_class is not important if prot_variant is an exact match
							#but if there's no info other than gene name, is an ambiguous match e.g. Grainne List
							elsif (grep /$variant_class/ , qw/NA misc/) {
								if (grep /p.$var$/, @aa_contexts) {
									$samevar = 1;
								}
								
								elsif ($var eq 'NA' and ($vars->{$gene}{variants}{$varpair}{mutation_type} =~ /\w{3,}/)) {
									unless ($vars->{$gene}{variants}{$varpair}{mutation_type} =~ /altered promoter/) {
										$samevar = 0.5;
									}
								}
							}
						}

						#ambiguous match with only gene name may happen when variant is associated with no aa_context
						elsif (grep /$variant_class/, qw/NA misc/) {
							if ($var eq 'NA' and ($vars->{$gene}{variants}{$varpair}{mutation_type} =~ /^\w{3,}/)
								 and (not ($vars->{$gene}{variants}{$varpair}{mutation_type} =~ /altered promoter/))) {
								$samevar = 0.5;
							}
						}
						
						#then match by variant_class for mutation types without exact protein changes
						#13. stopgain
						if ($variant_class =~ /stopgain/) {
							if ($vars->{$gene}{variants}{$varpair}{mutation_type} =~ /stopgain/) {
									$samevar = 1;
							 }
						}
						elsif ($variant_class =~ /homozygous deletion/) {
							 if ($vars->{$gene}{variants}{$varpair}{mutation_type} =~ /homozygous deletion/) {
									$samevar = 1;
							 }
						}
						elsif ($variant_class =~ /amplification/) {
							if ($vars->{$gene}{variants}{$varpair}{mutation_type} =~ /amplification/) {
									$samevar = 1;
							 }
						}
						elsif ($variant_class =~ /fusion/) {
							if (exists $vars->{$gene}{variants}{$varpair}{fusion_genes}) {
								if ($var =~ /^$gene-([A-Z0-9]+)$/ and (grep /$1/, $vars->{$gene}{variants}{$varpair}{fusion_genes})) {
									$samevar = 1;
								}
								elsif ($var =~ /^([A-Z0-9]+)-$gene$/ and (grep /$1/, $vars->{$gene}{variants}{$varpair}{fusion_genes})) {
									$samevar = 1;
								}
							}
						}
	
						#18. nonspecfic fusions with another gene are ambiguous
						elsif ($var =~ /^Fusions$/) {
							if (exists $vars->{$gene}{variants}{$varpair}{fusion_genes}) {
								unless ($vars->{$gene}{variants}{$varpair}{fusion_genes} eq 'NA' or
										($vars->{$gene}{variants}{$varpair}{fusion_genes} eq '.') or
										($vars->{$gene}{variants}{$varpair}{fusion_genes} =~ /^\s*$/)) {
										 $samevar = 0.5;
								}
							}
						}
						#each variant might consistent of multiple concurrent variants e.g. BRAF V600E amplification, and these would be individually matched
						#match for each scenario is stored
						push (@samevar, $samevar);
						
					}
					#if any of the multiple concurrent variants are not matched, then it is not the same variant
					
					if (any{$_ eq '0'} @samevar) {
						$match = 0;
					}
					#if all variants matched, but at least one is ambiguous, then it is an ambiguous match
					elsif (any{$_ eq '0.5'} @samevar) {
						$match = 0.5;
					}
					else {
						$match = 1;
					}
			
		}

					#if there is a final match, ambiguous or exact, add basic variant id and ab_counts information to finalTargetHash
					#also add all info from vars to finalTargetHash as well
					if ($match > 0) {
						
						my @idheaders = qw/gene variant_class prot_variant dbsnp germline database_name var_id match/;
	
						for my $header (@idheaders) {
							
							if ($header eq 'gene') {
								$finalTargetHash{"$gene\t$varpair\t$match"}{$header} = $gene;
							}
							if ($header eq 'var_id') {
								$finalTargetHash{"$gene\t$varpair\t$match"}{$header} = $varpair;
							}
							elsif ($header eq 'match') {
								if ($match == 0.5) {
									$finalTargetHash{"$gene\t$varpair\t$match"}{$header} = "ambiguous";
								}
								elsif ($match == 1) {
									$finalTargetHash{"$gene\t$varpair\t$match"}{$header} = "exact";
								}
							}
							else {
								$finalTargetHash{"$gene\t$varpair\t$match"}{$header} = $info[$columns->{$header}];
							}
						}
						
						for my $header (keys %{$vars->{$gene}{variants}{$varpair}}) {
							$finalTargetHash{"$gene\t$varpair\t$match"}{$header} = $vars->{$gene}{variants}{$varpair}{$header};
						}
						$finalTargetHash{"$gene\t$varpair\t$match"}{ab_counts} = $vars->{$gene}{ab_counts};
						push (@matchedvars, "$gene\t$varpair\t$match");
					}
				}
			}
				 
	}	
	return (\@matchedvars);
}


#reformats and condenses variant annotation from finaldrughash into more compact outputhash containing only output headers of interest
#separates outputhash into confidently matched and ambiguously matched for printing out HTML tables
sub printOncoHTML {
	my $htmlFile = shift;
	my $vars = shift;
	
	#$outputhash{$gene\t$var_id\$match}{no, match (exact, ambiguous), gene, origin (somatic, germline or NA), mutation_type, position,
	#	prot_variant, copies (ab_count), consequence = G/LOF, oncogenicity
	#	OncoKBDrugs{$drug}{$cancer}{oncokb_levels}, CivicDrugs{$drug}{$cancer}{@civic_summary, @_support, @_outcome, @_rating, @_summary},
	#	GrainneDrugs{$drug}{$cancer}{trials_id, trials_drugclass, trials_status, trials_sponsor}
	my @outputheaders = qw/no match gene origin mutation_type position prot_variant copies consequence oncogenicity OncoKBDrugs GrainneDrugs CivicDrugs/;
	
	my $count = 1;
	my %alias = ( "consequence"=> "oncokb_consequence",
					 "oncogenicity" => "oncokb_oncogenicity",
					 "OncoKBDrugs" => "drugs-oncokb",
					 "GrainneDrugs" => "drugs-grainne",
					 "CivicDrugs" => "drugs-civic" ,);
	
	for my $var (sort keys %finalTargetHash) {
		my $gene = $finalTargetHash{$var}{gene};
		my $varpair = $finalTargetHash{$var}{var_id};
		my (@nuc_context, @aa_context, @currdrugs, @cancers, @sponsors, @trialdrugs, @trialdrugclass, @trialstatus, @trialsid);
		my (@trials_germline, @trials_cancer);

		my @trials;
		#this information is needed for HTML output
		unless (exists $outputhash{$var}{position}) {
			$outputhash{$var}{position} = $finalTargetHash{$var}{position};
			$outputhash{$var}{match} = $finalTargetHash{$var}{match};
			$outputhash{$var}{var_id} = $finalTargetHash{$var}{var_id};
		}
		
		for my $header (@outputheaders) {
			if ($header eq 'no') {
				$outputhash{$var}{no} = $count;
			}
			
			#prot_variant comes from the actual sample, rather than the prot_variant matching criteria from the database
			elsif ($header eq 'prot_variant') {
				if (exists $finalTargetHash{$var}{aa_context}) {
					for my $aacontext (split /\|/, $finalTargetHash{$var}{aa_context}) {
						if ($aacontext =~ /p.(.*)$/) {
							unless (grep /$1/, @aa_context) {
								push (@aa_context, $1);
							}
						}
						else {
							push (@aa_context, "NA");
						}
					}
				}
				else {
					push (@aa_context, "NA");
				}
				$outputhash{$var}{$header} = \@aa_context;
			}
			elsif ($header eq 'origin') {
				if ($finalTargetHash{$var}{mutation_class} =~ /somatic/) {
					$outputhash{$var}{$header} = "somatic";  
				}
				elsif ($finalTargetHash{$var}{mutation_class} =~ /germline/) {
					$outputhash{$var}{$header} = "germline"; 
				}
				else {
					$outputhash{$var}{$header} = "NA"; 
				}
			}
			elsif ($header eq 'copies') {
				if (defined $finalTargetHash{$var}{ab_counts}) {
					$outputhash{$var}{$header} = $finalTargetHash{$var}{ab_counts};
				}
				else {
					$outputhash{$var}{$header} = "NA";
				}
			}
			#some renaming from finalTargetHash to outputhash to make headers more concise
			else {
				if (exists $alias{$header}) {
					if (defined $finalTargetHash{$var}{$alias{$header}}) {
						if (scalar $finalTargetHash{$var}{$alias{$header}} =~ /HASH/) {
							$outputhash{$var}{$header} = \%{$finalTargetHash{$var}{$alias{$header}}};
						}
						else {
						$outputhash{$var}{$header} = $finalTargetHash{$var}{$alias{$header}};
						}
					}
				}
				else {
					$outputhash{$var}{$header} = $finalTargetHash{$var}{$header};
				}
			}
		}
		$count++;
	}
	
	#print HTML style
	open (HTML, ">", $htmlFile) or die;
	print HTML "<!DOCTYPE html>
   <html>
   <head>
	<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
	<link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css\">
   <script src=\"https://ajax.googleapis.com/ajax/libs/jquery/3.2.1/jquery.min.js\"></script>
   <script src=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js\"></script>
   <style type=\"text/css\">
   table {
   border-collapse: collapse !important;
   }
   table.solid{
      border-style: solid;
   }
   /* Sortable tables */
   table.sortable thead {
           background-color:#FFFFFF;
           color:#000000;
           font-weight: bold;
           cursor: default;
   }
	\@media screen {
	#main {
		max-width: 1400px;
		margin: auto;
	}
	#header {
		max-width: 1500px;
		margin: auto;
	}
	#footer {
		max-width: 1500px;
		margin: auto;
	}
	}
	th {
		border-bottom: 1px solid grey; border-top: 2px solid grey; text-align:left;
		padding-right: 0.2em;
	}
	/* striped table */
	tr:nth-child(even) {background-color:#f2f2f2}
	tr {
		vertical-align:top;
	}
	tr.spaceUnder>td {
		padding-bottom: 0.2em;
		padding-right: 0.3em;
		font-size: 85%;
		word-wrap: break-all;
	}

	h1 {
      color: #FFFFFF;
      background-color: #444444;
      text-align: center;
      line-height: 1.8;
      border-radius: 10px;
      page-break-after: avoid;
   }
   
   h2 {
      color: #FFFFFF;
      background-color: #444444;
      text-align: center;
      line-height: 1.5;
      border-radius: 10px;
      page-break-after: avoid;
   }
   
   h3 {
      color: #000000;
      background-color: #DDDDDD;
      text-align: center;
      line-height: 1.5;
      border-radius: 10px;
      page-break-after: avoid;
   }
	a {
		color: #00000;
	}
	.button {
    border: none;
	 background-color:Transparent;
	 margin: 0;
	 padding: 0px 2px 0px 0px;
	 display:block;
	}
	</style>
	</head>\n";
	
	#print HTML for confidently and ambiguously matched tables separately
	print HTML "<div id = \"main\">
	<h2>$sample Hypothetically Actionable Genomic Variants</h2>
	<h3>Confidently Matched Variants</h3>
	</div>\n";

	my (%exactmatched, %ambiguousmatched);
	
	my $no = 1;
	for my $var (sort keys %outputhash) {
		if ($outputhash{$var}{match} eq 'exact') {
			$exactmatched{$var} = $outputhash{$var};
			$exactmatched{$var}{no} = $no;
			$no++;
		}
	}
	
	$no = 1;
	for my $var (sort keys %outputhash) {
		if ($outputhash{$var}{match} eq 'ambiguous') {
			$ambiguousmatched{$var} = $outputhash{$var};
			$ambiguousmatched{$var}{no} = $no;
			$no++
		}
	}
	close HTML;
	
	printOncoHTMLTable(\%exactmatched, $htmlFile);
	open (HTML, ">>", $htmlFile);
	print HTML "<br>\n";
	print HTML "<div id = \"main\">
	<h3>Ambiguously Matched Variants</h3>
	</div>\n";
	close HTML;
	printOncoHTMLTable(\%ambiguousmatched, $htmlFile);
	
	open (HTML, ">>", $htmlFile);
	
	print HTML "</html>\n";
	
	close HTML;
	
}

#print HTML tables for confidently matched and ambiguous matched variants
#tables expand on click with more info on cancers drug has been tested in, and evidence(s) found in each of those cancers
sub printOncoHTMLTable {
	my $suboutputhash = shift;
	my $htmlfile = shift;
	my %outputhash = %{$suboutputhash};
	#set column widths for each header type
	my %colwidth = ('no' => '3%',
						 'gene' => '6%',
						 'origin' => '6%',
						 'mutation_type' => '9.5%',
						 'consequence' => '10%',
						 'position' => '7%',
						 'prot_variant' => '6%',
						 'copies' => '6%',
						 'oncogenicity' => '10%',
						 'OncoKBDrugs' => '24.5%',
						 'GrainneDrugs' => '15%',
						 'CivicDrugs' => '20%');
	
	
	my @outputheaders = qw/no gene origin mutation_type consequence position prot_variant copies oncogenicity OncoKBDrugs CivicDrugs GrainneDrugs /;
	
	open (HTML, ">>", $htmlfile) or die ("can't open html output file");
	print HTML "<table style = \"margin: 1em auto; width: 1250px; float:center; table-layout:fixed; word-wrap:break-word;\">\n";
	print HTML "<tr padding-bottom:1em>\n";
	#further condensing of headers, inelegant but works
	for my $header (@outputheaders) {
		if ($header =~ /prot_variant/) {
			print HTML "<th width=\"$colwidth{$header}\">protein</th>\n";
		}
		elsif ($header =~ /mutation_type/) {
			print HTML "<th width=\"$colwidth{$header}\">mutation</th>\n";
		}
		else {
			print HTML "<th width=\"$colwidth{$header}\">$header</th>\n";
		}
	}
	print HTML "</tr>\n";
	
	#get the max number of rows needed for each variant i.e. list each prot_variant (if multiple isoforms) and drug on a separate line
	for my $var (sort keys %outputhash) {
		my @drugs;
		my $maxrowspan = 0;
		for my $header (@outputheaders) {
			if (exists $outputhash{$var}{$header}) {
				if (defined $outputhash{$var}{$header}) {
					if (scalar $outputhash{$var}{$header} =~ /HASH/) {
					my $rowspan = scalar keys %{$outputhash{$var}{$header}};
						if ($rowspan > $maxrowspan) {
							$maxrowspan = $rowspan;
						}
					}
					elsif (scalar $outputhash{$var}{$header} =~ /ARRAY/) {
						my $rowspan = scalar @{$outputhash{$var}{$header}};
						if ($rowspan > $maxrowspan) {
							$maxrowspan = $rowspan;
						}
					}
				}
			}
			else {
				$outputhash{$var}{$header} = "NA";
			}
		}
		
		#simple variant annotations like name, origin, gene are scalars, prot_variants are in arrays (splice variants)
		#drug info is stored in multi layered hashes, to allow for clickable expansions on its corresponding cancer models, and evidences
		print HTML "<tr class = \"spaceUnder\" rowspan = \"$maxrowspan\">\n";
		for my $header (@outputheaders) {
			my $width = $colwidth{$header};
			if (exists $outputhash{$var}{$header}) {
				if (defined $outputhash{$var}{$header}) {
					#for prot_variants in arrays, if there are multiple isoforms, print each as separate element on diff lines
					if (scalar $outputhash{$var}{$header} =~ /ARRAY/) {
						print HTML "<td>\n";
						
						for (my $i = 0; $i < scalar @{$outputhash{$var}{$header}}; $i++) {
							my $element = $outputhash{$var}{$header}[$i];
							
							if ($i == scalar @{$outputhash{$var}{$header}} - 1) {
								print HTML "$element";
							}
							else {
								print HTML "$element<br>";
							}
						}
						print HTML "</td>\n";
					}
					#cancer associated info are stored in hashes
					elsif (scalar $outputhash{$var}{$header} =~ /HASH/) {		
						print HTML "<td width=\"$width\">\n";
						
						#print hyperlink associated with the variant in the database (does not exist for grainne list)
						if ($header =~ /(.*)Drugs$/i) {
							my $database = lc($1);
							@drugs = sort (keys %{$outputhash{$var}{$header}});
							unless ($database =~ /grainne/i) {
								print HTML "<a href=$finalTargetHash{$var}{$database.\"_link\"}> link </a><br>";
							}
						}
						#for each drug, create an expandable div for cancers this drug has been tested in
						#then for each cancer, create an expandable div for cancer related info relevant to the database
						for (my $i = 0; $i < scalar @drugs; $i++) {
							my $drug = $drugs[$i];
							if ($drug eq 'NA') {
								print HTML "NA";
							}
							else {
								my $drugcode = "$outputhash{$var}{match}"."$outputhash{$var}{no}".$header."drug$i";
								my @cancers = sort (keys %{$outputhash{$var}{$header}{$drug}});
								my $cancerrowspan = scalar @cancers;
								
								print HTML "<a href =\"#$drugcode\" data-toggle=\"collapse\"> <font color=\"000000\"><u>$drug</u></font></a><br>\n";
							
								print HTML "<div id =\"$drugcode\" class=\"collapse\" rowspan=\"$cancerrowspan\">\n";
								
								for (my $j = 0; $j < scalar @cancers; $j++) {
									my $cancer = $cancers[$j];
									my $cancercode = $outputhash{$var}{match}.$outputhash{$var}{no}.$header."drug$i"."cancer$j";
									if ($header eq 'CivicDrugs') {
										print HTML "
										<a href = \"#$cancercode\" data-toggle=\"collapse\"><i><font color=\"000000\">$cancer</font></i></a><br>
										<div id = \"$cancercode\" style=\"margin-left: 1em\" class=\"collapse\">";
										if (defined $outputhash{$var}{$header}{$drug}{$cancer}{civic_summary}) {
											for (my $i = 0; $i < scalar @{$outputhash{$var}{$header}{$drug}{$cancer}{civic_summary}}; $i++) {
												print HTML "Evidence: $outputhash{$var}{$header}{$drug}{$cancer}{civic_evidence}[$i]<br>
												Prediction: $outputhash{$var}{$header}{$drug}{$cancer}{civic_support}[$i]
												$outputhash{$var}{$header}{$drug}{$cancer}{civic_outcome}[$i]<br>
												Rating: $outputhash{$var}{$header}{$drug}{$cancer}{civic_rating}[$i]<br>
												Summary: $outputhash{$var}{$header}{$drug}{$cancer}{civic_summary}[$i]<br>
													";
											}
										}
										else {
											print HTML "NA";
										}
										print HTML "</div>\n";
									}
									elsif ($header eq 'OncoKBDrugs') {
										print HTML "
										<a href = \"#$cancercode\" data-toggle=\"collapse\"><i><font color=\"000000\">$cancer</font></i></a>
										<div id = \"$cancercode\" style=\"text-indent: 1em\" class=\"collapse\">
										Evidence Level: @{$outputhash{$var}{$header}{$drug}{$cancer}{oncokb_levels}}
										</div>\n";
									}
									elsif ($header eq 'GrainneDrugs') {
										print HTML "
										<a href = \"#$cancercode\" data-toggle=\"collapse\"><i><font color=\"000000\">$cancer</font></i></a>
										<div id = \"$cancercode\" class=\"collapse\" style=\"margin-left: 1em;\">
										";
										
										if (defined $outputhash{$var}{$header}{$drug}{$cancer}{trials_id}) {
											print HTML "Trial ID: @{$outputhash{$var}{$header}{$drug}{$cancer}{trials_id}}<br>
											Drug Class: @{$outputhash{$var}{$header}{$drug}{$cancer}{trials_drugclass}}<br>
											Status: @{$outputhash{$var}{$header}{$drug}{$cancer}{trials_status}}<br>
											Sponsor: @{$outputhash{$var}{$header}{$drug}{$cancer}{trials_sponsor}}
											</div>\n";
										}
										else {
											print HTML "Drug Class: @{$outputhash{$var}{$header}{$drug}{$cancer}{trials_drugclass}}
											</div>\n";
										}
									}
								}
							
								print HTML "</div>\n";	
							}
						}
						print HTML "</td>\n";
					}
					else {
						if ($header eq 'origin') {
							print HTML "<td width=\"$width\">$outputhash{$var}{$header}</td>\n";
						}
						elsif ($header eq 'position') {
							my $position = (split (/\|/, $outputhash{$var}{$header}))[0];
							print HTML "<td width=\"$width\">$position</td>\n";
						}
						elsif (exists $outputhash{$var}{$header}) {
							print HTML "<td width=\"$width\">$outputhash{$var}{$header}</td>\n";
						}
						else {
							print HTML "<td width=\"$width\">NA</td>\n";
						}
					}
				}	
			
			}
		}
		print HTML "</tr>\n";
	}
	print HTML "</table>\n";
	close HTML;
}


sub addDrugBank {
	my $file = shift;
	my $vars = shift;
	my $database = shift;
	
	my %columns;
	my @allheaders;
	
	open (FILE, "<", $file) or die;
	while (my $l = <FILE>) {
		chomp $l;
		if ($l =~ /variant_class/) {
			@allheaders = split /\t/, $l;
			
			for (my $i = 0; $i < scalar @allheaders; $i++) {
				my $header = $allheaders[$i];
				$columns{"$header"} = $i;
			}
		}
		else {
			my @info = split /\t/, $l;
			my $gene = $info[0];
			
			if (exists $vars->{$gene}) {
				for my $varpair (sort keys %{$vars->{$gene}{variants}}) {
					if ($vars->{$gene}{variants}{$varpair}{rarity} ne 'common' and
						 ($vars->{$gene}{variants}{$varpair}{mutation_type} !~ /altered promoter/)) {
						push (@{$DrugBankMatches{$gene}{variants}}, $varpair);
						push (@{$DrugBankMatches{$gene}{nuc_context}}, $vars->{$gene}{variants}{$varpair}{nuc_context});
						push (@{$DrugBankMatches{$gene}{aa_context}}, $vars->{$gene}{variants}{$varpair}{aa_context});
						push (@{$DrugBankMatches{$gene}{mutation_type}}, $vars->{$gene}{variants}{$varpair}{mutation_type});
						push (@{$DrugBankMatches{$gene}{ab_counts}}, $vars->{$gene}{variants}{$varpair}{ab_counts});
						$DrugBankMatches{$gene}{copy_number} = $vars->{$gene}{variants}{$varpair}{copy_number};
						
						
						my @drugs = split /\|/, $info[$columns{drugbank_id}];
						my @status = split /\|/, $info[$columns{drugbank_status}];
						my @name = split /\|/, $info[$columns{drugbank_name}];
						#print Dumper @name;
						for (my $j = 0; $j < scalar @drugs; $j++) {
							my $drug_id = $drugs[$j];
							my $url = "https://www.drugbank.ca/drugs/".$drug_id;
							my $status = $status[$j];
							my $name = $name[$j];
							
							$DrugBankMatches{$gene}{drugs}{$drug_id}{URL} = $url;
							$DrugBankMatches{$gene}{drugs}{$drug_id}{drugbank_status} = $status;
							$DrugBankMatches{$gene}{drugs}{$drug_id}{drugbank_name} = $name;
						}
					}
				}
			}
		}
	}
}


sub printTrialsHTML {
		
	my $htmlFile = shift;
	my $datahash = shift;
	my $order = shift; #is it first (1) or not (2)
	my $label = shift; #Clinical Trial or Drug Bank
	
	if (`grep '</html>' $htmlFile`) {
		`sed -i 's%</html>%%g' $htmlFile`;
	}
	
	#print HTML doctype and CSS only if this is the first database to be printed out
	if ($order == 1) {
		open (HTML, ">", $htmlFile) or die;
		print HTML "<!DOCTYPE html>
		<html>
		<head>
		<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
		<link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css\">
		<script src=\"https://ajax.googleapis.com/ajax/libs/jquery/3.2.1/jquery.min.js\"></script>
		<script src=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js\"></script>
		<style type=\"text/css\">
		table {
		border-collapse: collapse !important;
		}
		table.solid{
			border-style: solid;
		}
		/* Sortable tables */
		table.sortable thead {
				  background-color:#FFFFFF;
				  color:#000000;
				  font-weight: bold;
				  cursor: default;
		}
		\@media screen {
		#main {
			max-width: 1400px;
			margin: auto;
		}
		#header {
			max-width: 1500px;
			margin: auto;
		}
		#footer {
			max-width: 1500px;
			margin: auto;
		}
		}
		tr.whole {
			border-bottom: 1px solid grey; border-top: 2px solid grey; text-align:left;
			padding-right: 0.2em;
		}
		/* striped table */
		tr.whole:nth-child(even) {background-color:#f2f2f2}
		tr {
			vertical-align:top;
		}
		tr.spaceUnder>td {
			padding-bottom: 0.2em;
			padding-right: 0.3em;
			font-size: 85%;
			word-wrap: break-all;
		}
	
		h1 {
			color: #FFFFFF;
			background-color: #444444;
			text-align: center;
			line-height: 1.8;
			border-radius: 10px;
			page-break-after: avoid;
		}
		
		h2 {
			color: #FFFFFF;
			background-color: #444444;
			text-align: center;
			line-height: 1.5;
			border-radius: 10px;
			page-break-after: avoid;
		}
		
		h3 {
			color: #000000;
			background-color: #DDDDDD;
			text-align: center;
			line-height: 1.5;
			border-radius: 10px;
			page-break-after: avoid;
		}
		a {
			color: #00000;
		}
		.button {
		 border: none;
		 background-color:Transparent;
		 margin: 0;
		 padding: 0px 2px 0px 0px;
		 display:block;
		}
		</style>
		</head>\n";
		
		print HTML "<div id = \"main\">
		<h2>$sample Hypothetically Actionable Genomic Variants</h2></div>\n";
		
		close HTML;
	}
	
	
	
	
	
	
	open (HTML, ">>", $htmlFile) or die;
	my $targets = scalar (keys %{$datahash});
	
	#print HTML for confidently and ambiguously matched tables separately
	
	print HTML "<h3>$targets Hypothetical Gene Targets in $label</h3><br>\n";
	
	my $no = 1;
	for my $gene (sort keys %{$datahash}) {
		$datahash->{$gene}{no} = $no;
		$no++;
	}
	
	print HTML "<table style = \"margin: 1em auto; width: 1250px; float:center; table-layout:fixed; word-wrap:break-word;\">\n";
	
	#print headers########################################################
	print HTML "<tr padding-bottom:1em class=\"whole\">\n";
	#define headers and column widths
	my %colwidths = ("No" => "8%",
						 "Gene" => "10%",
						 "Gene Variants in Sample" => "25%",
						 "Clinical Trial Info" => "40%",
						 "Drug Bank Info" => "40%");
	my @headers = ("No", "Gene", "Gene Variants in Sample", "$label Info");
	for (my $i = 0; $i < scalar @headers; $i++) {
		my $header = $headers[$i];
		my $width = $colwidths{"$header"};
		print HTML "<th width=$width>$header</th>\n";
	}
	print HTML "</tr>\n";
		
	#print table##########################################################
	for my $gene (sort keys %{$datahash}) {
		my $numvariants = scalar @{$datahash->{$gene}{variants}};
		my $rowsofvar = $numvariants + 1;
		my $numdatabase;
		
		if ($label =~ /Trial/) {
			$numdatabase = scalar (keys %{$datahash->{$gene}{trials}});
		}
		elsif ($label =~ /Drug/) {
			$numdatabase = scalar (keys %{$datahash->{$gene}{drugs}});
		}
	
		print HTML "<tr class = \"spaceUnder\">\n";
		print HTML "<td>$datahash->{$gene}{no}</td>\n";
		print HTML "<td>$gene</td>\n";
		
		#print variants###################################################
		#variants for each gene is a collapsible div class
		print HTML "<td><a href=\"#$gene-Variants\" data-toggle=\"collapse\">$numvariants Variants</a>\n";
		print HTML "<div id=\"$gene-Variants\" class=\"collapse\" rowspan=\"$rowsofvar\">\n";
		#the collapsible div class has an embedded table where each row is a variant
		print HTML "<table border=\"0\">\n";
		print HTML "<tr padding-bottom:1em>\n";
		print HTML "<th width=\"70%\">Variant Info</th>\n";
		print HTML "<th width=\"10%\">Copies</th>\n";
		print HTML "<th width=\"10%\">Alleles</th>\n</tr>\n";
		#print out condensed, unique nuc and aa_contexts
		for (my $i = 0; $i < $numvariants; $i++) {
			$datahash->{$gene}{nuc_context}[$i] =~ s/NM_[0-9]*://g;
			#my @nuc_context = uniq (split /\|/, $datahash->{$gene}{nuc_context}[$i]);
			my $nuc_context = (uniq (split /\|/, $datahash->{$gene}{nuc_context}[$i]))[0];
			
			$datahash->{$gene}{aa_context}[$i] =~ s/NM_[0-9]*://g;
			#my @aa_context = uniq (split /\|/, $datahash->{$gene}{aa_context}[$i]);
			my $aa_context = (uniq (split /\|/, $datahash->{$gene}{aa_context}[$i]))[0];
			
			print HTML "<tr>\n";
			print HTML "<td>$datahash->{$gene}{mutation_type}[$i] ($nuc_context; $aa_context)</td>\n";
			print HTML "<td>$datahash->{$gene}{copy_number}</td>\n";
			print HTML "<td>$datahash->{$gene}{ab_counts}[$i]</td>\n";
			print HTML "</tr>\n";
		}
		print HTML "</table>\n";
		print HTML "</div>\n";
		print HTML "</td>\n";
		
		#print trials#####################################################
		
		if ($label =~ /Trial/) {
			print HTML "<td width=\"$colwidths{\"$label Info\"}\"><a href=\"#$gene-Trials\" data-toggle=\"collapse\">$numdatabase Trials</a>\n";
			print HTML "<div id=\"$gene-Trials\" class=\"collapse\" rowspan=\"$numdatabase\">\n";
			
			for my $dataid (sort keys %{$datahash->{$gene}{trials}}) {
				print HTML "<a href=$datahash->{$gene}{trials}{$dataid}{URL}>$dataid</a>&emsp;
				$datahash->{$gene}{trials}{$dataid}{trial_status}\n<br>";
			}
		}
		elsif ($label =~ /Drug/) {
			print HTML "<td width=\"$colwidths{\"$label Info\"}\"><a href=\"#$gene-Drugs\" data-toggle=\"collapse\">$numdatabase Drugs</a>\n";
			print HTML "<div id=\"$gene-Drugs\" class=\"collapse\" rowspan=\"$numdatabase\">\n";
			
			#drug bank info has 3 columns, ID, status, name
			for my $dataid (sort keys %{$datahash->{$gene}{drugs}}) {
				print HTML "<a href=$datahash->{$gene}{drugs}{$dataid}{URL}>$dataid</a>
				&emsp;$datahash->{$gene}{drugs}{$dataid}{drugbank_status}
				&emsp;$datahash->{$gene}{drugs}{$dataid}{drugbank_name}
				&emsp;<a href=$datahash->{$gene}{drugs}{$dataid}{URL}#clinical-trials>trials</a><br>\n";
			}
		}
			
			print HTML "</div>\n";
			print HTML "</td>\n";
			
	}
	
	print HTML "</table>\n";
	
	print HTML "</html>\n";
	
	close HTML;
}
	
sub addClinicalTrials {
	my $file = shift;
	my $vars = shift;
	my $database = shift;

	my %alias = ("NCT Number" => "trial_id",
					 "Title" => "title",
					 "Recruitment" => "trial_status",
					 "Conditions" => "cancers",
					 "Interventions" => "drugs",
					 "Sponsor/Collaborators" => "sponsors",
					 "Summary" => "trial_summary",
					 "Full_Description" => "trial_description",
					 "URL" => "URL");
	
	my %columns;
	
	my @allheaders;
	
	open (FILE, "<", $file) or die;
	while (my $l = <FILE>) {

		chomp $l;
		
		#get all headers and the columns they're in
		if ($l =~ /NCT Number/) {
			@allheaders = split /\t/, $l;
			for (my $i = 0; $i < scalar @allheaders; $i++) {
				my $header = $allheaders[$i];
				$columns{$header} = $i;
			}
		}
	}
	
	#the gene name can be framed by multiple delimiters, should search for all these potential variations on the gene name
	my @delim = (',', ':', ' ', '-', ';');
	
	for my $gene (sort keys %vars) {
		for my $varpair (sort keys %{$vars->{$gene}{variants}}) {
			if ($vars->{$gene}{variants}{$varpair}{rarity} ne 'common' and
				 ($vars->{$gene}{variants}{$varpair}{mutation_type} !~ /altered promoter/)) {
				
				my @patterns;

				for my $delim1 (@delim) {
					for my $delim2 (@delim) {
					 my $pattern = $delim1.$gene.$delim2;
					 push (@patterns, $pattern);
					}
				}
				
				my $grepword = join ('\|', @patterns);

				my @results = `grep '$grepword' $ClinicalTrialsFormattedFile`;
				
				unless (scalar @results == 0) {
				my $match = 0.5;
				my ($aa_context, $nuc_context) = "NA";
					
						unless (grep /$varpair/, @{$ClinicalTrials{$gene}{variants}}) {
							push (@{$ClinicalTrials{$gene}{variants}}, $varpair);
							if (defined $vars->{$gene}{variants}{$varpair}{nuc_context}) {
								$nuc_context = $vars->{$gene}{variants}{$varpair}{nuc_context};	
							}
							if (defined $vars->{$gene}{variants}{$varpair}{aa_context}) {
								$aa_context = $vars->{$gene}{variants}{$varpair}{aa_context};	
							}
							push (@{$ClinicalTrials{$gene}{nuc_context}}, $nuc_context);
							push (@{$ClinicalTrials{$gene}{aa_context}}, $aa_context);
							push (@{$ClinicalTrials{$gene}{mutation_type}}, $vars->{$gene}{variants}{$varpair}{mutation_type});
							push (@{$ClinicalTrials{$gene}{ab_counts}}, $vars->{$gene}{variants}{$varpair}{ab_counts});
							$ClinicalTrials{$gene}{copy_number} = $vars->{$gene}{variants}{$varpair}{copy_number};
						
						
						for my $result (@results) {
							my @info = split /\t/, $result;
							my $trialid = $info[0];
							unless (grep /$varpair/, @{$ClinicalTrials{$gene}{variants}}) {
								push (@{$ClinicalTrials{$gene}{variants}}, $varpair);
							}
							
							for my $header (sort keys %alias) {
								my $outheader = $alias{$header};
								my @simpleheaders = qw/title trial_summary trial_status trial_description URL/;
								my @compoundheaders = qw/drugs cancers sponsors/;
								if (grep /$outheader/, @simpleheaders) {
									$ClinicalTrials{$gene}{trials}{$trialid}{$outheader} = $info[$columns{$header}];
								}
								elsif (grep /$outheader/, @compoundheaders) {
									@{$ClinicalTrials{$gene}{trials}{$trialid}{$outheader}} = split /\|/, $info[$columns{$header}];
								}
							}
							
							#what if there are multiple variants for a single gene; should create 1 entry per gene by compiling all variants into 1 entry
							#what if there are multiple trials for a single gene; should add by unique trial IDs
							#should do the same for DrugBank ID; rewrite from scratch since OncoKB is mostly taken care of by cBioPortal etc.
							
						}
					}
				
				}

			}
		}
	
	}
}


sub parseVariantsFile {
	
	my $file = shift;
	my $vars = shift;
	
	my @headers;
	my %columns;
	
	open (FILE, "<", $file) or die;
	while (my $l = <FILE>) {
		chomp $l;
		
		if ($l =~ /^donor/) {
			@headers = split /,/, $l;
			
			for (my $i = 0; $i < scalar @headers; $i++) {
				my $header = $headers[$i];
				$columns{$header} = $i;
			}
		}
		else {
			my @info = split /,/, $l;
			my $mut_type = $info[$columns{mutation_type}];
			my $mut_class = $info[$columns{mutation_class}];
			my $gene = $info[$columns{gene}];
			my $base_change = $info[$columns{base_change}];
			my $position = $info[$columns{position}];
			my $var = "$mut_class\t$mut_type\t$position\t$base_change";
			
			my @headersofinterest = qw/copy_number ab_counts mutation_class mutation_type position fusion_genes base_change tumour_freq/;
			push (@headersofinterest, qw/normal_freq nuc_context aa_context dbsnp cosmic rarity clinvar cosmic_census_data/);
			
			unless ($mut_type eq 'NA') {
				for my $hoi (@headersofinterest) {
					$vars->{$gene}{variants}{$var}{$hoi} = $info[$columns{$hoi}];
				}
			}
		}
		
	}
	#print Dumper (%columns);
	
	
}
