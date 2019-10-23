#! /usr/bin/env bash
#SBATCH --cpus=1
#SBATCH --mem=2gb
#SBATCH --time=7-00:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

############################################# README

# You need to create a file like rdm_models_info given in this repository
# You need to create a file like conf_pats_ids not given in this repository

############################################# README



source ~soft_bio_267/initializes/init_R
module load ruby/2.4.1
currDir=`pwd`
export PATH=$currDir'/scripts':$PATH
export PATH=$current_dir'/sys_bio_lab_scripts':$PATH


# SELECT NETWORK
original_file=processed_data/patient_data.txt
results_source=$SCRATCH'/analysed_Networks/'
source_networks=$SCRATCH'/build_Networks'
original_rels=$source_networks/get_network_nodes.rb_0000/net.txt
raw_original_rels=$source_networks/get_network_nodes.rb_0000/net_notE.txt
outf="results"

# CONGIF VARIABLES
th_enrichments=1e-03
th_rdm_filters=0.01
rdm_models=100
diff_rdms=3


# OUTPUT FOLD
mkdir results

#### UNIFY DATA

# Bring net info (2pats filter)
cp $raw_original_rels $outf'/source_net'
cp $original_rels $outf'/net'

# Obtain HPO of used networks
awk '$1 ~ /HP/ {print $1} ' $outf'/net' | uniq > $outf'/net_hpos'

# Obtain networks stats
hpo_stats.R -s $outf'/source_net' -n $outf'/net' -o $outf'/net_stats'


# Bring LOCIS info
# Filtered (after HyI)
cp $source_networks'/NetAnalyzer.rb_0000/loci2phen_coords_filtered.txt' $outf'/raw_sors'
awk '{print $1 "\t" $2 "\t" $3 "\t" (($3 - $2)) "\t" $4 "\t" $5 "\t" $6}' $outf'/raw_sors' > $outf'/real_sors'

# Unfiltered (source patients file, NO FILTERS)
awk '{print $2 "\t" $3 "\t" $4 "\t" (($4 - $3))}' $original_file > $outf'/source_cnvs'

# Unfiltered (validated by 2 patients)
cp $source_networks'/NetAnalyzer.rb_0000/loci2phen_coords.txt' $outf'/raw_sors'
awk '{print $1 "\t" $2 "\t" $3 "\t" (($3 - $2)) "\t" $4 "\t" $5 "\t" $6}' $outf'/raw_sors' > $outf'/twopats_validated_sors'
rm $outf'/raw_sors'

# ETL over all GO files
topGO2clusterProfiler.R -p $results_source -f 'enrichment_table' -o $outf'/go_enrichments' -P $th_enrichments

# Concat all HPO-Loci_Genes networks
grep -h tag $source_networks'/merge_locis_hpos_genes.R_0000/rdm_3layers_networks/'* | head -n 1 > $outf'/full_networks'
cat $source_networks'/merge_locis_hpos_genes.R_0000/rdm_3layers_networks/'* | grep -v tag >> $outf'/full_networks' 

# Concat all KEGG enrichments
grep -h Model $results_source/*'/enrich_by_onto.R_0000/kegg_enrichment' | head -n 1 > $outf'/kegg_enrichments'
cat $results_source/*'/enrich_by_onto.R_0000/kegg_enrichment' | grep -v Model >> $outf'/kegg_enrichments'

# Concat all REACTOME enrichments
grep -h Model $results_source/*'/enrich_by_onto.R_0000/reactome_enrichment' | head -n 1 > $outf'/reactome_enrichments_raw'
cat $results_source/*'/enrich_by_onto.R_0000/reactome_enrichment' | grep -v Model >> $outf'/reactome_enrichments_raw'
genesName2Entrez.R -e $outf'/reactome_enrichments_raw' -g '/' -i 10 -o $outf'/reactome_enrichments'

# Obtain affected genes after HyI filter
awk '{if( $1 == "dec" ) print $4}' $outf'/full_networks' | tr ':' '\n' | sort -u | tail -n +2 > $outf'/hyi_genes' 


############################################# 
# Calculate ratios (diference between RDM signal and RAW signals)
calc_ratios.R -e $outf'/go_enrichments' -s GO -r rdm_models_info -m dec -i "1,2,4,9" -o $outf'/deviance_ratios'
calc_ratios.R -e $outf'/kegg_enrichments' -s KEGG -r rdm_models_info -m dec -a -o $outf'/deviance_ratios'
calc_ratios.R -e $outf'/reactome_enrichments' -s REACTOME -r rdm_models_info -m dec -a -o $outf'/deviance_ratios'



### CALCULATE METADATA AND OTHER PROCESSED DATA

# Obtain metadata
metadata_calculator.R -l $outf'/full_networks' -g $outf'/go_enrichments' -k $outf'/kegg_enrichments' -r $outf'/reactome_enrichments' -o $outf'/meta'

## Apply special filtering 
tail -n +2 rdm_models_info | cut -f 2 > rdm_info
# KEGG
get_significant_terms.rb -r dec -t $th_rdm_filters -m rdm_info -n $rdm_models -c Model,Set,ID -i $outf'/kegg_enrichments' > $outf'/filtered_kegg_signal'
# GO 
get_significant_terms.rb -r dec -t $th_rdm_filters -m rdm_info -n $rdm_models -c Model,Term,GO -i $outf'/go_enrichments' > $outf'/filtered_go_signal'
# REACTOME 
get_significant_terms.rb -r dec -t $th_rdm_filters -m rdm_info -n $rdm_models -c Model,Set,ID -i $outf'/reactome_enrichments' > $outf'/filtered_reactome_signal'
# rm rdm_info


# Filter GO parentals into final filtered signal
remove_parentals.R -f $outf'/filtered_go_signal' -i "3,1,2" -o $outf'/filtered_go_cleaned'


# GO: Obtain extended metadata info for filtered signal
# node_freqs.R -i $original_rels -t 1 -o $outf'/hpo_freqs'
node_freqs.R -i $original_rels -t 1 -o $outf'/hpo_freqs' -a
extend_filtered_signal.R -s $outf'/filtered_go_cleaned' -r $outf'/deviance_ratios' -T GO -f $outf'/hpo_freqs' -o $outf'/filtered_go_extended' -t $th_enrichments -m rdm_l_
extend_filtered_signal.R -s $outf'/filtered_go_cleaned' -r $outf'/deviance_ratios' -T GO -f $outf'/hpo_freqs' -o $outf'/filtered_go_extended' -t $th_enrichments -m rdm_vl_ -a
extend_filtered_signal.R -s $outf'/filtered_go_cleaned' -r $outf'/deviance_ratios' -T GO -f $outf'/hpo_freqs' -o $outf'/filtered_go_extended' -t $th_enrichments -m rdm_g_ -a

# Include Real model into GO
awk '{ if ($1 == "dec") {print $1 "\t" $2 "\t" $4 "\t" $9} }' $outf'/go_enrichments' > $outf'/go_enrichments_fFormat'
tail -n +2 $outf/go_enrichments_fFormat > $outf/go_enrichments_filteredFormat
rm $outf/go_enrichments_fFormat
extend_filtered_signal.R -s $outf'/go_enrichments_filteredFormat' -T GO -f $outf'/hpo_freqs' -o $outf'/filtered_go_extended' -t $th_enrichments -m dec -a


# KEGG: Obtain extended metadata info for filtered signal
extend_filtered_signal.R -s $outf'/filtered_kegg_signal' -r $outf'/deviance_ratios' -T KEGG -f $outf'/hpo_freqs' -o $outf'/filtered_kegg_extended' -t $th_enrichments -m rdm_l_
extend_filtered_signal.R -s $outf'/filtered_kegg_signal' -r $outf'/deviance_ratios' -T KEGG -f $outf'/hpo_freqs' -o $outf'/filtered_kegg_extended' -t $th_enrichments -m rdm_vl_ -a
extend_filtered_signal.R -s $outf'/filtered_kegg_signal' -r $outf'/deviance_ratios' -T KEGG -f $outf'/hpo_freqs' -o $outf'/filtered_kegg_extended' -t $th_enrichments -m rdm_g_ -a

# REACTOME: Obtain extended metadata info for filtered signal
extend_filtered_signal.R -s $outf'/filtered_reactome_signal' -r $outf'/deviance_ratios' -T REACTOME -f $outf'/hpo_freqs' -o $outf'/filtered_reactome_extended' -t $th_enrichments -m rdm_l_
extend_filtered_signal.R -s $outf'/filtered_reactome_signal' -r $outf'/deviance_ratios' -T REACTOME -f $outf'/hpo_freqs' -o $outf'/filtered_reactome_extended' -t $th_enrichments -m rdm_vl_ -a
extend_filtered_signal.R -s $outf'/filtered_reactome_signal' -r $outf'/deviance_ratios' -T REACTOME -f $outf'/hpo_freqs' -o $outf'/filtered_reactome_extended' -t $th_enrichments -m rdm_g_ -a

# Obtain unified filtered signals
unify_filtered_signal.R -e $outf'/filtered_go_extended' -t $diff_rdms -o $outf'/filtered_go_unified' -r rdm_models_info
unify_filtered_signal.R -e $outf'/filtered_kegg_extended' -t $diff_rdms -o $outf'/filtered_kegg_unified' -r rdm_models_info
unify_filtered_signal.R -e $outf'/filtered_reactome_extended' -t $diff_rdms -o $outf'/filtered_reactome_unified' -r rdm_models_info


# Filter final set by redundant Pheno-FunSys pairs by ontology hierarchy
filter_pairs_by_most_specific.R -e $outf'/filtered_go_unified' -o $outf'/filtered_go_unified_trimmed'


# # Calculate PACKAGE EFFECT
sor_package_effect.R -e $outf'/meta_go' -n $outf'/full_networks' -t 1e-03 -p 2 -o $outf'/pack_go'
sor_package_effect.R -e $outf'/meta_kegg' -n $outf'/full_networks' -t 1e-03 -p 2 -o $outf'/pack_kegg'
sor_package_effect.R -e $outf'/meta_reac' -n $outf'/full_networks' -t 1e-03 -p 2 -o $outf'/pack_reac'


# Obtain experiment results patient stats and triplets
patients_final_results.R -t $raw_original_rels -n $outf'/full_networks' -e $outf'/go_enrichments' -i "1,2,4,10" -g ':' -f $outf'/filtered_go_unified' -o $outf'/patients_go_results' -P $outf'/patients_go_triplets'
patients_final_results.R -t $raw_original_rels -n $outf'/full_networks' -e $outf'/kegg_enrichments' -i "1,2,3,10" -g '/' -f $outf'/filtered_kegg_unified' -o $outf'/patients_kegg_results' -P $outf'/patients_kegg_triplets'
patients_final_results.R -t $raw_original_rels -n $outf'/full_networks' -e $outf'/reactome_enrichments' -i "1,2,3,10" -g '/' -f $outf'/filtered_reactome_unified' -o $outf'/patients_reactome_results' -P $outf'/patients_reactome_triplets'

# Obtain enrichment dictionaries
funsys_dict_generator.R -a -g $outf'/filtered_go_extended' -k $outf'/filtered_kegg_unified' -r $outf'/filtered_reactome_unified' -o $outf'/ures'


## MAIN REPORT
# Main report files
files_report=`echo -e "
>	external_data/hpo_db_phen2gene.txt,    # External data #
	$outf/full_networks,                 # Enrichments #
	$outf/source_cnvs,                   # Original CNVs #
>	$outf/go_enrichments,
>	$outf/reactome_enrichments,                
>	$outf/kegg_enrichments,
	$outf/meta_loci,                     # Metadata #
	$outf/meta_go,
	$outf/meta_kegg,
	$outf/meta_reac,
	rdm_models_info,                       # Models info #
>	$outf/filtered_go_cleaned,           # Filtered signal #
>	$outf/filtered_kegg_signal,
>	$outf/filtered_reactome_signal,
	$outf/filtered_go_extended,          # Special filtered signal #
	$outf/filtered_kegg_extended,
	$outf/filtered_reactome_extended,
	$outf/filtered_go_unified,
	$outf/filtered_kegg_unified,
	$outf/filtered_reactome_unified,                          
	$outf/patients_go_results,              # Patients results #
	$outf/patients_go_triplets,
>	experiment/patients_results,              # Patients results #
>	experiment/patients_triplets,
	$outf/deviance_ratios,               # Ratios info #
>	$outf/deviance_ratios_u,
	$outf/real_sors,                    # Real SORs #
	$outf/twopats_validated_sors,       # Source SORs #
	$outf/source_net,                   # Real net #
	$outf/net,                          # Used net #
>	$outf/pack_go,                      # Package effect metadata #
	$outf/ures_go_dict                  # GO frequencies #	
"| awk "{gsub(\"#.*#\",\"\"); print}" | 
awk "{gsub(\">.*(,|$)\",\"\"); print}" |   
tr -d [:space:] |
awk "{gsub(\",$\",\"\"); print}"`


files_report_headers=`echo -e "
>	f,  # hpo_db_phen2gene #
	t,  # full_networks #
	f,  # CNVs #
>	t,  # go_enrichments #
>	t,  # kegg_enrichments #
>	t,  # reactome_enrichments #
	t,  # meta_loci #
	t,  # meta_go #
	t,  # meta_kegg #
	t,  # meta_reac #
	t,  # rdm_models_info #
>	f,  # filtered_go_cleaned #
>	f,  # filtered_kegg_signal #
>	f,  # filtered_reactome_signal # 
	t,  # filtered_go_extended #
	t,  # filtered_kegg_extended #
	t,  # filtered_reactome_extended #
	t,  # filtered_go_unified #
	t,  # filtered_kegg_unified #
	t,  # filtered_reactome_unified #
	t,  # patients_go_results #
	t,  # patients_go_triplets #
	t,  # deviance_ratios #
>	t,  # deviance_ratios_u #
	f,  # real_sors #
	f,  # twopats_validated_sors #
>	t,  # pack_go #
	f,  # source_net #
	f,  # net #
	t   # ures_go_dict #
"| awk "{gsub(\"#.*#\",\"\"); print}" | 
awk "{gsub(\">.*(,|$)\",\"\"); print}" |   
tr -d [:space:] |
awk "{gsub(\",$\",\"\"); print}"` 

# # Generate report
create_report.R -t reports/general_report.Rmd -o $outf'/generalreport.html' -d $files_report -H $files_report_headers 





#############################################################  HPO REPORT
# To activate CONFIDENTIAL MODE remove ">" from patients_ids_confidential entrances bellow and include patient's IDs into named file
# WARNING!: confidential file name must be "conf_pats_ids" without extension
# touch conf_pats_ids
files_report=`echo -e "
	$outf/filtered_go_unified,
	$outf/ures_go_dict,
	$outf/diseases_profiles,
	$outf/patients_go_triplets,
	conf_pats_ids
"| awk "{gsub(\"#.*#\",\"\"); print}" | 
awk "{gsub(\">.*(,|$)\",\"\"); print}" |   
tr -d [:space:] |
awk "{gsub(\",$\",\"\"); print}"`


files_report_headers=`echo -e "
	t, # filtered_go_unified #
	t, # ures_*_dict #
	f, # diseases_profiles #
	t,  # patients_triplets #
	f  # Confidential patients #
"| awk "{gsub(\"#.*#\",\"\"); print}" | 
awk "{gsub(\">.*(,|$)\",\"\"); print}" |   
tr -d [:space:] |
awk "{gsub(\",$\",\"\"); print}"` 


# # Obtain diseases info
diseases_profiles.R -f external_data/phenotype_annotation.tab -i "1,2,3,6,5" -o $outf'/diseases_profiles'


# Generate HPO_GO report
create_report.R -t reports/phenotype_report.Rmd -o $outf'/pheno_go_report.html' -d $files_report -H $files_report_headers

# Generate HPO_KEGG report
files_report=`echo $files_report | awk "{gsub(\"go\",\"kegg\"); print}"`
create_report.R -t reports/phenotype_report.Rmd -o $outf'/pheno_kegg_report.html' -d $files_report -H $files_report_headers

# Generate HPO_KEGG report
files_report=`echo $files_report | awk "{gsub(\"kegg\",\"reactome\"); print}"`
create_report.R -t reports/phenotype_report.Rmd -o $outf'/pheno_reactome_report.html' -d $files_report -H $files_report_headers





#############################################################  PATIENTS REPORT

files_report=`echo -e "
	$outf/patients_go_triplets,
	$outf/patients_go_results,
	$outf/ures_go_dict,
	conf_pats_ids,
>	$outf/hpo_freqs
"| awk "{gsub(\"#.*#\",\"\"); print}" | 
awk "{gsub(\">.*(,|$)\",\"\"); print}" |   
tr -d [:space:] |
awk "{gsub(\",$\",\"\"); print}"`


files_report_headers=`echo -e "
	t, # patients_triplets #
	t, # patients_results #
	t, # ures_dict #
	f, # confidential patients #
>	f  # hpo_freqs #
"| awk "{gsub(\"#.*#\",\"\"); print}" | 
awk "{gsub(\">.*(,|$)\",\"\"); print}" |   
tr -d [:space:] |
awk "{gsub(\",$\",\"\"); print}"` 

# Generate GO Patient report
create_report.R -t reports/patient_results_report.Rmd -o $outf'/patient_go_report.html' -d $files_report -H $files_report_headers

# Generate KEGG Patient report
files_report=`echo $files_report | awk "{gsub(\"go\",\"kegg\"); print}"`
create_report.R -t reports/patient_results_report.Rmd -o $outf'/patient_kegg_report.html' -d $files_report -H $files_report_headers

# Generate Reactome Patient report
files_report=`echo $files_report | awk "{gsub(\"kegg\",\"reactome\"); print}"`
create_report.R -t reports/patient_results_report.Rmd -o $outf'/patient_reactome_report.html' -d $files_report -H $files_report_headers