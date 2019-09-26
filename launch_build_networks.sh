#! /usr/bin/env bash

# @author Pedro Seoane Zonjic, Fernando Moreno Jabato

# Initialize DEPENDENCIES
#	> Autoflow
source ~soft_bio_267/initializes/init_autoflow

# Add necessary scripts
current_dir=`pwd`
export PATH=$current_dir'/scripts':$PATH


#establish the variables we need in the workflow
mkdir external_data

wget 'ftp://ftp.ncbi.nih.gov/genomes/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/GFF/ref_GRCh37.p13_top_level.gff3.gz' -O external_data/genome.gz
gunzip -d  external_data/genome.gz
wget 'http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt' -O external_data/hpo_db.txt
tail -n +2 external_data/hpo_db.txt | cut -f 1,3 > external_data/hpo_db_phen2gene.txt  
wget http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/misc/phenotype_annotation.tab -O external_data/phenotype_annotation.tab
wget -O external_data/hp.obo http://purl.obolibrary.org/obo/hp.obo --no-check-certificate


# Patients info
mkdir processed_data

# Convert YOUR DATA format to our processing format (if it's necessary)
# [1] : Patient [2]: Chr [3]: Start [4]: End [5]: HPO_Name/Code
# ## ## ## ## ## ## 

# Create HPOs dictionary
parse_hpo_file.rb external_data/hp.obo > processed_data/hpo2name.txt


# About $hpo_enrichment
#	'-r' => no enrichment, '' => enrichment
# Prepare variables
variables=`echo -e "
	\\$patients_file=$current_dir'/processed_data/patient_data.txt',
	\\$hpo_dict=$current_dir'/processed_data/hpo2name.txt',
	\\$genome_annotation=$current_dir'/external_data/genome',
	\\$number_of_random_models=2,
	\\$association_thresold=2,
	\\$official_hpo=$current_dir/external_data/hpo_db_phen2gene.txt,
	\\$hpo_enrichment='',
	\\$hpo_ontology=$current_dir/external_data/hp.obo
" | tr -d [:space:]`

#Lauch autoflow 
AutoFlow -w build_networks.af -o $SCRATCH'/build_Networks' -V $variables -m 30gb $1
