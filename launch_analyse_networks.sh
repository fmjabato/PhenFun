#! /usr/bin/env bash

# @author Pedro Seoane Zonjic, Fernando Moreno Jabato


source ~soft_bio_267/initializes/init_autoflow
currDir=`pwd`
export PATH=$currDir'/scripts':$PATH


# SELECT NETWORK
networks_source=$SCRATCH/build_Networks/ln_0000/working_networks
mkdir $SCRATCH'/analysed_Networks'
ls $networks_source > working_nets


while read NETWORK
do
	variables=`echo -e "
		\\$working_net=$networks_source/$NETWORK,
		\\$n_frags=4,
		\\$clprof_enrichments=kegg:reactome,
		\\$pval_thresold=0.001
	" | tr -d [:space:]`
	AutoFlow -w analyse_networks.af -o $SCRATCH'/analysed_Networks/'$NETWORK -V $variables $1 -m 2gb -t '7-00:00:00'



done < working_nets
