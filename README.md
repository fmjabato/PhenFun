# PhenFun
Workflow to obtain Phenotype-(GO/KEGG/Reactome) relationships starting at a Phenotype-Patient-Mutation tripartite network.

## Usage
1. Clone this repository using 'git clone' and special option '--recurse-submodules' to obtain referenced scripts from our general lab repository (https://github.com/seoanezonjic/sys_bio_lab_scripts.git)
2. Configure and execute `launch_build_networks.sh` script to obtain your original and random networks in the correct format. Configure also `build_networks.af` if it's necessary.
3. Change `working_nets` generated file if you don't want to execute step 4, for all networks.
4. Configure and execute `launch_analyse_networks.sh` script to obtain functional enrichment of all networks. Configure also `analyse_networks.af` if it's necessary.
5. Conigure and execute `get_comparative_results.sh` script to obtain final results table and reports

## Important
* Remember give execution permission to scripts
* Remember clone the repository as Usage 1 indicates
