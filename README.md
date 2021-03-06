# PhenFun
Workflow to obtain Phenotype-(GO/KEGG/Reactome) relationships starting at a Phenotype-Patient-Mutation tripartite network.

## Installation
1. Clone this repository using 'git clone' and special option '--recurse-submodules' to obtain referenced scripts from our general lab repository (https://github.com/seoanezonjic/sys_bio_lab_scripts.git). **Code:** `git clone https://github.com/fmjabato/PhenFun.git --recurse-submodules`
2. Install R and R packages dependencies stored at `r_dependencies.txt` file. We recommend R version 3.6 or greater.
3. Install ruby. We recommend to use the RVM manager: RVM
4. Install AutoFlow. Code: `gem install autoflow`.
5. Make PATH accesible all the installed software.
6. Make PATH accesible the folder scripts.

## Usage
These workflow have three main steps. First one consists on **compute Phenotype-Patient-Mutation network** and necessary random models, including a quality process to remove not significant relationships. Second step performs a **functional analysis** over all Phenotype-Mutations relationships obtained. Third and last step is used to join all generated results and **generate reports** with visual information.

To launch these steps you must:

- Build tripartite (Phenotype-Patient-Genotype) networks, randomize and remove not significant relationships:
	1. Configure `launch_build_networks.sh` script to obtain your original and random networks in the correct format. This script will download necessary files, which makes some variables optional. Configurable variables are:
		- `patients_file`: input network file
		- `hpo_dict`: (OPTIONAL) HPO dictionary. Will be automatically downloaded and generated by this script.
		- `genome_annotation`: (OPTIONAL) Genome to be used. Genome GRCh37 will be automatically downloaded and generated by this script.
		- `number_of_random_models`: number of random models, of each type, to be generated. Default: 2.
		- `association_thresold`: threshold applied over HyI significance value between tripartite links. Default: 2.
		- `official_hpo`: (OPTIONAL) Human Phenotype Ontology official Phenotype-Gene file. Will be automatically downloaded by this script.
		- `hpo_enrichment`: flag to include parental phenotypes into patient profiles. Use `-r` to avoid enrichment or set empty to enrich. 
		- `hpo_ontology`: (OPTIONAL) Human Phenotype Ontology OBO file. Will be automatically downloaded by this script.
		- `regions_blacklist`: (OPTIONAL) A file with regions which is used to filter patients file and tripartite network.
		- `regions_filter`: (OPTIONAL) type of filter to be applied if regions blacklist file is provided. Allowed: remove (r) or cut in fragments (c) if overlaps with a blacklist region.
	2. Optionally you can configure `build_networks.af` to perform special experiments.
	3. Execute `launch_build_network.sh` script. 
- Perform functional analysis:
	1. Optionally, change `working_nets` generated file if you don't want to execute step 2, for all networks.
	2. Configure `launch_analyse_networks.sh` script. Configurable variables are:
		- `networks_source`: folder with generated networks. By default: *BuildStepFolder/ln_0000/working_netowrks*
		- `n_frags`: number of fragmets to split network. Used to speed up analysis. Default: 4. 
		- `clprof_enrichments`: terms sets to be used, separated by colon (:). Allowed: `go`, `kegg`, `reactome`.
		- `pval_thresold`: to be applied. Default: 0.001.
	3. Optionally, configure `analyse_networks.af` to perform special experiments.
	4. Execute `launch_analyse_networks.sh`.
- Generate result reports:
	1. Configure `get_comparative_results.sh`. Configurable variables are:
		- `original_file`: build step input file.
		- `results_source`: folder with analyse step results.
		- `source_networks`: folder with build step results.
		- `original_rels`: (optional) original calculated tripartite network. Default: source_networks/get_network_nodes.rb_0000/net.txt.
		- `raw_original_rels`: (optional) original calculated tripartite network without enrichment. Default source_networks/get_network_nodes.rb_0000/net_notE.txt
		- `outf`: (optional) folder name where results will be stored. Default: results.
		- `th_enrichments`: (optional) functional analysis threshold to be applied. Must be equal or greater than used at functional analysis step. Default: 1e-03
		- `th_rdm_filters`: (optional) threshold used to filter functional enrichments using random probability. Random probability is calculated using termInstancesAmongRDModels/numRDModels. Default : 0.01
		- `rdm_models`: check value to identify if any problem has occur during simulation. Must be equal to number of rdm models generated at build step.
		- `diff_rdms`: (optional) number of random models which must validate a Phenotype-FunctionalTerm tuple to be validated. Default: 3 
	2. Execute `get_comparative_results.sh`.

## License
This software is under MIT License.

Copyright (c) 2020 Fernando Moreno Jabato and collaborators.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
