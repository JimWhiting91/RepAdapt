<a href="https://yeamanlab.weebly.com"><img src="https://img.shields.io/badge/yeaman-lab-blue"></a>
<a href="https://github.com/yeamanlab/template/LICENSE"><img src="https://img.shields.io/badge/license-MIT-green" alt="License Badge"></a>

### Author Contact Information
james.whiting@ucalgary.ca

sam.yeaman@ucalgary.ca

### Usage and license information
If you use or are inspired by code in this repository please cite the following work or contact me about how to cite. Please also see [license information](LICENSE).

Whiting et al. (20XX) *in prep* [doi_link]()

---
![Global Mean Temp](./figs/repadapt_repo_fig.png?raw=true "Global Mean Temp")
# RepAdapt

The population genomics of repeated local adaptation to climate across plant species. The aims of this project are to explore the repeatability of genes associated with various facets of climate adaptation, and assess the contingencies and sources of variation driving differences in repeatabilty among diverse species.

---

# Environment setup

The scripts in this repository were written for a SLURM-based HPC and a local 48-CPU server.

Each dataset is comprised of a single VCF or SNPTabke, these are available at FigShare (doi.XXX) and should be placed into the `data/VCFs` subdirectory.

Similarly, links to reference genomes used are available in `metadata/fastq_accessions_and_doi.csv`. Each reference genome should have its own subdirectory within `data/reference_genomes`, that includes the relevant annotation in `.gff` format.

---

# Metadata
### Sampling Data Information
A metadata file is required for each dataset that is used to assign individuals to populations and provide geospatial information.
This must be formatted and included in the same directory as the VCF/SNPTable as `sampling_data.csv`:

| Taxon_name    | Pop_ID    | Sample_name   | Country   | Pool_Size | Lat   | Long  | Alt   | Author    | Data_type |
| :---:         | :---:     | :---:         | :---:     | :---:     | :---: | :---: | :---: | :---:     | :---:     |
| Arabidopsis lyrata | Grovudalen | Alyr_Norway_Grovudalen_T3-6 | Norway | 62.44 | 8.9 | 900 | Savolainen | WGS |
| Arabidopsis lyrata | Grovudalen | Alyr_Norway_Grovudalen_T3-5 | Norway | 62.44 | 8.9 | 900 | Savolainen | WGS |
| Arabidopsis lyrata | Grovudalen | Alyr_Norway_Grovudalen_T3-3	| Norway | 62.44 | 8.9 | 900 | Savolainen | WGS |

This metadata file contains the per-dataset information to build dataset names as `Taxon_name_Author_Data_type`.
The `Lat` `Long` columns are used to group individuals into populations to calculate site-level allele frequencies.

### Reference Genome and Annotation Information
A number of additional metadata files are described here:
 * `vcf_paths.txt` - List of relative paths to VCF files that can be provided to loop over for GEA analyses
 * `vcf_genome_gff_map.txt` - Table matching VCFs to genomes and annotations, for e.g.

| #VCF                                                                                          | GENOME                                                             | GFF                                                              | IS_POOL | GENOME_SPECIES | READY | RUN _OUTPUT | OUTPUT_DIR                               |
|-----------------------------------------------------------------------------------------------|--------------------------------------------------------------------|------------------------------------------------------------------|---------|----------------|-------|-------------|------------------------------------------|
| data/VCFs/02_Alyrata_Willi_pool/Alyr_full_concatened.filtered_NoMiss0.25_maf0.05.SNPTable.txt | data/reference_genomes/A.lyrata_reference/Alyr_v.1.0_genomic.fasta | data/reference_genomes/A.lyrata_reference/Alyr_v.1.0_genomic.gff | TRUE    | Alyrata        | YES   | YES         | Arabidopsis_lyrata_Willi_PoolSeq         |
| data/VCFs/03_Alyrata_Savolainen_wgs/savolainen_Alyrata_full_concatened_ScandScot_maf05.vcf.gz | data/reference_genomes/A.lyrata_reference/Alyr_v.1.0_genomic.fasta | data/reference_genomes/A.lyrata_reference/Alyr_v.1.0_genomic.gff | FALSE   | Alyrata        | YES   | YES         | Arabidopsis_lyrata_Savolainen_Individual |

---

# Usage 

The pipeline is split into several steps:
Details of individual scripts are given as READMEs within each directory.

### SNP Calling
There are two separate sets of scripts for SNP calling, depending on whether we're handling poolseq or wgs/capture data.

The mpileup scripts for wgs/capture data are found in `snp_calling_pipeline/mpileup_pipeline/`. A README there describes it in more detail.

The poolseq pipeline is maintained elsewhere and is available at `https://github.com/CoAdapTree/varscan_pipeline`

### Orthology
Using the reference genomes used for SNP-calling. Orthology is performed using OrthoFinder2.
This is split into two steps (although `run_orthofinder.sh` has loops to perform both tasks):
 * `scripts/prepare_proteome.sh` - Takes as input a reference genome and its associated gff annotation and exports a `.faa` proteome
 * `scripts/run_orthofinder.sh` - Runs OrthoFinder2 on a directory of proteomes

### Climate Data
Here we use global climate data from the worldclim database (https://www.worldclim.org/data/bioclim.html).
We use data collected at two resolutions: 2.5 and 10 minutes.

Current estimates for 2.5 min used in GEA analyses are pulled during the analysis.

Historical estimates to estimate site-level climate change were downloaded at source and stored in `data/worldclim/`

### General Setup
Prior to running GEAs etc. there are a few scripts to run to prepare some metadata and intermediate files
These scripts are found in `R/00_setup/`
 * Create liftover maps of genome co-ordinates to orthofinder gene IDs (`create_liftover_mapping_files_for_OF2_ids_and_gffs.R`)
 * Calculate climate change variables (`calculating_climatechange_clines.R`)
 * General processing of recombination maps where available, including smoothing (`process_recombination_maps.R`)
 * If required, filter SNP tables for maf and quality (`filter_snptable_add_maf.R`)

### Genotype x Environment Associations (GEA)
GEAs are split into 4 steps:
These scripts are found in `R/01_GEA_pipeline/`
 * Quantification of population structure and dataset summary (`quantify_population_structure.R`)
 * GEA through the weighted-Z analysis (WZA) (`perform_GEA.R`)
 * Collation of per-SNP GEA to gene-level summaries (`collate_GEA_to_gene_WZA.R`)
 * Variance partitioning of genetic variance by space and climate with RDA (`partition_variance_with_RDA.R`)
 
Each of these scripts is written as a piece of software taking cmd inputs for each dataset.
A loop for running them across all datasets is provided at: `scripts/GEA_pipeline_runner.sh`
More information is given on each of these scripts in `R/01_GEA_pipeline/README`

### Quantifying repeatability
To calculate orthogroup-level repeatability, this part is split into two ordered steps:
These scripts are found in `R/02_repeatability_calculating`
 * Process the per-gene WZA scores to per-orthogroup pvals (`process_GEA_WZA_to_corrected_OG_pvals_withOFcodes.R`)
 * Calculate orthogroup-level repeatability with PicMin (`calculate_orthogroup_GEA_repeatability_picmin.R`)
 * Also summarise the orthogroups and those tested (`Orthogroup_summary_stats.R`)

### Analysis of PicMin results
Following PicMin runs, outputs can be analysed using the scripts in: `R/03_picmin_results_analysis/`
 * General summary of PicMin results and plotting of Figure 2 (`summarise_picmin_convergence_results_fig2.R`)
 * STRING analysis of RAO interactions, GO enrichments and Figure 3 (`analyse_picmin_RAO_through_string_networks_fig3.R`)
 * Calculate Arabidopsis and Medicago coexpression network stats (`calculate_arabidopsis_coexpression_network_statistics.R`)
 * Analysis of pleiotropy enrichment in RAO and strongest PicMin orthogroups (`functional_pleiotropy_analysis_of_RAO.R`)
 * Analysis of associations between gene duplication and RAO and Figure 4 (`analysis_of_duplications_and_picmin_res_fig4.R`)
 * Analysis of niche-breadth and RDA results with PicMin (`niche_breadth_and_RDA_analysis_of_picmin_variation.R`)
 * Comparison between recombination rate and WZA variance (`recombination_rate_and_wza_association.R`)

### General Plots
Some assorted scripts for plots are also available in `R/04_plot_figures/`
 * Figure 1 - summary of study design and orthogroup assignment (`Figure1_plot.R`)
 * Agreement between orthogroup tree and timetree (`FigureSX_phylo_agreement.R`)
 * Network centrality statistics demo (`plot_network_centrality_demo.R`)

### Miscellaneous
There are two files with general functions that are used in various scripts:
 * `PicMin.R` - Contains all functions necessary for running PicMin analysis in vectorised form
 * `repadapt_functions.R` - Library of general functions used across scripts

# License
[LICENSE](LICENSE) - our license information

# Notes and updates
```
 * 24/03/03 - Repo updated with additional scripts and edits following NEE revisions
 * 23/04/21 - Repository updated with scripts up to final draft of manuscript
```

