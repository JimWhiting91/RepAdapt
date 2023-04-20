# RepAdapt GEA Analysis
This repository includes 4 R scripts that developed as command-line software.

---

## Metadata and Setup
In the same directory as the genotype data, there must be a csv file including sampling information called `sampling_data.csv`
This looks like:
```
| Taxon_name    | Pop_ID    | Sample_name   | Country   | Pool_Size | Lat   | Long  | Alt   | Author    | Data_type |
| :---:         | :---:     | :---:         | :---:     | :---:     | :---: | :---: | :---: | :---:     | :---:     |
| Arabidopsis lyrata | Grovudalen | Alyr_Norway_Grovudalen_T3-6 | Norway | 62.44 | 8.9 | 900 | Savolainen | WGS |
| Arabidopsis lyrata | Grovudalen | Alyr_Norway_Grovudalen_T3-5 | Norway | 62.44 | 8.9 | 900 | Savolainen | WGS |
| Arabidopsis lyrata | Grovudalen | Alyr_Norway_Grovudalen_T3-3	| Norway | 62.44 | 8.9 | 900 | Savolainen | WGS |
```

In order to perform GEA over climate change variables, these should also be calculated prior to running and deposited in the same directory as the genotype data and named `climate_change_env.txt`. These can be calculated using the script `R/00_setup/calculating_climatechange_clines.R`.

We also need to make sure that the genome that the SNP data has been mapped onto has a liftover map which includes the OrthoFinder2 gene ID that each gene is assigned.
Having made proteomes from each genome and each gff annotation, these liftovers are made using `R/00_setup/create_liftover_mapping_files_for_OF2_ids_and_gffs.R`.
These are essential in order to convert SNP results to gene WZA results and eventually map these back to Orthogroups.

---

## GEA Pipeline
### Step 1 - Quantify and summarise population structure
Script = `quantify_population_structure.R`

The purpose of this script is to quantify, summarise and visualise population structure and spatial FST among populations.
This information is useful for generally debugging issues with SNP datasets and highlighting any outlier populations.
The associations between Fst and physical distance among sites also gives an approximate sense of isolation by distance, although this is also covered by the later RDA analyses.

The summary of the workflow is that SNP data is randomly downsampled to a user-defined extent to remove linkage.
PCAs are then performed to describe pop structure using SNPRelate.
We also calculate FST among 'populations' defined as individuals with the same Lat/Long co-ordinates.
FST is compared against physical distance to approximate isolation by distance.

#### Parameters
```
--vcf = path to the vcf or snptable
--sub_SNP = number of SNPs to downsample to (default of 10,000)
--pool = boolean TRUE/FALSE as to whether the dataset is poolsequencing data or not
--dataset_dir = path to the working directory for all outputs and temporary files
--n_cores = number of computing cores for parallel operations
```

### Step 2 - Genotype x Environment Associations (GEA)
Script = `perform_GEA.R`

This script performs GEA analysis at a per-SNP level across all 19 BIOCLIM variables and the 2 climate change variables.
GEA in this case are Kendall's Tau correlations between per-site allele frequencies and per-site climatic variation.
The key output here is the per-SNP p-value of the correlation, which is later used to yield gene-level scores.
These are saved per-variable in the working directory, for e.g. `mean_temp_GEA.rds`

This script also produces a number of other outputs that are used elsewhere and are stored in the working directory, including:
 * Per-population allele frequencies (`pop_allele_frqs.rds`)
 * Per-population climate variation (`climate_cline.tsv`)

#### Parameters
```
--vcf = path to the vcf or snptable
--pool = boolean TRUE/FALSE as to whether the dataset is poolsequencing data or not
--dataset_dir = path to the working directory for all outputs and temporary files
--n_cores = number of computing cores for parallel operations
```

### Step 3 - Collate SNP-level GEA to gene-level scores
Script = `collate_GEA_to_gene_WZA_empirical_pvals.R`

This script takes the outputs from GEA, along with information on the boundaries of genes described in the annotation, and produces gene-level scores.

Gene-level scores are based on the weighted-Z analysis (WZA) (see - https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13768?af=R).
We also calculate as a comparison the 'Top Candidate' score (see - https://www.science.org/doi/10.1126/science.aaf7812).

The mapping between SNP co-ordinates and gene co-ordinates is done using data.table()

The script requires pointing towards where the gff annotation is, and assumes there is a liftover map with OrthoFinder2 ID named `*_proteome_to_OF_id_map.txt`.
This takes the form of:
```
| gene_ID     | OF_ID | seqid | start    | end      | gea_gene               |
|-------------|-------|-------|----------|----------|------------------------|
| Aa_G10.h1   | 9723  | chr3  | 16039538 | 16041808 | chr3:16039538-16041808 |
| Aa_G100.h1  | 9714  | chr3  | 15925024 | 15928903 | chr3:15925024-15928903 |
| Aa_G1000.h1 | 9013  | chr3  | 9651994  | 9653039  | chr3:9651994-9653039   |
```
The following script is therefore able to use the gene boundaries to group together SNPs in the same gene, and the OF_ID to map these genes back to orthogroups.

Note: Orthofinder renames genes as `GenomeN_ProteinN`, and the OF_ID column here represents the `ProteinN`, i.e. the order each gene appears in the proteome given to Orthofinder. Later on, we have to assign the `GenomeN` based on whichever set of orthogroups we are using, i.e. the `GenomeN` will depend on which genomes were actually included in the Orthofinder run.

Note: The script also adds flanking regions to all of these genes, which are user-defined parameters (default of 500bp)

#### Parameters
```
--vcf = path to the vcf or snptable
--gff = path to the gff (more specifically we want the directory with the orthofinder liftover map)
--tc_threshold = The binomial quantile cutoff used for the top-candidate approach (default of 99%)
--snp_per_gene = Either quantile or integer for the number of SNPs to downsample above (e.g. a value of 0.75 will downsample genes with SNPs more than the 75% quantile) (default of 20)
--gene_flank_size = Size of flanking region added to gene co-ordinates (default of 500bp)
--dataset_dir = path to the working directory for all outputs and temporary files
--n_cores = number of computing cores for parallel operations
```

### Step 4 - Variance partitioning using redundancy analysis (RDA)
This step of the analysis is optional, i.e. it doesn't produce anything that's strictly necessary for the WZA > PicMin > Repeatability analyses.
It uses RDA to partition genetic variation into spatial and climate predictors.

Models are built using per-site allele frequencies, so the script uses the allele frequencies from step 2.

Datasets are downsampled to a random subset of N (user-defined) SNPs, with preference given to genes with minimal missing data.
Only SNPs with at most 10% missing data across all populations are considered for models.
SNPs with some missing data (<10%) are imputed using the global mean allele frequency.

#### Parameters
```
--vcf = path to the vcf or snptable
--dataset_dir = path to the working directory for all outputs and temporary files
--snp_downsample = number of SNPs to downsample to (default of 10,000)
--n_cores = number of computing cores for parallel operations
```