<a href="https://yeamanlab.weebly.com"><img src="https://img.shields.io/badge/yeaman-lab-blue"></a>
<a href="https://github.com/yeamanlab/template/LICENSE"><img src="https://img.shields.io/badge/license-MIT-green" alt="License Badge"></a>

### Author Contact Information
james.whiting@ucalgary.ca
sam.yeaman@ucalgary.ca

### Usage and license information
If you use or are inspired by code in this repository please cite the following work or contact me about how to cite. Please also see [license information](LICENSE).

Whiting et al. (20XX) *in prep* [doi_link]()

---

# RepAdapt

The population genomics of convergent adaptation to climate across plant species. The aims of this project are to explore the repeatability of genes associated with various facets of climate adaptation, and assess the contingencies and sources of variation driving differences in repeatabilty among diverse species.


---

# Environment setup

The scripts in this repository were written for a SLURM-based HPC and a local 48-CPU server. The `bin/repadapt_env.yaml` file can be used to setup the environment.

The following can be used to setup the necessary directory structure for directories not included in the repo:
```sh
# Set up directories and subdirectories
mkdir data data/VCFs data/reference_genomes metadata figs scripts/logs outputs
```

Each dataset is comprised of a single VCF, these are available at FigShare (doi.XXX) and should be placed into the `data/VCFs` subdirectory.

Similarly, links to reference genomes used are available in `metadata/fastq_accessions_and_doi.csv`. Each reference genome should have its own subdirectory within `data/reference_genomes`, that includes the relevant annotation in `.gff` format.

---

# Metadata

Several scripts depend on a single metadata file formatted as below. Individuals are grouped into populations based on `Lat` and `Long` columns, and datasets are named according to {Taxon_name}_{Author}_{Data_type}. The metadata csv used is available at `metadata/XXX.csv`.

|Taxon_name             |Sample_name        |Lat  |Long  |Author|Data_type |VCF                                                                                                 |
|-----------------------|-------------------|-----|------|------|----------|----------------------------------------------------------------------------------------------------|
|Amaranthus tuberculatus|Atuberculatus_15_10|42.62|-82.49|Wright|Individual|/lu213/james.whiting/RepAdapt/data/VCFs/09_Atuberculatus_Wright/Atuberculatus_full_concatened.vcf.gz|
|Amaranthus tuberculatus|Atuberculatus_15_11|42.62|-82.49|Wright|Individual|/lu213/james.whiting/RepAdapt/data/VCFs/09_Atuberculatus_Wright/Atuberculatus_full_concatened.vcf.gz|
|Amaranthus tuberculatus|Atuberculatus_15_3 |42.62|-82.49|Wright|Individual|/lu213/james.whiting/RepAdapt/data/VCFs/09_Atuberculatus_Wright/Atuberculatus_full_concatened.vcf.gz|
|Amaranthus tuberculatus|Atuberculatus_15_5 |42.62|-82.49|Wright|Individual|/lu213/james.whiting/RepAdapt/data/VCFs/09_Atuberculatus_Wright/Atuberculatus_full_concatened.vcf.gz|

---

# Usage 

The pipeline is split into several steps:

### SNP Calling

The SNP calling is described in detail in its own README, available at XXX.

### Genotype x Environment Associations (GEA)

GEA's are split into 3 steps: 
 * Quantification of population structure and dataset summary (`R/quantify_population_structure.R`)
 * GEA (`R/perform_GEA.R`)
 * Collation of per-SNP GEA to gene-level summaries (`R/collate_GEA_to_gene_WZA.R`)
 
These three scripts can be run using the wrapper: `scripts/GEA_pipeline_runner.sh`

### Orthology

TBC

### Quantifying contingencies

The following scripts assess predicted contingencies for repeatability:
 * Comparison of climates (`analyse_cline_similarity.R`)
 * Gene network analysis (XXX)
 * TBC

# License
[LICENSE](LICENSE) - our license information

# Additional notes/scripts
The following scripts are also included that are used for various bits and bobs.

