![Ursa logo](https://user-images.githubusercontent.com/5945741/165857896-912bfe07-f290-483c-bb96-d5ff21db1ab6.png)

# Ursa: an automated multi-omics package for single-cell analysis

__Ursa__ is an R package consisting of seven single-cell omics automated analysis workflows. One-liner command for each omics to run a full post-quantification analysis for the omics.

Six single-cell (sc) omics and one bulk omics include:

1. scRNA-sequencing (sc)
2. scATAC-sequencing (sc)
3. scImmune profiling (sc)
4. scCNV (sc)
5. CyTOF (sc)
6. Flow cytometry (sc)
7. Spatial transcriptomics (bulk)

## Installation

Ursa can be installed in R via the command:
```sh
devtools::install_github('eudoraleer/Ursa')
```
Please download the example meta files from the Github directory (__https://github.com/eudoraleer/Ursa/tree/master/DB/Examples/Meta_Files__) with respect to the omics you will be running.

## Running single-cell analysis with Ursa
##### (1) scRNA-sequencing* (download example dataset from [__10X__]())
```sh
scRNASEQPip(project_name = 'My_scRNASeq', input_dir = '/home/input/', output_dir = '/home/output/', pheno_file = ‘/home/input/meta.txt’)
```
##### (2) scATAC-sequencing* (download example dataset from [__10X__](https://www.10xgenomics.com/resources/datasets/10k-human-pbmcs-atac-v2-chromium-controller-2-standard))
###### For this omics, running this workflow on a computer with memory >=16GB is recommended due to large input file size
###### The following input files are needed in the input directory before running the analysis:
###### 1) fragment file and its index file (i.e. fragments_filtered.tsv.gz, fragments_filtered.tsv.gz.tbi)
###### 2) single cell file (i.e. singlecell.csv)
###### 3) filtered peak matrix .h5 file (i.e. filtered_peak_bc_matrix.h5)

```sh
scATACPip(project_name = 'My_scATAC', input_dir = '/home/input/', output_dir = '/ home/output/', pheno_file = ‘/home/input/meta.txt’)
```
##### (3) scImmune profiling* (download example dataset from [__10X__]())
```sh
scImmunePip(project_name = 'My_scImmune', input_dir = '/home/input/', output_dir = '/home/output/', pheno_file = ‘/home/input/meta.txt’)
```
##### (4) scCNV* (download example dataset from [__10X__]())
```sh
scCNVPip(project_name = 'My_scCNV', input_dir = '/home/input/', output_dir = '/ home/output/', pheno_file = ‘/home/input/meta.txt’)
```
##### (5) CyTOF (download example dataset from [__Nowicka, M., et al. (2017)__](http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow/PBMC8_fcs_files.zip))
```sh
CyTOFPip(project_name = 'My_CyTOF', input_dir = '/home/input/', output_dir = '/ home/output/', pheno_file = ‘/home/input/meta.txt’)
```
##### (6) Flow (download example dataset from [__Dillon Hammill,2021__](https://github.com/DillonHammill/CytoExploreRData/tree/master/inst/extdata/Activation))
```sh
FlowPip(project_name = 'My_Flow', input_dir = '/home/input/', output_dir = '/home/ output/', pheno_file = ‘/home/input/meta.txt’)
```
##### (7) Spatial transcriptomics* (download example dataset from [__10X__]())
```sh
SpatialPip(project_name = 'My_Spatial', input_dir = '/home/input/', output_dir = '/ home/output/', pheno_file = ‘/home/input/meta.txt’)
```

###### *Registration is needed for downloading the data for the first time on 10X website. Subsequent download would no longer require any registration.
