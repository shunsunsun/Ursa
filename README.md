![Ursa logo](https://user-images.githubusercontent.com/5945741/165857896-912bfe07-f290-483c-bb96-d5ff21db1ab6.png)

# Ursa: an automated multi-omics package for single-cell analysis

__Ursa__ is an R package consisting of seven single-cell omics automated analysis workflows. One-liner command for each omics to run a full post-quantification analysis for the omics.

Seven single-cell omics include:

1. scRNA-sequencing
2. scATAC-sequencing
3. scImmune profiling
4. Spatial transcriptomics
5. scCNV
6. CyTOF
7. Flow cytometry

## Installation

Ursa can be installed in R via the command:
```sh
devtools::install_github('eudoraleer/Ursa')
```

## Running single-cell analysis with Ursa
##### (1) scRNA-sequencing
```sh
scRNASEQPip(project_name = 'My_scRNASeq', input_dir = '/home/input/', output_dir = '/home/output/', pheno_file = ‘/home/input/meta.txt’)
```
##### (2) scATAC-sequencing*
```sh
scATACPip(project_name = 'My_scATAC', input_dir = '/home/input/', output_dir = '/ home/output/', pheno_file = ‘/home/input/meta.txt’)
```
##### (3) scImmune profiling
```sh
scImmunePip(project_name = 'My_scImmune', input_dir = '/home/input/', output_dir = '/home/output/', pheno_file = ‘/home/input/meta.txt’)
```
##### (4) Spatial transcriptomics
```sh
SpatialPip(project_name = 'My_Spatial', input_dir = '/home/input/', output_dir = '/ home/output/', pheno_file = ‘/home/input/meta.txt’)
```
##### (5) scCNV
```sh
scCNVPip(project_name = 'My_scCNV', input_dir = '/home/input/', output_dir = '/ home/output/', pheno_file = ‘/home/input/meta.txt’)
```
##### (6) CyTOF
```sh
CyTOFPip(project_name = 'My_CyTOF', input_dir = '/home/input/', output_dir = '/ home/output/', pheno_file = ‘/home/input/meta.txt’)
```
##### (7) Flow
```sh
FlowPip(project_name = 'My_Flow', input_dir = '/home/input/', output_dir = '/home/ output/', pheno_file = ‘/home/input/meta.txt’)
```

###### *Strong recommend running this workflow on a computer with memory >=16GB due to the large size of this omics data
