% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Alioth.R
\name{scRNASEQPip}
\alias{scRNASEQPip}
\title{Run a post-quantification 10X technology-format scRNASEQ pipeline}
\usage{
scRNASEQPip(
  project_name = "Ursa_scRNASEQ",
  input_dir = "./",
  output_dir = "./",
  pheno_file = "study_meta.csv",
  num_genes_lower_bound = 200,
  num_genes_upper_bound = 25000,
  mito_cutoff = 5,
  integration_method = "Seurat",
  cc_regression = 0,
  pc_selection_method = "none",
  num_pcs = 30,
  to_impute = "NO",
  find_doublet = "NO",
  run_unbias_vis = "NO"
)
}
\arguments{
\item{project_name}{Project name. 'Ursa_scRNASEQ' by default.}

\item{input_dir}{Directory to all input files. Default directory
is the current working directory.}

\item{output_dir}{Output directory. Default directory is the
current working directory.
A new folder with the given project name with time stamp as suffix
will be created under the specified output directory.}

\item{pheno_file}{Meta data file directory. Accept only .csv/.txt
format files.}

\item{num_genes_lower_bound}{Lower bound to select for cells with
number of genes more than this value. Default is 200.}

\item{num_genes_upper_bound}{Upper bound to select for cells with
number of genes less than or equals to this value Default is 25000.}

\item{mito_cutoff}{Threshold percentage to select for cells with
percentage of mitochondria genes less than this value.
This is remove cells which are debris or artefacts and doublets.
Default is 5.}

\item{integration_method}{Integration method.
Accepts 'Seurat' or "Harmony' integration. Default is Seurat.}

\item{cc_regression}{Cell cycle regression method.
Accepts 0 for no regression,
1 for regression with two-phases, 2 for regression with phase
difference (See Seurat). Default is 0.}

\item{pc_selection_method}{Select a method to run an automatic selection
of the number of principal components (PCs) or Harmonys (if harmony was
chosen for integration). Methods include 'none', 'all', 'piecewise linear
model', 'first derivative', 'second derivative', 'preceding residual',
'perpendicular line', and 'k-means clustering'. Default is set to 'none',
which will ignore this option and use the stated number of PCs in the
num_pcs option. If 'all' is used, all methods will be assessed and the
final PC number will be determined via the mean of the output of all
methods. For more information, please visit:
https://github.com/haotian-zhuang/findPC
(Zhuang et al., Bioinformatics, 2022)}

\item{num_pcs}{Number of PCs or Harmonys (if harmony was chosen for
integration) to be used for the analysis. Default to 30. Please
state a reasonable number which is no more than the total number
of cells submitted to avoid running into errors. This option will
only be considered if option pc_selection_method is set to 'none'.}

\item{to_impute}{Pass 'YES' to perform imputation on individual
samples or 'NO' to skip imputation process. Default is set to 'NO'.}

\item{find_doublet}{Pass 'YES' to perform doublets removal using
 DoubletFinder (Christopher S. M., Cell Systems, 2019). Default
is set to 'NO'. If set to 'YES', doublets will be removed using
the true doublet rate (TDR).}

\item{run_unbias_vis}{Plot additional unbias visualizations for top
features using SCUBI (Wenpin H. & Zhicheng J., Cell Reports Methods,
2021) for dimensionally reduced plots such as UMAP. Pass 'YES' to
plot, default to 'NO'.}
}
\description{
This function will run a bioinformatics analysis of post-quantification scRNASEQ
pipeline. Supports multiple samples analysis.
}
