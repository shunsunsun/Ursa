% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Merak.R
\name{FlowPip}
\alias{FlowPip}
\title{Run a post-quantification 10X technology-format Flow Cytometry pipeline}
\usage{
FlowPip(
  project_name = "Ursa_Flow",
  input_dir = "./",
  output_dir = "./",
  pheno_file
)
}
\arguments{
\item{project_name}{Project name. 'Ursa_Flow' by default.}

\item{input_dir}{Directory to all input files. Current working directory by
default.}

\item{output_dir}{Output directory. Current working directory by default.
A new folder with the given project name with time stamp as suffix will be
created under the specified output directory.}

\item{pheno_file}{Meta data file directory. Accept only .csv/.txt format
files.}
}
\description{
This function will run a bioinformatics analysis of post-quantification Flow
pipeline. Supports multiple samples analysis.
}
