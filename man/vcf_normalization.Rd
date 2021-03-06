% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf.R
\name{vcf_normalization}
\alias{vcf_normalization}
\title{Normalize VCF file to bi-allelic variants}
\usage{
vcf_normalization(vcf.file, ref.file, output.file,
  bcftools.exec = "bcftools", num.threads)
}
\arguments{
\item{vcf.file}{[\code{string}]\cr
The input VCF file path.}

\item{ref.file}{[\code{string}]\cr
A human reference genome \code{fasta} file to normalize indels against.}

\item{output.file}{[\code{string}]\cr
The output VCF file path.}

\item{bcftools.exec}{[\code{string}]\cr
Path of bcftools executable.}

\item{num.threads}{[\code{int}]\cr
Number of CPUs usable by bcftools
Default is determined by SLURM environment variables and at least 1.}
}
\value{
Captured system output as \code{character} vector.
}
\description{
Converting VCF files to plink format has never been easier. However, there are a few issues related to some intrinsic limitations of the plink format. The first is related to the fact that variants in a plink file are bi-allelic only, while variants in a VCF file can be multi-allelic. The second is related to an intrinsic limitation of plink which makes indel definitions ambiguous. Here is an example: is the following variant an insertion or a deletion compared to the GRCh37 reference?
}
\details{
20 31022441 A AG

There is no way to tell, as the plink format does not record this information.

Keeping this in mind, we are going to split mulit-allelic variants into bi-allelic ones, left-normalize indels, and assign unique idetifiers.
}
