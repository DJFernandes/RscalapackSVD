% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scala.svd.finish.R
\name{scala.svd.finish}
\alias{scala.svd.finish}
\title{Converts ScaLAPACK SVD output to R data}
\usage{
scala.svd.finish(proc.dir, BIGarr = TRUE, write.function = NULL,
  extension = "", results.dir)
}
\arguments{
\item{proc.dir}{directory where all the ScaLAPACK SVD computations were run.}

\item{BIGarr}{LOGICAL. If true, output is folders "Left_Singular_Vectors" and "Right_Singular_Vectors_transpose", inside which are many files. Each file corresponds to columns of the Left Singular and Right Singular (transpose) matrices. If false, Left_Singular_Vectors.RData and Right_Singular_Vectors_transpose.RData files are made containing the Left Singular and Right Singular (transpose) matrices.}

\item{write.function}{CHARACTER. String of function name that writes vectors (needed if BIGarr=TRUE) to files. Function can only have two arguments. First is for the vector, Second is for the filename.}

\item{extension}{CHARACTER Extension of filenames created}

\item{results.dir}{Directory to output results}
}
\description{
Reads directory with successful ScaLAPACK SVD run and creates directory with results:
   1) Singular Values
   2) Left Singular Vectors
   3) Right Singular Vectors transpose
}
\examples{
\dontrun{
#### Example to read directory with successful ScaLAPACK SVD run ("gene_svd")
####  and output results to directory ("gene_svd_results").
####  With BIGarr=TRUE

# Create Write Function
simple.write.function=function(vec,fl){
   cat(vec,sep='\\n',file=fl)
}

# call scala.svd.finish
scala.svd.finish(proc.dir="gene_svd",
                 BIGarr=TRUE,write.function="simple.write.function",
                 extension=".txt",results.dir="gene_svd_results")

#### With BIGarr=FALSE
# call scala.svd.finish
scala.svd.finish(proc.dir="gene_svd",BIGarr=FALSE,results.dir="gene_svd_results")
}
}

