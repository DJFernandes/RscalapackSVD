#' @title Converts ScaLAPACK SVD output to R data
#' @description Reads directory with successful ScaLAPACK SVD run and creates directory with results:
#'    1) Singular Values
#'    2) Left Singular Vectors
#'    3) Right Singular Vectors transpose
#'
#' @param proc.dir directory where all the ScaLAPACK SVD computations were run.
#' @param BIGarr LOGICAL. If true, output is folders "Left_Singular_Vectors" and "Right_Singular_Vectors_transpose", inside which are many files. Each file corresponds to columns of the Left Singular and Right Singular (transpose) matrices. If false, Left_Singular_Vectors.RData and Right_Singular_Vectors_transpose.RData files are made containing the Left Singular and Right Singular (transpose) matrices.
#' @param write.function CHARACTER. String of function name that writes vectors (needed if BIGarr=TRUE) to files. Function can only have two arguments. First is for the vector, Second is for the filename.
#' @param extension CHARACTER Extension of filenames created
#' @param results.dir Directory to output results
#'
#' @examples
#' \dontrun{
#' #### Example to read directory with successful ScaLAPACK SVD run ("gene_svd")
#' ####  and output results to directory ("gene_svd_results").
#' ####  With BIGarr=TRUE
#' 
#' # Create Write Function
#' simple.write.function=function(vec,fl){
#'    cat(vec,sep='\n',file=fl)
#' }
#' 
#' # call scala.svd.finish
#' scala.svd.finish(proc.dir="gene_svd",
#'                  BIGarr=TRUE,write.function="simple.write.function",
#'                  extension=".txt",results.dir="gene_svd_results")
#'
#' #### With BIGarr=FALSE
#' # call scala.svd.finish
#' scala.svd.finish(proc.dir="gene_svd",BIGarr=FALSE,results.dir="gene_svd_results")
#' }
#' @export

scala.svd.finish <- function(       proc.dir,
                                    BIGarr=TRUE,
                                    write.function=NULL,
                                    extension="",
                                    results.dir) {

       # Check if arguments are valid
       if (is(tryCatch(
         finish.defence(BIGarr,write.function)
           ,error = function(e) e),"error")) {stop("Setup Failed")}

       # Make Directories for results
       make_finish_directories(results.dir,BIGarr)

       #parse the parameter file
       parameterfile=file.path(proc.dir,"SVDparameters.dat")
       params=extract_parameter_file(parameterfile)
       matdims=params[1:2]
       blockdims=params[3:4]
       procdims=params[5:6]

       #aggregate Left Singular Vectors block cyclic data
       LSVinDIR=file.path(proc.dir,"left_singular_vectors")
        if (BIGarr) { 
          LSVoutDIR=file.path(results.dir,"Left_Singular_Vectors")
          block.cyc.aggregate.big(LSVinDIR,LSVoutDIR,
                                   c("mode",extension),write.function,
                                   c(matdims[1],min(matdims)),blockdims,procdims)
        } else { 
          LSVoutFLE=file.path(results.dir,"Left_Singular_Vectors.RData")
          block.cyc.aggregate.small(LSVinDIR,LSVoutFLE,c(matdims[1],min(matdims)),blockdims,procdims)
        }

        #aggregate Right Singular Vectors block cyclic data
       RSVinDIR=file.path(proc.dir,"right_singular_vectors_transpose")
        if (BIGarr) { 
          RSVoutDIR=file.path(results.dir,"Right_Singular_Vectors_transpose")
          block.cyc.aggregate.big(RSVinDIR,RSVoutDIR,
                                   c("loading",extension),write.function,
                                   c(min(matdims),matdims[2]),blockdims,procdims)
        } else { 
          RSVoutFLE=file.path(results.dir,"Right_Singular_Vectors_transpose.RData")
          block.cyc.aggregate.small(RSVinDIR,RSVoutFLE,c(min(matdims),matdims[2]),blockdims,procdims)
        }

        #copy singular values into results
        sv=file.copy(file.path(proc.dir,"singularvalues.dat"),file.path(results.dir,"Singular_values.txt"),overwrite = TRUE)
}

