#' @title Initializes Calls to ScaLAPACK SVD
#' @description
#' Creates directory with everything needed to run scalapack
#'    1) Fortran program with call to scalapack
#'    2) Parameters needed for fortran program
#'    3) Block Cyclic Distribution of input array
#'    4) Folders for output
#'    5) Script to compile (NOTE: designed to work on SCINET, but should be easily modifiable to your preferred cluster)
#'    6) Script to execute (NOTE: same note as 5)
#'
#' @param input 2D numeric matrix (if BIGarr=FALSE) or vector of filenames (if BIGarr=TRUE), where filenames contain numeric vectors corresponding to columns of the input 2D matrix
#' @param BIGarr LOGICAL. If FALSE, input must be a 2D numeric matrix. If TRUE, input must be a vector of filenames.
#' @param read.function CHARACTER. String of function name that reads input filenames (needed if BIGarr=TRUE). Function can only have one argument and that MUST be the filename.
#' @param matdims dimensions of input matrix, for which SVD will be computed
#' @param blockdims dimensions of blocks for block-cyclic distribution (see http://netlib.org/scalapack/slug/node75.html for details on block-cyclic distribution)
#' @param procdims dimensions of processor array (Example procdims=c(8,7) means processor array has 8 rows and 7 columns, thus 56 processors in all).
#' @param proc.dir directory where all the ScaLAPACK SVD computations will be run.
#' @param walltime Max time for jobs on a PBS queuing system. You can ignore this is you have a different queuing system.
#' @param progbar progress bar when reading files (BIGarr=TRUE). Requires 'pbapply' package.
#'
#' @examples
#' \dontrun{
#' #### Example where the input array is huge (BIGarr=TRUE)
#' # Read csv with filenames
#' # Make read.function for the filenames. 
#' #    Read function requires functions from package 'ABIgeneRMINC'.
#' #    First mask the data, then interpolate and scale (mean 0, sd 1), and output the resulting vector
#' # Call scala.svd.setup to set up directory for ScaLAPACK SVD
#' 
#' # Read CSV
#' library(ABIgeneRMINC)
#' df=read.csv("/projects/egerek/matthijs/2015-07-Allen-Brain/Allen_Gene_Expression/gene_ids_stats2.csv")
#' files=as.character(subset(df,cutting_plane=='coronal')$raw_data)
#'
#' # Create Function That Reads Filename
#' mask=read.raw.gene("/projects/egerek/matthijs/2015-07-Allen-Brain/Allen_Gene_Expression/labels/allen_grid_labels.raw",labels=TRUE)>0
#' read.gene.int.scale=function(filename){
#'      genevec=read.raw.gene(filename)
#'      intgene=interpolate.gene(genevec,mask)
#'      sgene=scale(intgene[mask],center=TRUE,scale=TRUE)
#'      return(scale(intgene[mask],center=TRUE,scale=TRUE))
#' }
#' 
#' # call scala.svd.setup
#' scala.svd.setup(files,read.function="read.gene.int.scale",matdims=c(62529,4345),blockdims=c(100,100),procdims=c(8,8),proc.dir="gene_svd",walltime="01:00:00")
#'
#' 
#' #### Example where the input array is small (BIGarr=FALSE)
#' #Read files into columns of input matrix
#' inptmatrix=do.call('cbind',lapply(files,read.gene.int.scale))
#' 
#' # call scala.svd.setup
#' scala.svd.setup(inptmatrix,BIGarr=FALSE,matdims=dim(inptmatrix),blockdims=c(100,100),procdims=c(8,8),proc.dir="gene_svd",walltime="01:00:00")
#' }
#' @export


scala.svd.setup <- function(       input,
                                   BIGarr=TRUE,
                                   read.function=NULL,
                                   matdims,
                                   blockdims,
                                   procdims,
                                   proc.dir,
                                   walltime="01:00:00",
                                   progbar=TRUE) {
        
        # Check if arguments are valid
       if (is(tryCatch(
         setup.defence(input,BIGarr,read.function,matdims,blockdims,procdims,proc.dir,walltime)
           ,error = function(e) e),"error")) {stop("Setup Failed")}

	#make directories needed (and copy FORTRAN source code)
        make_setup_directories(proc.dir)
        file.copy(system.file('extdata/SVD.f',package="RscalapackSVD"),proc.dir)

        #distribut input data
        if (BIGarr) { 
          block.cyc.distribute.big(input,read.function,matdims,blockdims,procdims,proc.dir,progbar)
        } else { 
          block.cyc.distribute.small(input,matdims,blockdims,procdims,proc.dir)
        }

	#make parameter file
	parameterfile=file.path(proc.dir,"SVDparameters.dat")
        make_parameter_file(parameterfile,
                            matdims[1],matdims[2],
                            blockdims[1],blockdims[2],
                            procdims[1],procdims[2],
                            "'input'","'singularvalues.dat'",
                            "'left_singular_vectors'","'right_singular_vectors_transpose'")

	#make qsub script
	qsubscript=file.path(proc.dir,"qsub_script.sh")
        make_qsub_script(qsubscript,prod(procdims),walltime)

        #make compile and execute script
	execscript=file.path(proc.dir,"compile.sh")
        make_compile_script(execscript)
        
        #print instructions
        print_instructions(proc.dir)
}
