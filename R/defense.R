setup.defence <- function(       input,
                                 BIGarr,
                                 read.function,
                                 matdims,
                                 blockdims,
                                 procdims,
                                 proc.dir,
                                 walltime ) {
        if (is.null(input)) {stop("Must supply Input as either files (with bigarr=TRUE) or matrix")}
        if (BIGarr) { 
         if ( !(is.character(input) & is.null(dim(input))) )  {
             stop("If BIGarr=TRUE, input must be 1D vector of files ")}
         if (is.null(read.function)) { stop("If BIGarr=TRUE, You must have a function (as a string) to read the input files")}
         if ( length(input)!=matdims[2] )  {
             stop("number of files must equal number of columns (matdims[2])")}
        }
        if (!BIGarr) { 
         if ( !(is.numeric(input) & (length(dim(input))==2)) ) {
             stop("If BIGarr=FALSE, input must be 2D matrix of numbers ")}
         if ( prod(dim(input)==matdims)==0 )  {
             stop("matdims must be the dimensions of input array")}
        }
        
	if (blockdims[1]*procdims[1]>matdims[1]) { stop(paste("blockdims[1]*procdims[1]",blockdims[1]*procdims[1],"is greater than number of rows",matdims[1])) }
	if (blockdims[2]*procdims[2]>matdims[2]) { stop(paste("blockdims[2]*procdims[2]",blockdims[2]*procdims[2],"is greater than number of columns",matdims[2])) }
	#walltime must not be less than 30 minutes
	if (sum(as.numeric(strsplit(walltime,split=":")[[1]])*c(3600,60,1))<1800) {
		stop("walltime must NOT be less than 00:30:00")
		}
}

finish.defence <- function(BIGarr,write.function) {
    if (BIGarr) {
    if (is.null(write.function)) {stop("IF BIGarr is true, You must have a function (as a string) to write the output files")}
    }
}
