make.vecmask <- function(vecPROCi,vecMATsize,vecBLOCKsize,vecPROCsize) {
    
 if (vecPROCi>=vecPROCsize) {stop("vecPROCi must be less than or equal to vecPROCsize. Remember vecPROCi is the processor number starting from 0 (ex {0,1}) and vecPROCsize is the number of processors (ex {2}).")}

 if (vecBLOCKsize*vecPROCsize>vecMATsize) {stop("vecBLOCKsize*vecPROCsize must not be greater than vecMATsize")}

 # Starting (ROW/Column) Index processor (ROW/Column)
 SRTprocVECidx=1+(vecPROCi)*vecBLOCKsize
 # number of column blocks
 NUMvecBL=floor((vecMATsize-SRTprocVECidx)/(vecBLOCKsize*vecPROCsize))
 # Start Index of (ROW/Column) Blocks
 SRTblVECidx=SRTprocVECidx+(0:NUMvecBL)*vecBLOCKsize*vecPROCsize
 # End Index of (ROW/Column) Blocks
 ENDblVECidx=SRTblVECidx+vecBLOCKsize-1
 # If End Index of (ROW/Column) Blocks is larger than matrix, make the end index the size in that dimension
 if (ENDblVECidx[length(ENDblVECidx)]>vecMATsize) {
  ENDblVECidx[length(ENDblVECidx)]=vecMATsize}
 # Make (ROW/Column) Mask
 mask=rep(FALSE,vecMATsize)
 for (blk in 1:length(SRTblVECidx)) {mask[SRTblVECidx[blk]:ENDblVECidx[blk]]=TRUE}
 return(mask)

}

block.cyc.distribute.small <- function(
                                   input,
                                   matdims,
                                   blockdims,
                                   procdims,
                                   proc.dir) {
    #Cycle through columns and rows
    for (j in 0:(procdims[2]-1)) { 
     print(paste("setting processor column",j+1,"of",procdims[2]))
     # make column mask (which column in matrix belong to this processor column)
     colmask=make.vecmask(j,matdims[2],blockdims[2],procdims[2])
    for (i in 0:(procdims[1]-1)) {
     # make row mask (which rows in matrix belong to this processor row)
     rowmask=make.vecmask(i,matdims[1],blockdims[1],procdims[1])
     blockmat=matrix(NA,sum(rowmask),sum(colmask))
     blockmat=input[rowmask,colmask]
     cat(as.vector(blockmat),sep='\n',
       file=paste0(proc.dir,"/input/p_",sprintf("%02d",i),"-",sprintf("%02d",j)))
     }}
}

block.cyc.distribute.big <- function(
                                   colfiles,
                                   read.function,
                                   matdims,
                                   blockdims,
                                   procdims,
                                   proc.dir,
                                   progbar=TRUE) {
    #Cycle through columns
    for (j in 0:(procdims[2]-1)) {
    print(paste("setting processor column",j+1,"of",procdims[2]))
    # make column mask (which column in matrix belong to this processor column)
    colmask=make.vecmask(j,matdims[2],blockdims[2],procdims[2])
    # Subset files in this processor column
    cfiles=colfiles[colmask]
    # Read Files using a read function
    if (progbar) { colblock=do.call('cbind',pbapply::pblapply(cfiles,read.function)) 
     } else {      colblock=do.call('cbind',lapply(cfiles,read.function)) }
    # Check
    ncheck=nrow(colblock)
    if (ncheck!=matdims[1]) {
     stop(paste("number of rows",ncheck,"not equal to matdims[1]",matdims[1]))}

    # Cycle through Rows
    for (i in 0:(procdims[1]-1)) {
    # make row mask (which rows in matrix belong to this processor row)
     rowmask=make.vecmask(i,matdims[1],blockdims[1],procdims[1])
    # Create block this processor will process
      blockmat=matrix(NA,sum(rowmask),sum(colmask))
    # Take information and put in block
      blockmat=colblock[rowmask,]
    # Write information to file
      cat(as.vector(blockmat),sep='\n',
        file=paste0(proc.dir,"/input/p_",sprintf("%02d",i),"-",sprintf("%02d",j)))
   }}
}

block.cyc.aggregate.small <-function( proc.files.dir,
                                      aggregate.output.file,
                                      aggregatedims,
                                      blockdims,
                                      procdims) {
    #Initialize matrix to aggregate data
     agg=matrix(NA,aggregatedims[1],aggregatedims[2])
    #Cycle through columns and rows
    for (j in 0:(procdims[2]-1)) { 
     # make column mask (which column in matrix belong to this processor column)
     colmask=make.vecmask(j,aggregatedims[2],blockdims[2],procdims[2])
    for (i in 0:(procdims[1]-1)) {
     # make row mask (which rows in matrix belong to this processor row)
     rowmask=make.vecmask(i,aggregatedims[1],blockdims[1],procdims[1])
     #read data
     toread=file.path(proc.files.dir,paste0("p_",sprintf("%02d", i),"-",sprintf("%02d", j)))
     print(paste("reading",toread))
     blockmat=scan(toread)
     dim(blockmat)=c(sum(rowmask),sum(colmask))
     agg[rowmask,colmask]=blockmat
    }} 
    save(agg,file=aggregate.output.file)
}

block.cyc.aggregate.big <-function( proc.files.dir,
                                    aggregate.output.dir,
                                    aggregate.output.pre.ext=c('mode','.file'),
                                    write.function,
                                    aggregatedims,
                                    blockdims,
                                    procdims) {
    #Cycle through columns 
    PCcounter=0
    for (j in 0:(procdims[2]-1)) { 
     # make column mask (which column in matrix belong to this processor column)
     colmask=make.vecmask(j,aggregatedims[2],blockdims[2],procdims[2])
     # Initialize matrix to store all the data in this processor column
     colblock=matrix(NA,aggregatedims[1],sum(colmask))
     #column number of aggregated matrix each column in processor column belongs to
     PCnumbers=(1:aggregatedims[2])[colmask]
    #Cycle through rows
    for (i in 0:(procdims[1]-1)) {
     # make row mask (which rows in matrix belong to this processor row)
     rowmask=make.vecmask(i,aggregatedims[1],blockdims[1],procdims[1])
     #read data
     toread=file.path(proc.files.dir,paste0("p_",sprintf("%02d", i),"-",sprintf("%02d", j)))
     print(paste("reading",toread))
     blockmat=scan(toread)
     dim(blockmat)=c(sum(rowmask),sum(colmask))
     colblock[rowmask,]=blockmat
    }
    #Write Column by column
    for (j2 in 1:sum(colmask)) { 
        outfile=file.path(aggregate.output.dir,
                   paste0(aggregate.output.pre.ext[1],sprintf("%05d", PCnumbers[j2]),aggregate.output.pre.ext[2]))
        do.call(write.function,list(colblock[,j2],outfile)) }
    } 
}



