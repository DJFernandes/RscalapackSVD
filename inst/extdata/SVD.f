	Program SVD

!	PDGESVD computes the singular value decomposition (SVD) of an
!	M-by-N matrix A, optionally computing the left and/or right
!	singular vectors. The SVD is written as
!
!		A = U * S * VT
!
!	where SIGMA is an M-by-N matrix which is zero except for its
!	min(M,N) diagonal elements, U is an M-by-M orthogonal matrix, and
!	V is an N-by-N orthogonal matrix. The diagonal elements of SIGMA
!	are the singular values of A and the columns of U and V are the
!	corresponding right and left singular vectors, respectively. The
!	singular values are returned in array S in decreasing order and
!	only the first min(M,N) columns of U and rows of VT = V**T are
!	computed.

	USE mpi
	IMPLICIT NONE
!	integers for loops
	INTEGER		I

!	M and N are the number of rows and columns respectively in A
!	MB and NB are the number of rows and columns of the submatrixes of A
!	Z is min(M,N) (i.e. the number of singular values computed)
!	NBU is the number of rows of the submatrixes of U
!	MBVT  is the number of columns of the submatrixes of VT
	integer :: M,N,MB,NB,Z,NBU,MBVT

!	myXrows, myXcols is the size of the local subset of array X [X is A,U,VT]
	integer :: myArows, myAcols
	integer :: myUrows, myUcols
	integer :: myVTrows, myVTcols

!	submatrixes
	DOUBLE PRECISION, dimension(:), allocatable :: A,U,VT

!	work defines a workspace for processes to communicate with one another
!	worksize is the size of the work array
	DOUBLE PRECISION, dimension(:), allocatable :: work
	integer :: worksize

!	singular values
	DOUBLE PRECISION, dimension(:), allocatable :: S


	integer, external :: numroc   ! blacs routine

! 	me is an ID for the process
!	procs is the number of processes
!	icontxt specifies the BLACS context handle identifying the 
!	 created process grid.
!	prow, pcol are the number of rows and columns respectively of the process grid
!	myrow, mycol are the coordinates of the process grid associated with the current process
	integer :: me, procs, icontxt, prow, pcol, myrow, mycol

!	info is the scalapack return value
	integer :: info

!	scalapack array descriptors for array A,U,VT
	integer, dimension(9)   :: DESCA,DESCU,DESCVT

!	input/output data filename
	character(len=128) :: inarr,outleft,outright,outsingular,
     $				paramfile,oltpr,ortpr,inp

!	mpirank is an ID for the process
!	ierr is the error code associated with MPI
	integer :: mpirank,ierr

!	pdims is the dimension of the process grid
!	 pdims=[prow,pcol]
!	dims is the dimension of A
!	 dims = [M,N]
!	distribs tells MPI to distribute CYCLIC when reading file from array
!	dargs is the dimension of submatrices of A
!	 dargs = [MB,NB]
	integer, dimension(2) :: pdims, dims, distribs, dargs

!	various variables and arrays needed for MPI
	integer :: infile, mpistatus(MPI_STATUS_SIZE)
	integer :: darray,locsize,nelements
	integer(kind=MPI_ADDRESS_KIND) :: lb, locextent
	integer(kind=MPI_OFFSET_KIND) :: disp
	integer :: nargs

!	to check if parameter file exists
	LOGICAL :: THERE


!_________________________________________________________________________

!MAIN PROGRAM

!check if parameter file exist
	paramfile="SVDparameters.dat"
	INQUIRE( FILE=paramfile, EXIST=THERE ) 
	IF ( .not.THERE ) THEN
	 print*, "parameter file", paramfile, " does not exist"
	 GO TO 10
	END IF

! open and read parameter file
	OPEN( 11, FILE=paramfile, STATUS='OLD' )
	 READ( 11, FMT = * ) M
	 READ( 11, FMT = * ) N
	 READ( 11, FMT = * ) MB
	 READ( 11, FMT = * ) NB
	 READ( 11, FMT = * ) prow
	 READ( 11, FMT = * ) pcol
	 READ( 11, FMT = * ) inarr
	 READ( 11, FMT = * ) outsingular
	 READ( 11, FMT = * ) outleft
	 READ( 11, FMT = * ) outright
	close(11) 

! define other parameters
	Z=MIN(M,N)
	pdims=[prow,pcol]
	dims = [M,N]

	if (MB > (M/prow)) then
	 MB = (M/prow)
	 print*, "MB > M/prow, changing A grid row length to", MB
	end if
	if (NB > (N/pcol)) then
	 NB = (N/pcol)
	 print*, "NB > (N/pcol), changing A grid col length to", NB
	end if
	if (NB > (Z/pcol)) then
	 NBU = (Z/pcol)
	 print*, "NB > (Z/pcol), changing U grid col length to", NBU
	else
	 NBU=NB
	end if
	if (MB > (Z/prow)) then
	 MBVT = (Z/prow)
	 print*, "MB > M/prow, changing VT grid row length to", MBVT
	else
	 MBVT=MB
	end if	

! Initialize MPI (for MPI-IO)

	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_COMM_WORLD,procs,ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD,mpirank,ierr)

! check if the number of elements in processor grid equals the number of processors
	if (mpirank==0) then
	 if (procs.ne.(prow*pcol)) then
	  print*, "number of processors called", procs 
	  print*, "must equal number of processors on the grid"
	  print*,  prow,"*",pcol,"=",prow*pcol
	  GO TO 10	  
	 end if
	end if
	

! Use MPI-IO to read matrix A and distributed block-cyclically,
	distribs = [MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC]
	dargs = [MB,NB]

	call MPI_Type_create_darray(procs, mpirank, 2, dims,
     $		distribs, dargs, pdims, 
     $		MPI_ORDER_FORTRAN, MPI_DOUBLE,
     $		darray, ierr)

	call MPI_Type_commit(darray,ierr)
	call MPI_Type_size(darray, locsize, ierr)

! locsize is the space in bytes
! nelements is the number of elements allowed (the division by 8 is because of double precision)
	nelements = locsize/8
	call MPI_Type_get_extent(darray, lb, locextent, ierr)

! Initialize local arrays    
	allocate(A(nelements))
	allocate(U(nelements))
	allocate(VT(nelements))
	allocate(S(Z))

! Initialize blacs processor grid
	call blacs_pinfo	(me,procs)
	call blacs_get		(-1, 0, icontxt)
	call blacs_gridinit	(icontxt, 'R', prow, pcol)
	call blacs_gridinfo	(icontxt, prow, pcol, myrow, mycol)

! regurgitate input
	
	if (mpirank==0) then
         WRITE( 6, FMT = 9999 )
     $			'SVD of Matrix A using Scalapack'//
     $			' subroutine PDGESVD'
         WRITE( 6, FMT = * )
         WRITE( 6, FMT = 9999 )
     $			'An explanation of the input/output'//
     $			' parameters follows'
         WRITE( 6, FMT = 9999 )
     $			'M/N         : The number of rows/colums'//
     $			' of Matrix A'
         WRITE( 6, FMT = 9999 )
     $			'MB/NB       : The number of rows/colums'//
     $			' of submatrices of A'
         WRITE( 6, FMT = 9999 )
     $			'prow/pcol   : The number of rows/colums'//
     $			' of processor grid'
         WRITE( 6, FMT = 9999 )
     $			'Input File  : Binary file where'//
     $			' matrix A is read from'
         WRITE( 6, FMT = 9999 )
     $			'SV File     : File where'//
     $			' singular values of A are written'
         WRITE( 6, FMT = 9999 )
     $			'L/RT SV File : Directory'//
     $			' where left/right(Transpose) singular'//
     $			'  vectors are written'

	WRITE( 6, FMT = 9998 ) 'M            ', M
	WRITE( 6, FMT = 9998 ) 'N            ', N
	WRITE( 6, FMT = 9998 ) 'MB           ', MB
	WRITE( 6, FMT = 9998 ) 'NB           ', NB
	WRITE( 6, FMT = 9998 ) 'prow         ', prow
	WRITE( 6, FMT = 9998 ) 'pcol         ', pcol
	WRITE( 6, FMT = 9997 ) 'Input File   ', trim(inarr)
	WRITE( 6, FMT = 9997 ) 'SV File      ', trim(outsingular)
	WRITE( 6, FMT = 9997 ) 'L SV File    ', trim(outleft)
	WRITE( 6, FMT = 9997 ) 'RT SV File   ', trim(outright)
	end if

 9999 FORMAT( A )	
 9998 FORMAT( 2X, A, '   :        ', I6 )
 9997 FORMAT( 2X, A, '   :        ', A )

! find the number of local rows and columns of matrices A,U,VT
! ONLY THE FIRST Z=min(M,N) COLUMNS OF U AND ROWS OF VT ARE COMPUTED
	myArows 	= numroc(M,MB, myrow, 0, prow)
	myAcols 	= numroc(N,NB, mycol, 0, pcol)
	myUrows 	= myArows
	myUcols 	= numroc(Z,NBU, mycol, 0, pcol)
	myVTrows	= numroc(Z,MBVT, myrow, 0, prow)
	myVTcols	= myAcols
	
! read in the data
	write(inp,'("p","_",I2.2,"-",I2.2)') myrow,mycol
	inp=trim(inarr)//'/'//trim(inp)
	OPEN( 11, FILE=trim(inp), STATUS='OLD' )
	  DO 15, I = 1,(myArows*myAcols)
	   READ(11,*) A(I)
15	  CONTINUE
	 close(11)	
	close(11)

! Prepare array descriptors for ScaLAPACK 
	CALL DESCINIT( DESCA, M, N, MB,NB, 0, 0, icontxt,
     $               MAX(1,myArows), info )
	CALL DESCINIT( DESCU, M, Z, MB,NBU, 0, 0, icontxt,
     $               MAX(1,myUrows), info )
	CALL DESCINIT( DESCVT, Z, N, MBVT,NB, 0, 0, icontxt,
     $               MAX(1,myVTrows), info )

	
! Call ScaLAPACK library routine

	allocate(work(1))
	work(1)  = -1.
	
	Call PDGESVD('V','V',M,N,A,1,1,DESCA,S,U,1,1,
     $			DESCU,VT,1,1,DESCVT,work,-1,info)

	worksize  = int(work(1))

	if (mpirank == 0)
     $		print *, "size of the work array is ", worksize

	deallocate(work)
	allocate(work(worksize))
	Call PDGESVD('V','V',M,N,A,1,1,DESCA,S,U,1,1,
     $			DESCU,VT,1,1,DESCVT,work,worksize,info)

	if (info /= 0) then
         print *, 'Error: info = ', info
	end if
	if (mpirank == 0) then
	 OPEN( 11, FILE=outsingular, STATUS='REPLACE' )
	  DO 12, I = 1,Z
	   WRITE(11,*) S(I)
12	  CONTINUE
	 close(11)         
	endif

	write(oltpr,'("p","_",I2.2,"-",I2.2)') myrow,mycol
	oltpr=trim(outleft)//'/'//trim(oltpr)
	OPEN( 11, FILE=trim(oltpr), STATUS='REPLACE' )
	  DO 13, I = 1,(myUrows*myUcols)
	   WRITE(11,*) U(I)
13	  CONTINUE
	 close(11)	
	close(11)

	write(ortpr,'("p","_",I2.2,"-",I2.2)') myrow,mycol
	ortpr=trim(outright)//'/'//trim(ortpr)
	OPEN( 11, FILE=trim(ortpr), STATUS='REPLACE' )
	  DO 14, I = 1,(myVTrows*myVTcols)
	   WRITE(11,*) VT(I)
14	  CONTINUE
	close(11)

! Deallocate the local arrays

	deallocate(A,U,VT)
	deallocate(work)

! End blacs for processors that are used

	call blacs_gridexit(icontxt)

! Print results

	deallocate(S)

	call MPI_Finalize(ierr)

10	CONTINUE
	STOP
	END


