make_parameter_file=function(outfile,M,N,mb,nb,pr,pc,indir,s,lsv,rsvt) {
	write(file=outfile,append=FALSE,
              paste0(M,"\t\t\t\t\t\t","number of rows"))
	write(file=outfile,append=TRUE,
              paste0(N,"\t\t\t\t\t\t","number of columns"))
	write(file=outfile,append=TRUE,
              paste0(mb,"\t\t\t\t\t\t","number of rows in blocks"))
	write(file=outfile,append=TRUE,
              paste0(nb,"\t\t\t\t\t\t","number of columns in blocks"))
	write(file=outfile,append=TRUE,
              paste0(pr,"\t\t\t\t\t\t","number of processor rows"))
	write(file=outfile,append=TRUE,
              paste0(pc,"\t\t\t\t\t\t","number of processor columns"))
	write(file=outfile,append=TRUE,
              paste0(indir,"\t\t\t\t\t\t","input dir files to read"))
	write(file=outfile,append=TRUE,
              paste0(s,"\t\t\t\t","output file for singular values"))
	write(file=outfile,append=TRUE,
              paste0(lsv,"\t\t\t\t","output dir for left singular vectors"))
	write(file=outfile,append=TRUE,
              paste0(rsvt,"\t\t","output dir for right singular vectors(transpose)"))
}

extract_parameter_file=function(paramfile) {
       parsethis=readLines(paramfile)

       Ma=as.integer(stringi::stri_extract_first_regex(parsethis[1],"[0-9]+"))
       Na=as.integer(stringi::stri_extract_first_regex(parsethis[2],"[0-9]+"))
       mba=as.integer(stringi::stri_extract_first_regex(parsethis[3],"[0-9]+"))
       nba=as.integer(stringi::stri_extract_first_regex(parsethis[4],"[0-9]+"))
       pr=as.integer(stringi::stri_extract_first_regex(parsethis[5],"[0-9]+"))
       pc=as.integer(stringi::stri_extract_first_regex(parsethis[6],"[0-9]+")) 
       return(c(Ma,Na,mba,nba,pr,pc))
}

make_qsub_script=function(outfile,num.processors,walltime) {
	write(file=outfile,append=FALSE,
             "#!/bin/bash")	
	write(file=outfile,append=TRUE,
             paste0("#PBS -l nodes=",ceiling(num.processors/8),":ppn=",8))
	write(file=outfile,append=TRUE,
             paste0("#PBS -l walltime=",walltime))
	write(file=outfile,append=TRUE,
             "#PBS -N SVD")
	write(file=outfile,append=TRUE,
             "module purge")
	write(file=outfile,append=TRUE,
             "module load intel/15.0.2 intelmpi/5.0.3.048 scalapack/2.0.1-intel-intelmpi")
	write(file=outfile,append=TRUE,
             "cd !!!.MUST.BE.CHANGED.!!!")
	write(file=outfile,append=TRUE,
             paste("mpirun -np ",num.processors," ./SVD"))
}

make_compile_script=function(outfile) {

       	write(file=outfile,append=FALSE,
              "#!/bin/bash")
        #load modules
	write(file=outfile,append=TRUE,
              "module purge;module load intel/15.0.2 intelmpi/5.0.3.048 scalapack/2.0.1-intel-intelmpi")
	#compile SVD program
	write(file=outfile,append=TRUE,
              "mpif77 -openmp -I/scinet/gpc/intel/ics/composer_xe_2015.2.164/mkl/lib/intel64/include SVD.f -L/scinet/gpc/intel/ics/composer_xe_2015.2.164/mkl/lib/intel64/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm -o SVD")
	#modify script to contain execution directory
	write(file=outfile,append=TRUE,
              'la="s+!!!.MUST.BE.CHANGED.!!!+$( pwd )+g"')
	write(file=outfile,append=TRUE,
              'sed -i "$la" qsub_script.sh')
        write('echo "Compilation complete"',file=outfile,append=TRUE)
        write('echo "If you are happy with it,"',file=outfile,append=TRUE)
        write('echo "Run the command below to run the computation"',file=outfile,append=TRUE)
	write('echo qsub ./qsub_script.sh',file=outfile,append=TRUE)
}

make_setup_directories=function(proc.dir) {
	#make proc.dir
	dir.create(proc.dir)

	#make directories that will store the input data
        dir.create(file.path(proc.dir,"input"))

	#make directories that will store the return data
        dir.create(file.path(proc.dir,"left_singular_vectors"))
        dir.create(file.path(proc.dir,"right_singular_vectors_transpose"))
}

print_instructions=function(proc.dir) {
        print("Setup is done. To run, follow the steps below:")
	print(paste("1) copy",proc.dir,"to scinet (or whatever supercomputer centre you are affiliated with)"))
	print(paste("2) cd into",proc.dir))
        print("You may have to change things in the compile script and qsub script, based on the cluster you are running on.")
        print("They are designed to work on scinet")
	print("3) execute the compile script")
	print("4) qsub ./qsub_script.sh")
        print(paste("5) Once job finishes, copy",proc.dir,"back to local (if you want) and run scala.svd.finish function"))
}

make_finish_directories=function(results.dir,BIGarr) {
        dir.create(results.dir)
        if (BIGarr) {
         dir.create(file.path(results.dir,"Left_Singular_Vectors")) 
         dir.create(file.path(results.dir,"Right_Singular_Vectors_transpose"))
        }
}
