# Introduction #

Below is an outline of the steps necessary to install the cufflink library for use in OpenFOAM-extend.  In general, if you are able to compile OpenFOAM ext solvers on your own; compile and run CUDA; run CUSP examples; and use Open MPI to run parallel cases, then you will be able to compile and run cufflink.  Read the details below for more information.


# Details #

## System Requirements ##
> Due to subtle differences in the forks of [OpenFOAM®](http://openfoam.org/), this library is only compatible with [OpenFOAM-extend](http://www.extend-project.de).  cufflink was developed and tested on [Linux Ubuntu 10.04 amd64 LTS](http://www.ubuntu.com/) with [CUDA 4.2](http://developer.nvidia.com/cuda-downloads), Thrust v1.5 (included in CUDA 4.2), and [Cusp 0.3.0](http://code.google.com/p/cusp-library/). To check your version of Thrust, Cusp, and cufflink, please compile the version.cu file in the cufflink install directory with
```
     nvcc version.cu -o version
```

and run in the `./version` command in the terminal to give the output

```
        Thrust v1.5
        Cusp   v0.3.0
        Cufflink   v0.1.0
```

> Check your version of nvcc with the command nvcc --version and the output should be similar to:

```
        nvcc: NVIDIA (R) Cuda compiler driver
        Copyright (c) 2005-2012 NVIDIA Corporation
        Built on Thu_Apr__5_00:24:31_PDT_2012
        Cuda compilation tools, release 4.2, V0.2.1221
```

> Also, ensure that your open mpi installation is installed correctly so that the runParallel test script will work and the multi-GPU solvers can be run.  To see whether you OpenFOAM install will work in parallel, test Allrun script located in

```
tutorials/multiphase/interFoam/laminar
```

folder by executing

```
./Allrun
```

Without openMPI working correctly, the multi-gpu solvers will not work.

## Installation ##

There are several keys to a successful compilation of cufflink:

  1. Ensure a [working](http://code.google.com/p/cufflink-library/wiki/InstallOpenFOAM) of OpenFOAM-1.6-ext and that one can compile a simple OpenFOAM solver such as [myIcoFoam](http://openfoamwiki.net/index.php/How_to_add_temperature_to_icoFoam) (Specifically, part 2 of the tutorial of "How to add temperature to icoFoam").
  1. Ensure that CUDA 4.2, Thrust v1.5, and Cusp 0.3.0 are installed and running examples.  A simple cusp example is provided in this directory as testcg.cu and is compiled with
```
nvcc -o testcusp testcg.cu 
```

> and run with the command

```
./testcusp  
```

> The output from this test solver (with the provided A.mtx file from Cusp) should be similar to:

```
	Solver will continue until residual norm 0.01 or reaching 100 iterations 
	  Iteration Number  | Residual Norm
	                0       1.000000e+01
	                1       1.414214e+01
	                2       1.093707e+01
        	        3       8.949321e+00
        	        4       6.190057e+00
        	        5       3.835191e+00
        	        6       1.745482e+00
        	        7       5.963550e-01
        	        8       2.371136e-01
        	        9       1.152524e-01
        	       10       3.134469e-02
        	       11       1.144416e-02
        	       12       1.824177e-03
	Successfully converged after 12 iterations.
```

> 3. To compile cufflink, move to the installation directory of cufflink and use the nvccWmakeAll script:

```
	./nvccWmakeAll <arch> > make.log 2>&1
```

> Where `<arch>` is sm\_10, sm\_13, or sm\_20 depending on your device architecture.  The standard output stream and error messages should be printed to a file called make.log for you.  **Many warnings are given dealing with white spaces and tuple.inl functions and can be ignored with this version of cufflink.**

## Multi-GPU Implementation **REQUIRED** ##
### Changes in `lduInterface.H` ###
> The multi-gpu implementation **REQUIRES the file `lduInterface.H`** in cufflink to replace the OpenFOAM core file located in:
```
	OpenFOAM/OpenFOAM-1.6-ext/src/OpenFOAM/matrices/lduMatrix/lduAddressing/lduInterface 
```

> This change requires a complete recompilation of openfoam.  The reason the new file is required is that the parallel code needs to know information about neighbor processors for each processor patch.  New virtual functions were added to `lduInterface.H`.

### Possible changes in compute mode ###
> Each time the multi-gpu solver is used the first iteration is usually the longest when it bumps into the MPI functions for the first time.  This is most likely caused by having the card in compute-mode 0 which allows multiple parallel cpu threads to access the same gpu.  In order to reduce this overhead associated with cpu-gpu thread initialization, change the compute-mode to 1 by running:

```
	sudo /usr/bin/nvidia-smi -c 1 -g <gpu to change>
```

> where the `<gpu to change>` is an integer labeling the device number.  To find the device number, use the deviceQuery program from the Nvidia SDK to see before and after your call of nvidia-smi.  To see supported devices that allow compute mode changes, simply call nvidia-smi without flags or arguements.  See the following thread on the cusp-users group for more information ( http://goo.gl/ojQc9 and also http://goo.gl/kc1YU ).

## Running in OpenFOAM extend ##
> Once the cufflink library has been compiled, in order to use the code in openfoam one needs to include the line

```
	libs ("libCufflink.so");
```

> in your controlDict file of the case you are running.  In addition, a solver must be chosen in the fvSolution file under solver key word similar to:

```
    p
    {
	solver		cufflink_CG;
        preconditioner  none;
        tolerance        1e-10;
        //relTol           1e-08;
	maxIter		 10000;
	storage		    1;//COO=0 CSR=1 DIA=2 ELL=3 HYB=4 all other numbers use default CSR
	gpusPerMachine	2;//for multi gpu version on a machine with 2 gpus per machine node
	AinvType	;
	dropTolerance	;
	linStrategy	;
    }

```

This particular setup uses an un-preconditioned conjugate gradient solver on a single GPU; compressed sparse row (CSR)matrix storage; 1e-8 absolute tolerance and 0 relative tolerance; with a maximum number of inner iterations of 10,000.  For more examples, see the cufflink test cases provided in the distribution.

### Testing the installation ###
> To ensure that your original test directory is not altered and/or accidentally uploaded to the svn server, make a copy of the cufflinkTest folder and place this in the $FOAM\_RUN directory.

> A group of test cases has been provided in order to ensure that the install is working.  They are located in the cufflinkTest/testCases directory and require that the cufflinkTest/testCufflinkFoam application be compiled.  Once you are ready, go the the test/testCases folder and run the runSerialTests bash script to compare cufflink and OpenFOAM.

> Running the script getTimes will grep all the log files and extract time information (this script needs improvement).  Clean the test cases with the Allclean script to remove data since a full test will take >11 gb of harddrive space.  To test the parallel cases, use the runParallelTests script.  It is currently designed to run up to 6 processors and cufflinkTest/decomposeParDict file must be changed accordingly if you workstation has less than 6 processors

> For additional information to be printed to the screen such as the storage method, normalization factor, etc. go to the openfoam file:

```
		$FOAM_INST_DIR/etc/controlDict
```

> and change lduMatrix to 2 This will print additional information to the screen to aid in debugging.