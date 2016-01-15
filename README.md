* TITLE: *cufflink README for version 0.0.1*
* +AUTHOR:   Daniel P. Combest.
* DATE:     7 Sept 2011
* LINK:     To Be Determined
* OPTIONS:  Daniel Combest
* COPYRIGHT Held by Daniel Combest

# Preamble
      ______  __    __   _______  _______  __       __   __   __   __  ___     
     /      ||  |  |  | |   ____||   ____||  |     |  | |  \ |  | |  |/  /     
    |  ,----'|  |  |  | |  |__   |  |__   |  |     |  | |   \|  | |  '  /  
    |  |     |  |  |  | |   __|  |   __|  |  |     |  | |  . `  | |    <   
    |  `----.|  `--'  | |  |     |  |     |  `----.|  | |  |\   | |  .  \
     \______| \______/  |__|     |__|     |_______||__| |__| \__| |__|\__\

Cuda For FOAM Link

cufflink is a library for linking numerical methods based on Nvidia's 
Compute Unified Device Architecture (CUDA™) C/C++ programming language
and OpenFOAM®.

Please note that cufflink is not approved or endorsed by OpenCFD® 
Limited, the owner of the OpenFOAM® and OpenCFD® trademarks and 
producer of OpenFOAM® software.

The official web-site of OpenCFD® Limited is www.openfoam.com .

# License
cufflink is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

cufflink is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with cufflink.  If not, see <http://www.gnu.org/licenses/>. 

# System requirements
cufflink was developed and tested on Linux Ubuntu 10.04 amd64 with CUDA 4.0, 
Thrust v1.4, and Cusp 0.2.0. To check your version of Thrust, Cusp, and cufflink, 
please compile the version.cu file and run in the command line to give the output

    Thrust     v1.4
    Cusp       v0.2.0
    Cufflink   v0.0.1

Check your version of nvcc with the command nvcc --version and the output should be similar to: 
    
    nvcc: NVIDIA (R) Cuda compiler driver
    Copyright (c) 2005-2011 NVIDIA Corporation
    Built on Thu_May_12_11:09:45_PDT_2011
    Cuda compilation tools, release 4.0, V0.2.1221

Also, ensure that your mpi installation is correct so that the runParallel test script can be run.  
To see whether your parallel version of openfoam is running, test the interFoam/damBreak tutorial

# Installation
Ensure a working OpenFOAM-1.6-ext version and that one can compile a simple 
OpenFOAM solver.  Ensure that CUDA 4.0, Thrust v1.4, and Cusp 0.2.0 
are installed and running examples.  A simple cusp example is provided in 
this directory as testcg.cu and is compiled with nvcc -o testcusp testcg.cu 
and run with the command ./testcusp  The output from this solver should be similar to:

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

To compile the solvers to be used in OpenFOAM, use the script:  

    ./nvccWmakeAll <arch> > make.log 2>&1

Where <arch> is sm_10, sm_13, or sm_20 depending on your device architecture.  Many warnings are given dealing 
with white spaces and tuple.inl functions and can be ignored with this version of cufflink.

# Multi-GPU Implementation *****REQUIRED*****
The multi-gpu implementation REQUIRES the file lduInterface.H in cufflink to replace the OpenFOAM core file 
located in:

    OpenFOAM/OpenFOAM-1.6-ext/src/OpenFOAM/matrices/lduMatrix/lduAddressing/lduInterface 

This change requires a complete recompilation of openfoam.  Also, each time the multi-gpu solver is used the 
first iteration is usually the longest.  This is most likely caused by having the card in compute-mode 0 
which allows multiple parallel cpu threads to access the same gpu.  In order to reduce this overhead 
associated with cpu-gpu thread initialization, change the compute-mode to 1 by running:
    
    sudo /usr/bin/nvidia-smi -c 1 -g <gpu to change>

where the <gpu to change> is an integer labeling the device number.  To find the device number, use the 
deviceQuery program from the Nvidia SDK to see before and after your call of nvidia-smi.  To see supported 
devices that allow compute mode changes, simply call nvidia-smi without flags or arguements.  See the 
following thread on the cusp-users group for more information ( http://goo.gl/ojQc9 and also http://goo.gl/kc1YU ).

# Running in OpenFoam
Once the cufflink library has been compiled, in order to use the code in openfoam one needs to include the line

    libs ("libCufflink.so");

in your controlDict file of the case you are running.  In addition, a solver must be chosen in the fvSolution 
file under solver key word similar to:

    solver cufflink_CG;
    preconditioner  none;
    tolerance        1e-10;
    maxIter      10000;
    storage         1; //COO=0 CSR=1 DIA=2 ELL=3 HYB=4 all other numbers use default CSR
    gpusPerMachine  2; //for multi gpu version on a machine with 2 gpus per machine node
    AinvType ;
    dropTolerance ;
    linStrategym ;

# Testing the installation
To ensure that your original test directory is not altered and/or accidentally uploaded to the svn server, 
make a copy of the cufflinkTest folder and place this in the $FOAM_RUN directory.

A group of test cases has been provided in order to ensure that the install is working.  They are located in 
the cufflinkTest/testCases directory and require that the cufflinkTest/testCufflinkFoam application be compiled.  Once you are ready,go the test/testCases folder and run the runSerialTests bash script to compare cufflink and OpenFOAM.  
Running the script getTimes will grep all the log files and extract time information (this script needs improvement).  
Clean the test cases with the Allclean script to remove data since a full test will take over 11 GB of hard drive space.
To test the parallel cases, use the runParallelTests script.  It is currently designed to run up to 6 processors and 
cufflinkTest/decomposeParDict file must be changed accordingly if you workstation has less than 6 processors

For additional information to be printed to the screen such as the storage method, normalization factor, etc. go 
to the openfoam file:

        $FOAM_INST_DIR/etc/controlDict

and change lduMatrix to 2.  This will print additional information to the screen to aid in debugging.


# List of Contributors
* Daniel P. Combest
* Jeremy Day


    


