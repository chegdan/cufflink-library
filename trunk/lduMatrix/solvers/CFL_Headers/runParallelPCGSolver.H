/**********************************************************************\
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

------------------------------------------------------------------------
This file is part of cufflink.

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

    Author
        Daniel P. Combest.  All rights reserved.

    Description
        Generic header file to pass ax=b system to the cuda solvers        
                                                         
\**********************************************************************/

// y <- Ax
cusp::multiply(A, X, y);

// r <- b - A*x
blas::axpby(B, y, r, ValueType(1), ValueType(-1));

//build the normfactor for residual calcualtions
#include "buildGlobalNormfactor.h"

//initial residual
OFSP->iRes=gpuSumMag(r)/normFactor;

MPI_Bcast(&OFSP->iRes, 1, MPI_SCALAR,0, MPI_COMM_WORLD);//send it out to all the nodes

if(OFSP->debugCusp){
	//create a verbose monitor to monitor solution and report iteration count and residual at each iteration
	cusp::verbose_par_openFOAM_monitor<ValueType> monitor(B, OFSP->maxIter, OFSP->relTol,OFSP->tol, normFactor, ValueType(OFSP->iRes));

	// solve the linear system A * x = b with the preconditioned Conjugate Gradient method
	cusp::krylov::par_of_cg(A, X, B, monitor,M);

	//number of iterations
	OFSP->nIterations=monitor.iteration_count();
}else{
	//create a monitor to monitor solution
	cusp::par_openFOAM_monitor<ValueType> monitor(B, OFSP->maxIter, OFSP->relTol,OFSP->tol, normFactor, ValueType(OFSP->iRes));

	// solve the linear system A * x = b with the preconditioned Conjugate Gradient method
	cusp::krylov::par_of_cg(A, X, B, monitor,M);

	//number of iterations
	OFSP->nIterations=monitor.iteration_count();
}

// y <- Ax
cusp::multiply(A, X, y);

// r <- b - A*x
blas::axpby(B, y, r, ValueType(1), ValueType(-1));

//final residual
OFSP->fRes=gpuSumMag(r)/normFactor;
MPI_Bcast(&OFSP->fRes, 1, MPI_SCALAR,0, MPI_COMM_WORLD);//send it out to all the nodes


