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
        diagonal preconditioned conjugate gradient 
	solver for symmetric Matrices using a CUSP CUDA™ based solver.
                                                             
\**********************************************************************/

extern "C" void CFL_CG_Parallel(cusp_equation_system *CES,  OFSolverPerformance *OFSP, const cpuInterfaces *OFInterfaces)
{
	//convert the interface matrices from OpenFOAM to device memory interfaces.
	gpuInterfaces CFLInterfaces(OFInterfaces);

	ValueType small = 1e-20;// used to prevent floating point exception

	// Populate A
	#include "../CFL_Headers/fillCOOMatrix.H"
	
	#include "../CFL_Headers/setParallelGPUStorage.H"

	const size_t N = A.num_rows;

	// allocate workspace
	cusp::array1d<ValueType,MemorySpace> y(N,0);
	cusp::array1d<ValueType,MemorySpace> r(N,0);

    	// y <- Ax
	cusp::multiply(A, X, y);

	#include "AXLoop.H"

	// r <- b - A*x
	cusp::blas::axpby(B, y, r, ValueType(1), ValueType(-1));

	ValueType normFactor = 1.0;
	
	#include "../CFL_Headers/buildGlobalNormFactor.H"

	//start the krylov solver

	assert(A.num_rows == A.num_cols);        // sanity check
	
	// allocate workspace
	cusp::array1d<ValueType,MemorySpace> Ap(N);
	cusp::array1d<ValueType,MemorySpace> rold(N);
	cusp::array1d<ValueType,MemorySpace> p(N);

	cusp::array1d< ValueType, hostMemorySpace > ph(p.size(),0);//used in ApjLoop.H
	cusp::array1d<ValueType,hostMemorySpace> *pjh;
	pjh = new cusp::array1d<ValueType,hostMemorySpace>[CFLInterfaces.nParInterfaces];
	for(int j = 0;j<CFLInterfaces.nParInterfaces;j++){ pjh[j] = cusp::array1d<ValueType,hostMemorySpace> (OFInterfaces->nColsInterface[j]);	}

	//cusp::array1d< ValueType, hostMemorySpace > pjh(CFLInterfaces.nColsInterface[j],0);



	ValueType alpha = 1;
	ValueType beta = 1;
	        
  	// r0 -> p0
  	cusp::blas::copy(r, p);

    ValueType normR = gpuSumMag(r)/normFactor;
    ValueType normR0 = normR;//initial residual
    OFSP->iRes	= normR0;
    int count = 0;

 	if(0 == CFLInterfaces.myThreadNumber && OFSP->debugCusp){
 		std::cout << "   Iteration "<<count<<" residual = "<< std::setw(10) << normR << std::endl;
 	}
	
    while ( normR > (OFSP->tol) && count<= (OFSP->maxIter) && normR/normR0 >= (OFSP->relTol))
    {
        // Apj <- A*p
        cusp::multiply(A, p, Ap);

		#include "ApjLoop.H"
      
        // alpha <- <r,z>/<Ap,p>
    	alpha =  gpuSumProd(r,r) / gpuSumProd(Ap, p);

        // x <- x + alpha * p
        cusp::blas::axpy(p, X, alpha);

        //copy rold<-r
        cusp::blas::copy(r, rold);

        // r <- r - alpha * y		
        cusp::blas::axpy(Ap, r, -alpha);

        // beta <- <r,r>/<r_old,r_old> 
        beta = gpuSumProd(r, r) / gpuSumProd(rold, rold);
		
        // p <- r + beta*p should be p <- z + beta*p
        cusp::blas::axpby(r, p, p, ValueType(1), beta);

        normR = gpuSumMag(r)/normFactor;

        count++;

        if(0 == CFLInterfaces.myThreadNumber && OFSP->debugCusp){
        	std::cout << "   Iteration "<<count<<" residual = "<< std::setw(10) << normR << std::endl;
        }
    }

 	//end the krylov solver

	//final residual
	OFSP->fRes = normR;
	OFSP->nIterations=count;

	//converged?
	if(OFSP->fRes<=OFSP->tol || OFSP->fRes/OFSP->iRes<=OFSP->relTol)
		OFSP->converged=true;
	else
		OFSP->converged=false;

	//pass the solution vector back	
	CES->X = X;

	//delete the dynamically allocated arrays
	delete [] pjh; //dcombest Feb 13, 2012

}
