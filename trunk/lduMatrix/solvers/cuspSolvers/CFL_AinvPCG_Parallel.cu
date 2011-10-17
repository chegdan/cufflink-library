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
        Algebraic Multigrid (AMG) based on Smoothed Aggregation
	preconditioned conjugate gradient solver for symmetric 
	matrices using a CUSP CUDA™ based solver.
                                                             
\**********************************************************************/

extern "C" void CFL_AinvPCG_Parallel(cusp_equation_system *CES,  OFSolverPerformance *OFSP, const cpuInterfaces *OFInterfaces)
{
	//convert the interface matrices from OpenFOAM to device memory interfaces.
	gpuInterfaces CFLInterfaces(OFInterfaces);

	ValueType small = 1e-20;// used to prevent floating point exception

	// Populate A
	#include "../CFL_Headers/fillCOOMatrix.H"
	
	#include "../CFL_Headers/setParallelGPUStorage.H"

	if(0 == CFLInterfaces.myThreadNumber && OFSP->debugCusp){std::cout<<"   Using Cusp_Ainv preconditioner - Bridson Ainv with novel dropping strategy and linStrategy = "<< OFSP->linStrategy<<"\n";}

	cusp::precond::bridson_ainv<ValueType, MemorySpace> M(A, 0, -1, true, int(OFSP->linStrategy));

	//start the krylov solver
    assert(A.num_rows == A.num_cols);        // sanity check

    const size_t N = A.num_rows;

    // allocate workspace
    cusp::array1d<ValueType,MemorySpace> y(N);
    cusp::array1d<ValueType,MemorySpace> z(N);
    cusp::array1d<ValueType,MemorySpace> r(N);
    cusp::array1d<ValueType,MemorySpace> p(N);
    cusp::array1d<ValueType,MemorySpace> Ap(N);
        
    // y <- Ax
    cusp::multiply(A, X, y);

    //take into account the interfaces
    #include "AXLoop.H"

    //define the normalization factor
    ValueType normFactor = 1.0;

    #include "../CFL_Headers/buildGlobalNormFactor.H"

    // r <- b - A*x
    cusp::blas::axpby(B, y, r, ValueType(1), ValueType(-1));
   
    // z <- M*r
    cusp::multiply(M, r, z);

    // p <- z
    cusp::blas::copy(z, p);
		
    // rz = <r^H, z>
    ValueType rz = gpuSumProd(r, z);

    ValueType normR = gpuSumMag(r)/normFactor;
    ValueType normR0 = normR;//initial residual
    OFSP->iRes	= normR0;
    int count = 0;

 	if(0 == CFLInterfaces.myThreadNumber && OFSP->debugCusp){
 		std::cout << "   Iteration "<<count<<" residual = "<< std::setw(10) << normR << std::endl;
 	}
	
    while ( normR > (OFSP->tol) && count<= (OFSP->maxIter) && normR/normR0 >= (OFSP->relTol))
    {
        // Ap <- A*p
        cusp::multiply(A, p, Ap);

        //take into account the interfaces
		#include "ApjLoop.H"
        
        // alpha <- <r,z>/<Ap,p>
        ValueType alpha =  rz / gpuSumProd(Ap, p);

        // x <- x + alpha * p
        cusp::blas::axpy(p, X, alpha);

        // r <- r - alpha * Ap		
        cusp::blas::axpy(Ap, r, -alpha);

        // z <- M*r
        cusp::multiply(M, r, z);
		
        ValueType rz_old = rz;

        // rz = <r^H, z>
        rz = gpuSumProd(r, z);

        // beta <- <r_{i+1},r_{i+1}>/<r,r> 
        ValueType beta = rz / rz_old;
		
        // p <- r + beta*p should be p <- z + beta*p
        cusp::blas::axpby(z, p, p, ValueType(1), beta);

        normR = gpuSumMag(r)/normFactor;;

        count++;

        if(0 == CFLInterfaces.myThreadNumber && OFSP->debugCusp){
        	std::cout << "   Iteration "<<count<<" residual = "<< std::setw(10) << normR << std::endl;
        }
    } //end the krylov solver

	//final residual
	OFSP->fRes = normR;
	OFSP->nIterations = count;

	//converged?
	if(OFSP->fRes<=OFSP->tol || OFSP->fRes/OFSP->iRes<=OFSP->relTol)
		OFSP->converged=true;
	else
		OFSP->converged=false;

	//pass the solution vector back	
	CES->X = X;
}
