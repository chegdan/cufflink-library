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
        diagonal preconditioned conjugate gradient stabilized 
	solver for asymmetric Matrices using a CUSP CUDA™ based solver.
                                                             
\**********************************************************************/

extern "C" void CFL_DiagPBiCGStab_Parallel(cusp_equation_system *CES,  OFSolverPerformance *OFSP, const cpuInterfaces *OFInterfaces){

	//convert the interface matrices from OpenFOAM to device memory interfaces.
	gpuInterfaces CFLInterfaces(OFInterfaces);

	ValueType small = 1e-20;// used to prevent floating point exception

	cusp::coo_matrix<IndexType, ValueType, MemorySpace> A(CES->A);
	cusp::array1d<ValueType, MemorySpace> X(CES->X);
	cusp::array1d<ValueType, MemorySpace> B(CES->B);

	//fill in the rest of the diag (rows and col)
	thrust::sequence(A.row_indices.begin(),A.row_indices.begin()+OFSP->nCells);//determine row indices of diagonal values and fill A COO matrix
	thrust::sequence(A.column_indices.begin(),A.column_indices.begin()+OFSP->nCells);//determine column indices of diagonal values and fill A COO matrix

//	A.sort_by_row();//coo must be sorted by row.  Improves perfomance by 2X normal with this line
	A.sort_by_row_and_column();//sorted coo by row and column. speeds code up a little bit more

	//choose and convert to correct matrix format	
	#include "../CFL_Headers/setParallelGPUStorage.H"

	if(0 == CFLInterfaces.myThreadNumber && OFSP->debugCusp){ printf("   Using Cusp_Diagonal preconditioner\n");}

 			cusp::precond::diagonal<ValueType, MemorySpace> M(A);	

//start Krylov solver
    assert(A.num_rows == A.num_cols);        // sanity check

    const size_t N = A.num_rows;

    // allocate workspace
    cusp::array1d<ValueType,MemorySpace> y(N);

    cusp::array1d<ValueType,MemorySpace>   p(N);
    cusp::array1d<ValueType,MemorySpace>   r(N);
    cusp::array1d<ValueType,MemorySpace>   r_star(N);
    cusp::array1d<ValueType,MemorySpace>   s(N);
    cusp::array1d<ValueType,MemorySpace>  Mp(N);
    cusp::array1d<ValueType,MemorySpace> AMp(N);
    cusp::array1d<ValueType,MemorySpace>  Ms(N);
    cusp::array1d<ValueType,MemorySpace> AMs(N);

	cusp::array1d< ValueType, hostMemorySpace > Msh(Ms.size(),0);
	cusp::array1d< ValueType, hostMemorySpace > Mph(Mp.size(),0);
	cusp::array1d< ValueType, hostMemorySpace > *Mpjh;
	cusp::array1d< ValueType, hostMemorySpace > *Msjh;
	Mpjh = new cusp::array1d<ValueType,hostMemorySpace>[CFLInterfaces.nParInterfaces];
	Msjh = new cusp::array1d<ValueType,hostMemorySpace>[CFLInterfaces.nParInterfaces];

	for(int j = 0;j<CFLInterfaces.nParInterfaces;j++){ 
		Mpjh[j] = cusp::array1d<ValueType,hostMemorySpace> (OFInterfaces->nColsInterface[j]);	
		Msjh[j] = cusp::array1d<ValueType,hostMemorySpace> (OFInterfaces->nColsInterface[j]);	
	}

    // y <- Ax
    cusp::multiply(A, X, y);

    #include "AXLoop.H"

    //define the normalization factor
    ValueType normFactor = 1.0;

    #include "../CFL_Headers/buildGlobalNormFactor.H"

    // r <- b - A*x
    cusp::blas::axpby(B, y, r, ValueType(1), ValueType(-1));

    // p <- r
    cusp::blas::copy(r, p);

    // r_star <- r
    cusp::blas::copy(r, r_star);

    ValueType r_r_star_old = gpuSumProd(r_star, r);


    ValueType normR = gpuSumMag(r)/normFactor;
    ValueType normR0 = normR;//initial residual
    OFSP->iRes	= normR0;
    int count = 0;

 	if(0 == CFLInterfaces.myThreadNumber && OFSP->debugCusp){std::cout << "   Iteration "<<count<<" residual = "<< std::setw(10) << normR << std::endl;}
	
    while ( normR > (OFSP->tol) && count<= (OFSP->maxIter) && normR/normR0 >= (OFSP->relTol))
    {
        // Mp = M*p
        cusp::multiply(M, p, Mp);

        // AMp = A*Mp
        cusp::multiply(A, Mp, AMp);

    	#include "AMpLoop.H"

        // alpha = (r_j, r_star) / (A*M*p, r_star)
        ValueType alpha = r_r_star_old / gpuSumProd(r_star, AMp);
        
        // s_j = r_j - alpha * AMp
        cusp::blas::axpby(r, AMp, s, ValueType(1), ValueType(-alpha));

	ValueType normS = gpuSumMag(s)/normFactor;

	if (!( normS > (OFSP->tol) && normS/normR0 >= (OFSP->relTol))){//is this right?
	  // x += alpha*M*p_j
	  cusp::blas::axpby(X, Mp, X, ValueType(1), ValueType(alpha));
          
          // y <- Ax
          cusp::multiply(A, X, y);

	  #include "AX2Loop.H"

          // r <- b - A*x
          cusp::blas::axpby(B, y, r, ValueType(1), ValueType(-1));

	  normR = gpuSumMag(r)/normFactor;

	  count++;

	  if(0 == CFLInterfaces.myThreadNumber && OFSP->debugCusp){std::cout << "   Iteration "<<count<<" residual = "<< std::setw(10) << normR << std::endl;}
	  break;
	}

        // Ms = M*s_j
        cusp::multiply(M, s, Ms);
        
        // AMs = A*Ms
        cusp::multiply(A, Ms, AMs);

    	#include "AMsLoop.H"

        // omega = (AMs, s) / (AMs, AMs)
        ValueType omega = gpuSumProd(AMs, s) / gpuSumProd(AMs, AMs);
        
        // x_{j+1} = x_j + alpha*M*p_j + omega*M*s_j
        cusp::blas::axpbypcz(X, Mp, Ms, X, ValueType(1), alpha, omega);

        // r_{j+1} = s_j - omega*A*M*s
        cusp::blas::axpby(s, AMs, r, ValueType(1), -omega);

        // beta_j = (r_{j+1}, r_star) / (r_j, r_star) * (alpha/omega)
        ValueType r_r_star_new = gpuSumProd(r_star, r);
        ValueType beta = (r_r_star_new / r_r_star_old) * (alpha / omega);
        r_r_star_old = r_r_star_new;

        // p_{j+1} = r_{j+1} + beta*(p_j - omega*A*M*p)
        cusp::blas::axpbypcz(r, p, AMp, p, ValueType(1), beta, -beta*omega);

	normR = gpuSumMag(r)/normFactor;

	count++;

	if(0 == CFLInterfaces.myThreadNumber && OFSP->debugCusp){std::cout << "   Iteration "<<count<<" residual = "<< std::setw(10) << normR << std::endl;}
    }
//end Krylov Solver

	//final residual
	OFSP->fRes = gpuSumMag(r)/normFactor;
	OFSP->nIterations = count;

	//converged?
	if(OFSP->fRes<=OFSP->tol || OFSP->fRes/OFSP->iRes<=OFSP->relTol)
		OFSP->converged=true;
	else
		OFSP->converged=false;

	//pass the solution vector back	
	CES->X = X;


}
