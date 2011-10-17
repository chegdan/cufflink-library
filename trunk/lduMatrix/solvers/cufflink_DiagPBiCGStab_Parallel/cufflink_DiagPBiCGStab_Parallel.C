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
                                                         
\**********************************************************************/

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "cufflink_DiagPBiCGStab_Parallel.H"

//for COO and vector on host
#include <cusp/coo_matrix.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/copy.h> 
#include "../CFL_Headers/cuspTypeDefs.H"
//end headers for COO and vectors

#include <iomanip>
#include <iostream>
#include <string>
#include <stdio.h>
//using namespace std;

#include "../CFL_Headers/OFSolverPerformance.H"
#include "../CFL_Headers/cusp_equation_system.H"
#include "../CFL_Headers/cpuInterfaces.H"

extern "C" void CFL_DiagPBiCGStab_Parallel(cusp_equation_system *,  OFSolverPerformance *, const cpuInterfaces *);//extern function in CUDA that solves the equation system

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cufflink_DiagPBiCGStab_Parallel, 0);

    lduSolver::addasymMatrixConstructorToTable<cufflink_DiagPBiCGStab_Parallel>
        addcufflink_DiagPBiCGStab_ParallelAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cufflink_DiagPBiCGStab_Parallel::cufflink_DiagPBiCGStab_Parallel
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& dict
)
:
    lduSolver
    (
        fieldName,
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces,
        dict
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::lduSolverPerformance Foam::cufflink_DiagPBiCGStab_Parallel::solve
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{

if (!Pstream::parRun())
    {
	FatalErrorIn("cufflink:\nMulti-GPU Solver cannot be run serially.  Choose the serial version of the solver")<< "Fatal Error"<< abort(FatalError);
    }

OFSolverPerformance OFSP;//container for solver performance 

//get the sparse matrix storage method from fvSolution 
#include "../CFL_Headers/getGPUStorage.H"

    // --- Setup class containing solver performance data
    lduSolverPerformance solverPerf
    (
        "cufflink_DiagPBiCGStab_Parallel",
        fieldName()
    );


//insert my solver code here
    register label nCells = x.size();
    register label NFaces = matrix().lower().size();
    OFSP.nCells = nCells;
    OFSP.nFaces = NFaces;
    OFSP.debugCusp = false;
    
    //turn on the debug switch to print some important information
    if (lduMatrix::debug >= 2)
    {
         OFSP.debugCusp = true;
    }

//interface gathering

	//used when launching kernals in machines with multiple gpus
	int gpusPerMachine = readInt(dict().lookup("gpusPerMachine"));

	//get all the interface information
	#include "../CFL_Headers/getInterfaces.H"

	OFInterfaces.gpusPerMachine = gpusPerMachine;

	if(!OFInterfaces.checkGPUCount(OFInterfaces.gpusPerMachine))
		FatalErrorIn("checkGPUCount:\nThe number of GPUs per machine is not correct.\nChange the gpusPerMachine variable in fvSolution")<< "Fatal Error"<< abort(FatalError);

//end interface gathering

	cusp_equation_system CES;

//declare the COO, X and b and copy the OpenFOAM values here
	CES.A = cusp::coo_matrix<IndexType, ValueType, hostMemorySpace>(nCells,nCells,nCells+2*NFaces);
	CES.X = cusp::array1d< ValueType, hostMemorySpace>(nCells);
	CES.B = cusp::array1d< ValueType, hostMemorySpace>(nCells);

	//copy values from the lduMatrix to our equation system
	thrust::copy(matrix().diag().begin(),matrix().diag().end(),CES.A.values.begin());//copy values of lduMatrix diag to A COO matrix
	thrust::copy(matrix().lower().begin(),matrix().lower().end(),CES.A.values.begin()+nCells);//copy values of lduMatrix lower to A COO matrix
	thrust::copy(matrix().upper().begin(),matrix().upper().end(),CES.A.values.begin()+nCells+NFaces);//copy values of lduMatrix upper to A COO matrix

	//copy row and column indices of lower and upper to our equations system
	//do not initialize the row and column values of diag to save time-> performed on GPU
	thrust::copy(matrix().lduAddr().upperAddr().begin(),matrix().lduAddr().upperAddr().end(),CES.A.row_indices.begin()+nCells);//copy row indices of lower into A COO matrix
	thrust::copy(matrix().lduAddr().lowerAddr().begin(),matrix().lduAddr().lowerAddr().end(),CES.A.column_indices.begin()+nCells);//copy column indices of lower into A COO matrix
	thrust::copy(matrix().lduAddr().lowerAddr().begin(),matrix().lduAddr().lowerAddr().end(),CES.A.row_indices.begin()+nCells+NFaces);//copy row indices of upper into A COO matrix
	thrust::copy(matrix().lduAddr().upperAddr().begin(),matrix().lduAddr().upperAddr().end(),CES.A.column_indices.begin()+nCells+NFaces);//copy column indices of upper into A COO matrix

	thrust::copy(x.begin(),x.end(),CES.X.begin());//copy x of lower into x vector
	thrust::copy(b.begin(),b.end(),CES.B.begin());//copy b of lower into b vector
//end COO, X, b initialize and fill


	OFSP.iRes = -1;//initialization value for initial residual
	OFSP.fRes = -1;//initialization value for final residual
	OFSP.nIterations = 0;//initialization value for number of iterations
	OFSP.minIter = minIter();//minimum iterations
	OFSP.maxIter = maxIter();//maximum iterations
	OFSP.relTol = relTolerance();//relative tolerance
	OFSP.tol = tolerance();//tolerance

	//set the device
	cudaError_t cudareturn; 
	cudareturn = cudaSetDevice(Pstream::myProcNo() % (gpusPerMachine));
	if (cudareturn == cudaErrorInvalidDevice){ perror("cudaSetDevice returned  cudaErrorInvalidDevice");  } 

	CFL_DiagPBiCGStab_Parallel(&CES, &OFSP, &OFInterfaces);//call the CUDA code to solve the system
	
	thrust::copy(CES.X.begin(),CES.X.end(),x.begin());//copy the x vector back to Openfoam
	solverPerf.initialResidual() = OFSP.iRes;//return initial residual
	solverPerf.finalResidual() = OFSP.fRes;//return final residual
	solverPerf.nIterations() = OFSP.nIterations;//return the number of iterations
	//solverPerf.checkConvergence(tolerance(), relTolerance()) ; what should be passed here and does this waste time?
//	solverPerf.converged_=OFSP.converged;//private, cannot access
//	solverPerf.singular_=OFSP.singular;//private, cannot access

    return solverPerf;
}


// ************************************************************************* //
