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
        A header file to compile the extern "C" functions containing the 
	CUSP CUDA™ based solvers.  This avoids problems with the multiple 
	definitions of some functions.
                                                             
\**********************************************************************/
#include "cuda.h"

//System Includes
#include <mpi.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <stdio.h>
#include <limits>

//CUSP Includes
#include <cusp/detail/config.h>
#include <cusp/verify.h>
#include <cusp/precond/ainv.h>
#include <cusp/precond/diagonal.h>
#include <cusp/precond/smoothed_aggregation.h>
#include <cusp/csr_matrix.h>
#include <cusp/hyb_matrix.h>
#include <cusp/coo_matrix.h>
#include <cusp/ell_matrix.h>
#include <cusp/dia_matrix.h>
#include <cusp/blas.h>
#include <cusp/multiply.h>

//THRUST Includes
#include <thrust/reduce.h>//for the summation of components of a vector
#include <thrust/sequence.h>//for diagonal construction
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/copy.h> 

//CUFFLINK Includes
#include "CFL_Headers/cuspTypeDefs.H"//change this file ( along with passing the arch flag to sm_10 in nvccwmake script ) if you want single precision

#include "CFL_Headers/OFSolverPerformance.H"
#include "CFL_Headers/cusp_equation_system.H"
#include "CFL_Headers/cpuInterfaces.C"
#include "CFL_Headers/gpuInterfaces.C"
#include "CFL_Headers/globalOps.H"

//CUFFLINK extern function definitions
#include "cuspSolvers/CFL_AinvPCG.cu"
#include "cuspSolvers/CFL_AinvPCG_Parallel.cu"

#include "cuspSolvers/CFL_CG.cu"
#include "cuspSolvers/CFL_CG_Parallel.cu"

#include "cuspSolvers/CFL_DiagPCG.cu"
#include "cuspSolvers/CFL_DiagPCG_Parallel.cu"

#include "cuspSolvers/CFL_SmAPCG.cu"
#include "cuspSolvers/CFL_SmAPCG_Parallel.cu"

#include "cuspSolvers/CFL_DiagPBiCGStab.cu"
#include "cuspSolvers/CFL_DiagPBiCGStab_Parallel.cu"

#include "cuspSolvers/CFL_AinvPBiCGStab.cu"
#include "cuspSolvers/CFL_AinvPBiCGStab_Parallel.cu"


