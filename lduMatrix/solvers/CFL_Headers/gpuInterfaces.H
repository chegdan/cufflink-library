#ifndef GPUINTERFACES_H
#define GPUINTERFACES_H

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

    Class
	gpuInterfaces

    Author
        Daniel P. Combest.  All rights reserved.

    Description
        Collection of dynamic arrays to hold interfaces from cpu to be 
	passed to the GPU for use in a parallel implementation of the 
	Krylov space solvers. 

    Source File
	gpuInterfaces.C
                                                            
\**********************************************************************/

#include <cusp/coo_matrix.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <vector>
#include "cuspTypeDefs.H"
#include "cpuInterfaces.H"

class gpuInterfaces{

public:
	int gpusPerMachine;	//how many gpus per machine are present in the cluster?
	int myThreadNumber;	//holds the thread number of the cpu the gpuInterfaces originated
	std::vector<int> neighbProcNo;	//array to hold the neighbor processor number of each thread
	int nParInterfaces;	//holds the number of parallel interface for a particular processor
	std::vector<int> interfaceAddr;	//addressing array of patch id for parallel patches only
	int nRowsInterface;	//number of rows in the interface matrix.  always equal to Aii nCells
	std::vector<int> nColsInterface;	//array holding number of columns for each parallel interface
	std::vector<int> nnz;		//array holding the number of non-zeros that an interface possesses
	cusp::coo_matrix<IndexType, ValueType, MemorySpace> *Aij;//array of sparse matrices holding the interface matrices themselves
	cusp::array1d< ValueType, hostMemorySpace > *Xjh;//array of host vectors from neighbor processors

	gpuInterfaces(const cpuInterfaces *);//convert the cpuInterfaces
	gpuInterfaces();//do nothing constructor
	~gpuInterfaces();
	void printShortInfo();//used in testing


};
#endif //GPUINTERFACES_H
