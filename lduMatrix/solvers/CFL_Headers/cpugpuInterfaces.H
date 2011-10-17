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
        Collection of dynamic arrays to hold interfaces from cpu to be 
	passed to the GPU for use in a parallel implementation of the 
	Krylov space solvers. 
                                                            
\**********************************************************************/
#ifndef CPUGPUINTERFACES_H
#define CPUGPUINTERFACES_H


#include <cusp/coo_matrix.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <vector>
#include "cuspTypeDefs.H"

class cpuInterfaces{

public:
	int myThreadNumber;	//holds the thread number of the cpu the gpuInterfaces originated
	std::vector<int> neighbProcNo;	//array to hold the neighbor processor number of each thread
	int nParInterfaces;	//holds the number of parallel interface for a particular processor
	std::vector<int> interfaceAddr;	//addressing array of patch id for parallel patches only
	int nRowsInterface;	//number of rows in the interface matrix.  always equal to Aii nCells
	std::vector<int> nColsInterface;	//array holding number of columns for each parallel interface
	std::vector<int> nnz;		//array holding the number of non-zeros that an interface possesses
	std::vector< cusp::coo_matrix<IndexType, ValueType, hostMemorySpace> > Aij;//array of sparse matrices holding the interface matrices themselves
	std::vector< cusp::array1d< ValueType, hostMemorySpace> > Xj;//array holding the interface variable vectors that will be moved around in the CUDA code
	std::vector< cusp::array1d< ValueType, hostMemorySpace> > Bj;//array holding the interface source vectors that will be moved around in the CUDA code
	
	cpuInterfaces();
	~cpuInterfaces();


};

cpuInterfaces::cpuInterfaces(){

}

cpuInterfaces::~cpuInterfaces(){

}

class gpuInterfaces{

public:

	int myThreadNumber;	//holds the thread number of the cpu the gpuInterfaces originated
	std::vector<int> neighbProcNo;	//array to hold the neighbor processor number of each thread
	int nParInterfaces;	//holds the number of parallel interface for a particular processor
	std::vector<int> interfaceAddr;	//addressing array of patch id for parallel patches only
	int nRowsInterface;	//number of rows in the interface matrix.  always equal to Aii nCells
	std::vector<int> nColsInterface;	//array holding number of columns for each parallel interface
	std::vector<int> nnz;		//array holding the number of non-zeros that an interface possesses
	thrust::device_vector< cusp::coo_matrix<IndexType, ValueType, MemorySpace> > Aij;//array of sparse matrices holding the interface matrices themselves
	thrust::device_vector< cusp::array1d< ValueType, MemorySpace> > Xj;//array holding the interface variable vectors that will be moved around in the CUDA code
	thrust::device_vector< cusp::array1d< ValueType, MemorySpace> > Bj;//array holding the interface source vectors that will be moved around in the CUDA code
	//cusp::coo_matrix<IndexType, ValueType, MemorySpace> *Aij;//array of sparse matrices holding the interface matrices themselves
	//cusp::array1d< ValueType, MemorySpace> *Xj;//array holding the interface variable vectors that will be moved around in the CUDA code
	//cusp::array1d< ValueType, MemorySpace> *Bj;//array holding the interface source vectors that will be moved around in the CUDA code

	gpuInterfaces();
	~gpuInterfaces();


};

gpuInterfaces::gpuInterfaces(){

}

gpuInterfaces::~gpuInterfaces(){

}

#endif //CPUGPUINTERFACES_H
