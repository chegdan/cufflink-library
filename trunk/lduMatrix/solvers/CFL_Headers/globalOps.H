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
        Collection of methods to perform global sum, global dot product
	and global vector average compute.  Necessary for the 
	calculation of the normalization factro used in OpenFOAM®
                                                            
\**********************************************************************/
#ifndef GLOBALOPS_H
#define GLOBALOPS_H

//#include "cuspTypeDefs.H"
//#include <cusp/blas.h>
//#include <thrust/reduce.h>//for the summation of components of a vector

//include the mpi libraries
//#include <mpi.h>


#define MPI_CHECK(call) \
    if((call) != MPI_SUCCESS) { \
        std::cerr << "MPI error calling \""#call"\"\n"; \
        std::cout << "Test FAILED\n"; \
	MPI_Abort(MPI_COMM_WORLD, -1); }


ValueType gpuSumProd(cusp::array1d< ValueType, MemorySpace> &a, cusp::array1d< ValueType, MemorySpace> &b){

	ValueType dotNode = 0;
	ValueType SumProdValue = 0;

	// dotNode = <a,b>
        dotNode = cusp::blas::dotc(a, b);

	MPI_CHECK(MPI_Reduce(&dotNode, &SumProdValue, 1, MPI_SCALAR, MPI_SUM, 0, MPI_COMM_WORLD));
	MPI_CHECK(MPI_Bcast(&SumProdValue, 1, MPI_SCALAR, 0, MPI_COMM_WORLD));//send it out to all the nodes

	return SumProdValue;

}


ValueType gpuSumMag(cusp::array1d< ValueType, MemorySpace> &a){

	ValueType magNode = 0;
	ValueType SumMagValue = 0;

	// magNode = ||a||
        magNode = cusp::blas::nrm2(a);

	MPI_CHECK(MPI_Reduce(&magNode, &SumMagValue, 1, MPI_SCALAR, MPI_SUM, 0, MPI_COMM_WORLD));
	MPI_CHECK(MPI_Bcast(&SumMagValue, 1, MPI_SCALAR, 0, MPI_COMM_WORLD));//send it out to all the nodes

	return SumMagValue;

}

ValueType gpuAverage(cusp::array1d< ValueType, MemorySpace> &a){

	ValueType sumNode = 0;
	int sizeTotal = 1;
	ValueType SumValue = 1;

	int sizeNode = a.size();

	// sumNode = sum(a)
        sumNode = thrust::reduce(a.begin(),a.end());

	MPI_CHECK(MPI_Reduce(&sumNode, &SumValue, 1, MPI_SCALAR, MPI_SUM, 0, MPI_COMM_WORLD));

	MPI_CHECK(MPI_Reduce(&sizeNode, &sizeTotal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD));///not sure if this is right
	
	ValueType temp = 0;
	
	temp = SumValue/ValueType(sizeTotal);

	MPI_CHECK(MPI_Bcast(&temp, 1, MPI_SCALAR, 0, MPI_COMM_WORLD));//send it out to all the nodes

	return temp;

}

#endif //GLOBALOPS_H
