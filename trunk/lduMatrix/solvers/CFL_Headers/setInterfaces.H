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
        Naive code to convert openfoam interfaces to interfaces in device memory
                                                            
\**********************************************************************/

//#include "cpugpuInterfaces.H"


gpuInterfaces CFLInterfaces;
CFLInterfaces.myThreadNumber 	= OFInterfaces->myThreadNumber;
CFLInterfaces.nParInterfaces 	= OFInterfaces->nParInterfaces;
CFLInterfaces.nRowsInterface 	= OFInterfaces->nRowsInterface;
//CFLInterfaces.neighbProcNo 	= new int[CFLInterfaces.nParInterfaces];
//CFLInterfaces.nColsInterface 	= new int[CFLInterfaces.nParInterfaces];
//CFLInterfaces.nnz 		= new int[CFLInterfaces.nParInterfaces];
CFLInterfaces.Aij 		= new cusp::coo_matrix<IndexType, ValueType, MemorySpace>[CFLInterfaces.nParInterfaces];
CFLInterfaces.Xj		= new cusp::array1d< ValueType, MemorySpace>[CFLInterfaces.nParInterfaces];
CFLInterfaces.Bj		= new cusp::array1d< ValueType, MemorySpace>[CFLInterfaces.nParInterfaces];


//fill the CFLInterfaces Aij matrices
int j;
for( j=0;j<CFLInterfaces.nParInterfaces;j++){
	//fill the CFLInterfaces Aij matrices
	//CFLInterfaces.Aij.push_back(cusp::coo_matrix<IndexType, ValueType, MemorySpace>(OFInterfaces->Aij[j]));//errors on this line. submitted to thrust users board

	CFLInterfaces.Aij[j] = cusp::coo_matrix<IndexType, ValueType, MemorySpace>(OFInterfaces->Aij[j]);//dynamic array using new

	CFLInterfaces.neighbProcNo.push_back(OFInterfaces->neighbProcNo[j]);//vector class
	CFLInterfaces.nColsInterface.push_back(OFInterfaces->nColsInterface[j]);//vector class
	CFLInterfaces.nnz.push_back(OFInterfaces->nnz[j]);//vector class

	//fill the CFLInterfaces Xj vectors
	//CFLInterfaces.Xj[j].push_back(cusp::array1d< ValueType, MemorySpace>(OFInterfaces->Xj[j]));

	//fill the CFLInterfaces Bj vectors
	//CFLInterfaces.Bj[j].push_back(cusp::array1d< ValueType, MemorySpace>(OFInterfaces->Bj[j]));
}










