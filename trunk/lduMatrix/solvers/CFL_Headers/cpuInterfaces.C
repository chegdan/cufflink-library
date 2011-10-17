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
	cpuInterfaces

    Author
        Daniel P. Combest.  All rights reserved.

    Description
        Collection of dynamic arrays to hold interfaces from cpu to be 
	passed to the GPU for use in a parallel implementation of the 
	Krylov space solvers. 

    Source File
	cpuInterfaces.C
                                                            
\**********************************************************************/

#include "cpuInterfaces.H"

cpuInterfaces::cpuInterfaces(){

}

cpuInterfaces::~cpuInterfaces(){

}

void cpuInterfaces::printShortInfo(){
	int j;

	std::cout<<"From the CPU associated with rank "<<myThreadNumber<<" with "<<nParInterfaces<<" interface(s) \n";

	for(j = 0; j<nParInterfaces;j++){

		std::cout<<"node (i) = "<< myThreadNumber <<" interfacing with node (j) = "
   		    << neighbProcNo[j] <<" nRows = " 
                    << nRowsInterface <<" nCols = "<< nColsInterface[j] <<" nnz = "<<nnz[j]
		    <<" first three values "<< Aij[j].values[0] <<" "<< Aij[j].values[1] <<" "<< Aij[j].values[2]
		    <<" first three row values "<<Aij[j].row_indices[0]<<" "<<Aij[j].row_indices[1]<<" "<<Aij[j].row_indices[2]
		    <<" first three col values "<<Aij[j].column_indices[0]<<" "<<Aij[j].column_indices[1]<<" "<<Aij[j].column_indices[2]<< "\n";
		}
}//end printShortInfo()

bool cpuInterfaces::checkGPUCount(int gpusPerMachine){

	int GPU_N;
	cudaGetDeviceCount(&GPU_N);

	//if the gpu count on this node does not equal what was given in fvSolution, return false
	return ((gpusPerMachine>GPU_N || gpusPerMachine < 1)? false : true);

}

