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
        File to convert the given COO matrix to a user-defined matrix type 
	COO: Coordinate format set storage in fvsolution to 1 
	CSR: Compressed Sparse Row format set storage in fvsolution to 2
	DIA: Diagonal format set storage in fvsolution to 3
	ELL: ELLPACK/ITPACK format set storage in fvsolution to 4 
	HYB: Hybrid format set storage in fvsolution to 5       
                                                         
\**********************************************************************/

//this was turned into a header

if(OFSP->GPUStorage == 2){//COO=1 CSR=2 DIA=3 ELL=4 HYB=5
	cusp::csr_matrix<IndexType, ValueType, MemorySpace> temp(A);
	A = temp;
	if(OFSP->debugCusp){
		std::cout<<"   Using CSR sparse matrix storage\n";
	}
}else if(OFSP->GPUStorage == 3){
	cusp::dia_matrix<IndexType, ValueType, MemorySpace> temp(A);
	A = temp;
	if(OFSP->debugCusp){
		std::cout<<"   Using DIA sparse matrix storage\n";
	}
}else if(OFSP->GPUStorage == 4){
	cusp::ell_matrix<IndexType, ValueType, MemorySpace> temp(A);
	A = temp;
	if(OFSP->debugCusp){
		std::cout<<"   Using ELL sparse matrix storage\n";
	}
}else if(OFSP->GPUStorage == 5){
	cusp::hyb_matrix<IndexType, ValueType, MemorySpace> temp(A);
	A = temp;
	if(OFSP->debugCusp){
		std::cout<<"   Using HYB sparse matrix storage\n";
	}
}else {
	if(OFSP->debugCusp)
		std::cout<<"   Using COO sparse matrix storage\n";
}
