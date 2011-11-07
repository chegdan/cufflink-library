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
        used inside while loop of PBiCGStab
                                                            
\**********************************************************************/

//BEGIN INTERFACES
//taking into account the interface influence

	//cusp::array1d< ValueType, hostMemorySpace > Xh(X.size(),0); //not needed in redefinition of AXLoop

	thrust::copy(X.begin(),X.end(),Xh.begin());//here to reduce redendant operations to save some time


for(int j = 0;j<CFLInterfaces.nParInterfaces;j++)
{
	MPI_Status sts;

	cusp::array1d< ValueType, MemorySpace > Xj(CFLInterfaces.nColsInterface[j],0);//must remain here since device memory is a lot smaller.
	//cusp::array1d< ValueType, hostMemorySpace > CFLInterfaces.Xjh[j](CFLInterfaces.nColsInterface[j],0);
	//cusp::array1d< ValueType, hostMemorySpace > Xh(X.size(),0);

	if(CFLInterfaces.neighbProcNo[j]>CFLInterfaces.myThreadNumber){

		//thrust::copy(X.begin(),X.end(),Xh.begin());

		//send my x vector to my neighbor
		MPI_CHECK(MPI_Send(&Xh[0],Xh.size(),MPI_SCALAR,CFLInterfaces.neighbProcNo[j],0,MPI_COMM_WORLD));

		//recieve my neighbors x vector
		MPI_CHECK(MPI_Recv(&CFLInterfaces.Xjh[j][0],CFLInterfaces.nColsInterface[j],MPI_SCALAR,CFLInterfaces.neighbProcNo[j],0,MPI_COMM_WORLD,&sts));

		//CFLInterfaces.Xj[j] = CFLInterfaces.Xjh[j];
		thrust::copy(CFLInterfaces.Xjh[j].begin(),CFLInterfaces.Xjh[j].end(),Xj.begin());

	}else{
		//recieve my neighbors x vector
		MPI_CHECK(MPI_Recv(&CFLInterfaces.Xjh[j][0],CFLInterfaces.nColsInterface[j],MPI_SCALAR,CFLInterfaces.neighbProcNo[j],0,MPI_COMM_WORLD,&sts));
		//CFLInterfaces.Xj[j] = CFLInterfaces.Xjh[j];
		thrust::copy(CFLInterfaces.Xjh[j].begin(),CFLInterfaces.Xjh[j].end(),Xj.begin());

		//thrust::copy(X.begin(),X.end(),Xh.begin());

		//send my x vector to my neighbor
		MPI_CHECK(MPI_Send(&Xh[0],Xh.size(),MPI_SCALAR,CFLInterfaces.neighbProcNo[j],0,MPI_COMM_WORLD));
	}

	//need a yTemp
	cusp::array1d<ValueType,MemorySpace> yTemp(CFLInterfaces.nRowsInterface,0);

	//perform yTemp = Aij[j]*Xj[j]
	cusp::multiply(CFLInterfaces.Aij[j], Xj,yTemp);

	//y = y - yTemp
	cusp::blas::axpy(yTemp, y, ValueType(-1));
}
//done with interface influence
//END INTERFACES
