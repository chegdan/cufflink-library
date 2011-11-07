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
        The first loop involving the interfaces for PBiCGStab
                                                            
\**********************************************************************/

	thrust::copy(Mp.begin(),Mp.end(),Mph.begin());//out here to cut down on redundant processes

//BEGIN INTERFACES
//taking into account the interface influence
	for(int j = 0;j<CFLInterfaces.nParInterfaces;j++){

		MPI_Status sts;

		cusp::array1d< ValueType, MemorySpace > Mpj(CFLInterfaces.nColsInterface[j],0);
		//cusp::array1d< ValueType, hostMemorySpace > Mpjh[j](CFLInterfaces.nColsInterface[j],0);
		//cusp::array1d< ValueType, hostMemorySpace > Mph(Mp.size(),0);

		if(CFLInterfaces.neighbProcNo[j]>CFLInterfaces.myThreadNumber){

			//thrust::copy(Mp.begin(),Mp.end(),Mph.begin());
	
			//send my Mp vector to my neighbor
			MPI_CHECK(MPI_Send(&Mph[0],Mph.size(),MPI_SCALAR,CFLInterfaces.neighbProcNo[j],0,MPI_COMM_WORLD));

			//recieve my neighbors Mp vector
			MPI_CHECK(MPI_Recv(&Mpjh[j][0],CFLInterfaces.nColsInterface[j],MPI_SCALAR,CFLInterfaces.neighbProcNo[j],0,MPI_COMM_WORLD,&sts));

			//CFLInterfaces.Mpj[j] = Mpjh[j];
			thrust::copy(Mpjh[j].begin(),Mpjh[j].end(),Mpj.begin());

		}else{
			//recieve my neighbors Mp vector
			MPI_CHECK(MPI_Recv(&Mpjh[j][0],CFLInterfaces.nColsInterface[j],MPI_SCALAR,CFLInterfaces.neighbProcNo[j],0,MPI_COMM_WORLD,&sts));
			//CFLInterfaces.Mpj[j] = Mpjh[j];
			thrust::copy(Mpjh[j].begin(),Mpjh[j].end(),Mpj.begin());

			//thrust::copy(Mp.begin(),Mp.end(),Mph.begin());

			//send my Mp vector to my neighbor
			MPI_CHECK(MPI_Send(&Mph[0],Mph.size(),MPI_SCALAR,CFLInterfaces.neighbProcNo[j],0,MPI_COMM_WORLD));
		}

		//need a yTemp
		cusp::array1d<ValueType,MemorySpace> yTemp(CFLInterfaces.nRowsInterface,0);
	
		//perform yTemp = Aij[j]*Mpj[j]
		cusp::multiply(CFLInterfaces.Aij[j], Mpj,yTemp);

		//y = y - yTemp
		cusp::blas::axpy(yTemp, AMp, ValueType(-1));
	}
//done with interface influence
//END INTERFACES
