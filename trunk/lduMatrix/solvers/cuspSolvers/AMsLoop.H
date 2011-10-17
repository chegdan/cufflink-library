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

//BEGIN INTERFACES
//taking into account the interface influence
	for(int j = 0;j<CFLInterfaces.nParInterfaces;j++){

		MPI_Status sts;

		cusp::array1d< ValueType, MemorySpace > Msj(CFLInterfaces.nColsInterface[j],0);
		cusp::array1d< ValueType, hostMemorySpace > Msjh(CFLInterfaces.nColsInterface[j],0);
		cusp::array1d< ValueType, hostMemorySpace > Msh(Ms.size(),0);

		if(CFLInterfaces.neighbProcNo[j]>CFLInterfaces.myThreadNumber){

			thrust::copy(Ms.begin(),Ms.end(),Msh.begin());
	
			//send my x vector to my neighbor
			MPI_CHECK(MPI_Send(&Msh[0],Msh.size(),MPI_SCALAR,CFLInterfaces.neighbProcNo[j],0,MPI_COMM_WORLD));

			//recieve my neighbors x vector
			MPI_CHECK(MPI_Recv(&Msjh[0],CFLInterfaces.nColsInterface[j],MPI_SCALAR,CFLInterfaces.neighbProcNo[j],0,MPI_COMM_WORLD,&sts));

			//CFLInterfaces.Msj[j] = Msjh;
			thrust::copy(Msjh.begin(),Msjh.end(),Msj.begin());

		}else{
			//recieve my neighbors x vector
			MPI_CHECK(MPI_Recv(&Msjh[0],CFLInterfaces.nColsInterface[j],MPI_SCALAR,CFLInterfaces.neighbProcNo[j],0,MPI_COMM_WORLD,&sts));
			//CFLInterfaces.Msj[j] = Msjh;
			thrust::copy(Msjh.begin(),Msjh.end(),Msj.begin());

			thrust::copy(Ms.begin(),Ms.end(),Msh.begin());

			//send my x vector to my neighbor
			MPI_CHECK(MPI_Send(&Msh[0],Msh.size(),MPI_SCALAR,CFLInterfaces.neighbProcNo[j],0,MPI_COMM_WORLD));
		}

		//need a yTemp
		cusp::array1d<ValueType,MemorySpace> yTemp(CFLInterfaces.nRowsInterface,0);
	
		//perform yTemp = Aij[j]*Msj[j]
		cusp::multiply(CFLInterfaces.Aij[j], Msj,yTemp);

		//y = y - yTemp
		cusp::blas::axpy(yTemp, AMs, ValueType(-1));
	}
//done with interface influence
//END INTERFACES
