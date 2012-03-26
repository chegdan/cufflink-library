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
        Naive code to grab openfoam interfaces
                                                            
\**********************************************************************/

cpuInterfaces OFInterfaces;//create container for OpenFOAM interfaces

OFInterfaces.myThreadNumber = Pstream::myProcNo();//get my thread number before I forget.

    OFInterfaces.nParInterfaces = 0;//number of parallel interfaces on this processor initializer

    forAll(interfaces(), patchI)
    {
	if (interfaces().set(patchI)){
	//count the processor interfaces on this processor
       if(interfaces()[patchI].interface().myProcNo()!=(-1)){ OFInterfaces.nParInterfaces++; }
	}
    }
   	int interfaceCounter = 0;

    forAll(interfaces(), patchI)
    {
	if (interfaces().set(patchI)){
		//get patch id of the interface and store in addressing array.  only store parallel interfaces
       		if(interfaces()[patchI].interface().myProcNo()!=(-1)){ 
			OFInterfaces.interfaceAddr.push_back(patchI);
			OFInterfaces.neighbProcNo.push_back(interfaces()[patchI].interface().neighbProcNo());//array to hold neighbor processor ids
			interfaceCounter++; 
			}
	}
    }

	
	OFInterfaces.nRowsInterface = nCells;//must always be this size to produce a vector the length of the x and b on this processor
	

    for(int j = 0; j<OFInterfaces.nParInterfaces;j++){
            const unallocLabelList& pa = matrix().lduAddr().patchAddr(OFInterfaces.interfaceAddr[j]);//row addressing of the interface matrix
	    OFInterfaces.nnz.push_back(pa.size());//store the number of non-zeros of the interface

		//send the size of the x vector here. example in combineGatherScatter.C
		OPstream::write
                (
                    Pstream::blocking,//Pstream::scheduled,//only works correctly with scheduled
                    interfaces()[OFInterfaces.interfaceAddr[j]].interface().neighbProcNo(),
                    reinterpret_cast<const char*>(&OFInterfaces.nRowsInterface),
                    sizeof(OFInterfaces.nRowsInterface)
                );
		
		int tempCols;
		//read the value from neighbors. example in combineGatherScatter.C
                IPstream::read
                (
                    Pstream::blocking,//Pstream::scheduled,//only works correctly with scheduled
                    interfaces()[OFInterfaces.interfaceAddr[j]].interface().neighbProcNo(),
                    reinterpret_cast<char*>(&tempCols),
                    sizeof(int)
                );
		
		//add the value of the columns for interface j
		OFInterfaces.nColsInterface.push_back(tempCols);

	    	//declare COO matrix now that all sizes are known	
	        OFInterfaces.Aij.push_back(cusp::coo_matrix<IndexType, ValueType, hostMemorySpace>(OFInterfaces.nRowsInterface,OFInterfaces.nColsInterface[j],OFInterfaces.nnz[j]));

	    //copy the row indices to the interface COO matrix
		thrust::copy(pa.begin(),pa.end(), OFInterfaces.Aij[j].row_indices.begin());

    	   //send the row indices to the neighbor that are used as neighbor column indices. example in combineGatherScatter.C
			OPstream::write
         	       (
       		        Pstream::blocking,//Pstream::scheduled,
       		        interfaces()[OFInterfaces.interfaceAddr[j]].interface().neighbProcNo(),
       		        reinterpret_cast<const char*>(&pa[0]),
       		        OFInterfaces.nnz[j]*sizeof(int)
       		       );

	  //read the column indices from neighbors. example in combineGatherScatter.C
               	       IPstream::read
                	(
                    	Pstream::blocking,//Pstream::scheduled,
                   	interfaces()[OFInterfaces.interfaceAddr[j]].interface().neighbProcNo(),
                    	reinterpret_cast<char*>(&OFInterfaces.Aij[j].column_indices[0]),
                   	OFInterfaces.nnz[j]*sizeof(int)
                	);


	  //get a reference to the correct interface matrix values
            const scalarField& pCoeffs = coupleBouCoeffs_[OFInterfaces.interfaceAddr[j]] ;

	 //copy the values of the interface matrix to the interface matrix itself
		thrust::copy(pCoeffs.begin(),pCoeffs.end(), OFInterfaces.Aij[j].values.begin());
   
		 if (lduMatrix::debug >= 3)
    		{
			OFInterfaces.printShortInfo();
		}

    }//end for loop
//parallel interfaces - End

