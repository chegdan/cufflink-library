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
        exampleClass

    Author
        Daniel P. Combest.  All rights reserved.

    Description
        Method to get cuspPreconditioner name from fvSolution file and 
	pass to solver.  This file is obsolete and not used by the 
	current solvers since it uses string comparison.
                                                             
\**********************************************************************/


word GPUPreconditionerType = dict().lookup("cuspPreconditioner");
OFSP.preconditionerType = 4;//default

     if ("Cusp_Ainv" == GPUPreconditionerType){
	OFSP.preconditionerType = 1;
	//OFSP.AinvType = dict().lookupOrDefault<int>("AinvType", 0);
	//OFSP.dropTolerance = dict().lookupOrDefault<double>("dropTolerance", 0.1);
	//OFSP.nonZeros = dict().lookupOrDefault<int>("nonZeros", 10);
	OFSP.linStrategy = dict().lookupOrDefault<int>("linStrategy", 2);

}else if("Cusp_Diagonal" == GPUPreconditionerType)
	OFSP.preconditionerType = 2;
else if ("Cusp_None" == GPUPreconditionerType)
	OFSP.preconditionerType = 3;
else if("Cusp_SmA" == GPUPreconditionerType)
	OFSP.preconditionerType = 4;
else{
FatalErrorIn("cuspPreconditioner not defined correctly in fvSolution.  No such preconditioner. choose from \n4\n(\n Cusp_Ainv\n Cusp_Diagonal\n Cusp_None\n Cusp_SmA\n)")<< "Fatal Error"<< abort(FatalError);
}
