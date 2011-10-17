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
        Struct containing solver performance information that is passed 
	back to OpenFOAM.
                                                             
\**********************************************************************/


struct OFSolverPerformance{

	int nCells;
	int nFaces;
	
	int minIter;
	int maxIter;
	double relTol;
	double tol;
	
	//return values for openfoam solverPerformance object
	int nIterations;
	double iRes;
	double fRes;
	bool converged;
	bool singular;
	bool debugCusp;
	
	//some values to return
	double aTol;
	double aRelTol;

	int GPUStorage;
	int deviceNumber;
	int preconditionerType;

	//for Ainv preconditioner
	int AinvType; //0 (standard) (defualt),1 (static),2 (novel)
	double dropTolerance;
	int nonZeros;
	int linStrategy;
};


