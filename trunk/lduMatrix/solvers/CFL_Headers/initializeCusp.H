cusp_equation_system CES;
CES.A = cusp::coo_matrix<IndexType, ValueType, hostMemorySpace>(nCells,nCells,nCells+2*NFaces);
CES.X = cusp::array1d< ValueType, hostMemorySpace>(nCells);
CES.B = cusp::array1d< ValueType, hostMemorySpace>(nCells);

//copy values of lduMatrix diag to A COO matrix
thrust::copy(matrix().diag().begin(),
			 matrix().diag().end(),
			 CES.A.values.begin());
//copy values of lduMatrix lower to A COO matrix
thrust::copy(matrix().lower().begin(),
			 matrix().lower().end(),
			 CES.A.values.begin()+nCells);
//since matrix is symmetric, do not copy upper values to save time-> performed on GPU
//copy row indices of lower into A COO matrix
thrust::copy(matrix().lduAddr().upperAddr().begin(),
			 matrix().lduAddr().upperAddr().end(),
			 CES.A.row_indices.begin()+nCells);
//copy column indices of lower into A COO matrix
thrust::copy(matrix().lduAddr().lowerAddr().begin(),
			 matrix().lduAddr().lowerAddr().end(),
			 CES.A.column_indices.begin()+nCells);

//do not initialize the row and column values of diag to save time-> performed on GPU
thrust::copy(x.begin(),x.end(),CES.X.begin());//copy x of lower into x vector
thrust::copy(b.begin(),b.end(),CES.B.begin());//copy b of lower into b vector
//end COO, X, b initialize and fill

OFSP.iRes = -1;//initialization value for initial residual
OFSP.fRes = -1;//initialization value for final residual
OFSP.nIterations = 0;//initialization value for number of iterations
OFSP.minIter = minIter();//minimum iterations
OFSP.maxIter = maxIter();//maximum iterations
OFSP.relTol = relTolerance();//relative tolerance
OFSP.tol = tolerance();//tolerance
