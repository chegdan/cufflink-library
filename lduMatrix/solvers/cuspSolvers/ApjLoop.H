//BEGIN INTERFACES
//taking into account the interface influence
	for(int j = 0;j<CFLInterfaces.nParInterfaces;j++){
		MPI_Status sts;

		cusp::array1d< ValueType, MemorySpace > pj(CFLInterfaces.nColsInterface[j],0);
		cusp::array1d< ValueType, hostMemorySpace > pjh(CFLInterfaces.nColsInterface[j],0);
		cusp::array1d< ValueType, hostMemorySpace > ph(p.size(),0);

		if(CFLInterfaces.neighbProcNo[j]>CFLInterfaces.myThreadNumber){

			thrust::copy(p.begin(),p.end(),ph.begin());
	
			//send my x vector to my neighbor
			MPI_CHECK(MPI_Send(&ph[0],ph.size(),MPI_SCALAR,CFLInterfaces.neighbProcNo[j],0,MPI_COMM_WORLD));

			//recieve my neighbors x vector
			MPI_CHECK(MPI_Recv(&pjh[0],CFLInterfaces.nColsInterface[j],MPI_SCALAR,CFLInterfaces.neighbProcNo[j],0,MPI_COMM_WORLD,&sts));

			//CFLInterfaces.Xj[j] = Xjh;
			thrust::copy(pjh.begin(),pjh.end(),pj.begin());

		}else{
			//recieve my neighbors x vector
			MPI_CHECK(MPI_Recv(&pjh[0],CFLInterfaces.nColsInterface[j],MPI_SCALAR,CFLInterfaces.neighbProcNo[j],0,MPI_COMM_WORLD,&sts));
			//CFLInterfaces.Xj[j] = Xjh;
			thrust::copy(pjh.begin(),pjh.end(),pj.begin());

			thrust::copy(p.begin(),p.end(),ph.begin());

			//send my x vector to my neighbor
			MPI_CHECK(MPI_Send(&ph[0],ph.size(),MPI_SCALAR,CFLInterfaces.neighbProcNo[j],0,MPI_COMM_WORLD));
		}

		//need a yTemp
		cusp::array1d<ValueType,MemorySpace> yTemp(CFLInterfaces.nRowsInterface,0);
	
		//perform yTemp = Aij[j]*Xj[j]
		cusp::multiply(CFLInterfaces.Aij[j], pj,yTemp);

		//Apj = Apj - yTemp
		cusp::blas::axpy(yTemp, Ap, ValueType(-1));
	}
//done with interface influence
//END INTERFACES
