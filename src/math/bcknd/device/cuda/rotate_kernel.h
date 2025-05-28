#ifndef __ROTATE_KERNEL_H__
#define __ROTATE_KERNEL_H__
template< typename T, const int LX, const int CHUNKS >
__global__ void rotate_cyc_kernel(
			   T * __restrict__ rx,
			   T * __restrict__ ry,
			   T * __restrict__ rz,
			   const T * __restrict__ x,
			   const T * __restrict__ y,
			   const T * __restrict__ z,
			   const T * __restrict__ cyclic_angle,
			   const int idir) 
	{ 

	int i,j,k;
	//int j, k;
	const int e = blockIdx.x;
	const int iii = threadIdx.x;
	const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;
	const int ele = e*LX*LX*LX;
	
	__shared__ T shu[LX * LX * LX];
	__shared__ T shv[LX * LX * LX];
	//__shared__ T shw[LX * LX * LX];
	__shared__ T sht[LX * LX * LX];

	j = iii;
	while(j < (LX * LX * LX)) {
		shu[j] = rx[j + e * LX * LX * LX];
		shv[j] = ry[j + e * LX * LX * LX];
		//shw[j] = rz[j + e * LX * LX * LX];
		sht[j] = cyclic_angle[j + e * LX * LX * LX];
		j = j + CHUNKS;
	}
  
  __syncthreads();	


	for (int n = 0; n < nchunks; n++) {
		const int ijk = iii + n * CHUNKS;
		const int jk = ijk / LX;
		i = ijk - jk * LX;
		k = jk / LX;
		j = jk - k * LX;

		T theta = 0.0;
		T cost = 0.0;
		T sint = 0.0;
		T rnor = 0.0;
		T rtn1 = 0.0;
		T dotprod = 0.0;
		T tol = 1e-15;
		theta = sht[ijk];
        //can copy u, v to __shared__ shu ...
		if (fabs(theta)>tol){ // warp divergence here?
			cost = cosf(theta);
			sint = sinf(theta);
			dotprod = y[ijk + ele] * cost - x[ijk + ele] * sint;

			if (idir==1){
				rnor = shu[ijk] * cost + shv[ijk] * sint;
				rtn1 =-shu[ijk] * sint + shv[ijk] * cost;
			}	
			else if (idir==0){
				rnor = shu[ijk] * cost - shv[ijk] * sint;
				rtn1 = shu[ijk] * sint + shv[ijk] * cost;
			}
			//can avoid if block (warp divergence) with idir=1 and -1 for forward and backward ops
			// rnor = shu[ijk] * cost + shv[ijk] * sint*idir;
			// rtn1 = shu[ijk] * sint*idir + shv[ijk] * cost;
			if (dotprod>=0){
				rx[ijk + ele] = rnor;
				ry[ijk + ele] = rtn1;
			}
			else{
				rx[ijk + ele] = -rnor;
				ry[ijk + ele] = -rtn1;
			}
		}
	  } 
	}

#endif // __ROTATE_KERNEL_H__
