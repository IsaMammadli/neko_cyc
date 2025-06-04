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

	j = iii;
	while( j < (LX * LX * LX)) {
		const int i = j + e * LX * LX * LX;
		const T theta = cyclic_angle[i];
		T cost=0;
		T sint=0;
		T dotprod=0;
		T rnor = 0;
		T rtn1 = 0;
		const T tol = 1e-16;
		//if (fabs(theta)>0.0){
			cost = cos(theta);
	 		sint = sin(theta);
			dotprod = y[i] * cost - x[i] * sint;

			if (idir==1){
				rnor = rx[i] * cost + ry[i] * sint;
				rtn1 =-rx[i] * sint + ry[i] * cost;
			}	
			else if (idir==0){
				rnor = rx[i] * cost - ry[i] * sint;
				rtn1 = rx[i] * sint + ry[i] * cost;
			}
			if (dotprod>=0){
				rx[i] = rnor;
				ry[i] = rtn1;
			}
			else{
				rx[i] = -rnor;
				ry[i] = -rtn1;
			}
		//}
		j = j + CHUNKS;
		}
		
	//__syncthreads(); //needed for shared memory only
	}
#endif // __ROTATE_KERNEL_H__

	
// 	__shared__ T shu[LX * LX * LX];
// 	__shared__ T shv[LX * LX * LX];
// 	//__shared__ T shw[LX * LX * LX];
// 	__shared__ T sht[LX * LX * LX];

// 	j = iii;
// 	while(j < (LX * LX * LX)) {
// 		shu[j] = rx[j + e * LX * LX * LX];
// 		shv[j] = ry[j + e * LX * LX * LX];
// 		//shw[j] = rz[j + e * LX * LX * LX];
// 		sht[j] = cyclic_angle[j + e * LX * LX * LX];
// 		j = j + CHUNKS;
// 	}
  
//   __syncthreads();	






// for (int n = 0; n < nchunks; n++) {
// 	const int ijk = iii + n * CHUNKS;
// 	const int jk = ijk / LX;
// 	i = ijk - jk * LX;
// 	k = jk / LX;
// 	j = jk - k * LX;

// 	T theta = 0.0;
// 	T cost = 0.0;
// 	T sint = 0.0;
// 	T rnor = 0.0;
// 	T rtn1 = 0.0;
// 	T dotprod = 0.0;
// 	T tol = 1e-16;
// 	theta = cyclic_angle[ijk + e * LX * LX * LX];
// 	//can copy u, v to __shared__ shu ...
// 	if (fabs(theta)>tol){ // warp divergence here?
// 		cost = cos(theta);
// 		sint = sin(theta);
// 		dotprod = y[ijk + e * LX * LX * LX] * cost - x[ijk + e * LX * LX * LX] * sint;

// 		if (idir==1){
// 			rnor = rx[ijk + e * LX * LX * LX] * cost + ry[ijk + e * LX * LX * LX] * sint;
// 			rtn1 =-rx[ijk + e * LX * LX * LX] * sint + ry[ijk + e * LX * LX * LX] * cost;
// 		}	
// 		else if (idir==0){
// 			rnor = rx[ijk + e * LX * LX * LX] * cost - ry[ijk + e * LX * LX * LX] * sint;
// 			rtn1 = rx[ijk + e * LX * LX * LX] * sint + ry[ijk + e * LX * LX * LX] * cost;
// 		}
// 		if (dotprod>=0){
// 			rx[ijk + e * LX * LX * LX] = rnor;
// 			ry[ijk + e * LX * LX * LX] = rtn1;
// 		}
// 		else{
// 			rx[ijk + e * LX * LX * LX] = -rnor;
// 			ry[ijk + e * LX * LX * LX] = -rtn1;
// 		}
// 	}
//   } 