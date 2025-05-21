#ifndef __ROTATE_KERNEL_H__
#define __ROTATE_KERNEL_H__
template< typename T, const int LX, const int CHUNKS >
__global__ void rotate_cyc_kernel(
			   T * __restrict__ u,
			   T * __restrict__ v,
			   T * __restrict__ w,
			   const T * __restrict__ x,
			   const T * __restrict__ y,
			   const T * __restrict__ z,
			   const T * __restrict__ cyclic_angle,
			   const T * __restrict__ idir) 
	{ 

	int i,j,k;

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
		shu[j] = u[j + e * LX * LX * LX];
		shv[j] = v[j + e * LX * LX * LX];
		//shw[j] = w[j + e * LX * LX * LX];
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

		const int idir_val = *idir;
		T theta = 0.0;
		T cost = 0.0;
		T sint = 0.0;
		T rnor = 0.0;
		T rtn1 = 0.0;
		T dotprod = 0.0;
		T tol = 1e-15;
		theta = sht[ijk];
        //can copy u, v to __shared__ shu ...
		if (fabs(theta)>tol){
			cost = cosf(theta);
			sint = sinf(theta);
			dotprod = y[ijk + ele] * cost - x[ijk + ele] * sint;

			if (idir_val==1){
				rnor = shu[ijk] * cost + shv[ijk] * sint;
				rtn1 =-shu[ijk] * sint + shv[ijk] * cost;
			}	
			else if (idir_val==0){
				rnor = shu[ijk] * cost - shv[ijk] * sint;
				rtn1 = shu[ijk] * sint + shv[ijk] * cost;
			}
			//can avoid if block (warp divergence) with idir=1 and -1 for forward and backward ops
			// rnor = shu[ijk] * cost + shv[ijk] * sint*idir;
			// rtn1 = shu[ijk] * sint*idir + shv[ijk] * cost;
			if (dotprod>=0){
				u[ijk + ele] = rnor;
				v[ijk + ele] = rtn1;
			}
			else{
				u[ijk + ele] = -rnor;
				v[ijk + ele] = -rtn1;
			}
		}
	  } 
	}

#endif // __ROTATE_KERNEL_H__
			
// 	  do i = 1, ntot
// 	  theta = coef%cyc_angle(i, 1, 1, 1)
// 	  if(abs(theta) .gt. tol) then
// 		   x = coef%dof%x(i, 1, 1, 1)
// 		   y = coef%dof%y(i, 1, 1, 1)
// 		   z = coef%dof%z(i, 1, 1, 1)
// 		   cost = cos(theta)
// 		   sint = sin(theta)
// 		   dotprod = cost*y-sint*x
// 		   !can avoid this if block using idir=1 and -1 mutliplied by sint term
// 		   if (idir.eq.1) then
// 				rnor = rx(i, 1, 1, 1)*cost+ry(i, 1, 1, 1)*sint
// 				rtn1 =-rx(i, 1, 1, 1)*sint+ry(i, 1, 1, 1)*cost
// 		   else if(idir.eq.0) then
// 				rnor = rx(i, 1, 1, 1)*cost-ry(i, 1, 1, 1)*sint
// 				rtn1 = rx(i, 1, 1, 1)*sint+ry(i, 1, 1, 1)*cost
// 		   else
// 				write(*, *) "idir must be 0 or 1 in rotate_cyc"
// 				!stop !must be removed if do concurrent is used.
// 		   end if

// 		   !invert direction if vector comes out of plane (i.e. dotprod<0)
// 		   if (dotprod.ge.0) then
// 				rx(i, 1, 1, 1) = rnor
// 				ry(i, 1, 1, 1) = rtn1
// 		   else
// 				rx(i, 1, 1, 1) =-rnor
// 				ry(i, 1, 1, 1) =-rtn1
// 		   end if
// 	  end if 
//  end do