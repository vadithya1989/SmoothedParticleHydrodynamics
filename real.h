#ifndef __REAL_H__
#define	__REAL_H_

// - Single Precision
#ifdef	__OPENCL_VERSION__
typedef float real;
typedef float4 real4;
#else
typedef cl_float cl_real;
typedef cl_float4 cl_real4;
#endif

/*// - Double Precision
#ifdef	__OPENCL_VERSION__
#pragma	OPENCL EXTENSION cl_khr_fp64 : enable
typedef double real;
typedef double4 real4;
#else
typedef cl_double cl_real;
typedef cl_double4 cl_real4;
#endif
*/

#endif
