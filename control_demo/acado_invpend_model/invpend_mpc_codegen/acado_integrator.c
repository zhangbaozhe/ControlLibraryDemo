/*
 *    This file was auto-generated using the ACADO Toolkit.
 *    
 *    While ACADO Toolkit is free software released under the terms of
 *    the GNU Lesser General Public License (LGPL), the generated code
 *    as such remains the property of the user who used ACADO Toolkit
 *    to generate this code. In particular, user dependent data of the code
 *    do not inherit the GNU LGPL license. On the other hand, parts of the
 *    generated code that are a direct copy of source code from the
 *    ACADO Toolkit or the software tools it is based on, remain, as derived
 *    work, automatically covered by the LGPL license.
 *    
 *    ACADO Toolkit is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *    
 */


#include "acado_common.h"


real_t rk_dim8_swap;

/** Column vector of size: 8 */
real_t rk_dim8_bPerm[ 8 ];

/** Column vector of size: 2 */
real_t auxVar[ 2 ];

real_t rk_ttt;

/** Row vector of size: 5 */
real_t rk_xxx[ 5 ];

/** Matrix of size: 4 x 2 (row major format) */
real_t rk_kkk[ 8 ];

/** Matrix of size: 8 x 8 (row major format) */
real_t rk_A[ 64 ];

/** Column vector of size: 8 */
real_t rk_b[ 8 ];

/** Row vector of size: 8 */
int rk_dim8_perm[ 8 ];

/** Column vector of size: 4 */
real_t rk_rhsTemp[ 4 ];

/** Matrix of size: 2 x 20 (row major format) */
real_t rk_diffsTemp2[ 40 ];

/** Matrix of size: 4 x 2 (row major format) */
real_t rk_diffK[ 8 ];

/** Matrix of size: 4 x 5 (row major format) */
real_t rk_diffsNew2[ 20 ];

#pragma omp threadprivate( auxVar, rk_ttt, rk_xxx, rk_kkk, rk_diffK, rk_rhsTemp, rk_dim8_perm, rk_A, rk_b, rk_diffsNew2, rk_diffsTemp2, rk_dim8_swap, rk_dim8_bPerm )

void acado_rhs(const real_t* in, real_t* out)
{
const real_t* xd = in;
const real_t* u = in + 4;

/* Compute outputs: */
out[0] = xd[1];
out[1] = ((real_t)(1.8181818181818183e+00)*((u[0]-((real_t)(1.0000000000000001e-01)*xd[1]))-((real_t)(3.5280000000000006e-02)*xd[2])));
out[2] = xd[3];
out[3] = ((((real_t)(5.8800000000000008e-01)*xd[2])-((real_t)(5.9999999999999998e-02)*((real_t)(1.8181818181818183e+00)*((u[0]-((real_t)(1.0000000000000001e-01)*xd[1]))-((real_t)(3.5280000000000006e-02)*xd[2])))))/(real_t)(2.4000000000000000e-02));
}



void acado_diffs(const real_t* in, real_t* out)
{
const real_t* xd = in;
/* Vector of auxiliary variables; number of elements: 2. */
real_t* a = auxVar;

/* Compute intermediate quantities: */
a[0] = ((real_t)(1.0000000000000000e+00)/(real_t)(2.4000000000000000e-02));
a[1] = ((real_t)(1.0000000000000000e+00)/(real_t)(2.4000000000000000e-02));

/* Compute outputs: */
out[0] = (real_t)(0.0000000000000000e+00);
out[1] = (real_t)(1.0000000000000000e+00);
out[2] = (real_t)(0.0000000000000000e+00);
out[3] = (real_t)(0.0000000000000000e+00);
out[4] = (real_t)(0.0000000000000000e+00);
out[5] = (real_t)(0.0000000000000000e+00);
out[6] = ((real_t)(-1.8181818181818185e-01));
out[7] = ((real_t)(-6.4145454545454561e-02));
out[8] = (real_t)(0.0000000000000000e+00);
out[9] = (real_t)(1.8181818181818183e+00);
out[10] = (real_t)(0.0000000000000000e+00);
out[11] = (real_t)(0.0000000000000000e+00);
out[12] = (real_t)(0.0000000000000000e+00);
out[13] = (real_t)(1.0000000000000000e+00);
out[14] = (real_t)(0.0000000000000000e+00);
out[15] = (real_t)(0.0000000000000000e+00);
out[16] = (((real_t)(0.0000000000000000e+00)-((real_t)(-1.0909090909090910e-02)))*a[0]);
out[17] = (((real_t)(5.8800000000000008e-01)-((real_t)(-3.8487272727272736e-03)))*a[0]);
out[18] = (real_t)(0.0000000000000000e+00);
out[19] = (((real_t)(0.0000000000000000e+00)-((real_t)(1.0909090909090910e-01)))*a[1]);
}



void acado_solve_dim8_triangular( real_t* const A, real_t* const b )
{

b[7] = b[7]/A[63];
b[6] -= + A[55]*b[7];
b[6] = b[6]/A[54];
b[5] -= + A[47]*b[7];
b[5] -= + A[46]*b[6];
b[5] = b[5]/A[45];
b[4] -= + A[39]*b[7];
b[4] -= + A[38]*b[6];
b[4] -= + A[37]*b[5];
b[4] = b[4]/A[36];
b[3] -= + A[31]*b[7];
b[3] -= + A[30]*b[6];
b[3] -= + A[29]*b[5];
b[3] -= + A[28]*b[4];
b[3] = b[3]/A[27];
b[2] -= + A[23]*b[7];
b[2] -= + A[22]*b[6];
b[2] -= + A[21]*b[5];
b[2] -= + A[20]*b[4];
b[2] -= + A[19]*b[3];
b[2] = b[2]/A[18];
b[1] -= + A[15]*b[7];
b[1] -= + A[14]*b[6];
b[1] -= + A[13]*b[5];
b[1] -= + A[12]*b[4];
b[1] -= + A[11]*b[3];
b[1] -= + A[10]*b[2];
b[1] = b[1]/A[9];
b[0] -= + A[7]*b[7];
b[0] -= + A[6]*b[6];
b[0] -= + A[5]*b[5];
b[0] -= + A[4]*b[4];
b[0] -= + A[3]*b[3];
b[0] -= + A[2]*b[2];
b[0] -= + A[1]*b[1];
b[0] = b[0]/A[0];
}

real_t acado_solve_dim8_system( real_t* const A, real_t* const b, int* const rk_perm )
{
real_t det;

int i;
int j;
int k;

int indexMax;

int intSwap;

real_t valueMax;

real_t temp;

for (i = 0; i < 8; ++i)
{
rk_perm[i] = i;
}
det = 1.0000000000000000e+00;
for( i=0; i < (7); i++ ) {
	indexMax = i;
	valueMax = fabs(A[i*8+i]);
	for( j=(i+1); j < 8; j++ ) {
		temp = fabs(A[j*8+i]);
		if( temp > valueMax ) {
			indexMax = j;
			valueMax = temp;
		}
	}
	if( indexMax > i ) {
for (k = 0; k < 8; ++k)
{
	rk_dim8_swap = A[i*8+k];
	A[i*8+k] = A[indexMax*8+k];
	A[indexMax*8+k] = rk_dim8_swap;
}
	rk_dim8_swap = b[i];
	b[i] = b[indexMax];
	b[indexMax] = rk_dim8_swap;
	intSwap = rk_perm[i];
	rk_perm[i] = rk_perm[indexMax];
	rk_perm[indexMax] = intSwap;
	}
	det *= A[i*8+i];
	for( j=i+1; j < 8; j++ ) {
		A[j*8+i] = -A[j*8+i]/A[i*8+i];
		for( k=i+1; k < 8; k++ ) {
			A[j*8+k] += A[j*8+i] * A[i*8+k];
		}
		b[j] += A[j*8+i] * b[i];
	}
}
det *= A[63];
det = fabs(det);
acado_solve_dim8_triangular( A, b );
return det;
}

void acado_solve_dim8_system_reuse( real_t* const A, real_t* const b, int* const rk_perm )
{

rk_dim8_bPerm[0] = b[rk_perm[0]];
rk_dim8_bPerm[1] = b[rk_perm[1]];
rk_dim8_bPerm[2] = b[rk_perm[2]];
rk_dim8_bPerm[3] = b[rk_perm[3]];
rk_dim8_bPerm[4] = b[rk_perm[4]];
rk_dim8_bPerm[5] = b[rk_perm[5]];
rk_dim8_bPerm[6] = b[rk_perm[6]];
rk_dim8_bPerm[7] = b[rk_perm[7]];
rk_dim8_bPerm[1] += A[8]*rk_dim8_bPerm[0];

rk_dim8_bPerm[2] += A[16]*rk_dim8_bPerm[0];
rk_dim8_bPerm[2] += A[17]*rk_dim8_bPerm[1];

rk_dim8_bPerm[3] += A[24]*rk_dim8_bPerm[0];
rk_dim8_bPerm[3] += A[25]*rk_dim8_bPerm[1];
rk_dim8_bPerm[3] += A[26]*rk_dim8_bPerm[2];

rk_dim8_bPerm[4] += A[32]*rk_dim8_bPerm[0];
rk_dim8_bPerm[4] += A[33]*rk_dim8_bPerm[1];
rk_dim8_bPerm[4] += A[34]*rk_dim8_bPerm[2];
rk_dim8_bPerm[4] += A[35]*rk_dim8_bPerm[3];

rk_dim8_bPerm[5] += A[40]*rk_dim8_bPerm[0];
rk_dim8_bPerm[5] += A[41]*rk_dim8_bPerm[1];
rk_dim8_bPerm[5] += A[42]*rk_dim8_bPerm[2];
rk_dim8_bPerm[5] += A[43]*rk_dim8_bPerm[3];
rk_dim8_bPerm[5] += A[44]*rk_dim8_bPerm[4];

rk_dim8_bPerm[6] += A[48]*rk_dim8_bPerm[0];
rk_dim8_bPerm[6] += A[49]*rk_dim8_bPerm[1];
rk_dim8_bPerm[6] += A[50]*rk_dim8_bPerm[2];
rk_dim8_bPerm[6] += A[51]*rk_dim8_bPerm[3];
rk_dim8_bPerm[6] += A[52]*rk_dim8_bPerm[4];
rk_dim8_bPerm[6] += A[53]*rk_dim8_bPerm[5];

rk_dim8_bPerm[7] += A[56]*rk_dim8_bPerm[0];
rk_dim8_bPerm[7] += A[57]*rk_dim8_bPerm[1];
rk_dim8_bPerm[7] += A[58]*rk_dim8_bPerm[2];
rk_dim8_bPerm[7] += A[59]*rk_dim8_bPerm[3];
rk_dim8_bPerm[7] += A[60]*rk_dim8_bPerm[4];
rk_dim8_bPerm[7] += A[61]*rk_dim8_bPerm[5];
rk_dim8_bPerm[7] += A[62]*rk_dim8_bPerm[6];


acado_solve_dim8_triangular( A, rk_dim8_bPerm );
b[0] = rk_dim8_bPerm[0];
b[1] = rk_dim8_bPerm[1];
b[2] = rk_dim8_bPerm[2];
b[3] = rk_dim8_bPerm[3];
b[4] = rk_dim8_bPerm[4];
b[5] = rk_dim8_bPerm[5];
b[6] = rk_dim8_bPerm[6];
b[7] = rk_dim8_bPerm[7];
}



/** Matrix of size: 2 x 2 (row major format) */
static const real_t acado_Ah_mat[ 4 ] = 
{ 2.5000000000000001e-02, 5.3867513459481292e-02, 
-3.8675134594812867e-03, 2.5000000000000001e-02 };


/* Fixed step size:0.1 */
int acado_integrate( real_t* const rk_eta, int resetIntegrator )
{
int error;

int i;
int j;
int k;
int run;
int run1;
int tmp_index1;
int tmp_index2;

real_t det;

rk_ttt = 0.0000000000000000e+00;
rk_xxx[4] = rk_eta[24];

for (run = 0; run < 1; ++run)
{
if( resetIntegrator ) {
for (i = 0; i < 1; ++i)
{
for (run1 = 0; run1 < 2; ++run1)
{
for (j = 0; j < 4; ++j)
{
rk_xxx[j] = rk_eta[j];
tmp_index1 = j;
rk_xxx[j] += + acado_Ah_mat[run1 * 2]*rk_kkk[tmp_index1 * 2];
rk_xxx[j] += + acado_Ah_mat[run1 * 2 + 1]*rk_kkk[tmp_index1 * 2 + 1];
}
acado_diffs( rk_xxx, &(rk_diffsTemp2[ run1 * 20 ]) );
for (j = 0; j < 4; ++j)
{
tmp_index1 = (run1 * 4) + (j);
rk_A[tmp_index1 * 8] = + acado_Ah_mat[run1 * 2]*rk_diffsTemp2[(run1 * 20) + (j * 5)];
rk_A[tmp_index1 * 8 + 1] = + acado_Ah_mat[run1 * 2]*rk_diffsTemp2[(run1 * 20) + (j * 5 + 1)];
rk_A[tmp_index1 * 8 + 2] = + acado_Ah_mat[run1 * 2]*rk_diffsTemp2[(run1 * 20) + (j * 5 + 2)];
rk_A[tmp_index1 * 8 + 3] = + acado_Ah_mat[run1 * 2]*rk_diffsTemp2[(run1 * 20) + (j * 5 + 3)];
if( 0 == run1 ) rk_A[(tmp_index1 * 8) + (j)] -= 1.0000000000000000e+00;
rk_A[tmp_index1 * 8 + 4] = + acado_Ah_mat[run1 * 2 + 1]*rk_diffsTemp2[(run1 * 20) + (j * 5)];
rk_A[tmp_index1 * 8 + 5] = + acado_Ah_mat[run1 * 2 + 1]*rk_diffsTemp2[(run1 * 20) + (j * 5 + 1)];
rk_A[tmp_index1 * 8 + 6] = + acado_Ah_mat[run1 * 2 + 1]*rk_diffsTemp2[(run1 * 20) + (j * 5 + 2)];
rk_A[tmp_index1 * 8 + 7] = + acado_Ah_mat[run1 * 2 + 1]*rk_diffsTemp2[(run1 * 20) + (j * 5 + 3)];
if( 1 == run1 ) rk_A[(tmp_index1 * 8) + (j + 4)] -= 1.0000000000000000e+00;
}
acado_rhs( rk_xxx, rk_rhsTemp );
rk_b[run1 * 4] = rk_kkk[run1] - rk_rhsTemp[0];
rk_b[run1 * 4 + 1] = rk_kkk[run1 + 2] - rk_rhsTemp[1];
rk_b[run1 * 4 + 2] = rk_kkk[run1 + 4] - rk_rhsTemp[2];
rk_b[run1 * 4 + 3] = rk_kkk[run1 + 6] - rk_rhsTemp[3];
}
det = acado_solve_dim8_system( rk_A, rk_b, rk_dim8_perm );
for (j = 0; j < 2; ++j)
{
rk_kkk[j] += rk_b[j * 4];
rk_kkk[j + 2] += rk_b[j * 4 + 1];
rk_kkk[j + 4] += rk_b[j * 4 + 2];
rk_kkk[j + 6] += rk_b[j * 4 + 3];
}
}
}
for (i = 0; i < 5; ++i)
{
for (run1 = 0; run1 < 2; ++run1)
{
for (j = 0; j < 4; ++j)
{
rk_xxx[j] = rk_eta[j];
tmp_index1 = j;
rk_xxx[j] += + acado_Ah_mat[run1 * 2]*rk_kkk[tmp_index1 * 2];
rk_xxx[j] += + acado_Ah_mat[run1 * 2 + 1]*rk_kkk[tmp_index1 * 2 + 1];
}
acado_rhs( rk_xxx, rk_rhsTemp );
rk_b[run1 * 4] = rk_kkk[run1] - rk_rhsTemp[0];
rk_b[run1 * 4 + 1] = rk_kkk[run1 + 2] - rk_rhsTemp[1];
rk_b[run1 * 4 + 2] = rk_kkk[run1 + 4] - rk_rhsTemp[2];
rk_b[run1 * 4 + 3] = rk_kkk[run1 + 6] - rk_rhsTemp[3];
}
acado_solve_dim8_system_reuse( rk_A, rk_b, rk_dim8_perm );
for (j = 0; j < 2; ++j)
{
rk_kkk[j] += rk_b[j * 4];
rk_kkk[j + 2] += rk_b[j * 4 + 1];
rk_kkk[j + 4] += rk_b[j * 4 + 2];
rk_kkk[j + 6] += rk_b[j * 4 + 3];
}
}
for (run1 = 0; run1 < 2; ++run1)
{
for (j = 0; j < 4; ++j)
{
rk_xxx[j] = rk_eta[j];
tmp_index1 = j;
rk_xxx[j] += + acado_Ah_mat[run1 * 2]*rk_kkk[tmp_index1 * 2];
rk_xxx[j] += + acado_Ah_mat[run1 * 2 + 1]*rk_kkk[tmp_index1 * 2 + 1];
}
acado_diffs( rk_xxx, &(rk_diffsTemp2[ run1 * 20 ]) );
for (j = 0; j < 4; ++j)
{
tmp_index1 = (run1 * 4) + (j);
rk_A[tmp_index1 * 8] = + acado_Ah_mat[run1 * 2]*rk_diffsTemp2[(run1 * 20) + (j * 5)];
rk_A[tmp_index1 * 8 + 1] = + acado_Ah_mat[run1 * 2]*rk_diffsTemp2[(run1 * 20) + (j * 5 + 1)];
rk_A[tmp_index1 * 8 + 2] = + acado_Ah_mat[run1 * 2]*rk_diffsTemp2[(run1 * 20) + (j * 5 + 2)];
rk_A[tmp_index1 * 8 + 3] = + acado_Ah_mat[run1 * 2]*rk_diffsTemp2[(run1 * 20) + (j * 5 + 3)];
if( 0 == run1 ) rk_A[(tmp_index1 * 8) + (j)] -= 1.0000000000000000e+00;
rk_A[tmp_index1 * 8 + 4] = + acado_Ah_mat[run1 * 2 + 1]*rk_diffsTemp2[(run1 * 20) + (j * 5)];
rk_A[tmp_index1 * 8 + 5] = + acado_Ah_mat[run1 * 2 + 1]*rk_diffsTemp2[(run1 * 20) + (j * 5 + 1)];
rk_A[tmp_index1 * 8 + 6] = + acado_Ah_mat[run1 * 2 + 1]*rk_diffsTemp2[(run1 * 20) + (j * 5 + 2)];
rk_A[tmp_index1 * 8 + 7] = + acado_Ah_mat[run1 * 2 + 1]*rk_diffsTemp2[(run1 * 20) + (j * 5 + 3)];
if( 1 == run1 ) rk_A[(tmp_index1 * 8) + (j + 4)] -= 1.0000000000000000e+00;
}
}
for (run1 = 0; run1 < 4; ++run1)
{
for (i = 0; i < 2; ++i)
{
rk_b[i * 4] = - rk_diffsTemp2[(i * 20) + (run1)];
rk_b[i * 4 + 1] = - rk_diffsTemp2[(i * 20) + (run1 + 5)];
rk_b[i * 4 + 2] = - rk_diffsTemp2[(i * 20) + (run1 + 10)];
rk_b[i * 4 + 3] = - rk_diffsTemp2[(i * 20) + (run1 + 15)];
}
if( 0 == run1 ) {
det = acado_solve_dim8_system( rk_A, rk_b, rk_dim8_perm );
}
 else {
acado_solve_dim8_system_reuse( rk_A, rk_b, rk_dim8_perm );
}
for (i = 0; i < 2; ++i)
{
rk_diffK[i] = rk_b[i * 4];
rk_diffK[i + 2] = rk_b[i * 4 + 1];
rk_diffK[i + 4] = rk_b[i * 4 + 2];
rk_diffK[i + 6] = rk_b[i * 4 + 3];
}
for (i = 0; i < 4; ++i)
{
rk_diffsNew2[(i * 5) + (run1)] = (i == run1-0);
rk_diffsNew2[(i * 5) + (run1)] += + rk_diffK[i * 2]*(real_t)5.0000000000000003e-02 + rk_diffK[i * 2 + 1]*(real_t)5.0000000000000003e-02;
}
}
for (run1 = 0; run1 < 1; ++run1)
{
for (i = 0; i < 2; ++i)
{
for (j = 0; j < 4; ++j)
{
tmp_index1 = (i * 4) + (j);
tmp_index2 = (run1) + (j * 5);
rk_b[tmp_index1] = - rk_diffsTemp2[(i * 20) + (tmp_index2 + 4)];
}
}
acado_solve_dim8_system_reuse( rk_A, rk_b, rk_dim8_perm );
for (i = 0; i < 2; ++i)
{
rk_diffK[i] = rk_b[i * 4];
rk_diffK[i + 2] = rk_b[i * 4 + 1];
rk_diffK[i + 4] = rk_b[i * 4 + 2];
rk_diffK[i + 6] = rk_b[i * 4 + 3];
}
for (i = 0; i < 4; ++i)
{
rk_diffsNew2[(i * 5) + (run1 + 4)] = + rk_diffK[i * 2]*(real_t)5.0000000000000003e-02 + rk_diffK[i * 2 + 1]*(real_t)5.0000000000000003e-02;
}
}
rk_eta[0] += + rk_kkk[0]*(real_t)5.0000000000000003e-02 + rk_kkk[1]*(real_t)5.0000000000000003e-02;
rk_eta[1] += + rk_kkk[2]*(real_t)5.0000000000000003e-02 + rk_kkk[3]*(real_t)5.0000000000000003e-02;
rk_eta[2] += + rk_kkk[4]*(real_t)5.0000000000000003e-02 + rk_kkk[5]*(real_t)5.0000000000000003e-02;
rk_eta[3] += + rk_kkk[6]*(real_t)5.0000000000000003e-02 + rk_kkk[7]*(real_t)5.0000000000000003e-02;
for (i = 0; i < 4; ++i)
{
for (j = 0; j < 4; ++j)
{
tmp_index2 = (j) + (i * 4);
rk_eta[tmp_index2 + 4] = rk_diffsNew2[(i * 5) + (j)];
}
for (j = 0; j < 1; ++j)
{
tmp_index2 = (j) + (i);
rk_eta[tmp_index2 + 20] = rk_diffsNew2[(i * 5) + (j + 4)];
}
}
resetIntegrator = 0;
rk_ttt += 1.0000000000000000e+00;
}
for (i = 0; i < 4; ++i)
{
}
if( det < 1e-12 ) {
error = 2;
} else if( det < 1e-6 ) {
error = 1;
} else {
error = 0;
}
return error;
}



