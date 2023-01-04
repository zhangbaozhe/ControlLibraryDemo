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




/******************************************************************************/
/*                                                                            */
/* ACADO code generation                                                      */
/*                                                                            */
/******************************************************************************/


/** Row vector of size: 25 */
real_t state[ 25 ];

int acado_modelSimulation(  )
{
int ret;

int lRun1;
ret = 0;
#pragma omp parallel for private(lRun1, state) shared(acadoWorkspace, acadoVariables)
for (lRun1 = 0; lRun1 < 20; ++lRun1)
{
state[0] = acadoVariables.x[lRun1 * 4];
state[1] = acadoVariables.x[lRun1 * 4 + 1];
state[2] = acadoVariables.x[lRun1 * 4 + 2];
state[3] = acadoVariables.x[lRun1 * 4 + 3];

state[24] = acadoVariables.u[lRun1];

ret = acado_integrate(state, 1);

acadoWorkspace.d[lRun1 * 4] = state[0] - acadoVariables.x[lRun1 * 4 + 4];
acadoWorkspace.d[lRun1 * 4 + 1] = state[1] - acadoVariables.x[lRun1 * 4 + 5];
acadoWorkspace.d[lRun1 * 4 + 2] = state[2] - acadoVariables.x[lRun1 * 4 + 6];
acadoWorkspace.d[lRun1 * 4 + 3] = state[3] - acadoVariables.x[lRun1 * 4 + 7];

acadoWorkspace.evGx[lRun1 * 16] = state[4];
acadoWorkspace.evGx[lRun1 * 16 + 1] = state[5];
acadoWorkspace.evGx[lRun1 * 16 + 2] = state[6];
acadoWorkspace.evGx[lRun1 * 16 + 3] = state[7];
acadoWorkspace.evGx[lRun1 * 16 + 4] = state[8];
acadoWorkspace.evGx[lRun1 * 16 + 5] = state[9];
acadoWorkspace.evGx[lRun1 * 16 + 6] = state[10];
acadoWorkspace.evGx[lRun1 * 16 + 7] = state[11];
acadoWorkspace.evGx[lRun1 * 16 + 8] = state[12];
acadoWorkspace.evGx[lRun1 * 16 + 9] = state[13];
acadoWorkspace.evGx[lRun1 * 16 + 10] = state[14];
acadoWorkspace.evGx[lRun1 * 16 + 11] = state[15];
acadoWorkspace.evGx[lRun1 * 16 + 12] = state[16];
acadoWorkspace.evGx[lRun1 * 16 + 13] = state[17];
acadoWorkspace.evGx[lRun1 * 16 + 14] = state[18];
acadoWorkspace.evGx[lRun1 * 16 + 15] = state[19];

acadoWorkspace.evGu[lRun1 * 4] = state[20];
acadoWorkspace.evGu[lRun1 * 4 + 1] = state[21];
acadoWorkspace.evGu[lRun1 * 4 + 2] = state[22];
acadoWorkspace.evGu[lRun1 * 4 + 3] = state[23];
}
return ret;
}

void acado_evaluateLSQ(const real_t* in, real_t* out)
{
const real_t* xd = in;
const real_t* u = in + 4;

/* Compute outputs: */
out[0] = xd[0];
out[1] = xd[1];
out[2] = xd[2];
out[3] = xd[3];
out[4] = u[0];
}

void acado_evaluateLSQEndTerm(const real_t* in, real_t* out)
{
const real_t* xd = in;

/* Compute outputs: */
out[0] = xd[0];
out[1] = xd[1];
out[2] = xd[2];
out[3] = xd[3];
}

void acado_setObjQ1Q2( real_t* const tmpObjS, real_t* const tmpQ1, real_t* const tmpQ2 )
{
tmpQ2[0] = +tmpObjS[0];
tmpQ2[1] = +tmpObjS[1];
tmpQ2[2] = +tmpObjS[2];
tmpQ2[3] = +tmpObjS[3];
tmpQ2[4] = +tmpObjS[4];
tmpQ2[5] = +tmpObjS[5];
tmpQ2[6] = +tmpObjS[6];
tmpQ2[7] = +tmpObjS[7];
tmpQ2[8] = +tmpObjS[8];
tmpQ2[9] = +tmpObjS[9];
tmpQ2[10] = +tmpObjS[10];
tmpQ2[11] = +tmpObjS[11];
tmpQ2[12] = +tmpObjS[12];
tmpQ2[13] = +tmpObjS[13];
tmpQ2[14] = +tmpObjS[14];
tmpQ2[15] = +tmpObjS[15];
tmpQ2[16] = +tmpObjS[16];
tmpQ2[17] = +tmpObjS[17];
tmpQ2[18] = +tmpObjS[18];
tmpQ2[19] = +tmpObjS[19];
tmpQ1[0] = + tmpQ2[0];
tmpQ1[1] = + tmpQ2[1];
tmpQ1[2] = + tmpQ2[2];
tmpQ1[3] = + tmpQ2[3];
tmpQ1[4] = + tmpQ2[5];
tmpQ1[5] = + tmpQ2[6];
tmpQ1[6] = + tmpQ2[7];
tmpQ1[7] = + tmpQ2[8];
tmpQ1[8] = + tmpQ2[10];
tmpQ1[9] = + tmpQ2[11];
tmpQ1[10] = + tmpQ2[12];
tmpQ1[11] = + tmpQ2[13];
tmpQ1[12] = + tmpQ2[15];
tmpQ1[13] = + tmpQ2[16];
tmpQ1[14] = + tmpQ2[17];
tmpQ1[15] = + tmpQ2[18];
}

void acado_setObjR1R2( real_t* const tmpObjS, real_t* const tmpR1, real_t* const tmpR2 )
{
tmpR2[0] = +tmpObjS[20];
tmpR2[1] = +tmpObjS[21];
tmpR2[2] = +tmpObjS[22];
tmpR2[3] = +tmpObjS[23];
tmpR2[4] = +tmpObjS[24];
tmpR1[0] = + tmpR2[4];
}

void acado_setObjQN1QN2( real_t* const tmpObjSEndTerm, real_t* const tmpQN1, real_t* const tmpQN2 )
{
tmpQN2[0] = +tmpObjSEndTerm[0];
tmpQN2[1] = +tmpObjSEndTerm[1];
tmpQN2[2] = +tmpObjSEndTerm[2];
tmpQN2[3] = +tmpObjSEndTerm[3];
tmpQN2[4] = +tmpObjSEndTerm[4];
tmpQN2[5] = +tmpObjSEndTerm[5];
tmpQN2[6] = +tmpObjSEndTerm[6];
tmpQN2[7] = +tmpObjSEndTerm[7];
tmpQN2[8] = +tmpObjSEndTerm[8];
tmpQN2[9] = +tmpObjSEndTerm[9];
tmpQN2[10] = +tmpObjSEndTerm[10];
tmpQN2[11] = +tmpObjSEndTerm[11];
tmpQN2[12] = +tmpObjSEndTerm[12];
tmpQN2[13] = +tmpObjSEndTerm[13];
tmpQN2[14] = +tmpObjSEndTerm[14];
tmpQN2[15] = +tmpObjSEndTerm[15];
tmpQN1[0] = + tmpQN2[0];
tmpQN1[1] = + tmpQN2[1];
tmpQN1[2] = + tmpQN2[2];
tmpQN1[3] = + tmpQN2[3];
tmpQN1[4] = + tmpQN2[4];
tmpQN1[5] = + tmpQN2[5];
tmpQN1[6] = + tmpQN2[6];
tmpQN1[7] = + tmpQN2[7];
tmpQN1[8] = + tmpQN2[8];
tmpQN1[9] = + tmpQN2[9];
tmpQN1[10] = + tmpQN2[10];
tmpQN1[11] = + tmpQN2[11];
tmpQN1[12] = + tmpQN2[12];
tmpQN1[13] = + tmpQN2[13];
tmpQN1[14] = + tmpQN2[14];
tmpQN1[15] = + tmpQN2[15];
}

void acado_evaluateObjective(  )
{
int runObj;
for (runObj = 0; runObj < 20; ++runObj)
{
acadoWorkspace.objValueIn[0] = acadoVariables.x[runObj * 4];
acadoWorkspace.objValueIn[1] = acadoVariables.x[runObj * 4 + 1];
acadoWorkspace.objValueIn[2] = acadoVariables.x[runObj * 4 + 2];
acadoWorkspace.objValueIn[3] = acadoVariables.x[runObj * 4 + 3];
acadoWorkspace.objValueIn[4] = acadoVariables.u[runObj];

acado_evaluateLSQ( acadoWorkspace.objValueIn, acadoWorkspace.objValueOut );
acadoWorkspace.Dy[runObj * 5] = acadoWorkspace.objValueOut[0];
acadoWorkspace.Dy[runObj * 5 + 1] = acadoWorkspace.objValueOut[1];
acadoWorkspace.Dy[runObj * 5 + 2] = acadoWorkspace.objValueOut[2];
acadoWorkspace.Dy[runObj * 5 + 3] = acadoWorkspace.objValueOut[3];
acadoWorkspace.Dy[runObj * 5 + 4] = acadoWorkspace.objValueOut[4];

acado_setObjQ1Q2( &(acadoVariables.W[ runObj * 25 ]), &(acadoWorkspace.Q1[ runObj * 16 ]), &(acadoWorkspace.Q2[ runObj * 20 ]) );

acado_setObjR1R2( &(acadoVariables.W[ runObj * 25 ]), &(acadoWorkspace.R1[ runObj ]), &(acadoWorkspace.R2[ runObj * 5 ]) );

}
acadoWorkspace.objValueIn[0] = acadoVariables.x[80];
acadoWorkspace.objValueIn[1] = acadoVariables.x[81];
acadoWorkspace.objValueIn[2] = acadoVariables.x[82];
acadoWorkspace.objValueIn[3] = acadoVariables.x[83];
acado_evaluateLSQEndTerm( acadoWorkspace.objValueIn, acadoWorkspace.objValueOut );

acadoWorkspace.DyN[0] = acadoWorkspace.objValueOut[0];
acadoWorkspace.DyN[1] = acadoWorkspace.objValueOut[1];
acadoWorkspace.DyN[2] = acadoWorkspace.objValueOut[2];
acadoWorkspace.DyN[3] = acadoWorkspace.objValueOut[3];

acado_setObjQN1QN2( acadoVariables.WN, acadoWorkspace.QN1, acadoWorkspace.QN2 );

}

void acado_multGxGu( real_t* const Gx1, real_t* const Gu1, real_t* const Gu2 )
{
Gu2[0] = + Gx1[0]*Gu1[0] + Gx1[1]*Gu1[1] + Gx1[2]*Gu1[2] + Gx1[3]*Gu1[3];
Gu2[1] = + Gx1[4]*Gu1[0] + Gx1[5]*Gu1[1] + Gx1[6]*Gu1[2] + Gx1[7]*Gu1[3];
Gu2[2] = + Gx1[8]*Gu1[0] + Gx1[9]*Gu1[1] + Gx1[10]*Gu1[2] + Gx1[11]*Gu1[3];
Gu2[3] = + Gx1[12]*Gu1[0] + Gx1[13]*Gu1[1] + Gx1[14]*Gu1[2] + Gx1[15]*Gu1[3];
}

void acado_moveGuE( real_t* const Gu1, real_t* const Gu2 )
{
Gu2[0] = Gu1[0];
Gu2[1] = Gu1[1];
Gu2[2] = Gu1[2];
Gu2[3] = Gu1[3];
}

void acado_multBTW1( real_t* const Gu1, real_t* const Gu2, int iRow, int iCol )
{
acadoWorkspace.H[(iRow * 20) + (iCol)] = + Gu1[0]*Gu2[0] + Gu1[1]*Gu2[1] + Gu1[2]*Gu2[2] + Gu1[3]*Gu2[3];
}

void acado_multBTW1_R1( real_t* const R11, real_t* const Gu1, real_t* const Gu2, int iRow )
{
acadoWorkspace.H[iRow * 21] = + Gu1[0]*Gu2[0] + Gu1[1]*Gu2[1] + Gu1[2]*Gu2[2] + Gu1[3]*Gu2[3] + R11[0];
acadoWorkspace.H[iRow * 21] += 1.0000000000000000e+00;
}

void acado_multGxTGu( real_t* const Gx1, real_t* const Gu1, real_t* const Gu2 )
{
Gu2[0] = + Gx1[0]*Gu1[0] + Gx1[4]*Gu1[1] + Gx1[8]*Gu1[2] + Gx1[12]*Gu1[3];
Gu2[1] = + Gx1[1]*Gu1[0] + Gx1[5]*Gu1[1] + Gx1[9]*Gu1[2] + Gx1[13]*Gu1[3];
Gu2[2] = + Gx1[2]*Gu1[0] + Gx1[6]*Gu1[1] + Gx1[10]*Gu1[2] + Gx1[14]*Gu1[3];
Gu2[3] = + Gx1[3]*Gu1[0] + Gx1[7]*Gu1[1] + Gx1[11]*Gu1[2] + Gx1[15]*Gu1[3];
}

void acado_multQEW2( real_t* const Q11, real_t* const Gu1, real_t* const Gu2, real_t* const Gu3 )
{
Gu3[0] = + Q11[0]*Gu1[0] + Q11[1]*Gu1[1] + Q11[2]*Gu1[2] + Q11[3]*Gu1[3] + Gu2[0];
Gu3[1] = + Q11[4]*Gu1[0] + Q11[5]*Gu1[1] + Q11[6]*Gu1[2] + Q11[7]*Gu1[3] + Gu2[1];
Gu3[2] = + Q11[8]*Gu1[0] + Q11[9]*Gu1[1] + Q11[10]*Gu1[2] + Q11[11]*Gu1[3] + Gu2[2];
Gu3[3] = + Q11[12]*Gu1[0] + Q11[13]*Gu1[1] + Q11[14]*Gu1[2] + Q11[15]*Gu1[3] + Gu2[3];
}

void acado_macATw1QDy( real_t* const Gx1, real_t* const w11, real_t* const w12, real_t* const w13 )
{
w13[0] = + Gx1[0]*w11[0] + Gx1[4]*w11[1] + Gx1[8]*w11[2] + Gx1[12]*w11[3] + w12[0];
w13[1] = + Gx1[1]*w11[0] + Gx1[5]*w11[1] + Gx1[9]*w11[2] + Gx1[13]*w11[3] + w12[1];
w13[2] = + Gx1[2]*w11[0] + Gx1[6]*w11[1] + Gx1[10]*w11[2] + Gx1[14]*w11[3] + w12[2];
w13[3] = + Gx1[3]*w11[0] + Gx1[7]*w11[1] + Gx1[11]*w11[2] + Gx1[15]*w11[3] + w12[3];
}

void acado_macBTw1( real_t* const Gu1, real_t* const w11, real_t* const U1 )
{
U1[0] += + Gu1[0]*w11[0] + Gu1[1]*w11[1] + Gu1[2]*w11[2] + Gu1[3]*w11[3];
}

void acado_macQSbarW2( real_t* const Q11, real_t* const w11, real_t* const w12, real_t* const w13 )
{
w13[0] = + Q11[0]*w11[0] + Q11[1]*w11[1] + Q11[2]*w11[2] + Q11[3]*w11[3] + w12[0];
w13[1] = + Q11[4]*w11[0] + Q11[5]*w11[1] + Q11[6]*w11[2] + Q11[7]*w11[3] + w12[1];
w13[2] = + Q11[8]*w11[0] + Q11[9]*w11[1] + Q11[10]*w11[2] + Q11[11]*w11[3] + w12[2];
w13[3] = + Q11[12]*w11[0] + Q11[13]*w11[1] + Q11[14]*w11[2] + Q11[15]*w11[3] + w12[3];
}

void acado_macASbar( real_t* const Gx1, real_t* const w11, real_t* const w12 )
{
w12[0] += + Gx1[0]*w11[0] + Gx1[1]*w11[1] + Gx1[2]*w11[2] + Gx1[3]*w11[3];
w12[1] += + Gx1[4]*w11[0] + Gx1[5]*w11[1] + Gx1[6]*w11[2] + Gx1[7]*w11[3];
w12[2] += + Gx1[8]*w11[0] + Gx1[9]*w11[1] + Gx1[10]*w11[2] + Gx1[11]*w11[3];
w12[3] += + Gx1[12]*w11[0] + Gx1[13]*w11[1] + Gx1[14]*w11[2] + Gx1[15]*w11[3];
}

void acado_expansionStep( real_t* const Gx1, real_t* const Gu1, real_t* const U1, real_t* const w11, real_t* const w12 )
{
w12[0] += + Gx1[0]*w11[0] + Gx1[1]*w11[1] + Gx1[2]*w11[2] + Gx1[3]*w11[3];
w12[1] += + Gx1[4]*w11[0] + Gx1[5]*w11[1] + Gx1[6]*w11[2] + Gx1[7]*w11[3];
w12[2] += + Gx1[8]*w11[0] + Gx1[9]*w11[1] + Gx1[10]*w11[2] + Gx1[11]*w11[3];
w12[3] += + Gx1[12]*w11[0] + Gx1[13]*w11[1] + Gx1[14]*w11[2] + Gx1[15]*w11[3];
w12[0] += + Gu1[0]*U1[0];
w12[1] += + Gu1[1]*U1[0];
w12[2] += + Gu1[2]*U1[0];
w12[3] += + Gu1[3]*U1[0];
}

void acado_copyHTH( int iRow, int iCol )
{
acadoWorkspace.H[(iRow * 20) + (iCol)] = acadoWorkspace.H[(iCol * 20) + (iRow)];
}

void acado_multRDy( real_t* const R2, real_t* const Dy1, real_t* const RDy1 )
{
RDy1[0] = + R2[0]*Dy1[0] + R2[1]*Dy1[1] + R2[2]*Dy1[2] + R2[3]*Dy1[3] + R2[4]*Dy1[4];
}

void acado_multQDy( real_t* const Q2, real_t* const Dy1, real_t* const QDy1 )
{
QDy1[0] = + Q2[0]*Dy1[0] + Q2[1]*Dy1[1] + Q2[2]*Dy1[2] + Q2[3]*Dy1[3] + Q2[4]*Dy1[4];
QDy1[1] = + Q2[5]*Dy1[0] + Q2[6]*Dy1[1] + Q2[7]*Dy1[2] + Q2[8]*Dy1[3] + Q2[9]*Dy1[4];
QDy1[2] = + Q2[10]*Dy1[0] + Q2[11]*Dy1[1] + Q2[12]*Dy1[2] + Q2[13]*Dy1[3] + Q2[14]*Dy1[4];
QDy1[3] = + Q2[15]*Dy1[0] + Q2[16]*Dy1[1] + Q2[17]*Dy1[2] + Q2[18]*Dy1[3] + Q2[19]*Dy1[4];
}

void acado_condensePrep(  )
{
int lRun1;
int lRun2;
int lRun3;
for (lRun2 = 0; lRun2 < 20; ++lRun2)
{
lRun3 = ((lRun2) * (lRun2 * -1 + 41)) / (2);
acado_moveGuE( &(acadoWorkspace.evGu[ lRun2 * 4 ]), &(acadoWorkspace.E[ lRun3 * 4 ]) );
for (lRun1 = 1; lRun1 < lRun2 * -1 + 20; ++lRun1)
{
acado_multGxGu( &(acadoWorkspace.evGx[ ((((lRun2) + (lRun1)) * (4)) * (4)) + (0) ]), &(acadoWorkspace.E[ (((((lRun3) + (lRun1)) - (1)) * (4)) * (1)) + (0) ]), &(acadoWorkspace.E[ ((((lRun3) + (lRun1)) * (4)) * (1)) + (0) ]) );
}

acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ ((((((lRun3) - (lRun2)) + (20)) - (1)) * (4)) * (1)) + (0) ]), acadoWorkspace.W1 );
for (lRun1 = 19; lRun2 < lRun1; --lRun1)
{
acado_multBTW1( &(acadoWorkspace.evGu[ lRun1 * 4 ]), acadoWorkspace.W1, lRun1, lRun2 );
acado_multGxTGu( &(acadoWorkspace.evGx[ lRun1 * 16 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ lRun1 * 16 ]), &(acadoWorkspace.E[ ((((((lRun3) + (lRun1)) - (lRun2)) - (1)) * (4)) * (1)) + (0) ]), acadoWorkspace.W2, acadoWorkspace.W1 );
}
acado_multBTW1_R1( &(acadoWorkspace.R1[ lRun2 ]), &(acadoWorkspace.evGu[ lRun2 * 4 ]), acadoWorkspace.W1, lRun2 );
}

acado_copyHTH( 0, 1 );
acado_copyHTH( 0, 2 );
acado_copyHTH( 1, 2 );
acado_copyHTH( 0, 3 );
acado_copyHTH( 1, 3 );
acado_copyHTH( 2, 3 );
acado_copyHTH( 0, 4 );
acado_copyHTH( 1, 4 );
acado_copyHTH( 2, 4 );
acado_copyHTH( 3, 4 );
acado_copyHTH( 0, 5 );
acado_copyHTH( 1, 5 );
acado_copyHTH( 2, 5 );
acado_copyHTH( 3, 5 );
acado_copyHTH( 4, 5 );
acado_copyHTH( 0, 6 );
acado_copyHTH( 1, 6 );
acado_copyHTH( 2, 6 );
acado_copyHTH( 3, 6 );
acado_copyHTH( 4, 6 );
acado_copyHTH( 5, 6 );
acado_copyHTH( 0, 7 );
acado_copyHTH( 1, 7 );
acado_copyHTH( 2, 7 );
acado_copyHTH( 3, 7 );
acado_copyHTH( 4, 7 );
acado_copyHTH( 5, 7 );
acado_copyHTH( 6, 7 );
acado_copyHTH( 0, 8 );
acado_copyHTH( 1, 8 );
acado_copyHTH( 2, 8 );
acado_copyHTH( 3, 8 );
acado_copyHTH( 4, 8 );
acado_copyHTH( 5, 8 );
acado_copyHTH( 6, 8 );
acado_copyHTH( 7, 8 );
acado_copyHTH( 0, 9 );
acado_copyHTH( 1, 9 );
acado_copyHTH( 2, 9 );
acado_copyHTH( 3, 9 );
acado_copyHTH( 4, 9 );
acado_copyHTH( 5, 9 );
acado_copyHTH( 6, 9 );
acado_copyHTH( 7, 9 );
acado_copyHTH( 8, 9 );
acado_copyHTH( 0, 10 );
acado_copyHTH( 1, 10 );
acado_copyHTH( 2, 10 );
acado_copyHTH( 3, 10 );
acado_copyHTH( 4, 10 );
acado_copyHTH( 5, 10 );
acado_copyHTH( 6, 10 );
acado_copyHTH( 7, 10 );
acado_copyHTH( 8, 10 );
acado_copyHTH( 9, 10 );
acado_copyHTH( 0, 11 );
acado_copyHTH( 1, 11 );
acado_copyHTH( 2, 11 );
acado_copyHTH( 3, 11 );
acado_copyHTH( 4, 11 );
acado_copyHTH( 5, 11 );
acado_copyHTH( 6, 11 );
acado_copyHTH( 7, 11 );
acado_copyHTH( 8, 11 );
acado_copyHTH( 9, 11 );
acado_copyHTH( 10, 11 );
acado_copyHTH( 0, 12 );
acado_copyHTH( 1, 12 );
acado_copyHTH( 2, 12 );
acado_copyHTH( 3, 12 );
acado_copyHTH( 4, 12 );
acado_copyHTH( 5, 12 );
acado_copyHTH( 6, 12 );
acado_copyHTH( 7, 12 );
acado_copyHTH( 8, 12 );
acado_copyHTH( 9, 12 );
acado_copyHTH( 10, 12 );
acado_copyHTH( 11, 12 );
acado_copyHTH( 0, 13 );
acado_copyHTH( 1, 13 );
acado_copyHTH( 2, 13 );
acado_copyHTH( 3, 13 );
acado_copyHTH( 4, 13 );
acado_copyHTH( 5, 13 );
acado_copyHTH( 6, 13 );
acado_copyHTH( 7, 13 );
acado_copyHTH( 8, 13 );
acado_copyHTH( 9, 13 );
acado_copyHTH( 10, 13 );
acado_copyHTH( 11, 13 );
acado_copyHTH( 12, 13 );
acado_copyHTH( 0, 14 );
acado_copyHTH( 1, 14 );
acado_copyHTH( 2, 14 );
acado_copyHTH( 3, 14 );
acado_copyHTH( 4, 14 );
acado_copyHTH( 5, 14 );
acado_copyHTH( 6, 14 );
acado_copyHTH( 7, 14 );
acado_copyHTH( 8, 14 );
acado_copyHTH( 9, 14 );
acado_copyHTH( 10, 14 );
acado_copyHTH( 11, 14 );
acado_copyHTH( 12, 14 );
acado_copyHTH( 13, 14 );
acado_copyHTH( 0, 15 );
acado_copyHTH( 1, 15 );
acado_copyHTH( 2, 15 );
acado_copyHTH( 3, 15 );
acado_copyHTH( 4, 15 );
acado_copyHTH( 5, 15 );
acado_copyHTH( 6, 15 );
acado_copyHTH( 7, 15 );
acado_copyHTH( 8, 15 );
acado_copyHTH( 9, 15 );
acado_copyHTH( 10, 15 );
acado_copyHTH( 11, 15 );
acado_copyHTH( 12, 15 );
acado_copyHTH( 13, 15 );
acado_copyHTH( 14, 15 );
acado_copyHTH( 0, 16 );
acado_copyHTH( 1, 16 );
acado_copyHTH( 2, 16 );
acado_copyHTH( 3, 16 );
acado_copyHTH( 4, 16 );
acado_copyHTH( 5, 16 );
acado_copyHTH( 6, 16 );
acado_copyHTH( 7, 16 );
acado_copyHTH( 8, 16 );
acado_copyHTH( 9, 16 );
acado_copyHTH( 10, 16 );
acado_copyHTH( 11, 16 );
acado_copyHTH( 12, 16 );
acado_copyHTH( 13, 16 );
acado_copyHTH( 14, 16 );
acado_copyHTH( 15, 16 );
acado_copyHTH( 0, 17 );
acado_copyHTH( 1, 17 );
acado_copyHTH( 2, 17 );
acado_copyHTH( 3, 17 );
acado_copyHTH( 4, 17 );
acado_copyHTH( 5, 17 );
acado_copyHTH( 6, 17 );
acado_copyHTH( 7, 17 );
acado_copyHTH( 8, 17 );
acado_copyHTH( 9, 17 );
acado_copyHTH( 10, 17 );
acado_copyHTH( 11, 17 );
acado_copyHTH( 12, 17 );
acado_copyHTH( 13, 17 );
acado_copyHTH( 14, 17 );
acado_copyHTH( 15, 17 );
acado_copyHTH( 16, 17 );
acado_copyHTH( 0, 18 );
acado_copyHTH( 1, 18 );
acado_copyHTH( 2, 18 );
acado_copyHTH( 3, 18 );
acado_copyHTH( 4, 18 );
acado_copyHTH( 5, 18 );
acado_copyHTH( 6, 18 );
acado_copyHTH( 7, 18 );
acado_copyHTH( 8, 18 );
acado_copyHTH( 9, 18 );
acado_copyHTH( 10, 18 );
acado_copyHTH( 11, 18 );
acado_copyHTH( 12, 18 );
acado_copyHTH( 13, 18 );
acado_copyHTH( 14, 18 );
acado_copyHTH( 15, 18 );
acado_copyHTH( 16, 18 );
acado_copyHTH( 17, 18 );
acado_copyHTH( 0, 19 );
acado_copyHTH( 1, 19 );
acado_copyHTH( 2, 19 );
acado_copyHTH( 3, 19 );
acado_copyHTH( 4, 19 );
acado_copyHTH( 5, 19 );
acado_copyHTH( 6, 19 );
acado_copyHTH( 7, 19 );
acado_copyHTH( 8, 19 );
acado_copyHTH( 9, 19 );
acado_copyHTH( 10, 19 );
acado_copyHTH( 11, 19 );
acado_copyHTH( 12, 19 );
acado_copyHTH( 13, 19 );
acado_copyHTH( 14, 19 );
acado_copyHTH( 15, 19 );
acado_copyHTH( 16, 19 );
acado_copyHTH( 17, 19 );
acado_copyHTH( 18, 19 );

acadoWorkspace.sbar[4] = acadoWorkspace.d[0];
acadoWorkspace.sbar[5] = acadoWorkspace.d[1];
acadoWorkspace.sbar[6] = acadoWorkspace.d[2];
acadoWorkspace.sbar[7] = acadoWorkspace.d[3];
acadoWorkspace.sbar[8] = acadoWorkspace.d[4];
acadoWorkspace.sbar[9] = acadoWorkspace.d[5];
acadoWorkspace.sbar[10] = acadoWorkspace.d[6];
acadoWorkspace.sbar[11] = acadoWorkspace.d[7];
acadoWorkspace.sbar[12] = acadoWorkspace.d[8];
acadoWorkspace.sbar[13] = acadoWorkspace.d[9];
acadoWorkspace.sbar[14] = acadoWorkspace.d[10];
acadoWorkspace.sbar[15] = acadoWorkspace.d[11];
acadoWorkspace.sbar[16] = acadoWorkspace.d[12];
acadoWorkspace.sbar[17] = acadoWorkspace.d[13];
acadoWorkspace.sbar[18] = acadoWorkspace.d[14];
acadoWorkspace.sbar[19] = acadoWorkspace.d[15];
acadoWorkspace.sbar[20] = acadoWorkspace.d[16];
acadoWorkspace.sbar[21] = acadoWorkspace.d[17];
acadoWorkspace.sbar[22] = acadoWorkspace.d[18];
acadoWorkspace.sbar[23] = acadoWorkspace.d[19];
acadoWorkspace.sbar[24] = acadoWorkspace.d[20];
acadoWorkspace.sbar[25] = acadoWorkspace.d[21];
acadoWorkspace.sbar[26] = acadoWorkspace.d[22];
acadoWorkspace.sbar[27] = acadoWorkspace.d[23];
acadoWorkspace.sbar[28] = acadoWorkspace.d[24];
acadoWorkspace.sbar[29] = acadoWorkspace.d[25];
acadoWorkspace.sbar[30] = acadoWorkspace.d[26];
acadoWorkspace.sbar[31] = acadoWorkspace.d[27];
acadoWorkspace.sbar[32] = acadoWorkspace.d[28];
acadoWorkspace.sbar[33] = acadoWorkspace.d[29];
acadoWorkspace.sbar[34] = acadoWorkspace.d[30];
acadoWorkspace.sbar[35] = acadoWorkspace.d[31];
acadoWorkspace.sbar[36] = acadoWorkspace.d[32];
acadoWorkspace.sbar[37] = acadoWorkspace.d[33];
acadoWorkspace.sbar[38] = acadoWorkspace.d[34];
acadoWorkspace.sbar[39] = acadoWorkspace.d[35];
acadoWorkspace.sbar[40] = acadoWorkspace.d[36];
acadoWorkspace.sbar[41] = acadoWorkspace.d[37];
acadoWorkspace.sbar[42] = acadoWorkspace.d[38];
acadoWorkspace.sbar[43] = acadoWorkspace.d[39];
acadoWorkspace.sbar[44] = acadoWorkspace.d[40];
acadoWorkspace.sbar[45] = acadoWorkspace.d[41];
acadoWorkspace.sbar[46] = acadoWorkspace.d[42];
acadoWorkspace.sbar[47] = acadoWorkspace.d[43];
acadoWorkspace.sbar[48] = acadoWorkspace.d[44];
acadoWorkspace.sbar[49] = acadoWorkspace.d[45];
acadoWorkspace.sbar[50] = acadoWorkspace.d[46];
acadoWorkspace.sbar[51] = acadoWorkspace.d[47];
acadoWorkspace.sbar[52] = acadoWorkspace.d[48];
acadoWorkspace.sbar[53] = acadoWorkspace.d[49];
acadoWorkspace.sbar[54] = acadoWorkspace.d[50];
acadoWorkspace.sbar[55] = acadoWorkspace.d[51];
acadoWorkspace.sbar[56] = acadoWorkspace.d[52];
acadoWorkspace.sbar[57] = acadoWorkspace.d[53];
acadoWorkspace.sbar[58] = acadoWorkspace.d[54];
acadoWorkspace.sbar[59] = acadoWorkspace.d[55];
acadoWorkspace.sbar[60] = acadoWorkspace.d[56];
acadoWorkspace.sbar[61] = acadoWorkspace.d[57];
acadoWorkspace.sbar[62] = acadoWorkspace.d[58];
acadoWorkspace.sbar[63] = acadoWorkspace.d[59];
acadoWorkspace.sbar[64] = acadoWorkspace.d[60];
acadoWorkspace.sbar[65] = acadoWorkspace.d[61];
acadoWorkspace.sbar[66] = acadoWorkspace.d[62];
acadoWorkspace.sbar[67] = acadoWorkspace.d[63];
acadoWorkspace.sbar[68] = acadoWorkspace.d[64];
acadoWorkspace.sbar[69] = acadoWorkspace.d[65];
acadoWorkspace.sbar[70] = acadoWorkspace.d[66];
acadoWorkspace.sbar[71] = acadoWorkspace.d[67];
acadoWorkspace.sbar[72] = acadoWorkspace.d[68];
acadoWorkspace.sbar[73] = acadoWorkspace.d[69];
acadoWorkspace.sbar[74] = acadoWorkspace.d[70];
acadoWorkspace.sbar[75] = acadoWorkspace.d[71];
acadoWorkspace.sbar[76] = acadoWorkspace.d[72];
acadoWorkspace.sbar[77] = acadoWorkspace.d[73];
acadoWorkspace.sbar[78] = acadoWorkspace.d[74];
acadoWorkspace.sbar[79] = acadoWorkspace.d[75];
acadoWorkspace.sbar[80] = acadoWorkspace.d[76];
acadoWorkspace.sbar[81] = acadoWorkspace.d[77];
acadoWorkspace.sbar[82] = acadoWorkspace.d[78];
acadoWorkspace.sbar[83] = acadoWorkspace.d[79];

}

void acado_condenseFdb(  )
{
acadoWorkspace.Dx0[0] = acadoVariables.x0[0] - acadoVariables.x[0];
acadoWorkspace.Dx0[1] = acadoVariables.x0[1] - acadoVariables.x[1];
acadoWorkspace.Dx0[2] = acadoVariables.x0[2] - acadoVariables.x[2];
acadoWorkspace.Dx0[3] = acadoVariables.x0[3] - acadoVariables.x[3];
acadoWorkspace.Dy[0] -= acadoVariables.y[0];
acadoWorkspace.Dy[1] -= acadoVariables.y[1];
acadoWorkspace.Dy[2] -= acadoVariables.y[2];
acadoWorkspace.Dy[3] -= acadoVariables.y[3];
acadoWorkspace.Dy[4] -= acadoVariables.y[4];
acadoWorkspace.Dy[5] -= acadoVariables.y[5];
acadoWorkspace.Dy[6] -= acadoVariables.y[6];
acadoWorkspace.Dy[7] -= acadoVariables.y[7];
acadoWorkspace.Dy[8] -= acadoVariables.y[8];
acadoWorkspace.Dy[9] -= acadoVariables.y[9];
acadoWorkspace.Dy[10] -= acadoVariables.y[10];
acadoWorkspace.Dy[11] -= acadoVariables.y[11];
acadoWorkspace.Dy[12] -= acadoVariables.y[12];
acadoWorkspace.Dy[13] -= acadoVariables.y[13];
acadoWorkspace.Dy[14] -= acadoVariables.y[14];
acadoWorkspace.Dy[15] -= acadoVariables.y[15];
acadoWorkspace.Dy[16] -= acadoVariables.y[16];
acadoWorkspace.Dy[17] -= acadoVariables.y[17];
acadoWorkspace.Dy[18] -= acadoVariables.y[18];
acadoWorkspace.Dy[19] -= acadoVariables.y[19];
acadoWorkspace.Dy[20] -= acadoVariables.y[20];
acadoWorkspace.Dy[21] -= acadoVariables.y[21];
acadoWorkspace.Dy[22] -= acadoVariables.y[22];
acadoWorkspace.Dy[23] -= acadoVariables.y[23];
acadoWorkspace.Dy[24] -= acadoVariables.y[24];
acadoWorkspace.Dy[25] -= acadoVariables.y[25];
acadoWorkspace.Dy[26] -= acadoVariables.y[26];
acadoWorkspace.Dy[27] -= acadoVariables.y[27];
acadoWorkspace.Dy[28] -= acadoVariables.y[28];
acadoWorkspace.Dy[29] -= acadoVariables.y[29];
acadoWorkspace.Dy[30] -= acadoVariables.y[30];
acadoWorkspace.Dy[31] -= acadoVariables.y[31];
acadoWorkspace.Dy[32] -= acadoVariables.y[32];
acadoWorkspace.Dy[33] -= acadoVariables.y[33];
acadoWorkspace.Dy[34] -= acadoVariables.y[34];
acadoWorkspace.Dy[35] -= acadoVariables.y[35];
acadoWorkspace.Dy[36] -= acadoVariables.y[36];
acadoWorkspace.Dy[37] -= acadoVariables.y[37];
acadoWorkspace.Dy[38] -= acadoVariables.y[38];
acadoWorkspace.Dy[39] -= acadoVariables.y[39];
acadoWorkspace.Dy[40] -= acadoVariables.y[40];
acadoWorkspace.Dy[41] -= acadoVariables.y[41];
acadoWorkspace.Dy[42] -= acadoVariables.y[42];
acadoWorkspace.Dy[43] -= acadoVariables.y[43];
acadoWorkspace.Dy[44] -= acadoVariables.y[44];
acadoWorkspace.Dy[45] -= acadoVariables.y[45];
acadoWorkspace.Dy[46] -= acadoVariables.y[46];
acadoWorkspace.Dy[47] -= acadoVariables.y[47];
acadoWorkspace.Dy[48] -= acadoVariables.y[48];
acadoWorkspace.Dy[49] -= acadoVariables.y[49];
acadoWorkspace.Dy[50] -= acadoVariables.y[50];
acadoWorkspace.Dy[51] -= acadoVariables.y[51];
acadoWorkspace.Dy[52] -= acadoVariables.y[52];
acadoWorkspace.Dy[53] -= acadoVariables.y[53];
acadoWorkspace.Dy[54] -= acadoVariables.y[54];
acadoWorkspace.Dy[55] -= acadoVariables.y[55];
acadoWorkspace.Dy[56] -= acadoVariables.y[56];
acadoWorkspace.Dy[57] -= acadoVariables.y[57];
acadoWorkspace.Dy[58] -= acadoVariables.y[58];
acadoWorkspace.Dy[59] -= acadoVariables.y[59];
acadoWorkspace.Dy[60] -= acadoVariables.y[60];
acadoWorkspace.Dy[61] -= acadoVariables.y[61];
acadoWorkspace.Dy[62] -= acadoVariables.y[62];
acadoWorkspace.Dy[63] -= acadoVariables.y[63];
acadoWorkspace.Dy[64] -= acadoVariables.y[64];
acadoWorkspace.Dy[65] -= acadoVariables.y[65];
acadoWorkspace.Dy[66] -= acadoVariables.y[66];
acadoWorkspace.Dy[67] -= acadoVariables.y[67];
acadoWorkspace.Dy[68] -= acadoVariables.y[68];
acadoWorkspace.Dy[69] -= acadoVariables.y[69];
acadoWorkspace.Dy[70] -= acadoVariables.y[70];
acadoWorkspace.Dy[71] -= acadoVariables.y[71];
acadoWorkspace.Dy[72] -= acadoVariables.y[72];
acadoWorkspace.Dy[73] -= acadoVariables.y[73];
acadoWorkspace.Dy[74] -= acadoVariables.y[74];
acadoWorkspace.Dy[75] -= acadoVariables.y[75];
acadoWorkspace.Dy[76] -= acadoVariables.y[76];
acadoWorkspace.Dy[77] -= acadoVariables.y[77];
acadoWorkspace.Dy[78] -= acadoVariables.y[78];
acadoWorkspace.Dy[79] -= acadoVariables.y[79];
acadoWorkspace.Dy[80] -= acadoVariables.y[80];
acadoWorkspace.Dy[81] -= acadoVariables.y[81];
acadoWorkspace.Dy[82] -= acadoVariables.y[82];
acadoWorkspace.Dy[83] -= acadoVariables.y[83];
acadoWorkspace.Dy[84] -= acadoVariables.y[84];
acadoWorkspace.Dy[85] -= acadoVariables.y[85];
acadoWorkspace.Dy[86] -= acadoVariables.y[86];
acadoWorkspace.Dy[87] -= acadoVariables.y[87];
acadoWorkspace.Dy[88] -= acadoVariables.y[88];
acadoWorkspace.Dy[89] -= acadoVariables.y[89];
acadoWorkspace.Dy[90] -= acadoVariables.y[90];
acadoWorkspace.Dy[91] -= acadoVariables.y[91];
acadoWorkspace.Dy[92] -= acadoVariables.y[92];
acadoWorkspace.Dy[93] -= acadoVariables.y[93];
acadoWorkspace.Dy[94] -= acadoVariables.y[94];
acadoWorkspace.Dy[95] -= acadoVariables.y[95];
acadoWorkspace.Dy[96] -= acadoVariables.y[96];
acadoWorkspace.Dy[97] -= acadoVariables.y[97];
acadoWorkspace.Dy[98] -= acadoVariables.y[98];
acadoWorkspace.Dy[99] -= acadoVariables.y[99];
acadoWorkspace.DyN[0] -= acadoVariables.yN[0];
acadoWorkspace.DyN[1] -= acadoVariables.yN[1];
acadoWorkspace.DyN[2] -= acadoVariables.yN[2];
acadoWorkspace.DyN[3] -= acadoVariables.yN[3];

acado_multRDy( acadoWorkspace.R2, acadoWorkspace.Dy, acadoWorkspace.g );
acado_multRDy( &(acadoWorkspace.R2[ 5 ]), &(acadoWorkspace.Dy[ 5 ]), &(acadoWorkspace.g[ 1 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 10 ]), &(acadoWorkspace.Dy[ 10 ]), &(acadoWorkspace.g[ 2 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 15 ]), &(acadoWorkspace.Dy[ 15 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 20 ]), &(acadoWorkspace.Dy[ 20 ]), &(acadoWorkspace.g[ 4 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 25 ]), &(acadoWorkspace.Dy[ 25 ]), &(acadoWorkspace.g[ 5 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 30 ]), &(acadoWorkspace.Dy[ 30 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 35 ]), &(acadoWorkspace.Dy[ 35 ]), &(acadoWorkspace.g[ 7 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 40 ]), &(acadoWorkspace.Dy[ 40 ]), &(acadoWorkspace.g[ 8 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 45 ]), &(acadoWorkspace.Dy[ 45 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 50 ]), &(acadoWorkspace.Dy[ 50 ]), &(acadoWorkspace.g[ 10 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 55 ]), &(acadoWorkspace.Dy[ 55 ]), &(acadoWorkspace.g[ 11 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 60 ]), &(acadoWorkspace.Dy[ 60 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 65 ]), &(acadoWorkspace.Dy[ 65 ]), &(acadoWorkspace.g[ 13 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 70 ]), &(acadoWorkspace.Dy[ 70 ]), &(acadoWorkspace.g[ 14 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 75 ]), &(acadoWorkspace.Dy[ 75 ]), &(acadoWorkspace.g[ 15 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 80 ]), &(acadoWorkspace.Dy[ 80 ]), &(acadoWorkspace.g[ 16 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 85 ]), &(acadoWorkspace.Dy[ 85 ]), &(acadoWorkspace.g[ 17 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 90 ]), &(acadoWorkspace.Dy[ 90 ]), &(acadoWorkspace.g[ 18 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 95 ]), &(acadoWorkspace.Dy[ 95 ]), &(acadoWorkspace.g[ 19 ]) );

acado_multQDy( acadoWorkspace.Q2, acadoWorkspace.Dy, acadoWorkspace.QDy );
acado_multQDy( &(acadoWorkspace.Q2[ 20 ]), &(acadoWorkspace.Dy[ 5 ]), &(acadoWorkspace.QDy[ 4 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 40 ]), &(acadoWorkspace.Dy[ 10 ]), &(acadoWorkspace.QDy[ 8 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 60 ]), &(acadoWorkspace.Dy[ 15 ]), &(acadoWorkspace.QDy[ 12 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 80 ]), &(acadoWorkspace.Dy[ 20 ]), &(acadoWorkspace.QDy[ 16 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 100 ]), &(acadoWorkspace.Dy[ 25 ]), &(acadoWorkspace.QDy[ 20 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 120 ]), &(acadoWorkspace.Dy[ 30 ]), &(acadoWorkspace.QDy[ 24 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 140 ]), &(acadoWorkspace.Dy[ 35 ]), &(acadoWorkspace.QDy[ 28 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 160 ]), &(acadoWorkspace.Dy[ 40 ]), &(acadoWorkspace.QDy[ 32 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 180 ]), &(acadoWorkspace.Dy[ 45 ]), &(acadoWorkspace.QDy[ 36 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 200 ]), &(acadoWorkspace.Dy[ 50 ]), &(acadoWorkspace.QDy[ 40 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 220 ]), &(acadoWorkspace.Dy[ 55 ]), &(acadoWorkspace.QDy[ 44 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 240 ]), &(acadoWorkspace.Dy[ 60 ]), &(acadoWorkspace.QDy[ 48 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 260 ]), &(acadoWorkspace.Dy[ 65 ]), &(acadoWorkspace.QDy[ 52 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 280 ]), &(acadoWorkspace.Dy[ 70 ]), &(acadoWorkspace.QDy[ 56 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 300 ]), &(acadoWorkspace.Dy[ 75 ]), &(acadoWorkspace.QDy[ 60 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 320 ]), &(acadoWorkspace.Dy[ 80 ]), &(acadoWorkspace.QDy[ 64 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 340 ]), &(acadoWorkspace.Dy[ 85 ]), &(acadoWorkspace.QDy[ 68 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 360 ]), &(acadoWorkspace.Dy[ 90 ]), &(acadoWorkspace.QDy[ 72 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 380 ]), &(acadoWorkspace.Dy[ 95 ]), &(acadoWorkspace.QDy[ 76 ]) );

acadoWorkspace.QDy[80] = + acadoWorkspace.QN2[0]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[1]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[2]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[3]*acadoWorkspace.DyN[3];
acadoWorkspace.QDy[81] = + acadoWorkspace.QN2[4]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[5]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[6]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[7]*acadoWorkspace.DyN[3];
acadoWorkspace.QDy[82] = + acadoWorkspace.QN2[8]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[9]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[10]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[11]*acadoWorkspace.DyN[3];
acadoWorkspace.QDy[83] = + acadoWorkspace.QN2[12]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[13]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[14]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[15]*acadoWorkspace.DyN[3];

acadoWorkspace.sbar[0] = acadoWorkspace.Dx0[0];
acadoWorkspace.sbar[1] = acadoWorkspace.Dx0[1];
acadoWorkspace.sbar[2] = acadoWorkspace.Dx0[2];
acadoWorkspace.sbar[3] = acadoWorkspace.Dx0[3];
acado_macASbar( acadoWorkspace.evGx, acadoWorkspace.sbar, &(acadoWorkspace.sbar[ 4 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 16 ]), &(acadoWorkspace.sbar[ 4 ]), &(acadoWorkspace.sbar[ 8 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 32 ]), &(acadoWorkspace.sbar[ 8 ]), &(acadoWorkspace.sbar[ 12 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 48 ]), &(acadoWorkspace.sbar[ 12 ]), &(acadoWorkspace.sbar[ 16 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 64 ]), &(acadoWorkspace.sbar[ 16 ]), &(acadoWorkspace.sbar[ 20 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 80 ]), &(acadoWorkspace.sbar[ 20 ]), &(acadoWorkspace.sbar[ 24 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 96 ]), &(acadoWorkspace.sbar[ 24 ]), &(acadoWorkspace.sbar[ 28 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 112 ]), &(acadoWorkspace.sbar[ 28 ]), &(acadoWorkspace.sbar[ 32 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 128 ]), &(acadoWorkspace.sbar[ 32 ]), &(acadoWorkspace.sbar[ 36 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 144 ]), &(acadoWorkspace.sbar[ 36 ]), &(acadoWorkspace.sbar[ 40 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 160 ]), &(acadoWorkspace.sbar[ 40 ]), &(acadoWorkspace.sbar[ 44 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 176 ]), &(acadoWorkspace.sbar[ 44 ]), &(acadoWorkspace.sbar[ 48 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 192 ]), &(acadoWorkspace.sbar[ 48 ]), &(acadoWorkspace.sbar[ 52 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 208 ]), &(acadoWorkspace.sbar[ 52 ]), &(acadoWorkspace.sbar[ 56 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 224 ]), &(acadoWorkspace.sbar[ 56 ]), &(acadoWorkspace.sbar[ 60 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 240 ]), &(acadoWorkspace.sbar[ 60 ]), &(acadoWorkspace.sbar[ 64 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 256 ]), &(acadoWorkspace.sbar[ 64 ]), &(acadoWorkspace.sbar[ 68 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 272 ]), &(acadoWorkspace.sbar[ 68 ]), &(acadoWorkspace.sbar[ 72 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.sbar[ 72 ]), &(acadoWorkspace.sbar[ 76 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 304 ]), &(acadoWorkspace.sbar[ 76 ]), &(acadoWorkspace.sbar[ 80 ]) );

acadoWorkspace.w1[0] = + acadoWorkspace.QN1[0]*acadoWorkspace.sbar[80] + acadoWorkspace.QN1[1]*acadoWorkspace.sbar[81] + acadoWorkspace.QN1[2]*acadoWorkspace.sbar[82] + acadoWorkspace.QN1[3]*acadoWorkspace.sbar[83] + acadoWorkspace.QDy[80];
acadoWorkspace.w1[1] = + acadoWorkspace.QN1[4]*acadoWorkspace.sbar[80] + acadoWorkspace.QN1[5]*acadoWorkspace.sbar[81] + acadoWorkspace.QN1[6]*acadoWorkspace.sbar[82] + acadoWorkspace.QN1[7]*acadoWorkspace.sbar[83] + acadoWorkspace.QDy[81];
acadoWorkspace.w1[2] = + acadoWorkspace.QN1[8]*acadoWorkspace.sbar[80] + acadoWorkspace.QN1[9]*acadoWorkspace.sbar[81] + acadoWorkspace.QN1[10]*acadoWorkspace.sbar[82] + acadoWorkspace.QN1[11]*acadoWorkspace.sbar[83] + acadoWorkspace.QDy[82];
acadoWorkspace.w1[3] = + acadoWorkspace.QN1[12]*acadoWorkspace.sbar[80] + acadoWorkspace.QN1[13]*acadoWorkspace.sbar[81] + acadoWorkspace.QN1[14]*acadoWorkspace.sbar[82] + acadoWorkspace.QN1[15]*acadoWorkspace.sbar[83] + acadoWorkspace.QDy[83];
acado_macBTw1( &(acadoWorkspace.evGu[ 76 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 19 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 304 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 76 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 304 ]), &(acadoWorkspace.sbar[ 76 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 72 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 18 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 288 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 72 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.sbar[ 72 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 68 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 17 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 272 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 68 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 272 ]), &(acadoWorkspace.sbar[ 68 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 64 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 16 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 256 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 64 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 256 ]), &(acadoWorkspace.sbar[ 64 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 60 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 15 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 240 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 60 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 240 ]), &(acadoWorkspace.sbar[ 60 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 56 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 14 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 224 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 56 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 224 ]), &(acadoWorkspace.sbar[ 56 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 52 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 13 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 208 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 52 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 208 ]), &(acadoWorkspace.sbar[ 52 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 48 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 12 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 192 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 48 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 192 ]), &(acadoWorkspace.sbar[ 48 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 44 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 11 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 176 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 44 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 176 ]), &(acadoWorkspace.sbar[ 44 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 40 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 10 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 160 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 40 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 160 ]), &(acadoWorkspace.sbar[ 40 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 36 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 9 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 144 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 36 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 144 ]), &(acadoWorkspace.sbar[ 36 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 32 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 8 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 128 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 32 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 128 ]), &(acadoWorkspace.sbar[ 32 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 28 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 7 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 112 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 28 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 112 ]), &(acadoWorkspace.sbar[ 28 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 24 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 6 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 96 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 24 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 96 ]), &(acadoWorkspace.sbar[ 24 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 20 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 5 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 80 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 20 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 80 ]), &(acadoWorkspace.sbar[ 20 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 16 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 4 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 64 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 16 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 64 ]), &(acadoWorkspace.sbar[ 16 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 12 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 3 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 48 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 12 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 48 ]), &(acadoWorkspace.sbar[ 12 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 8 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 2 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 32 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 8 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 32 ]), &(acadoWorkspace.sbar[ 8 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 4 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 1 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 16 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 4 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 16 ]), &(acadoWorkspace.sbar[ 4 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( acadoWorkspace.evGu, acadoWorkspace.w1, acadoWorkspace.g );

acadoWorkspace.lb[0] = acadoVariables.lbValues[0] - acadoVariables.u[0];
acadoWorkspace.lb[1] = acadoVariables.lbValues[1] - acadoVariables.u[1];
acadoWorkspace.lb[2] = acadoVariables.lbValues[2] - acadoVariables.u[2];
acadoWorkspace.lb[3] = acadoVariables.lbValues[3] - acadoVariables.u[3];
acadoWorkspace.lb[4] = acadoVariables.lbValues[4] - acadoVariables.u[4];
acadoWorkspace.lb[5] = acadoVariables.lbValues[5] - acadoVariables.u[5];
acadoWorkspace.lb[6] = acadoVariables.lbValues[6] - acadoVariables.u[6];
acadoWorkspace.lb[7] = acadoVariables.lbValues[7] - acadoVariables.u[7];
acadoWorkspace.lb[8] = acadoVariables.lbValues[8] - acadoVariables.u[8];
acadoWorkspace.lb[9] = acadoVariables.lbValues[9] - acadoVariables.u[9];
acadoWorkspace.lb[10] = acadoVariables.lbValues[10] - acadoVariables.u[10];
acadoWorkspace.lb[11] = acadoVariables.lbValues[11] - acadoVariables.u[11];
acadoWorkspace.lb[12] = acadoVariables.lbValues[12] - acadoVariables.u[12];
acadoWorkspace.lb[13] = acadoVariables.lbValues[13] - acadoVariables.u[13];
acadoWorkspace.lb[14] = acadoVariables.lbValues[14] - acadoVariables.u[14];
acadoWorkspace.lb[15] = acadoVariables.lbValues[15] - acadoVariables.u[15];
acadoWorkspace.lb[16] = acadoVariables.lbValues[16] - acadoVariables.u[16];
acadoWorkspace.lb[17] = acadoVariables.lbValues[17] - acadoVariables.u[17];
acadoWorkspace.lb[18] = acadoVariables.lbValues[18] - acadoVariables.u[18];
acadoWorkspace.lb[19] = acadoVariables.lbValues[19] - acadoVariables.u[19];
acadoWorkspace.ub[0] = acadoVariables.ubValues[0] - acadoVariables.u[0];
acadoWorkspace.ub[1] = acadoVariables.ubValues[1] - acadoVariables.u[1];
acadoWorkspace.ub[2] = acadoVariables.ubValues[2] - acadoVariables.u[2];
acadoWorkspace.ub[3] = acadoVariables.ubValues[3] - acadoVariables.u[3];
acadoWorkspace.ub[4] = acadoVariables.ubValues[4] - acadoVariables.u[4];
acadoWorkspace.ub[5] = acadoVariables.ubValues[5] - acadoVariables.u[5];
acadoWorkspace.ub[6] = acadoVariables.ubValues[6] - acadoVariables.u[6];
acadoWorkspace.ub[7] = acadoVariables.ubValues[7] - acadoVariables.u[7];
acadoWorkspace.ub[8] = acadoVariables.ubValues[8] - acadoVariables.u[8];
acadoWorkspace.ub[9] = acadoVariables.ubValues[9] - acadoVariables.u[9];
acadoWorkspace.ub[10] = acadoVariables.ubValues[10] - acadoVariables.u[10];
acadoWorkspace.ub[11] = acadoVariables.ubValues[11] - acadoVariables.u[11];
acadoWorkspace.ub[12] = acadoVariables.ubValues[12] - acadoVariables.u[12];
acadoWorkspace.ub[13] = acadoVariables.ubValues[13] - acadoVariables.u[13];
acadoWorkspace.ub[14] = acadoVariables.ubValues[14] - acadoVariables.u[14];
acadoWorkspace.ub[15] = acadoVariables.ubValues[15] - acadoVariables.u[15];
acadoWorkspace.ub[16] = acadoVariables.ubValues[16] - acadoVariables.u[16];
acadoWorkspace.ub[17] = acadoVariables.ubValues[17] - acadoVariables.u[17];
acadoWorkspace.ub[18] = acadoVariables.ubValues[18] - acadoVariables.u[18];
acadoWorkspace.ub[19] = acadoVariables.ubValues[19] - acadoVariables.u[19];

}

void acado_expand(  )
{
acadoVariables.u[0] += acadoWorkspace.x[0];
acadoVariables.u[1] += acadoWorkspace.x[1];
acadoVariables.u[2] += acadoWorkspace.x[2];
acadoVariables.u[3] += acadoWorkspace.x[3];
acadoVariables.u[4] += acadoWorkspace.x[4];
acadoVariables.u[5] += acadoWorkspace.x[5];
acadoVariables.u[6] += acadoWorkspace.x[6];
acadoVariables.u[7] += acadoWorkspace.x[7];
acadoVariables.u[8] += acadoWorkspace.x[8];
acadoVariables.u[9] += acadoWorkspace.x[9];
acadoVariables.u[10] += acadoWorkspace.x[10];
acadoVariables.u[11] += acadoWorkspace.x[11];
acadoVariables.u[12] += acadoWorkspace.x[12];
acadoVariables.u[13] += acadoWorkspace.x[13];
acadoVariables.u[14] += acadoWorkspace.x[14];
acadoVariables.u[15] += acadoWorkspace.x[15];
acadoVariables.u[16] += acadoWorkspace.x[16];
acadoVariables.u[17] += acadoWorkspace.x[17];
acadoVariables.u[18] += acadoWorkspace.x[18];
acadoVariables.u[19] += acadoWorkspace.x[19];
acadoWorkspace.sbar[0] = acadoWorkspace.Dx0[0];
acadoWorkspace.sbar[1] = acadoWorkspace.Dx0[1];
acadoWorkspace.sbar[2] = acadoWorkspace.Dx0[2];
acadoWorkspace.sbar[3] = acadoWorkspace.Dx0[3];
acadoWorkspace.sbar[4] = acadoWorkspace.d[0];
acadoWorkspace.sbar[5] = acadoWorkspace.d[1];
acadoWorkspace.sbar[6] = acadoWorkspace.d[2];
acadoWorkspace.sbar[7] = acadoWorkspace.d[3];
acadoWorkspace.sbar[8] = acadoWorkspace.d[4];
acadoWorkspace.sbar[9] = acadoWorkspace.d[5];
acadoWorkspace.sbar[10] = acadoWorkspace.d[6];
acadoWorkspace.sbar[11] = acadoWorkspace.d[7];
acadoWorkspace.sbar[12] = acadoWorkspace.d[8];
acadoWorkspace.sbar[13] = acadoWorkspace.d[9];
acadoWorkspace.sbar[14] = acadoWorkspace.d[10];
acadoWorkspace.sbar[15] = acadoWorkspace.d[11];
acadoWorkspace.sbar[16] = acadoWorkspace.d[12];
acadoWorkspace.sbar[17] = acadoWorkspace.d[13];
acadoWorkspace.sbar[18] = acadoWorkspace.d[14];
acadoWorkspace.sbar[19] = acadoWorkspace.d[15];
acadoWorkspace.sbar[20] = acadoWorkspace.d[16];
acadoWorkspace.sbar[21] = acadoWorkspace.d[17];
acadoWorkspace.sbar[22] = acadoWorkspace.d[18];
acadoWorkspace.sbar[23] = acadoWorkspace.d[19];
acadoWorkspace.sbar[24] = acadoWorkspace.d[20];
acadoWorkspace.sbar[25] = acadoWorkspace.d[21];
acadoWorkspace.sbar[26] = acadoWorkspace.d[22];
acadoWorkspace.sbar[27] = acadoWorkspace.d[23];
acadoWorkspace.sbar[28] = acadoWorkspace.d[24];
acadoWorkspace.sbar[29] = acadoWorkspace.d[25];
acadoWorkspace.sbar[30] = acadoWorkspace.d[26];
acadoWorkspace.sbar[31] = acadoWorkspace.d[27];
acadoWorkspace.sbar[32] = acadoWorkspace.d[28];
acadoWorkspace.sbar[33] = acadoWorkspace.d[29];
acadoWorkspace.sbar[34] = acadoWorkspace.d[30];
acadoWorkspace.sbar[35] = acadoWorkspace.d[31];
acadoWorkspace.sbar[36] = acadoWorkspace.d[32];
acadoWorkspace.sbar[37] = acadoWorkspace.d[33];
acadoWorkspace.sbar[38] = acadoWorkspace.d[34];
acadoWorkspace.sbar[39] = acadoWorkspace.d[35];
acadoWorkspace.sbar[40] = acadoWorkspace.d[36];
acadoWorkspace.sbar[41] = acadoWorkspace.d[37];
acadoWorkspace.sbar[42] = acadoWorkspace.d[38];
acadoWorkspace.sbar[43] = acadoWorkspace.d[39];
acadoWorkspace.sbar[44] = acadoWorkspace.d[40];
acadoWorkspace.sbar[45] = acadoWorkspace.d[41];
acadoWorkspace.sbar[46] = acadoWorkspace.d[42];
acadoWorkspace.sbar[47] = acadoWorkspace.d[43];
acadoWorkspace.sbar[48] = acadoWorkspace.d[44];
acadoWorkspace.sbar[49] = acadoWorkspace.d[45];
acadoWorkspace.sbar[50] = acadoWorkspace.d[46];
acadoWorkspace.sbar[51] = acadoWorkspace.d[47];
acadoWorkspace.sbar[52] = acadoWorkspace.d[48];
acadoWorkspace.sbar[53] = acadoWorkspace.d[49];
acadoWorkspace.sbar[54] = acadoWorkspace.d[50];
acadoWorkspace.sbar[55] = acadoWorkspace.d[51];
acadoWorkspace.sbar[56] = acadoWorkspace.d[52];
acadoWorkspace.sbar[57] = acadoWorkspace.d[53];
acadoWorkspace.sbar[58] = acadoWorkspace.d[54];
acadoWorkspace.sbar[59] = acadoWorkspace.d[55];
acadoWorkspace.sbar[60] = acadoWorkspace.d[56];
acadoWorkspace.sbar[61] = acadoWorkspace.d[57];
acadoWorkspace.sbar[62] = acadoWorkspace.d[58];
acadoWorkspace.sbar[63] = acadoWorkspace.d[59];
acadoWorkspace.sbar[64] = acadoWorkspace.d[60];
acadoWorkspace.sbar[65] = acadoWorkspace.d[61];
acadoWorkspace.sbar[66] = acadoWorkspace.d[62];
acadoWorkspace.sbar[67] = acadoWorkspace.d[63];
acadoWorkspace.sbar[68] = acadoWorkspace.d[64];
acadoWorkspace.sbar[69] = acadoWorkspace.d[65];
acadoWorkspace.sbar[70] = acadoWorkspace.d[66];
acadoWorkspace.sbar[71] = acadoWorkspace.d[67];
acadoWorkspace.sbar[72] = acadoWorkspace.d[68];
acadoWorkspace.sbar[73] = acadoWorkspace.d[69];
acadoWorkspace.sbar[74] = acadoWorkspace.d[70];
acadoWorkspace.sbar[75] = acadoWorkspace.d[71];
acadoWorkspace.sbar[76] = acadoWorkspace.d[72];
acadoWorkspace.sbar[77] = acadoWorkspace.d[73];
acadoWorkspace.sbar[78] = acadoWorkspace.d[74];
acadoWorkspace.sbar[79] = acadoWorkspace.d[75];
acadoWorkspace.sbar[80] = acadoWorkspace.d[76];
acadoWorkspace.sbar[81] = acadoWorkspace.d[77];
acadoWorkspace.sbar[82] = acadoWorkspace.d[78];
acadoWorkspace.sbar[83] = acadoWorkspace.d[79];
acado_expansionStep( acadoWorkspace.evGx, acadoWorkspace.evGu, acadoWorkspace.x, acadoWorkspace.sbar, &(acadoWorkspace.sbar[ 4 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 16 ]), &(acadoWorkspace.evGu[ 4 ]), &(acadoWorkspace.x[ 1 ]), &(acadoWorkspace.sbar[ 4 ]), &(acadoWorkspace.sbar[ 8 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 32 ]), &(acadoWorkspace.evGu[ 8 ]), &(acadoWorkspace.x[ 2 ]), &(acadoWorkspace.sbar[ 8 ]), &(acadoWorkspace.sbar[ 12 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 48 ]), &(acadoWorkspace.evGu[ 12 ]), &(acadoWorkspace.x[ 3 ]), &(acadoWorkspace.sbar[ 12 ]), &(acadoWorkspace.sbar[ 16 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 64 ]), &(acadoWorkspace.evGu[ 16 ]), &(acadoWorkspace.x[ 4 ]), &(acadoWorkspace.sbar[ 16 ]), &(acadoWorkspace.sbar[ 20 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 80 ]), &(acadoWorkspace.evGu[ 20 ]), &(acadoWorkspace.x[ 5 ]), &(acadoWorkspace.sbar[ 20 ]), &(acadoWorkspace.sbar[ 24 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 96 ]), &(acadoWorkspace.evGu[ 24 ]), &(acadoWorkspace.x[ 6 ]), &(acadoWorkspace.sbar[ 24 ]), &(acadoWorkspace.sbar[ 28 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 112 ]), &(acadoWorkspace.evGu[ 28 ]), &(acadoWorkspace.x[ 7 ]), &(acadoWorkspace.sbar[ 28 ]), &(acadoWorkspace.sbar[ 32 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 128 ]), &(acadoWorkspace.evGu[ 32 ]), &(acadoWorkspace.x[ 8 ]), &(acadoWorkspace.sbar[ 32 ]), &(acadoWorkspace.sbar[ 36 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 144 ]), &(acadoWorkspace.evGu[ 36 ]), &(acadoWorkspace.x[ 9 ]), &(acadoWorkspace.sbar[ 36 ]), &(acadoWorkspace.sbar[ 40 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 160 ]), &(acadoWorkspace.evGu[ 40 ]), &(acadoWorkspace.x[ 10 ]), &(acadoWorkspace.sbar[ 40 ]), &(acadoWorkspace.sbar[ 44 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 176 ]), &(acadoWorkspace.evGu[ 44 ]), &(acadoWorkspace.x[ 11 ]), &(acadoWorkspace.sbar[ 44 ]), &(acadoWorkspace.sbar[ 48 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 192 ]), &(acadoWorkspace.evGu[ 48 ]), &(acadoWorkspace.x[ 12 ]), &(acadoWorkspace.sbar[ 48 ]), &(acadoWorkspace.sbar[ 52 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 208 ]), &(acadoWorkspace.evGu[ 52 ]), &(acadoWorkspace.x[ 13 ]), &(acadoWorkspace.sbar[ 52 ]), &(acadoWorkspace.sbar[ 56 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 224 ]), &(acadoWorkspace.evGu[ 56 ]), &(acadoWorkspace.x[ 14 ]), &(acadoWorkspace.sbar[ 56 ]), &(acadoWorkspace.sbar[ 60 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 240 ]), &(acadoWorkspace.evGu[ 60 ]), &(acadoWorkspace.x[ 15 ]), &(acadoWorkspace.sbar[ 60 ]), &(acadoWorkspace.sbar[ 64 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 256 ]), &(acadoWorkspace.evGu[ 64 ]), &(acadoWorkspace.x[ 16 ]), &(acadoWorkspace.sbar[ 64 ]), &(acadoWorkspace.sbar[ 68 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 272 ]), &(acadoWorkspace.evGu[ 68 ]), &(acadoWorkspace.x[ 17 ]), &(acadoWorkspace.sbar[ 68 ]), &(acadoWorkspace.sbar[ 72 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.evGu[ 72 ]), &(acadoWorkspace.x[ 18 ]), &(acadoWorkspace.sbar[ 72 ]), &(acadoWorkspace.sbar[ 76 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 304 ]), &(acadoWorkspace.evGu[ 76 ]), &(acadoWorkspace.x[ 19 ]), &(acadoWorkspace.sbar[ 76 ]), &(acadoWorkspace.sbar[ 80 ]) );
acadoVariables.x[0] += acadoWorkspace.sbar[0];
acadoVariables.x[1] += acadoWorkspace.sbar[1];
acadoVariables.x[2] += acadoWorkspace.sbar[2];
acadoVariables.x[3] += acadoWorkspace.sbar[3];
acadoVariables.x[4] += acadoWorkspace.sbar[4];
acadoVariables.x[5] += acadoWorkspace.sbar[5];
acadoVariables.x[6] += acadoWorkspace.sbar[6];
acadoVariables.x[7] += acadoWorkspace.sbar[7];
acadoVariables.x[8] += acadoWorkspace.sbar[8];
acadoVariables.x[9] += acadoWorkspace.sbar[9];
acadoVariables.x[10] += acadoWorkspace.sbar[10];
acadoVariables.x[11] += acadoWorkspace.sbar[11];
acadoVariables.x[12] += acadoWorkspace.sbar[12];
acadoVariables.x[13] += acadoWorkspace.sbar[13];
acadoVariables.x[14] += acadoWorkspace.sbar[14];
acadoVariables.x[15] += acadoWorkspace.sbar[15];
acadoVariables.x[16] += acadoWorkspace.sbar[16];
acadoVariables.x[17] += acadoWorkspace.sbar[17];
acadoVariables.x[18] += acadoWorkspace.sbar[18];
acadoVariables.x[19] += acadoWorkspace.sbar[19];
acadoVariables.x[20] += acadoWorkspace.sbar[20];
acadoVariables.x[21] += acadoWorkspace.sbar[21];
acadoVariables.x[22] += acadoWorkspace.sbar[22];
acadoVariables.x[23] += acadoWorkspace.sbar[23];
acadoVariables.x[24] += acadoWorkspace.sbar[24];
acadoVariables.x[25] += acadoWorkspace.sbar[25];
acadoVariables.x[26] += acadoWorkspace.sbar[26];
acadoVariables.x[27] += acadoWorkspace.sbar[27];
acadoVariables.x[28] += acadoWorkspace.sbar[28];
acadoVariables.x[29] += acadoWorkspace.sbar[29];
acadoVariables.x[30] += acadoWorkspace.sbar[30];
acadoVariables.x[31] += acadoWorkspace.sbar[31];
acadoVariables.x[32] += acadoWorkspace.sbar[32];
acadoVariables.x[33] += acadoWorkspace.sbar[33];
acadoVariables.x[34] += acadoWorkspace.sbar[34];
acadoVariables.x[35] += acadoWorkspace.sbar[35];
acadoVariables.x[36] += acadoWorkspace.sbar[36];
acadoVariables.x[37] += acadoWorkspace.sbar[37];
acadoVariables.x[38] += acadoWorkspace.sbar[38];
acadoVariables.x[39] += acadoWorkspace.sbar[39];
acadoVariables.x[40] += acadoWorkspace.sbar[40];
acadoVariables.x[41] += acadoWorkspace.sbar[41];
acadoVariables.x[42] += acadoWorkspace.sbar[42];
acadoVariables.x[43] += acadoWorkspace.sbar[43];
acadoVariables.x[44] += acadoWorkspace.sbar[44];
acadoVariables.x[45] += acadoWorkspace.sbar[45];
acadoVariables.x[46] += acadoWorkspace.sbar[46];
acadoVariables.x[47] += acadoWorkspace.sbar[47];
acadoVariables.x[48] += acadoWorkspace.sbar[48];
acadoVariables.x[49] += acadoWorkspace.sbar[49];
acadoVariables.x[50] += acadoWorkspace.sbar[50];
acadoVariables.x[51] += acadoWorkspace.sbar[51];
acadoVariables.x[52] += acadoWorkspace.sbar[52];
acadoVariables.x[53] += acadoWorkspace.sbar[53];
acadoVariables.x[54] += acadoWorkspace.sbar[54];
acadoVariables.x[55] += acadoWorkspace.sbar[55];
acadoVariables.x[56] += acadoWorkspace.sbar[56];
acadoVariables.x[57] += acadoWorkspace.sbar[57];
acadoVariables.x[58] += acadoWorkspace.sbar[58];
acadoVariables.x[59] += acadoWorkspace.sbar[59];
acadoVariables.x[60] += acadoWorkspace.sbar[60];
acadoVariables.x[61] += acadoWorkspace.sbar[61];
acadoVariables.x[62] += acadoWorkspace.sbar[62];
acadoVariables.x[63] += acadoWorkspace.sbar[63];
acadoVariables.x[64] += acadoWorkspace.sbar[64];
acadoVariables.x[65] += acadoWorkspace.sbar[65];
acadoVariables.x[66] += acadoWorkspace.sbar[66];
acadoVariables.x[67] += acadoWorkspace.sbar[67];
acadoVariables.x[68] += acadoWorkspace.sbar[68];
acadoVariables.x[69] += acadoWorkspace.sbar[69];
acadoVariables.x[70] += acadoWorkspace.sbar[70];
acadoVariables.x[71] += acadoWorkspace.sbar[71];
acadoVariables.x[72] += acadoWorkspace.sbar[72];
acadoVariables.x[73] += acadoWorkspace.sbar[73];
acadoVariables.x[74] += acadoWorkspace.sbar[74];
acadoVariables.x[75] += acadoWorkspace.sbar[75];
acadoVariables.x[76] += acadoWorkspace.sbar[76];
acadoVariables.x[77] += acadoWorkspace.sbar[77];
acadoVariables.x[78] += acadoWorkspace.sbar[78];
acadoVariables.x[79] += acadoWorkspace.sbar[79];
acadoVariables.x[80] += acadoWorkspace.sbar[80];
acadoVariables.x[81] += acadoWorkspace.sbar[81];
acadoVariables.x[82] += acadoWorkspace.sbar[82];
acadoVariables.x[83] += acadoWorkspace.sbar[83];
}

int acado_preparationStep(  )
{
int ret;

ret = acado_modelSimulation();
acado_evaluateObjective(  );
acado_condensePrep(  );
return ret;
}

int acado_feedbackStep(  )
{
int tmp;

acado_condenseFdb(  );

tmp = acado_solve( );

acado_expand(  );
return tmp;
}

int acado_initializeSolver(  )
{
int ret;

/* This is a function which must be called once before any other function call! */


ret = 0;

memset(&acadoWorkspace, 0, sizeof( acadoWorkspace ));
acadoVariables.lbValues[0] = -2.0000000000000000e+00;
acadoVariables.lbValues[1] = -2.0000000000000000e+00;
acadoVariables.lbValues[2] = -2.0000000000000000e+00;
acadoVariables.lbValues[3] = -2.0000000000000000e+00;
acadoVariables.lbValues[4] = -2.0000000000000000e+00;
acadoVariables.lbValues[5] = -2.0000000000000000e+00;
acadoVariables.lbValues[6] = -2.0000000000000000e+00;
acadoVariables.lbValues[7] = -2.0000000000000000e+00;
acadoVariables.lbValues[8] = -2.0000000000000000e+00;
acadoVariables.lbValues[9] = -2.0000000000000000e+00;
acadoVariables.lbValues[10] = -2.0000000000000000e+00;
acadoVariables.lbValues[11] = -2.0000000000000000e+00;
acadoVariables.lbValues[12] = -2.0000000000000000e+00;
acadoVariables.lbValues[13] = -2.0000000000000000e+00;
acadoVariables.lbValues[14] = -2.0000000000000000e+00;
acadoVariables.lbValues[15] = -2.0000000000000000e+00;
acadoVariables.lbValues[16] = -2.0000000000000000e+00;
acadoVariables.lbValues[17] = -2.0000000000000000e+00;
acadoVariables.lbValues[18] = -2.0000000000000000e+00;
acadoVariables.lbValues[19] = -2.0000000000000000e+00;
acadoVariables.ubValues[0] = 2.0000000000000000e+00;
acadoVariables.ubValues[1] = 2.0000000000000000e+00;
acadoVariables.ubValues[2] = 2.0000000000000000e+00;
acadoVariables.ubValues[3] = 2.0000000000000000e+00;
acadoVariables.ubValues[4] = 2.0000000000000000e+00;
acadoVariables.ubValues[5] = 2.0000000000000000e+00;
acadoVariables.ubValues[6] = 2.0000000000000000e+00;
acadoVariables.ubValues[7] = 2.0000000000000000e+00;
acadoVariables.ubValues[8] = 2.0000000000000000e+00;
acadoVariables.ubValues[9] = 2.0000000000000000e+00;
acadoVariables.ubValues[10] = 2.0000000000000000e+00;
acadoVariables.ubValues[11] = 2.0000000000000000e+00;
acadoVariables.ubValues[12] = 2.0000000000000000e+00;
acadoVariables.ubValues[13] = 2.0000000000000000e+00;
acadoVariables.ubValues[14] = 2.0000000000000000e+00;
acadoVariables.ubValues[15] = 2.0000000000000000e+00;
acadoVariables.ubValues[16] = 2.0000000000000000e+00;
acadoVariables.ubValues[17] = 2.0000000000000000e+00;
acadoVariables.ubValues[18] = 2.0000000000000000e+00;
acadoVariables.ubValues[19] = 2.0000000000000000e+00;
return ret;
}

void acado_initializeNodesByForwardSimulation(  )
{
int index;
for (index = 0; index < 20; ++index)
{
state[0] = acadoVariables.x[index * 4];
state[1] = acadoVariables.x[index * 4 + 1];
state[2] = acadoVariables.x[index * 4 + 2];
state[3] = acadoVariables.x[index * 4 + 3];
state[24] = acadoVariables.u[index];

acado_integrate(state, index == 0);

acadoVariables.x[index * 4 + 4] = state[0];
acadoVariables.x[index * 4 + 5] = state[1];
acadoVariables.x[index * 4 + 6] = state[2];
acadoVariables.x[index * 4 + 7] = state[3];
}
}

void acado_shiftStates( int strategy, real_t* const xEnd, real_t* const uEnd )
{
int index;
for (index = 0; index < 20; ++index)
{
acadoVariables.x[index * 4] = acadoVariables.x[index * 4 + 4];
acadoVariables.x[index * 4 + 1] = acadoVariables.x[index * 4 + 5];
acadoVariables.x[index * 4 + 2] = acadoVariables.x[index * 4 + 6];
acadoVariables.x[index * 4 + 3] = acadoVariables.x[index * 4 + 7];
}

if (strategy == 1 && xEnd != 0)
{
acadoVariables.x[80] = xEnd[0];
acadoVariables.x[81] = xEnd[1];
acadoVariables.x[82] = xEnd[2];
acadoVariables.x[83] = xEnd[3];
}
else if (strategy == 2) 
{
state[0] = acadoVariables.x[80];
state[1] = acadoVariables.x[81];
state[2] = acadoVariables.x[82];
state[3] = acadoVariables.x[83];
if (uEnd != 0)
{
state[24] = uEnd[0];
}
else
{
state[24] = acadoVariables.u[19];
}

acado_integrate(state, 1);

acadoVariables.x[80] = state[0];
acadoVariables.x[81] = state[1];
acadoVariables.x[82] = state[2];
acadoVariables.x[83] = state[3];
}
}

void acado_shiftControls( real_t* const uEnd )
{
int index;
for (index = 0; index < 19; ++index)
{
acadoVariables.u[index] = acadoVariables.u[index + 1];
}

if (uEnd != 0)
{
acadoVariables.u[19] = uEnd[0];
}
}

real_t acado_getKKT(  )
{
real_t kkt;

int index;
real_t prd;

kkt = + acadoWorkspace.g[0]*acadoWorkspace.x[0] + acadoWorkspace.g[1]*acadoWorkspace.x[1] + acadoWorkspace.g[2]*acadoWorkspace.x[2] + acadoWorkspace.g[3]*acadoWorkspace.x[3] + acadoWorkspace.g[4]*acadoWorkspace.x[4] + acadoWorkspace.g[5]*acadoWorkspace.x[5] + acadoWorkspace.g[6]*acadoWorkspace.x[6] + acadoWorkspace.g[7]*acadoWorkspace.x[7] + acadoWorkspace.g[8]*acadoWorkspace.x[8] + acadoWorkspace.g[9]*acadoWorkspace.x[9] + acadoWorkspace.g[10]*acadoWorkspace.x[10] + acadoWorkspace.g[11]*acadoWorkspace.x[11] + acadoWorkspace.g[12]*acadoWorkspace.x[12] + acadoWorkspace.g[13]*acadoWorkspace.x[13] + acadoWorkspace.g[14]*acadoWorkspace.x[14] + acadoWorkspace.g[15]*acadoWorkspace.x[15] + acadoWorkspace.g[16]*acadoWorkspace.x[16] + acadoWorkspace.g[17]*acadoWorkspace.x[17] + acadoWorkspace.g[18]*acadoWorkspace.x[18] + acadoWorkspace.g[19]*acadoWorkspace.x[19];
kkt = fabs( kkt );
for (index = 0; index < 20; ++index)
{
prd = acadoWorkspace.y[index];
if (prd > 1e-12)
kkt += fabs(acadoWorkspace.lb[index] * prd);
else if (prd < -1e-12)
kkt += fabs(acadoWorkspace.ub[index] * prd);
}
return kkt;
}

real_t acado_getObjective(  )
{
real_t objVal;

int lRun1;
/** Row vector of size: 5 */
real_t tmpDy[ 5 ];

/** Row vector of size: 4 */
real_t tmpDyN[ 4 ];

for (lRun1 = 0; lRun1 < 20; ++lRun1)
{
acadoWorkspace.objValueIn[0] = acadoVariables.x[lRun1 * 4];
acadoWorkspace.objValueIn[1] = acadoVariables.x[lRun1 * 4 + 1];
acadoWorkspace.objValueIn[2] = acadoVariables.x[lRun1 * 4 + 2];
acadoWorkspace.objValueIn[3] = acadoVariables.x[lRun1 * 4 + 3];
acadoWorkspace.objValueIn[4] = acadoVariables.u[lRun1];

acado_evaluateLSQ( acadoWorkspace.objValueIn, acadoWorkspace.objValueOut );
acadoWorkspace.Dy[lRun1 * 5] = acadoWorkspace.objValueOut[0] - acadoVariables.y[lRun1 * 5];
acadoWorkspace.Dy[lRun1 * 5 + 1] = acadoWorkspace.objValueOut[1] - acadoVariables.y[lRun1 * 5 + 1];
acadoWorkspace.Dy[lRun1 * 5 + 2] = acadoWorkspace.objValueOut[2] - acadoVariables.y[lRun1 * 5 + 2];
acadoWorkspace.Dy[lRun1 * 5 + 3] = acadoWorkspace.objValueOut[3] - acadoVariables.y[lRun1 * 5 + 3];
acadoWorkspace.Dy[lRun1 * 5 + 4] = acadoWorkspace.objValueOut[4] - acadoVariables.y[lRun1 * 5 + 4];
}
acadoWorkspace.objValueIn[0] = acadoVariables.x[80];
acadoWorkspace.objValueIn[1] = acadoVariables.x[81];
acadoWorkspace.objValueIn[2] = acadoVariables.x[82];
acadoWorkspace.objValueIn[3] = acadoVariables.x[83];
acado_evaluateLSQEndTerm( acadoWorkspace.objValueIn, acadoWorkspace.objValueOut );
acadoWorkspace.DyN[0] = acadoWorkspace.objValueOut[0] - acadoVariables.yN[0];
acadoWorkspace.DyN[1] = acadoWorkspace.objValueOut[1] - acadoVariables.yN[1];
acadoWorkspace.DyN[2] = acadoWorkspace.objValueOut[2] - acadoVariables.yN[2];
acadoWorkspace.DyN[3] = acadoWorkspace.objValueOut[3] - acadoVariables.yN[3];
objVal = 0.0000000000000000e+00;
for (lRun1 = 0; lRun1 < 20; ++lRun1)
{
tmpDy[0] = + acadoWorkspace.Dy[lRun1 * 5]*acadoVariables.W[lRun1 * 25] + acadoWorkspace.Dy[lRun1 * 5 + 1]*acadoVariables.W[lRun1 * 25 + 5] + acadoWorkspace.Dy[lRun1 * 5 + 2]*acadoVariables.W[lRun1 * 25 + 10] + acadoWorkspace.Dy[lRun1 * 5 + 3]*acadoVariables.W[lRun1 * 25 + 15] + acadoWorkspace.Dy[lRun1 * 5 + 4]*acadoVariables.W[lRun1 * 25 + 20];
tmpDy[1] = + acadoWorkspace.Dy[lRun1 * 5]*acadoVariables.W[lRun1 * 25 + 1] + acadoWorkspace.Dy[lRun1 * 5 + 1]*acadoVariables.W[lRun1 * 25 + 6] + acadoWorkspace.Dy[lRun1 * 5 + 2]*acadoVariables.W[lRun1 * 25 + 11] + acadoWorkspace.Dy[lRun1 * 5 + 3]*acadoVariables.W[lRun1 * 25 + 16] + acadoWorkspace.Dy[lRun1 * 5 + 4]*acadoVariables.W[lRun1 * 25 + 21];
tmpDy[2] = + acadoWorkspace.Dy[lRun1 * 5]*acadoVariables.W[lRun1 * 25 + 2] + acadoWorkspace.Dy[lRun1 * 5 + 1]*acadoVariables.W[lRun1 * 25 + 7] + acadoWorkspace.Dy[lRun1 * 5 + 2]*acadoVariables.W[lRun1 * 25 + 12] + acadoWorkspace.Dy[lRun1 * 5 + 3]*acadoVariables.W[lRun1 * 25 + 17] + acadoWorkspace.Dy[lRun1 * 5 + 4]*acadoVariables.W[lRun1 * 25 + 22];
tmpDy[3] = + acadoWorkspace.Dy[lRun1 * 5]*acadoVariables.W[lRun1 * 25 + 3] + acadoWorkspace.Dy[lRun1 * 5 + 1]*acadoVariables.W[lRun1 * 25 + 8] + acadoWorkspace.Dy[lRun1 * 5 + 2]*acadoVariables.W[lRun1 * 25 + 13] + acadoWorkspace.Dy[lRun1 * 5 + 3]*acadoVariables.W[lRun1 * 25 + 18] + acadoWorkspace.Dy[lRun1 * 5 + 4]*acadoVariables.W[lRun1 * 25 + 23];
tmpDy[4] = + acadoWorkspace.Dy[lRun1 * 5]*acadoVariables.W[lRun1 * 25 + 4] + acadoWorkspace.Dy[lRun1 * 5 + 1]*acadoVariables.W[lRun1 * 25 + 9] + acadoWorkspace.Dy[lRun1 * 5 + 2]*acadoVariables.W[lRun1 * 25 + 14] + acadoWorkspace.Dy[lRun1 * 5 + 3]*acadoVariables.W[lRun1 * 25 + 19] + acadoWorkspace.Dy[lRun1 * 5 + 4]*acadoVariables.W[lRun1 * 25 + 24];
objVal += + acadoWorkspace.Dy[lRun1 * 5]*tmpDy[0] + acadoWorkspace.Dy[lRun1 * 5 + 1]*tmpDy[1] + acadoWorkspace.Dy[lRun1 * 5 + 2]*tmpDy[2] + acadoWorkspace.Dy[lRun1 * 5 + 3]*tmpDy[3] + acadoWorkspace.Dy[lRun1 * 5 + 4]*tmpDy[4];
}

tmpDyN[0] = + acadoWorkspace.DyN[0]*acadoVariables.WN[0];
tmpDyN[1] = + acadoWorkspace.DyN[1]*acadoVariables.WN[5];
tmpDyN[2] = + acadoWorkspace.DyN[2]*acadoVariables.WN[10];
tmpDyN[3] = + acadoWorkspace.DyN[3]*acadoVariables.WN[15];
objVal += + acadoWorkspace.DyN[0]*tmpDyN[0] + acadoWorkspace.DyN[1]*tmpDyN[1] + acadoWorkspace.DyN[2]*tmpDyN[2] + acadoWorkspace.DyN[3]*tmpDyN[3];

objVal *= 0.5;
return objVal;
}

