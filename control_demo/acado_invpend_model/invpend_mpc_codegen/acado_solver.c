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
for (lRun1 = 0; lRun1 < 200; ++lRun1)
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
for (runObj = 0; runObj < 200; ++runObj)
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
acadoWorkspace.objValueIn[0] = acadoVariables.x[800];
acadoWorkspace.objValueIn[1] = acadoVariables.x[801];
acadoWorkspace.objValueIn[2] = acadoVariables.x[802];
acadoWorkspace.objValueIn[3] = acadoVariables.x[803];
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
acadoWorkspace.H[(iRow * 200) + (iCol)] = + Gu1[0]*Gu2[0] + Gu1[1]*Gu2[1] + Gu1[2]*Gu2[2] + Gu1[3]*Gu2[3];
}

void acado_multBTW1_R1( real_t* const R11, real_t* const Gu1, real_t* const Gu2, int iRow )
{
acadoWorkspace.H[iRow * 201] = + Gu1[0]*Gu2[0] + Gu1[1]*Gu2[1] + Gu1[2]*Gu2[2] + Gu1[3]*Gu2[3] + R11[0];
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
acadoWorkspace.H[(iRow * 200) + (iCol)] = acadoWorkspace.H[(iCol * 200) + (iRow)];
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
for (lRun2 = 0; lRun2 < 200; ++lRun2)
{
lRun3 = ((lRun2) * (lRun2 * -1 + 401)) / (2);
acado_moveGuE( &(acadoWorkspace.evGu[ lRun2 * 4 ]), &(acadoWorkspace.E[ lRun3 * 4 ]) );
for (lRun1 = 1; lRun1 < lRun2 * -1 + 200; ++lRun1)
{
acado_multGxGu( &(acadoWorkspace.evGx[ ((((lRun2) + (lRun1)) * (4)) * (4)) + (0) ]), &(acadoWorkspace.E[ (((((lRun3) + (lRun1)) - (1)) * (4)) * (1)) + (0) ]), &(acadoWorkspace.E[ ((((lRun3) + (lRun1)) * (4)) * (1)) + (0) ]) );
}

acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ ((((((lRun3) - (lRun2)) + (200)) - (1)) * (4)) * (1)) + (0) ]), acadoWorkspace.W1 );
for (lRun1 = 199; lRun2 < lRun1; --lRun1)
{
acado_multBTW1( &(acadoWorkspace.evGu[ lRun1 * 4 ]), acadoWorkspace.W1, lRun1, lRun2 );
acado_multGxTGu( &(acadoWorkspace.evGx[ lRun1 * 16 ]), acadoWorkspace.W1, acadoWorkspace.W2 );
acado_multQEW2( &(acadoWorkspace.Q1[ lRun1 * 16 ]), &(acadoWorkspace.E[ ((((((lRun3) + (lRun1)) - (lRun2)) - (1)) * (4)) * (1)) + (0) ]), acadoWorkspace.W2, acadoWorkspace.W1 );
}
acado_multBTW1_R1( &(acadoWorkspace.R1[ lRun2 ]), &(acadoWorkspace.evGu[ lRun2 * 4 ]), acadoWorkspace.W1, lRun2 );
}

for (lRun1 = 0; lRun1 < 200; ++lRun1)
{
for (lRun2 = 0; lRun2 < lRun1; ++lRun2)
{
acado_copyHTH( lRun2, lRun1 );
}
}

for (lRun1 = 0; lRun1 < 800; ++lRun1)
acadoWorkspace.sbar[lRun1 + 4] = acadoWorkspace.d[lRun1];


}

void acado_condenseFdb(  )
{
int lRun1;
acadoWorkspace.Dx0[0] = acadoVariables.x0[0] - acadoVariables.x[0];
acadoWorkspace.Dx0[1] = acadoVariables.x0[1] - acadoVariables.x[1];
acadoWorkspace.Dx0[2] = acadoVariables.x0[2] - acadoVariables.x[2];
acadoWorkspace.Dx0[3] = acadoVariables.x0[3] - acadoVariables.x[3];
for (lRun1 = 0; lRun1 < 1000; ++lRun1)
acadoWorkspace.Dy[lRun1] -= acadoVariables.y[lRun1];

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
acado_multRDy( &(acadoWorkspace.R2[ 100 ]), &(acadoWorkspace.Dy[ 100 ]), &(acadoWorkspace.g[ 20 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 105 ]), &(acadoWorkspace.Dy[ 105 ]), &(acadoWorkspace.g[ 21 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 110 ]), &(acadoWorkspace.Dy[ 110 ]), &(acadoWorkspace.g[ 22 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 115 ]), &(acadoWorkspace.Dy[ 115 ]), &(acadoWorkspace.g[ 23 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 120 ]), &(acadoWorkspace.Dy[ 120 ]), &(acadoWorkspace.g[ 24 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 125 ]), &(acadoWorkspace.Dy[ 125 ]), &(acadoWorkspace.g[ 25 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 130 ]), &(acadoWorkspace.Dy[ 130 ]), &(acadoWorkspace.g[ 26 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 135 ]), &(acadoWorkspace.Dy[ 135 ]), &(acadoWorkspace.g[ 27 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 140 ]), &(acadoWorkspace.Dy[ 140 ]), &(acadoWorkspace.g[ 28 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 145 ]), &(acadoWorkspace.Dy[ 145 ]), &(acadoWorkspace.g[ 29 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 150 ]), &(acadoWorkspace.Dy[ 150 ]), &(acadoWorkspace.g[ 30 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 155 ]), &(acadoWorkspace.Dy[ 155 ]), &(acadoWorkspace.g[ 31 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 160 ]), &(acadoWorkspace.Dy[ 160 ]), &(acadoWorkspace.g[ 32 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 165 ]), &(acadoWorkspace.Dy[ 165 ]), &(acadoWorkspace.g[ 33 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 170 ]), &(acadoWorkspace.Dy[ 170 ]), &(acadoWorkspace.g[ 34 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 175 ]), &(acadoWorkspace.Dy[ 175 ]), &(acadoWorkspace.g[ 35 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 180 ]), &(acadoWorkspace.Dy[ 180 ]), &(acadoWorkspace.g[ 36 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 185 ]), &(acadoWorkspace.Dy[ 185 ]), &(acadoWorkspace.g[ 37 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 190 ]), &(acadoWorkspace.Dy[ 190 ]), &(acadoWorkspace.g[ 38 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 195 ]), &(acadoWorkspace.Dy[ 195 ]), &(acadoWorkspace.g[ 39 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 200 ]), &(acadoWorkspace.Dy[ 200 ]), &(acadoWorkspace.g[ 40 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 205 ]), &(acadoWorkspace.Dy[ 205 ]), &(acadoWorkspace.g[ 41 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 210 ]), &(acadoWorkspace.Dy[ 210 ]), &(acadoWorkspace.g[ 42 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 215 ]), &(acadoWorkspace.Dy[ 215 ]), &(acadoWorkspace.g[ 43 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 220 ]), &(acadoWorkspace.Dy[ 220 ]), &(acadoWorkspace.g[ 44 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 225 ]), &(acadoWorkspace.Dy[ 225 ]), &(acadoWorkspace.g[ 45 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 230 ]), &(acadoWorkspace.Dy[ 230 ]), &(acadoWorkspace.g[ 46 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 235 ]), &(acadoWorkspace.Dy[ 235 ]), &(acadoWorkspace.g[ 47 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 240 ]), &(acadoWorkspace.Dy[ 240 ]), &(acadoWorkspace.g[ 48 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 245 ]), &(acadoWorkspace.Dy[ 245 ]), &(acadoWorkspace.g[ 49 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 250 ]), &(acadoWorkspace.Dy[ 250 ]), &(acadoWorkspace.g[ 50 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 255 ]), &(acadoWorkspace.Dy[ 255 ]), &(acadoWorkspace.g[ 51 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 260 ]), &(acadoWorkspace.Dy[ 260 ]), &(acadoWorkspace.g[ 52 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 265 ]), &(acadoWorkspace.Dy[ 265 ]), &(acadoWorkspace.g[ 53 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 270 ]), &(acadoWorkspace.Dy[ 270 ]), &(acadoWorkspace.g[ 54 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 275 ]), &(acadoWorkspace.Dy[ 275 ]), &(acadoWorkspace.g[ 55 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 280 ]), &(acadoWorkspace.Dy[ 280 ]), &(acadoWorkspace.g[ 56 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 285 ]), &(acadoWorkspace.Dy[ 285 ]), &(acadoWorkspace.g[ 57 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 290 ]), &(acadoWorkspace.Dy[ 290 ]), &(acadoWorkspace.g[ 58 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 295 ]), &(acadoWorkspace.Dy[ 295 ]), &(acadoWorkspace.g[ 59 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 300 ]), &(acadoWorkspace.Dy[ 300 ]), &(acadoWorkspace.g[ 60 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 305 ]), &(acadoWorkspace.Dy[ 305 ]), &(acadoWorkspace.g[ 61 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 310 ]), &(acadoWorkspace.Dy[ 310 ]), &(acadoWorkspace.g[ 62 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 315 ]), &(acadoWorkspace.Dy[ 315 ]), &(acadoWorkspace.g[ 63 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 320 ]), &(acadoWorkspace.Dy[ 320 ]), &(acadoWorkspace.g[ 64 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 325 ]), &(acadoWorkspace.Dy[ 325 ]), &(acadoWorkspace.g[ 65 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 330 ]), &(acadoWorkspace.Dy[ 330 ]), &(acadoWorkspace.g[ 66 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 335 ]), &(acadoWorkspace.Dy[ 335 ]), &(acadoWorkspace.g[ 67 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 340 ]), &(acadoWorkspace.Dy[ 340 ]), &(acadoWorkspace.g[ 68 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 345 ]), &(acadoWorkspace.Dy[ 345 ]), &(acadoWorkspace.g[ 69 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 350 ]), &(acadoWorkspace.Dy[ 350 ]), &(acadoWorkspace.g[ 70 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 355 ]), &(acadoWorkspace.Dy[ 355 ]), &(acadoWorkspace.g[ 71 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 360 ]), &(acadoWorkspace.Dy[ 360 ]), &(acadoWorkspace.g[ 72 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 365 ]), &(acadoWorkspace.Dy[ 365 ]), &(acadoWorkspace.g[ 73 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 370 ]), &(acadoWorkspace.Dy[ 370 ]), &(acadoWorkspace.g[ 74 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 375 ]), &(acadoWorkspace.Dy[ 375 ]), &(acadoWorkspace.g[ 75 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 380 ]), &(acadoWorkspace.Dy[ 380 ]), &(acadoWorkspace.g[ 76 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 385 ]), &(acadoWorkspace.Dy[ 385 ]), &(acadoWorkspace.g[ 77 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 390 ]), &(acadoWorkspace.Dy[ 390 ]), &(acadoWorkspace.g[ 78 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 395 ]), &(acadoWorkspace.Dy[ 395 ]), &(acadoWorkspace.g[ 79 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 400 ]), &(acadoWorkspace.Dy[ 400 ]), &(acadoWorkspace.g[ 80 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 405 ]), &(acadoWorkspace.Dy[ 405 ]), &(acadoWorkspace.g[ 81 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 410 ]), &(acadoWorkspace.Dy[ 410 ]), &(acadoWorkspace.g[ 82 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 415 ]), &(acadoWorkspace.Dy[ 415 ]), &(acadoWorkspace.g[ 83 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 420 ]), &(acadoWorkspace.Dy[ 420 ]), &(acadoWorkspace.g[ 84 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 425 ]), &(acadoWorkspace.Dy[ 425 ]), &(acadoWorkspace.g[ 85 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 430 ]), &(acadoWorkspace.Dy[ 430 ]), &(acadoWorkspace.g[ 86 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 435 ]), &(acadoWorkspace.Dy[ 435 ]), &(acadoWorkspace.g[ 87 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 440 ]), &(acadoWorkspace.Dy[ 440 ]), &(acadoWorkspace.g[ 88 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 445 ]), &(acadoWorkspace.Dy[ 445 ]), &(acadoWorkspace.g[ 89 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 450 ]), &(acadoWorkspace.Dy[ 450 ]), &(acadoWorkspace.g[ 90 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 455 ]), &(acadoWorkspace.Dy[ 455 ]), &(acadoWorkspace.g[ 91 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 460 ]), &(acadoWorkspace.Dy[ 460 ]), &(acadoWorkspace.g[ 92 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 465 ]), &(acadoWorkspace.Dy[ 465 ]), &(acadoWorkspace.g[ 93 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 470 ]), &(acadoWorkspace.Dy[ 470 ]), &(acadoWorkspace.g[ 94 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 475 ]), &(acadoWorkspace.Dy[ 475 ]), &(acadoWorkspace.g[ 95 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 480 ]), &(acadoWorkspace.Dy[ 480 ]), &(acadoWorkspace.g[ 96 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 485 ]), &(acadoWorkspace.Dy[ 485 ]), &(acadoWorkspace.g[ 97 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 490 ]), &(acadoWorkspace.Dy[ 490 ]), &(acadoWorkspace.g[ 98 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 495 ]), &(acadoWorkspace.Dy[ 495 ]), &(acadoWorkspace.g[ 99 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 500 ]), &(acadoWorkspace.Dy[ 500 ]), &(acadoWorkspace.g[ 100 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 505 ]), &(acadoWorkspace.Dy[ 505 ]), &(acadoWorkspace.g[ 101 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 510 ]), &(acadoWorkspace.Dy[ 510 ]), &(acadoWorkspace.g[ 102 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 515 ]), &(acadoWorkspace.Dy[ 515 ]), &(acadoWorkspace.g[ 103 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 520 ]), &(acadoWorkspace.Dy[ 520 ]), &(acadoWorkspace.g[ 104 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 525 ]), &(acadoWorkspace.Dy[ 525 ]), &(acadoWorkspace.g[ 105 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 530 ]), &(acadoWorkspace.Dy[ 530 ]), &(acadoWorkspace.g[ 106 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 535 ]), &(acadoWorkspace.Dy[ 535 ]), &(acadoWorkspace.g[ 107 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 540 ]), &(acadoWorkspace.Dy[ 540 ]), &(acadoWorkspace.g[ 108 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 545 ]), &(acadoWorkspace.Dy[ 545 ]), &(acadoWorkspace.g[ 109 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 550 ]), &(acadoWorkspace.Dy[ 550 ]), &(acadoWorkspace.g[ 110 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 555 ]), &(acadoWorkspace.Dy[ 555 ]), &(acadoWorkspace.g[ 111 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 560 ]), &(acadoWorkspace.Dy[ 560 ]), &(acadoWorkspace.g[ 112 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 565 ]), &(acadoWorkspace.Dy[ 565 ]), &(acadoWorkspace.g[ 113 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 570 ]), &(acadoWorkspace.Dy[ 570 ]), &(acadoWorkspace.g[ 114 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 575 ]), &(acadoWorkspace.Dy[ 575 ]), &(acadoWorkspace.g[ 115 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 580 ]), &(acadoWorkspace.Dy[ 580 ]), &(acadoWorkspace.g[ 116 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 585 ]), &(acadoWorkspace.Dy[ 585 ]), &(acadoWorkspace.g[ 117 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 590 ]), &(acadoWorkspace.Dy[ 590 ]), &(acadoWorkspace.g[ 118 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 595 ]), &(acadoWorkspace.Dy[ 595 ]), &(acadoWorkspace.g[ 119 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 600 ]), &(acadoWorkspace.Dy[ 600 ]), &(acadoWorkspace.g[ 120 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 605 ]), &(acadoWorkspace.Dy[ 605 ]), &(acadoWorkspace.g[ 121 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 610 ]), &(acadoWorkspace.Dy[ 610 ]), &(acadoWorkspace.g[ 122 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 615 ]), &(acadoWorkspace.Dy[ 615 ]), &(acadoWorkspace.g[ 123 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 620 ]), &(acadoWorkspace.Dy[ 620 ]), &(acadoWorkspace.g[ 124 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 625 ]), &(acadoWorkspace.Dy[ 625 ]), &(acadoWorkspace.g[ 125 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 630 ]), &(acadoWorkspace.Dy[ 630 ]), &(acadoWorkspace.g[ 126 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 635 ]), &(acadoWorkspace.Dy[ 635 ]), &(acadoWorkspace.g[ 127 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 640 ]), &(acadoWorkspace.Dy[ 640 ]), &(acadoWorkspace.g[ 128 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 645 ]), &(acadoWorkspace.Dy[ 645 ]), &(acadoWorkspace.g[ 129 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 650 ]), &(acadoWorkspace.Dy[ 650 ]), &(acadoWorkspace.g[ 130 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 655 ]), &(acadoWorkspace.Dy[ 655 ]), &(acadoWorkspace.g[ 131 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 660 ]), &(acadoWorkspace.Dy[ 660 ]), &(acadoWorkspace.g[ 132 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 665 ]), &(acadoWorkspace.Dy[ 665 ]), &(acadoWorkspace.g[ 133 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 670 ]), &(acadoWorkspace.Dy[ 670 ]), &(acadoWorkspace.g[ 134 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 675 ]), &(acadoWorkspace.Dy[ 675 ]), &(acadoWorkspace.g[ 135 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 680 ]), &(acadoWorkspace.Dy[ 680 ]), &(acadoWorkspace.g[ 136 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 685 ]), &(acadoWorkspace.Dy[ 685 ]), &(acadoWorkspace.g[ 137 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 690 ]), &(acadoWorkspace.Dy[ 690 ]), &(acadoWorkspace.g[ 138 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 695 ]), &(acadoWorkspace.Dy[ 695 ]), &(acadoWorkspace.g[ 139 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 700 ]), &(acadoWorkspace.Dy[ 700 ]), &(acadoWorkspace.g[ 140 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 705 ]), &(acadoWorkspace.Dy[ 705 ]), &(acadoWorkspace.g[ 141 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 710 ]), &(acadoWorkspace.Dy[ 710 ]), &(acadoWorkspace.g[ 142 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 715 ]), &(acadoWorkspace.Dy[ 715 ]), &(acadoWorkspace.g[ 143 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 720 ]), &(acadoWorkspace.Dy[ 720 ]), &(acadoWorkspace.g[ 144 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 725 ]), &(acadoWorkspace.Dy[ 725 ]), &(acadoWorkspace.g[ 145 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 730 ]), &(acadoWorkspace.Dy[ 730 ]), &(acadoWorkspace.g[ 146 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 735 ]), &(acadoWorkspace.Dy[ 735 ]), &(acadoWorkspace.g[ 147 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 740 ]), &(acadoWorkspace.Dy[ 740 ]), &(acadoWorkspace.g[ 148 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 745 ]), &(acadoWorkspace.Dy[ 745 ]), &(acadoWorkspace.g[ 149 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 750 ]), &(acadoWorkspace.Dy[ 750 ]), &(acadoWorkspace.g[ 150 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 755 ]), &(acadoWorkspace.Dy[ 755 ]), &(acadoWorkspace.g[ 151 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 760 ]), &(acadoWorkspace.Dy[ 760 ]), &(acadoWorkspace.g[ 152 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 765 ]), &(acadoWorkspace.Dy[ 765 ]), &(acadoWorkspace.g[ 153 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 770 ]), &(acadoWorkspace.Dy[ 770 ]), &(acadoWorkspace.g[ 154 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 775 ]), &(acadoWorkspace.Dy[ 775 ]), &(acadoWorkspace.g[ 155 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 780 ]), &(acadoWorkspace.Dy[ 780 ]), &(acadoWorkspace.g[ 156 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 785 ]), &(acadoWorkspace.Dy[ 785 ]), &(acadoWorkspace.g[ 157 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 790 ]), &(acadoWorkspace.Dy[ 790 ]), &(acadoWorkspace.g[ 158 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 795 ]), &(acadoWorkspace.Dy[ 795 ]), &(acadoWorkspace.g[ 159 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 800 ]), &(acadoWorkspace.Dy[ 800 ]), &(acadoWorkspace.g[ 160 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 805 ]), &(acadoWorkspace.Dy[ 805 ]), &(acadoWorkspace.g[ 161 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 810 ]), &(acadoWorkspace.Dy[ 810 ]), &(acadoWorkspace.g[ 162 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 815 ]), &(acadoWorkspace.Dy[ 815 ]), &(acadoWorkspace.g[ 163 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 820 ]), &(acadoWorkspace.Dy[ 820 ]), &(acadoWorkspace.g[ 164 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 825 ]), &(acadoWorkspace.Dy[ 825 ]), &(acadoWorkspace.g[ 165 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 830 ]), &(acadoWorkspace.Dy[ 830 ]), &(acadoWorkspace.g[ 166 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 835 ]), &(acadoWorkspace.Dy[ 835 ]), &(acadoWorkspace.g[ 167 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 840 ]), &(acadoWorkspace.Dy[ 840 ]), &(acadoWorkspace.g[ 168 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 845 ]), &(acadoWorkspace.Dy[ 845 ]), &(acadoWorkspace.g[ 169 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 850 ]), &(acadoWorkspace.Dy[ 850 ]), &(acadoWorkspace.g[ 170 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 855 ]), &(acadoWorkspace.Dy[ 855 ]), &(acadoWorkspace.g[ 171 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 860 ]), &(acadoWorkspace.Dy[ 860 ]), &(acadoWorkspace.g[ 172 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 865 ]), &(acadoWorkspace.Dy[ 865 ]), &(acadoWorkspace.g[ 173 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 870 ]), &(acadoWorkspace.Dy[ 870 ]), &(acadoWorkspace.g[ 174 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 875 ]), &(acadoWorkspace.Dy[ 875 ]), &(acadoWorkspace.g[ 175 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 880 ]), &(acadoWorkspace.Dy[ 880 ]), &(acadoWorkspace.g[ 176 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 885 ]), &(acadoWorkspace.Dy[ 885 ]), &(acadoWorkspace.g[ 177 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 890 ]), &(acadoWorkspace.Dy[ 890 ]), &(acadoWorkspace.g[ 178 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 895 ]), &(acadoWorkspace.Dy[ 895 ]), &(acadoWorkspace.g[ 179 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 900 ]), &(acadoWorkspace.Dy[ 900 ]), &(acadoWorkspace.g[ 180 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 905 ]), &(acadoWorkspace.Dy[ 905 ]), &(acadoWorkspace.g[ 181 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 910 ]), &(acadoWorkspace.Dy[ 910 ]), &(acadoWorkspace.g[ 182 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 915 ]), &(acadoWorkspace.Dy[ 915 ]), &(acadoWorkspace.g[ 183 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 920 ]), &(acadoWorkspace.Dy[ 920 ]), &(acadoWorkspace.g[ 184 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 925 ]), &(acadoWorkspace.Dy[ 925 ]), &(acadoWorkspace.g[ 185 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 930 ]), &(acadoWorkspace.Dy[ 930 ]), &(acadoWorkspace.g[ 186 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 935 ]), &(acadoWorkspace.Dy[ 935 ]), &(acadoWorkspace.g[ 187 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 940 ]), &(acadoWorkspace.Dy[ 940 ]), &(acadoWorkspace.g[ 188 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 945 ]), &(acadoWorkspace.Dy[ 945 ]), &(acadoWorkspace.g[ 189 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 950 ]), &(acadoWorkspace.Dy[ 950 ]), &(acadoWorkspace.g[ 190 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 955 ]), &(acadoWorkspace.Dy[ 955 ]), &(acadoWorkspace.g[ 191 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 960 ]), &(acadoWorkspace.Dy[ 960 ]), &(acadoWorkspace.g[ 192 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 965 ]), &(acadoWorkspace.Dy[ 965 ]), &(acadoWorkspace.g[ 193 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 970 ]), &(acadoWorkspace.Dy[ 970 ]), &(acadoWorkspace.g[ 194 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 975 ]), &(acadoWorkspace.Dy[ 975 ]), &(acadoWorkspace.g[ 195 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 980 ]), &(acadoWorkspace.Dy[ 980 ]), &(acadoWorkspace.g[ 196 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 985 ]), &(acadoWorkspace.Dy[ 985 ]), &(acadoWorkspace.g[ 197 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 990 ]), &(acadoWorkspace.Dy[ 990 ]), &(acadoWorkspace.g[ 198 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 995 ]), &(acadoWorkspace.Dy[ 995 ]), &(acadoWorkspace.g[ 199 ]) );

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
acado_multQDy( &(acadoWorkspace.Q2[ 400 ]), &(acadoWorkspace.Dy[ 100 ]), &(acadoWorkspace.QDy[ 80 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 420 ]), &(acadoWorkspace.Dy[ 105 ]), &(acadoWorkspace.QDy[ 84 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 440 ]), &(acadoWorkspace.Dy[ 110 ]), &(acadoWorkspace.QDy[ 88 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 460 ]), &(acadoWorkspace.Dy[ 115 ]), &(acadoWorkspace.QDy[ 92 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 480 ]), &(acadoWorkspace.Dy[ 120 ]), &(acadoWorkspace.QDy[ 96 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 500 ]), &(acadoWorkspace.Dy[ 125 ]), &(acadoWorkspace.QDy[ 100 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 520 ]), &(acadoWorkspace.Dy[ 130 ]), &(acadoWorkspace.QDy[ 104 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 540 ]), &(acadoWorkspace.Dy[ 135 ]), &(acadoWorkspace.QDy[ 108 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 560 ]), &(acadoWorkspace.Dy[ 140 ]), &(acadoWorkspace.QDy[ 112 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 580 ]), &(acadoWorkspace.Dy[ 145 ]), &(acadoWorkspace.QDy[ 116 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 600 ]), &(acadoWorkspace.Dy[ 150 ]), &(acadoWorkspace.QDy[ 120 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 620 ]), &(acadoWorkspace.Dy[ 155 ]), &(acadoWorkspace.QDy[ 124 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 640 ]), &(acadoWorkspace.Dy[ 160 ]), &(acadoWorkspace.QDy[ 128 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 660 ]), &(acadoWorkspace.Dy[ 165 ]), &(acadoWorkspace.QDy[ 132 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 680 ]), &(acadoWorkspace.Dy[ 170 ]), &(acadoWorkspace.QDy[ 136 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 700 ]), &(acadoWorkspace.Dy[ 175 ]), &(acadoWorkspace.QDy[ 140 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 720 ]), &(acadoWorkspace.Dy[ 180 ]), &(acadoWorkspace.QDy[ 144 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 740 ]), &(acadoWorkspace.Dy[ 185 ]), &(acadoWorkspace.QDy[ 148 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 760 ]), &(acadoWorkspace.Dy[ 190 ]), &(acadoWorkspace.QDy[ 152 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 780 ]), &(acadoWorkspace.Dy[ 195 ]), &(acadoWorkspace.QDy[ 156 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 800 ]), &(acadoWorkspace.Dy[ 200 ]), &(acadoWorkspace.QDy[ 160 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 820 ]), &(acadoWorkspace.Dy[ 205 ]), &(acadoWorkspace.QDy[ 164 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 840 ]), &(acadoWorkspace.Dy[ 210 ]), &(acadoWorkspace.QDy[ 168 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 860 ]), &(acadoWorkspace.Dy[ 215 ]), &(acadoWorkspace.QDy[ 172 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 880 ]), &(acadoWorkspace.Dy[ 220 ]), &(acadoWorkspace.QDy[ 176 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 900 ]), &(acadoWorkspace.Dy[ 225 ]), &(acadoWorkspace.QDy[ 180 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 920 ]), &(acadoWorkspace.Dy[ 230 ]), &(acadoWorkspace.QDy[ 184 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 940 ]), &(acadoWorkspace.Dy[ 235 ]), &(acadoWorkspace.QDy[ 188 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 960 ]), &(acadoWorkspace.Dy[ 240 ]), &(acadoWorkspace.QDy[ 192 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 980 ]), &(acadoWorkspace.Dy[ 245 ]), &(acadoWorkspace.QDy[ 196 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1000 ]), &(acadoWorkspace.Dy[ 250 ]), &(acadoWorkspace.QDy[ 200 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1020 ]), &(acadoWorkspace.Dy[ 255 ]), &(acadoWorkspace.QDy[ 204 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1040 ]), &(acadoWorkspace.Dy[ 260 ]), &(acadoWorkspace.QDy[ 208 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1060 ]), &(acadoWorkspace.Dy[ 265 ]), &(acadoWorkspace.QDy[ 212 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1080 ]), &(acadoWorkspace.Dy[ 270 ]), &(acadoWorkspace.QDy[ 216 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1100 ]), &(acadoWorkspace.Dy[ 275 ]), &(acadoWorkspace.QDy[ 220 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1120 ]), &(acadoWorkspace.Dy[ 280 ]), &(acadoWorkspace.QDy[ 224 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1140 ]), &(acadoWorkspace.Dy[ 285 ]), &(acadoWorkspace.QDy[ 228 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1160 ]), &(acadoWorkspace.Dy[ 290 ]), &(acadoWorkspace.QDy[ 232 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1180 ]), &(acadoWorkspace.Dy[ 295 ]), &(acadoWorkspace.QDy[ 236 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1200 ]), &(acadoWorkspace.Dy[ 300 ]), &(acadoWorkspace.QDy[ 240 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1220 ]), &(acadoWorkspace.Dy[ 305 ]), &(acadoWorkspace.QDy[ 244 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1240 ]), &(acadoWorkspace.Dy[ 310 ]), &(acadoWorkspace.QDy[ 248 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1260 ]), &(acadoWorkspace.Dy[ 315 ]), &(acadoWorkspace.QDy[ 252 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1280 ]), &(acadoWorkspace.Dy[ 320 ]), &(acadoWorkspace.QDy[ 256 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1300 ]), &(acadoWorkspace.Dy[ 325 ]), &(acadoWorkspace.QDy[ 260 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1320 ]), &(acadoWorkspace.Dy[ 330 ]), &(acadoWorkspace.QDy[ 264 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1340 ]), &(acadoWorkspace.Dy[ 335 ]), &(acadoWorkspace.QDy[ 268 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1360 ]), &(acadoWorkspace.Dy[ 340 ]), &(acadoWorkspace.QDy[ 272 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1380 ]), &(acadoWorkspace.Dy[ 345 ]), &(acadoWorkspace.QDy[ 276 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1400 ]), &(acadoWorkspace.Dy[ 350 ]), &(acadoWorkspace.QDy[ 280 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1420 ]), &(acadoWorkspace.Dy[ 355 ]), &(acadoWorkspace.QDy[ 284 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1440 ]), &(acadoWorkspace.Dy[ 360 ]), &(acadoWorkspace.QDy[ 288 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1460 ]), &(acadoWorkspace.Dy[ 365 ]), &(acadoWorkspace.QDy[ 292 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1480 ]), &(acadoWorkspace.Dy[ 370 ]), &(acadoWorkspace.QDy[ 296 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1500 ]), &(acadoWorkspace.Dy[ 375 ]), &(acadoWorkspace.QDy[ 300 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1520 ]), &(acadoWorkspace.Dy[ 380 ]), &(acadoWorkspace.QDy[ 304 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1540 ]), &(acadoWorkspace.Dy[ 385 ]), &(acadoWorkspace.QDy[ 308 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1560 ]), &(acadoWorkspace.Dy[ 390 ]), &(acadoWorkspace.QDy[ 312 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1580 ]), &(acadoWorkspace.Dy[ 395 ]), &(acadoWorkspace.QDy[ 316 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1600 ]), &(acadoWorkspace.Dy[ 400 ]), &(acadoWorkspace.QDy[ 320 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1620 ]), &(acadoWorkspace.Dy[ 405 ]), &(acadoWorkspace.QDy[ 324 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1640 ]), &(acadoWorkspace.Dy[ 410 ]), &(acadoWorkspace.QDy[ 328 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1660 ]), &(acadoWorkspace.Dy[ 415 ]), &(acadoWorkspace.QDy[ 332 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1680 ]), &(acadoWorkspace.Dy[ 420 ]), &(acadoWorkspace.QDy[ 336 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1700 ]), &(acadoWorkspace.Dy[ 425 ]), &(acadoWorkspace.QDy[ 340 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1720 ]), &(acadoWorkspace.Dy[ 430 ]), &(acadoWorkspace.QDy[ 344 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1740 ]), &(acadoWorkspace.Dy[ 435 ]), &(acadoWorkspace.QDy[ 348 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1760 ]), &(acadoWorkspace.Dy[ 440 ]), &(acadoWorkspace.QDy[ 352 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1780 ]), &(acadoWorkspace.Dy[ 445 ]), &(acadoWorkspace.QDy[ 356 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1800 ]), &(acadoWorkspace.Dy[ 450 ]), &(acadoWorkspace.QDy[ 360 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1820 ]), &(acadoWorkspace.Dy[ 455 ]), &(acadoWorkspace.QDy[ 364 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1840 ]), &(acadoWorkspace.Dy[ 460 ]), &(acadoWorkspace.QDy[ 368 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1860 ]), &(acadoWorkspace.Dy[ 465 ]), &(acadoWorkspace.QDy[ 372 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1880 ]), &(acadoWorkspace.Dy[ 470 ]), &(acadoWorkspace.QDy[ 376 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1900 ]), &(acadoWorkspace.Dy[ 475 ]), &(acadoWorkspace.QDy[ 380 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1920 ]), &(acadoWorkspace.Dy[ 480 ]), &(acadoWorkspace.QDy[ 384 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1940 ]), &(acadoWorkspace.Dy[ 485 ]), &(acadoWorkspace.QDy[ 388 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1960 ]), &(acadoWorkspace.Dy[ 490 ]), &(acadoWorkspace.QDy[ 392 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1980 ]), &(acadoWorkspace.Dy[ 495 ]), &(acadoWorkspace.QDy[ 396 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2000 ]), &(acadoWorkspace.Dy[ 500 ]), &(acadoWorkspace.QDy[ 400 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2020 ]), &(acadoWorkspace.Dy[ 505 ]), &(acadoWorkspace.QDy[ 404 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2040 ]), &(acadoWorkspace.Dy[ 510 ]), &(acadoWorkspace.QDy[ 408 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2060 ]), &(acadoWorkspace.Dy[ 515 ]), &(acadoWorkspace.QDy[ 412 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2080 ]), &(acadoWorkspace.Dy[ 520 ]), &(acadoWorkspace.QDy[ 416 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2100 ]), &(acadoWorkspace.Dy[ 525 ]), &(acadoWorkspace.QDy[ 420 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2120 ]), &(acadoWorkspace.Dy[ 530 ]), &(acadoWorkspace.QDy[ 424 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2140 ]), &(acadoWorkspace.Dy[ 535 ]), &(acadoWorkspace.QDy[ 428 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2160 ]), &(acadoWorkspace.Dy[ 540 ]), &(acadoWorkspace.QDy[ 432 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2180 ]), &(acadoWorkspace.Dy[ 545 ]), &(acadoWorkspace.QDy[ 436 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2200 ]), &(acadoWorkspace.Dy[ 550 ]), &(acadoWorkspace.QDy[ 440 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2220 ]), &(acadoWorkspace.Dy[ 555 ]), &(acadoWorkspace.QDy[ 444 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2240 ]), &(acadoWorkspace.Dy[ 560 ]), &(acadoWorkspace.QDy[ 448 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2260 ]), &(acadoWorkspace.Dy[ 565 ]), &(acadoWorkspace.QDy[ 452 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2280 ]), &(acadoWorkspace.Dy[ 570 ]), &(acadoWorkspace.QDy[ 456 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2300 ]), &(acadoWorkspace.Dy[ 575 ]), &(acadoWorkspace.QDy[ 460 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2320 ]), &(acadoWorkspace.Dy[ 580 ]), &(acadoWorkspace.QDy[ 464 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2340 ]), &(acadoWorkspace.Dy[ 585 ]), &(acadoWorkspace.QDy[ 468 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2360 ]), &(acadoWorkspace.Dy[ 590 ]), &(acadoWorkspace.QDy[ 472 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2380 ]), &(acadoWorkspace.Dy[ 595 ]), &(acadoWorkspace.QDy[ 476 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2400 ]), &(acadoWorkspace.Dy[ 600 ]), &(acadoWorkspace.QDy[ 480 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2420 ]), &(acadoWorkspace.Dy[ 605 ]), &(acadoWorkspace.QDy[ 484 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2440 ]), &(acadoWorkspace.Dy[ 610 ]), &(acadoWorkspace.QDy[ 488 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2460 ]), &(acadoWorkspace.Dy[ 615 ]), &(acadoWorkspace.QDy[ 492 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2480 ]), &(acadoWorkspace.Dy[ 620 ]), &(acadoWorkspace.QDy[ 496 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2500 ]), &(acadoWorkspace.Dy[ 625 ]), &(acadoWorkspace.QDy[ 500 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2520 ]), &(acadoWorkspace.Dy[ 630 ]), &(acadoWorkspace.QDy[ 504 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2540 ]), &(acadoWorkspace.Dy[ 635 ]), &(acadoWorkspace.QDy[ 508 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2560 ]), &(acadoWorkspace.Dy[ 640 ]), &(acadoWorkspace.QDy[ 512 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2580 ]), &(acadoWorkspace.Dy[ 645 ]), &(acadoWorkspace.QDy[ 516 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2600 ]), &(acadoWorkspace.Dy[ 650 ]), &(acadoWorkspace.QDy[ 520 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2620 ]), &(acadoWorkspace.Dy[ 655 ]), &(acadoWorkspace.QDy[ 524 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2640 ]), &(acadoWorkspace.Dy[ 660 ]), &(acadoWorkspace.QDy[ 528 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2660 ]), &(acadoWorkspace.Dy[ 665 ]), &(acadoWorkspace.QDy[ 532 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2680 ]), &(acadoWorkspace.Dy[ 670 ]), &(acadoWorkspace.QDy[ 536 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2700 ]), &(acadoWorkspace.Dy[ 675 ]), &(acadoWorkspace.QDy[ 540 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2720 ]), &(acadoWorkspace.Dy[ 680 ]), &(acadoWorkspace.QDy[ 544 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2740 ]), &(acadoWorkspace.Dy[ 685 ]), &(acadoWorkspace.QDy[ 548 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2760 ]), &(acadoWorkspace.Dy[ 690 ]), &(acadoWorkspace.QDy[ 552 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2780 ]), &(acadoWorkspace.Dy[ 695 ]), &(acadoWorkspace.QDy[ 556 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2800 ]), &(acadoWorkspace.Dy[ 700 ]), &(acadoWorkspace.QDy[ 560 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2820 ]), &(acadoWorkspace.Dy[ 705 ]), &(acadoWorkspace.QDy[ 564 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2840 ]), &(acadoWorkspace.Dy[ 710 ]), &(acadoWorkspace.QDy[ 568 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2860 ]), &(acadoWorkspace.Dy[ 715 ]), &(acadoWorkspace.QDy[ 572 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2880 ]), &(acadoWorkspace.Dy[ 720 ]), &(acadoWorkspace.QDy[ 576 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2900 ]), &(acadoWorkspace.Dy[ 725 ]), &(acadoWorkspace.QDy[ 580 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2920 ]), &(acadoWorkspace.Dy[ 730 ]), &(acadoWorkspace.QDy[ 584 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2940 ]), &(acadoWorkspace.Dy[ 735 ]), &(acadoWorkspace.QDy[ 588 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2960 ]), &(acadoWorkspace.Dy[ 740 ]), &(acadoWorkspace.QDy[ 592 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 2980 ]), &(acadoWorkspace.Dy[ 745 ]), &(acadoWorkspace.QDy[ 596 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3000 ]), &(acadoWorkspace.Dy[ 750 ]), &(acadoWorkspace.QDy[ 600 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3020 ]), &(acadoWorkspace.Dy[ 755 ]), &(acadoWorkspace.QDy[ 604 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3040 ]), &(acadoWorkspace.Dy[ 760 ]), &(acadoWorkspace.QDy[ 608 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3060 ]), &(acadoWorkspace.Dy[ 765 ]), &(acadoWorkspace.QDy[ 612 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3080 ]), &(acadoWorkspace.Dy[ 770 ]), &(acadoWorkspace.QDy[ 616 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3100 ]), &(acadoWorkspace.Dy[ 775 ]), &(acadoWorkspace.QDy[ 620 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3120 ]), &(acadoWorkspace.Dy[ 780 ]), &(acadoWorkspace.QDy[ 624 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3140 ]), &(acadoWorkspace.Dy[ 785 ]), &(acadoWorkspace.QDy[ 628 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3160 ]), &(acadoWorkspace.Dy[ 790 ]), &(acadoWorkspace.QDy[ 632 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3180 ]), &(acadoWorkspace.Dy[ 795 ]), &(acadoWorkspace.QDy[ 636 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3200 ]), &(acadoWorkspace.Dy[ 800 ]), &(acadoWorkspace.QDy[ 640 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3220 ]), &(acadoWorkspace.Dy[ 805 ]), &(acadoWorkspace.QDy[ 644 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3240 ]), &(acadoWorkspace.Dy[ 810 ]), &(acadoWorkspace.QDy[ 648 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3260 ]), &(acadoWorkspace.Dy[ 815 ]), &(acadoWorkspace.QDy[ 652 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3280 ]), &(acadoWorkspace.Dy[ 820 ]), &(acadoWorkspace.QDy[ 656 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3300 ]), &(acadoWorkspace.Dy[ 825 ]), &(acadoWorkspace.QDy[ 660 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3320 ]), &(acadoWorkspace.Dy[ 830 ]), &(acadoWorkspace.QDy[ 664 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3340 ]), &(acadoWorkspace.Dy[ 835 ]), &(acadoWorkspace.QDy[ 668 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3360 ]), &(acadoWorkspace.Dy[ 840 ]), &(acadoWorkspace.QDy[ 672 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3380 ]), &(acadoWorkspace.Dy[ 845 ]), &(acadoWorkspace.QDy[ 676 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3400 ]), &(acadoWorkspace.Dy[ 850 ]), &(acadoWorkspace.QDy[ 680 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3420 ]), &(acadoWorkspace.Dy[ 855 ]), &(acadoWorkspace.QDy[ 684 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3440 ]), &(acadoWorkspace.Dy[ 860 ]), &(acadoWorkspace.QDy[ 688 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3460 ]), &(acadoWorkspace.Dy[ 865 ]), &(acadoWorkspace.QDy[ 692 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3480 ]), &(acadoWorkspace.Dy[ 870 ]), &(acadoWorkspace.QDy[ 696 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3500 ]), &(acadoWorkspace.Dy[ 875 ]), &(acadoWorkspace.QDy[ 700 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3520 ]), &(acadoWorkspace.Dy[ 880 ]), &(acadoWorkspace.QDy[ 704 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3540 ]), &(acadoWorkspace.Dy[ 885 ]), &(acadoWorkspace.QDy[ 708 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3560 ]), &(acadoWorkspace.Dy[ 890 ]), &(acadoWorkspace.QDy[ 712 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3580 ]), &(acadoWorkspace.Dy[ 895 ]), &(acadoWorkspace.QDy[ 716 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3600 ]), &(acadoWorkspace.Dy[ 900 ]), &(acadoWorkspace.QDy[ 720 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3620 ]), &(acadoWorkspace.Dy[ 905 ]), &(acadoWorkspace.QDy[ 724 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3640 ]), &(acadoWorkspace.Dy[ 910 ]), &(acadoWorkspace.QDy[ 728 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3660 ]), &(acadoWorkspace.Dy[ 915 ]), &(acadoWorkspace.QDy[ 732 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3680 ]), &(acadoWorkspace.Dy[ 920 ]), &(acadoWorkspace.QDy[ 736 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3700 ]), &(acadoWorkspace.Dy[ 925 ]), &(acadoWorkspace.QDy[ 740 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3720 ]), &(acadoWorkspace.Dy[ 930 ]), &(acadoWorkspace.QDy[ 744 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3740 ]), &(acadoWorkspace.Dy[ 935 ]), &(acadoWorkspace.QDy[ 748 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3760 ]), &(acadoWorkspace.Dy[ 940 ]), &(acadoWorkspace.QDy[ 752 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3780 ]), &(acadoWorkspace.Dy[ 945 ]), &(acadoWorkspace.QDy[ 756 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3800 ]), &(acadoWorkspace.Dy[ 950 ]), &(acadoWorkspace.QDy[ 760 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3820 ]), &(acadoWorkspace.Dy[ 955 ]), &(acadoWorkspace.QDy[ 764 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3840 ]), &(acadoWorkspace.Dy[ 960 ]), &(acadoWorkspace.QDy[ 768 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3860 ]), &(acadoWorkspace.Dy[ 965 ]), &(acadoWorkspace.QDy[ 772 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3880 ]), &(acadoWorkspace.Dy[ 970 ]), &(acadoWorkspace.QDy[ 776 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3900 ]), &(acadoWorkspace.Dy[ 975 ]), &(acadoWorkspace.QDy[ 780 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3920 ]), &(acadoWorkspace.Dy[ 980 ]), &(acadoWorkspace.QDy[ 784 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3940 ]), &(acadoWorkspace.Dy[ 985 ]), &(acadoWorkspace.QDy[ 788 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3960 ]), &(acadoWorkspace.Dy[ 990 ]), &(acadoWorkspace.QDy[ 792 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 3980 ]), &(acadoWorkspace.Dy[ 995 ]), &(acadoWorkspace.QDy[ 796 ]) );

acadoWorkspace.QDy[800] = + acadoWorkspace.QN2[0]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[1]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[2]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[3]*acadoWorkspace.DyN[3];
acadoWorkspace.QDy[801] = + acadoWorkspace.QN2[4]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[5]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[6]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[7]*acadoWorkspace.DyN[3];
acadoWorkspace.QDy[802] = + acadoWorkspace.QN2[8]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[9]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[10]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[11]*acadoWorkspace.DyN[3];
acadoWorkspace.QDy[803] = + acadoWorkspace.QN2[12]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[13]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[14]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[15]*acadoWorkspace.DyN[3];

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
acado_macASbar( &(acadoWorkspace.evGx[ 320 ]), &(acadoWorkspace.sbar[ 80 ]), &(acadoWorkspace.sbar[ 84 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 336 ]), &(acadoWorkspace.sbar[ 84 ]), &(acadoWorkspace.sbar[ 88 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 352 ]), &(acadoWorkspace.sbar[ 88 ]), &(acadoWorkspace.sbar[ 92 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 368 ]), &(acadoWorkspace.sbar[ 92 ]), &(acadoWorkspace.sbar[ 96 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 384 ]), &(acadoWorkspace.sbar[ 96 ]), &(acadoWorkspace.sbar[ 100 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 400 ]), &(acadoWorkspace.sbar[ 100 ]), &(acadoWorkspace.sbar[ 104 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 416 ]), &(acadoWorkspace.sbar[ 104 ]), &(acadoWorkspace.sbar[ 108 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 432 ]), &(acadoWorkspace.sbar[ 108 ]), &(acadoWorkspace.sbar[ 112 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 448 ]), &(acadoWorkspace.sbar[ 112 ]), &(acadoWorkspace.sbar[ 116 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 464 ]), &(acadoWorkspace.sbar[ 116 ]), &(acadoWorkspace.sbar[ 120 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 480 ]), &(acadoWorkspace.sbar[ 120 ]), &(acadoWorkspace.sbar[ 124 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 496 ]), &(acadoWorkspace.sbar[ 124 ]), &(acadoWorkspace.sbar[ 128 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 512 ]), &(acadoWorkspace.sbar[ 128 ]), &(acadoWorkspace.sbar[ 132 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 528 ]), &(acadoWorkspace.sbar[ 132 ]), &(acadoWorkspace.sbar[ 136 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 544 ]), &(acadoWorkspace.sbar[ 136 ]), &(acadoWorkspace.sbar[ 140 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 560 ]), &(acadoWorkspace.sbar[ 140 ]), &(acadoWorkspace.sbar[ 144 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.sbar[ 144 ]), &(acadoWorkspace.sbar[ 148 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 592 ]), &(acadoWorkspace.sbar[ 148 ]), &(acadoWorkspace.sbar[ 152 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 608 ]), &(acadoWorkspace.sbar[ 152 ]), &(acadoWorkspace.sbar[ 156 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 624 ]), &(acadoWorkspace.sbar[ 156 ]), &(acadoWorkspace.sbar[ 160 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 640 ]), &(acadoWorkspace.sbar[ 160 ]), &(acadoWorkspace.sbar[ 164 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 656 ]), &(acadoWorkspace.sbar[ 164 ]), &(acadoWorkspace.sbar[ 168 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 672 ]), &(acadoWorkspace.sbar[ 168 ]), &(acadoWorkspace.sbar[ 172 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 688 ]), &(acadoWorkspace.sbar[ 172 ]), &(acadoWorkspace.sbar[ 176 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 704 ]), &(acadoWorkspace.sbar[ 176 ]), &(acadoWorkspace.sbar[ 180 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 720 ]), &(acadoWorkspace.sbar[ 180 ]), &(acadoWorkspace.sbar[ 184 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 736 ]), &(acadoWorkspace.sbar[ 184 ]), &(acadoWorkspace.sbar[ 188 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 752 ]), &(acadoWorkspace.sbar[ 188 ]), &(acadoWorkspace.sbar[ 192 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 768 ]), &(acadoWorkspace.sbar[ 192 ]), &(acadoWorkspace.sbar[ 196 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 784 ]), &(acadoWorkspace.sbar[ 196 ]), &(acadoWorkspace.sbar[ 200 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 800 ]), &(acadoWorkspace.sbar[ 200 ]), &(acadoWorkspace.sbar[ 204 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 816 ]), &(acadoWorkspace.sbar[ 204 ]), &(acadoWorkspace.sbar[ 208 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 832 ]), &(acadoWorkspace.sbar[ 208 ]), &(acadoWorkspace.sbar[ 212 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 848 ]), &(acadoWorkspace.sbar[ 212 ]), &(acadoWorkspace.sbar[ 216 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 864 ]), &(acadoWorkspace.sbar[ 216 ]), &(acadoWorkspace.sbar[ 220 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 880 ]), &(acadoWorkspace.sbar[ 220 ]), &(acadoWorkspace.sbar[ 224 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 896 ]), &(acadoWorkspace.sbar[ 224 ]), &(acadoWorkspace.sbar[ 228 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 912 ]), &(acadoWorkspace.sbar[ 228 ]), &(acadoWorkspace.sbar[ 232 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 928 ]), &(acadoWorkspace.sbar[ 232 ]), &(acadoWorkspace.sbar[ 236 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 944 ]), &(acadoWorkspace.sbar[ 236 ]), &(acadoWorkspace.sbar[ 240 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 960 ]), &(acadoWorkspace.sbar[ 240 ]), &(acadoWorkspace.sbar[ 244 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 976 ]), &(acadoWorkspace.sbar[ 244 ]), &(acadoWorkspace.sbar[ 248 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 992 ]), &(acadoWorkspace.sbar[ 248 ]), &(acadoWorkspace.sbar[ 252 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1008 ]), &(acadoWorkspace.sbar[ 252 ]), &(acadoWorkspace.sbar[ 256 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1024 ]), &(acadoWorkspace.sbar[ 256 ]), &(acadoWorkspace.sbar[ 260 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1040 ]), &(acadoWorkspace.sbar[ 260 ]), &(acadoWorkspace.sbar[ 264 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1056 ]), &(acadoWorkspace.sbar[ 264 ]), &(acadoWorkspace.sbar[ 268 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1072 ]), &(acadoWorkspace.sbar[ 268 ]), &(acadoWorkspace.sbar[ 272 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1088 ]), &(acadoWorkspace.sbar[ 272 ]), &(acadoWorkspace.sbar[ 276 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1104 ]), &(acadoWorkspace.sbar[ 276 ]), &(acadoWorkspace.sbar[ 280 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1120 ]), &(acadoWorkspace.sbar[ 280 ]), &(acadoWorkspace.sbar[ 284 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1136 ]), &(acadoWorkspace.sbar[ 284 ]), &(acadoWorkspace.sbar[ 288 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1152 ]), &(acadoWorkspace.sbar[ 288 ]), &(acadoWorkspace.sbar[ 292 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1168 ]), &(acadoWorkspace.sbar[ 292 ]), &(acadoWorkspace.sbar[ 296 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1184 ]), &(acadoWorkspace.sbar[ 296 ]), &(acadoWorkspace.sbar[ 300 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1200 ]), &(acadoWorkspace.sbar[ 300 ]), &(acadoWorkspace.sbar[ 304 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1216 ]), &(acadoWorkspace.sbar[ 304 ]), &(acadoWorkspace.sbar[ 308 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1232 ]), &(acadoWorkspace.sbar[ 308 ]), &(acadoWorkspace.sbar[ 312 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1248 ]), &(acadoWorkspace.sbar[ 312 ]), &(acadoWorkspace.sbar[ 316 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1264 ]), &(acadoWorkspace.sbar[ 316 ]), &(acadoWorkspace.sbar[ 320 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1280 ]), &(acadoWorkspace.sbar[ 320 ]), &(acadoWorkspace.sbar[ 324 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1296 ]), &(acadoWorkspace.sbar[ 324 ]), &(acadoWorkspace.sbar[ 328 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1312 ]), &(acadoWorkspace.sbar[ 328 ]), &(acadoWorkspace.sbar[ 332 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1328 ]), &(acadoWorkspace.sbar[ 332 ]), &(acadoWorkspace.sbar[ 336 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1344 ]), &(acadoWorkspace.sbar[ 336 ]), &(acadoWorkspace.sbar[ 340 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1360 ]), &(acadoWorkspace.sbar[ 340 ]), &(acadoWorkspace.sbar[ 344 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1376 ]), &(acadoWorkspace.sbar[ 344 ]), &(acadoWorkspace.sbar[ 348 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1392 ]), &(acadoWorkspace.sbar[ 348 ]), &(acadoWorkspace.sbar[ 352 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1408 ]), &(acadoWorkspace.sbar[ 352 ]), &(acadoWorkspace.sbar[ 356 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1424 ]), &(acadoWorkspace.sbar[ 356 ]), &(acadoWorkspace.sbar[ 360 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1440 ]), &(acadoWorkspace.sbar[ 360 ]), &(acadoWorkspace.sbar[ 364 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1456 ]), &(acadoWorkspace.sbar[ 364 ]), &(acadoWorkspace.sbar[ 368 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1472 ]), &(acadoWorkspace.sbar[ 368 ]), &(acadoWorkspace.sbar[ 372 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1488 ]), &(acadoWorkspace.sbar[ 372 ]), &(acadoWorkspace.sbar[ 376 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1504 ]), &(acadoWorkspace.sbar[ 376 ]), &(acadoWorkspace.sbar[ 380 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1520 ]), &(acadoWorkspace.sbar[ 380 ]), &(acadoWorkspace.sbar[ 384 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1536 ]), &(acadoWorkspace.sbar[ 384 ]), &(acadoWorkspace.sbar[ 388 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1552 ]), &(acadoWorkspace.sbar[ 388 ]), &(acadoWorkspace.sbar[ 392 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1568 ]), &(acadoWorkspace.sbar[ 392 ]), &(acadoWorkspace.sbar[ 396 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1584 ]), &(acadoWorkspace.sbar[ 396 ]), &(acadoWorkspace.sbar[ 400 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1600 ]), &(acadoWorkspace.sbar[ 400 ]), &(acadoWorkspace.sbar[ 404 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1616 ]), &(acadoWorkspace.sbar[ 404 ]), &(acadoWorkspace.sbar[ 408 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1632 ]), &(acadoWorkspace.sbar[ 408 ]), &(acadoWorkspace.sbar[ 412 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1648 ]), &(acadoWorkspace.sbar[ 412 ]), &(acadoWorkspace.sbar[ 416 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1664 ]), &(acadoWorkspace.sbar[ 416 ]), &(acadoWorkspace.sbar[ 420 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1680 ]), &(acadoWorkspace.sbar[ 420 ]), &(acadoWorkspace.sbar[ 424 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1696 ]), &(acadoWorkspace.sbar[ 424 ]), &(acadoWorkspace.sbar[ 428 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1712 ]), &(acadoWorkspace.sbar[ 428 ]), &(acadoWorkspace.sbar[ 432 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1728 ]), &(acadoWorkspace.sbar[ 432 ]), &(acadoWorkspace.sbar[ 436 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1744 ]), &(acadoWorkspace.sbar[ 436 ]), &(acadoWorkspace.sbar[ 440 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1760 ]), &(acadoWorkspace.sbar[ 440 ]), &(acadoWorkspace.sbar[ 444 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1776 ]), &(acadoWorkspace.sbar[ 444 ]), &(acadoWorkspace.sbar[ 448 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1792 ]), &(acadoWorkspace.sbar[ 448 ]), &(acadoWorkspace.sbar[ 452 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1808 ]), &(acadoWorkspace.sbar[ 452 ]), &(acadoWorkspace.sbar[ 456 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1824 ]), &(acadoWorkspace.sbar[ 456 ]), &(acadoWorkspace.sbar[ 460 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1840 ]), &(acadoWorkspace.sbar[ 460 ]), &(acadoWorkspace.sbar[ 464 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1856 ]), &(acadoWorkspace.sbar[ 464 ]), &(acadoWorkspace.sbar[ 468 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1872 ]), &(acadoWorkspace.sbar[ 468 ]), &(acadoWorkspace.sbar[ 472 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1888 ]), &(acadoWorkspace.sbar[ 472 ]), &(acadoWorkspace.sbar[ 476 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1904 ]), &(acadoWorkspace.sbar[ 476 ]), &(acadoWorkspace.sbar[ 480 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1920 ]), &(acadoWorkspace.sbar[ 480 ]), &(acadoWorkspace.sbar[ 484 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1936 ]), &(acadoWorkspace.sbar[ 484 ]), &(acadoWorkspace.sbar[ 488 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1952 ]), &(acadoWorkspace.sbar[ 488 ]), &(acadoWorkspace.sbar[ 492 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1968 ]), &(acadoWorkspace.sbar[ 492 ]), &(acadoWorkspace.sbar[ 496 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 1984 ]), &(acadoWorkspace.sbar[ 496 ]), &(acadoWorkspace.sbar[ 500 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2000 ]), &(acadoWorkspace.sbar[ 500 ]), &(acadoWorkspace.sbar[ 504 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2016 ]), &(acadoWorkspace.sbar[ 504 ]), &(acadoWorkspace.sbar[ 508 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2032 ]), &(acadoWorkspace.sbar[ 508 ]), &(acadoWorkspace.sbar[ 512 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2048 ]), &(acadoWorkspace.sbar[ 512 ]), &(acadoWorkspace.sbar[ 516 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2064 ]), &(acadoWorkspace.sbar[ 516 ]), &(acadoWorkspace.sbar[ 520 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2080 ]), &(acadoWorkspace.sbar[ 520 ]), &(acadoWorkspace.sbar[ 524 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2096 ]), &(acadoWorkspace.sbar[ 524 ]), &(acadoWorkspace.sbar[ 528 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2112 ]), &(acadoWorkspace.sbar[ 528 ]), &(acadoWorkspace.sbar[ 532 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2128 ]), &(acadoWorkspace.sbar[ 532 ]), &(acadoWorkspace.sbar[ 536 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2144 ]), &(acadoWorkspace.sbar[ 536 ]), &(acadoWorkspace.sbar[ 540 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2160 ]), &(acadoWorkspace.sbar[ 540 ]), &(acadoWorkspace.sbar[ 544 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2176 ]), &(acadoWorkspace.sbar[ 544 ]), &(acadoWorkspace.sbar[ 548 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2192 ]), &(acadoWorkspace.sbar[ 548 ]), &(acadoWorkspace.sbar[ 552 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2208 ]), &(acadoWorkspace.sbar[ 552 ]), &(acadoWorkspace.sbar[ 556 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2224 ]), &(acadoWorkspace.sbar[ 556 ]), &(acadoWorkspace.sbar[ 560 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2240 ]), &(acadoWorkspace.sbar[ 560 ]), &(acadoWorkspace.sbar[ 564 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2256 ]), &(acadoWorkspace.sbar[ 564 ]), &(acadoWorkspace.sbar[ 568 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2272 ]), &(acadoWorkspace.sbar[ 568 ]), &(acadoWorkspace.sbar[ 572 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2288 ]), &(acadoWorkspace.sbar[ 572 ]), &(acadoWorkspace.sbar[ 576 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2304 ]), &(acadoWorkspace.sbar[ 576 ]), &(acadoWorkspace.sbar[ 580 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2320 ]), &(acadoWorkspace.sbar[ 580 ]), &(acadoWorkspace.sbar[ 584 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2336 ]), &(acadoWorkspace.sbar[ 584 ]), &(acadoWorkspace.sbar[ 588 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2352 ]), &(acadoWorkspace.sbar[ 588 ]), &(acadoWorkspace.sbar[ 592 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2368 ]), &(acadoWorkspace.sbar[ 592 ]), &(acadoWorkspace.sbar[ 596 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2384 ]), &(acadoWorkspace.sbar[ 596 ]), &(acadoWorkspace.sbar[ 600 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2400 ]), &(acadoWorkspace.sbar[ 600 ]), &(acadoWorkspace.sbar[ 604 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2416 ]), &(acadoWorkspace.sbar[ 604 ]), &(acadoWorkspace.sbar[ 608 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2432 ]), &(acadoWorkspace.sbar[ 608 ]), &(acadoWorkspace.sbar[ 612 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2448 ]), &(acadoWorkspace.sbar[ 612 ]), &(acadoWorkspace.sbar[ 616 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2464 ]), &(acadoWorkspace.sbar[ 616 ]), &(acadoWorkspace.sbar[ 620 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2480 ]), &(acadoWorkspace.sbar[ 620 ]), &(acadoWorkspace.sbar[ 624 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2496 ]), &(acadoWorkspace.sbar[ 624 ]), &(acadoWorkspace.sbar[ 628 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2512 ]), &(acadoWorkspace.sbar[ 628 ]), &(acadoWorkspace.sbar[ 632 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2528 ]), &(acadoWorkspace.sbar[ 632 ]), &(acadoWorkspace.sbar[ 636 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2544 ]), &(acadoWorkspace.sbar[ 636 ]), &(acadoWorkspace.sbar[ 640 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2560 ]), &(acadoWorkspace.sbar[ 640 ]), &(acadoWorkspace.sbar[ 644 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2576 ]), &(acadoWorkspace.sbar[ 644 ]), &(acadoWorkspace.sbar[ 648 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2592 ]), &(acadoWorkspace.sbar[ 648 ]), &(acadoWorkspace.sbar[ 652 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2608 ]), &(acadoWorkspace.sbar[ 652 ]), &(acadoWorkspace.sbar[ 656 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2624 ]), &(acadoWorkspace.sbar[ 656 ]), &(acadoWorkspace.sbar[ 660 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2640 ]), &(acadoWorkspace.sbar[ 660 ]), &(acadoWorkspace.sbar[ 664 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2656 ]), &(acadoWorkspace.sbar[ 664 ]), &(acadoWorkspace.sbar[ 668 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2672 ]), &(acadoWorkspace.sbar[ 668 ]), &(acadoWorkspace.sbar[ 672 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2688 ]), &(acadoWorkspace.sbar[ 672 ]), &(acadoWorkspace.sbar[ 676 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2704 ]), &(acadoWorkspace.sbar[ 676 ]), &(acadoWorkspace.sbar[ 680 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2720 ]), &(acadoWorkspace.sbar[ 680 ]), &(acadoWorkspace.sbar[ 684 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2736 ]), &(acadoWorkspace.sbar[ 684 ]), &(acadoWorkspace.sbar[ 688 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2752 ]), &(acadoWorkspace.sbar[ 688 ]), &(acadoWorkspace.sbar[ 692 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2768 ]), &(acadoWorkspace.sbar[ 692 ]), &(acadoWorkspace.sbar[ 696 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2784 ]), &(acadoWorkspace.sbar[ 696 ]), &(acadoWorkspace.sbar[ 700 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2800 ]), &(acadoWorkspace.sbar[ 700 ]), &(acadoWorkspace.sbar[ 704 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2816 ]), &(acadoWorkspace.sbar[ 704 ]), &(acadoWorkspace.sbar[ 708 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2832 ]), &(acadoWorkspace.sbar[ 708 ]), &(acadoWorkspace.sbar[ 712 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2848 ]), &(acadoWorkspace.sbar[ 712 ]), &(acadoWorkspace.sbar[ 716 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2864 ]), &(acadoWorkspace.sbar[ 716 ]), &(acadoWorkspace.sbar[ 720 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2880 ]), &(acadoWorkspace.sbar[ 720 ]), &(acadoWorkspace.sbar[ 724 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2896 ]), &(acadoWorkspace.sbar[ 724 ]), &(acadoWorkspace.sbar[ 728 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2912 ]), &(acadoWorkspace.sbar[ 728 ]), &(acadoWorkspace.sbar[ 732 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2928 ]), &(acadoWorkspace.sbar[ 732 ]), &(acadoWorkspace.sbar[ 736 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2944 ]), &(acadoWorkspace.sbar[ 736 ]), &(acadoWorkspace.sbar[ 740 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2960 ]), &(acadoWorkspace.sbar[ 740 ]), &(acadoWorkspace.sbar[ 744 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2976 ]), &(acadoWorkspace.sbar[ 744 ]), &(acadoWorkspace.sbar[ 748 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 2992 ]), &(acadoWorkspace.sbar[ 748 ]), &(acadoWorkspace.sbar[ 752 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 3008 ]), &(acadoWorkspace.sbar[ 752 ]), &(acadoWorkspace.sbar[ 756 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 3024 ]), &(acadoWorkspace.sbar[ 756 ]), &(acadoWorkspace.sbar[ 760 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 3040 ]), &(acadoWorkspace.sbar[ 760 ]), &(acadoWorkspace.sbar[ 764 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 3056 ]), &(acadoWorkspace.sbar[ 764 ]), &(acadoWorkspace.sbar[ 768 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 3072 ]), &(acadoWorkspace.sbar[ 768 ]), &(acadoWorkspace.sbar[ 772 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 3088 ]), &(acadoWorkspace.sbar[ 772 ]), &(acadoWorkspace.sbar[ 776 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 3104 ]), &(acadoWorkspace.sbar[ 776 ]), &(acadoWorkspace.sbar[ 780 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 3120 ]), &(acadoWorkspace.sbar[ 780 ]), &(acadoWorkspace.sbar[ 784 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 3136 ]), &(acadoWorkspace.sbar[ 784 ]), &(acadoWorkspace.sbar[ 788 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 3152 ]), &(acadoWorkspace.sbar[ 788 ]), &(acadoWorkspace.sbar[ 792 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 3168 ]), &(acadoWorkspace.sbar[ 792 ]), &(acadoWorkspace.sbar[ 796 ]) );
acado_macASbar( &(acadoWorkspace.evGx[ 3184 ]), &(acadoWorkspace.sbar[ 796 ]), &(acadoWorkspace.sbar[ 800 ]) );

acadoWorkspace.w1[0] = + acadoWorkspace.QN1[0]*acadoWorkspace.sbar[800] + acadoWorkspace.QN1[1]*acadoWorkspace.sbar[801] + acadoWorkspace.QN1[2]*acadoWorkspace.sbar[802] + acadoWorkspace.QN1[3]*acadoWorkspace.sbar[803] + acadoWorkspace.QDy[800];
acadoWorkspace.w1[1] = + acadoWorkspace.QN1[4]*acadoWorkspace.sbar[800] + acadoWorkspace.QN1[5]*acadoWorkspace.sbar[801] + acadoWorkspace.QN1[6]*acadoWorkspace.sbar[802] + acadoWorkspace.QN1[7]*acadoWorkspace.sbar[803] + acadoWorkspace.QDy[801];
acadoWorkspace.w1[2] = + acadoWorkspace.QN1[8]*acadoWorkspace.sbar[800] + acadoWorkspace.QN1[9]*acadoWorkspace.sbar[801] + acadoWorkspace.QN1[10]*acadoWorkspace.sbar[802] + acadoWorkspace.QN1[11]*acadoWorkspace.sbar[803] + acadoWorkspace.QDy[802];
acadoWorkspace.w1[3] = + acadoWorkspace.QN1[12]*acadoWorkspace.sbar[800] + acadoWorkspace.QN1[13]*acadoWorkspace.sbar[801] + acadoWorkspace.QN1[14]*acadoWorkspace.sbar[802] + acadoWorkspace.QN1[15]*acadoWorkspace.sbar[803] + acadoWorkspace.QDy[803];
acado_macBTw1( &(acadoWorkspace.evGu[ 796 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 199 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 3184 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 796 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 3184 ]), &(acadoWorkspace.sbar[ 796 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 792 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 198 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 3168 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 792 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 3168 ]), &(acadoWorkspace.sbar[ 792 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 788 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 197 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 3152 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 788 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 3152 ]), &(acadoWorkspace.sbar[ 788 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 784 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 196 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 3136 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 784 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 3136 ]), &(acadoWorkspace.sbar[ 784 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 780 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 195 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 3120 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 780 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 3120 ]), &(acadoWorkspace.sbar[ 780 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 776 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 194 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 3104 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 776 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 3104 ]), &(acadoWorkspace.sbar[ 776 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 772 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 193 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 3088 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 772 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 3088 ]), &(acadoWorkspace.sbar[ 772 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 768 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 192 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 3072 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 768 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 3072 ]), &(acadoWorkspace.sbar[ 768 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 764 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 191 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 3056 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 764 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 3056 ]), &(acadoWorkspace.sbar[ 764 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 760 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 190 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 3040 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 760 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 3040 ]), &(acadoWorkspace.sbar[ 760 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 756 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 189 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 3024 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 756 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 3024 ]), &(acadoWorkspace.sbar[ 756 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 752 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 188 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 3008 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 752 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 3008 ]), &(acadoWorkspace.sbar[ 752 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 748 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 187 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2992 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 748 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2992 ]), &(acadoWorkspace.sbar[ 748 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 744 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 186 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2976 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 744 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2976 ]), &(acadoWorkspace.sbar[ 744 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 740 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 185 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2960 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 740 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2960 ]), &(acadoWorkspace.sbar[ 740 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 736 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 184 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2944 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 736 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2944 ]), &(acadoWorkspace.sbar[ 736 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 732 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 183 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2928 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 732 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2928 ]), &(acadoWorkspace.sbar[ 732 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 728 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 182 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2912 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 728 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2912 ]), &(acadoWorkspace.sbar[ 728 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 724 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 181 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2896 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 724 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2896 ]), &(acadoWorkspace.sbar[ 724 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 720 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 180 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2880 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 720 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2880 ]), &(acadoWorkspace.sbar[ 720 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 716 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 179 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2864 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 716 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2864 ]), &(acadoWorkspace.sbar[ 716 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 712 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 178 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2848 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 712 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2848 ]), &(acadoWorkspace.sbar[ 712 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 708 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 177 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2832 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 708 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2832 ]), &(acadoWorkspace.sbar[ 708 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 704 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 176 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2816 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 704 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2816 ]), &(acadoWorkspace.sbar[ 704 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 700 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 175 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2800 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 700 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2800 ]), &(acadoWorkspace.sbar[ 700 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 696 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 174 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2784 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 696 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2784 ]), &(acadoWorkspace.sbar[ 696 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 692 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 173 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2768 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 692 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2768 ]), &(acadoWorkspace.sbar[ 692 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 688 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 172 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2752 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 688 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2752 ]), &(acadoWorkspace.sbar[ 688 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 684 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 171 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2736 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 684 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2736 ]), &(acadoWorkspace.sbar[ 684 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 680 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 170 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2720 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 680 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2720 ]), &(acadoWorkspace.sbar[ 680 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 676 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 169 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2704 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 676 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2704 ]), &(acadoWorkspace.sbar[ 676 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 672 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 168 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2688 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 672 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2688 ]), &(acadoWorkspace.sbar[ 672 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 668 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 167 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2672 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 668 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2672 ]), &(acadoWorkspace.sbar[ 668 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 664 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 166 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2656 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 664 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2656 ]), &(acadoWorkspace.sbar[ 664 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 660 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 165 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2640 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 660 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2640 ]), &(acadoWorkspace.sbar[ 660 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 656 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 164 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2624 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 656 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2624 ]), &(acadoWorkspace.sbar[ 656 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 652 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 163 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2608 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 652 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2608 ]), &(acadoWorkspace.sbar[ 652 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 648 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 162 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2592 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 648 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2592 ]), &(acadoWorkspace.sbar[ 648 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 644 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 161 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2576 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 644 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2576 ]), &(acadoWorkspace.sbar[ 644 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 640 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 160 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2560 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 640 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2560 ]), &(acadoWorkspace.sbar[ 640 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 636 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 159 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2544 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 636 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2544 ]), &(acadoWorkspace.sbar[ 636 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 632 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 158 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2528 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 632 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2528 ]), &(acadoWorkspace.sbar[ 632 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 628 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 157 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2512 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 628 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2512 ]), &(acadoWorkspace.sbar[ 628 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 624 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 156 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2496 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 624 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2496 ]), &(acadoWorkspace.sbar[ 624 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 620 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 155 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2480 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 620 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2480 ]), &(acadoWorkspace.sbar[ 620 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 616 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 154 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2464 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 616 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2464 ]), &(acadoWorkspace.sbar[ 616 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 612 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 153 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2448 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 612 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2448 ]), &(acadoWorkspace.sbar[ 612 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 608 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 152 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2432 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 608 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2432 ]), &(acadoWorkspace.sbar[ 608 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 604 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 151 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2416 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 604 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2416 ]), &(acadoWorkspace.sbar[ 604 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 600 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 150 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2400 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 600 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2400 ]), &(acadoWorkspace.sbar[ 600 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 596 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 149 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2384 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 596 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2384 ]), &(acadoWorkspace.sbar[ 596 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 592 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 148 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2368 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 592 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2368 ]), &(acadoWorkspace.sbar[ 592 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 588 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 147 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2352 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 588 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2352 ]), &(acadoWorkspace.sbar[ 588 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 584 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 146 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2336 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 584 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2336 ]), &(acadoWorkspace.sbar[ 584 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 580 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 145 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2320 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 580 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2320 ]), &(acadoWorkspace.sbar[ 580 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 576 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 144 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2304 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 576 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2304 ]), &(acadoWorkspace.sbar[ 576 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 572 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 143 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2288 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 572 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2288 ]), &(acadoWorkspace.sbar[ 572 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 568 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 142 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2272 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 568 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2272 ]), &(acadoWorkspace.sbar[ 568 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 564 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 141 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2256 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 564 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2256 ]), &(acadoWorkspace.sbar[ 564 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 560 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 140 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2240 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 560 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2240 ]), &(acadoWorkspace.sbar[ 560 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 556 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 139 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2224 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 556 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2224 ]), &(acadoWorkspace.sbar[ 556 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 552 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 138 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2208 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 552 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2208 ]), &(acadoWorkspace.sbar[ 552 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 548 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 137 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2192 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 548 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2192 ]), &(acadoWorkspace.sbar[ 548 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 544 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 136 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2176 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 544 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2176 ]), &(acadoWorkspace.sbar[ 544 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 540 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 135 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2160 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 540 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2160 ]), &(acadoWorkspace.sbar[ 540 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 536 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 134 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2144 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 536 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2144 ]), &(acadoWorkspace.sbar[ 536 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 532 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 133 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2128 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 532 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2128 ]), &(acadoWorkspace.sbar[ 532 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 528 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 132 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2112 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 528 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2112 ]), &(acadoWorkspace.sbar[ 528 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 524 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 131 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2096 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 524 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2096 ]), &(acadoWorkspace.sbar[ 524 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 520 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 130 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2080 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 520 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2080 ]), &(acadoWorkspace.sbar[ 520 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 516 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 129 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2064 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 516 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2064 ]), &(acadoWorkspace.sbar[ 516 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 512 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 128 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2048 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 512 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2048 ]), &(acadoWorkspace.sbar[ 512 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 508 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 127 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2032 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 508 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2032 ]), &(acadoWorkspace.sbar[ 508 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 504 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 126 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2016 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 504 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2016 ]), &(acadoWorkspace.sbar[ 504 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 500 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 125 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 2000 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 500 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 2000 ]), &(acadoWorkspace.sbar[ 500 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 496 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 124 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1984 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 496 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1984 ]), &(acadoWorkspace.sbar[ 496 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 492 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 123 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1968 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 492 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1968 ]), &(acadoWorkspace.sbar[ 492 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 488 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 122 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1952 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 488 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1952 ]), &(acadoWorkspace.sbar[ 488 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 484 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 121 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1936 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 484 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1936 ]), &(acadoWorkspace.sbar[ 484 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 480 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 120 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1920 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 480 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1920 ]), &(acadoWorkspace.sbar[ 480 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 476 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 119 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1904 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 476 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1904 ]), &(acadoWorkspace.sbar[ 476 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 472 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 118 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1888 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 472 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1888 ]), &(acadoWorkspace.sbar[ 472 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 468 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 117 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1872 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 468 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1872 ]), &(acadoWorkspace.sbar[ 468 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 464 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 116 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1856 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 464 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1856 ]), &(acadoWorkspace.sbar[ 464 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 460 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 115 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1840 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 460 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1840 ]), &(acadoWorkspace.sbar[ 460 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 456 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 114 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1824 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 456 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1824 ]), &(acadoWorkspace.sbar[ 456 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 452 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 113 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1808 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 452 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1808 ]), &(acadoWorkspace.sbar[ 452 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 448 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 112 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1792 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 448 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1792 ]), &(acadoWorkspace.sbar[ 448 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 444 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 111 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1776 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 444 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1776 ]), &(acadoWorkspace.sbar[ 444 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 440 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 110 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1760 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 440 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1760 ]), &(acadoWorkspace.sbar[ 440 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 436 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 109 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1744 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 436 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1744 ]), &(acadoWorkspace.sbar[ 436 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 432 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 108 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1728 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 432 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1728 ]), &(acadoWorkspace.sbar[ 432 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 428 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 107 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1712 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 428 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1712 ]), &(acadoWorkspace.sbar[ 428 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 424 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 106 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1696 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 424 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1696 ]), &(acadoWorkspace.sbar[ 424 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 420 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 105 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1680 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 420 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1680 ]), &(acadoWorkspace.sbar[ 420 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 416 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 104 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1664 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 416 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1664 ]), &(acadoWorkspace.sbar[ 416 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 412 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 103 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1648 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 412 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1648 ]), &(acadoWorkspace.sbar[ 412 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 408 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 102 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1632 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 408 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1632 ]), &(acadoWorkspace.sbar[ 408 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 404 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 101 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1616 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 404 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1616 ]), &(acadoWorkspace.sbar[ 404 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 400 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 100 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1600 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 400 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1600 ]), &(acadoWorkspace.sbar[ 400 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 396 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 99 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1584 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 396 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1584 ]), &(acadoWorkspace.sbar[ 396 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 392 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 98 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1568 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 392 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1568 ]), &(acadoWorkspace.sbar[ 392 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 388 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 97 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1552 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 388 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1552 ]), &(acadoWorkspace.sbar[ 388 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 384 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 96 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1536 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 384 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1536 ]), &(acadoWorkspace.sbar[ 384 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 380 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 95 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1520 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 380 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1520 ]), &(acadoWorkspace.sbar[ 380 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 376 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 94 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1504 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 376 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1504 ]), &(acadoWorkspace.sbar[ 376 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 372 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 93 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1488 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 372 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1488 ]), &(acadoWorkspace.sbar[ 372 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 368 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 92 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1472 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 368 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1472 ]), &(acadoWorkspace.sbar[ 368 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 364 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 91 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1456 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 364 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1456 ]), &(acadoWorkspace.sbar[ 364 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 360 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 90 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1440 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 360 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1440 ]), &(acadoWorkspace.sbar[ 360 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 356 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 89 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1424 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 356 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1424 ]), &(acadoWorkspace.sbar[ 356 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 352 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 88 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1408 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 352 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1408 ]), &(acadoWorkspace.sbar[ 352 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 348 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 87 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1392 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 348 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1392 ]), &(acadoWorkspace.sbar[ 348 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 344 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 86 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1376 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 344 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1376 ]), &(acadoWorkspace.sbar[ 344 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 340 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 85 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1360 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 340 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1360 ]), &(acadoWorkspace.sbar[ 340 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 336 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 84 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1344 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 336 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1344 ]), &(acadoWorkspace.sbar[ 336 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 332 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 83 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1328 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 332 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1328 ]), &(acadoWorkspace.sbar[ 332 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 328 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 82 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1312 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 328 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1312 ]), &(acadoWorkspace.sbar[ 328 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 324 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 81 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1296 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 324 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1296 ]), &(acadoWorkspace.sbar[ 324 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 320 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 80 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1280 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 320 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1280 ]), &(acadoWorkspace.sbar[ 320 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 316 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 79 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1264 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 316 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1264 ]), &(acadoWorkspace.sbar[ 316 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 312 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 78 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1248 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 312 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1248 ]), &(acadoWorkspace.sbar[ 312 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 308 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 77 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1232 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 308 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1232 ]), &(acadoWorkspace.sbar[ 308 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 304 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 76 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1216 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 304 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1216 ]), &(acadoWorkspace.sbar[ 304 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 300 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 75 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1200 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 300 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1200 ]), &(acadoWorkspace.sbar[ 300 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 296 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 74 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1184 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 296 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1184 ]), &(acadoWorkspace.sbar[ 296 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 292 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 73 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1168 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 292 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1168 ]), &(acadoWorkspace.sbar[ 292 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 288 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 72 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1152 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 288 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1152 ]), &(acadoWorkspace.sbar[ 288 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 284 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 71 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1136 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 284 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1136 ]), &(acadoWorkspace.sbar[ 284 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 280 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 70 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1120 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 280 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1120 ]), &(acadoWorkspace.sbar[ 280 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 276 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 69 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1104 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 276 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1104 ]), &(acadoWorkspace.sbar[ 276 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 272 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 68 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1088 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 272 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1088 ]), &(acadoWorkspace.sbar[ 272 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 268 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 67 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1072 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 268 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1072 ]), &(acadoWorkspace.sbar[ 268 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 264 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 66 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1056 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 264 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1056 ]), &(acadoWorkspace.sbar[ 264 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 260 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 65 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1040 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 260 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1040 ]), &(acadoWorkspace.sbar[ 260 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 256 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 64 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1024 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 256 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1024 ]), &(acadoWorkspace.sbar[ 256 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 252 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 63 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 1008 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 252 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 1008 ]), &(acadoWorkspace.sbar[ 252 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 248 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 62 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 992 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 248 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 992 ]), &(acadoWorkspace.sbar[ 248 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 244 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 61 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 976 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 244 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 976 ]), &(acadoWorkspace.sbar[ 244 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 240 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 60 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 960 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 240 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 960 ]), &(acadoWorkspace.sbar[ 240 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 236 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 59 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 944 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 236 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 944 ]), &(acadoWorkspace.sbar[ 236 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 232 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 58 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 928 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 232 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 928 ]), &(acadoWorkspace.sbar[ 232 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 228 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 57 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 912 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 228 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 912 ]), &(acadoWorkspace.sbar[ 228 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 224 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 56 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 896 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 224 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 896 ]), &(acadoWorkspace.sbar[ 224 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 220 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 55 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 880 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 220 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 880 ]), &(acadoWorkspace.sbar[ 220 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 216 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 54 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 864 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 216 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 864 ]), &(acadoWorkspace.sbar[ 216 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 212 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 53 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 848 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 212 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 848 ]), &(acadoWorkspace.sbar[ 212 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 208 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 52 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 832 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 208 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 832 ]), &(acadoWorkspace.sbar[ 208 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 204 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 51 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 816 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 204 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 816 ]), &(acadoWorkspace.sbar[ 204 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 200 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 50 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 800 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 200 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 800 ]), &(acadoWorkspace.sbar[ 200 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 196 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 49 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 784 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 196 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 784 ]), &(acadoWorkspace.sbar[ 196 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 192 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 48 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 768 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 192 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 768 ]), &(acadoWorkspace.sbar[ 192 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 188 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 47 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 752 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 188 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 752 ]), &(acadoWorkspace.sbar[ 188 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 184 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 46 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 736 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 184 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 736 ]), &(acadoWorkspace.sbar[ 184 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 180 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 45 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 720 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 180 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 720 ]), &(acadoWorkspace.sbar[ 180 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 176 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 44 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 704 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 176 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 704 ]), &(acadoWorkspace.sbar[ 176 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 172 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 43 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 688 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 172 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 688 ]), &(acadoWorkspace.sbar[ 172 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 168 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 42 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 672 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 168 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 672 ]), &(acadoWorkspace.sbar[ 168 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 164 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 41 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 656 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 164 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 656 ]), &(acadoWorkspace.sbar[ 164 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 160 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 40 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 640 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 160 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 640 ]), &(acadoWorkspace.sbar[ 160 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 156 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 39 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 624 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 156 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 624 ]), &(acadoWorkspace.sbar[ 156 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 152 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 38 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 608 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 152 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 608 ]), &(acadoWorkspace.sbar[ 152 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 148 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 37 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 592 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 148 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 592 ]), &(acadoWorkspace.sbar[ 148 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 144 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 36 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 576 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 144 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 576 ]), &(acadoWorkspace.sbar[ 144 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 140 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 35 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 560 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 140 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 560 ]), &(acadoWorkspace.sbar[ 140 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 136 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 34 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 544 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 136 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 544 ]), &(acadoWorkspace.sbar[ 136 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 132 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 33 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 528 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 132 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 528 ]), &(acadoWorkspace.sbar[ 132 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 128 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 32 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 512 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 128 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 512 ]), &(acadoWorkspace.sbar[ 128 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 124 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 31 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 496 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 124 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 496 ]), &(acadoWorkspace.sbar[ 124 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 120 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 30 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 480 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 120 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 480 ]), &(acadoWorkspace.sbar[ 120 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 116 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 29 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 464 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 116 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 464 ]), &(acadoWorkspace.sbar[ 116 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 112 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 28 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 448 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 112 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 448 ]), &(acadoWorkspace.sbar[ 112 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 108 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 27 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 432 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 108 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 432 ]), &(acadoWorkspace.sbar[ 108 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 104 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 26 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 416 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 104 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 416 ]), &(acadoWorkspace.sbar[ 104 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 100 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 25 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 400 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 100 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 400 ]), &(acadoWorkspace.sbar[ 100 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 96 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 24 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 384 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 96 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 384 ]), &(acadoWorkspace.sbar[ 96 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 92 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 23 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 368 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 92 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 368 ]), &(acadoWorkspace.sbar[ 92 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 88 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 22 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 352 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 88 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 352 ]), &(acadoWorkspace.sbar[ 88 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 84 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 21 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 336 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 84 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 336 ]), &(acadoWorkspace.sbar[ 84 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
acado_macBTw1( &(acadoWorkspace.evGu[ 80 ]), acadoWorkspace.w1, &(acadoWorkspace.g[ 20 ]) );
acado_macATw1QDy( &(acadoWorkspace.evGx[ 320 ]), acadoWorkspace.w1, &(acadoWorkspace.QDy[ 80 ]), acadoWorkspace.w2 );
acado_macQSbarW2( &(acadoWorkspace.Q1[ 320 ]), &(acadoWorkspace.sbar[ 80 ]), acadoWorkspace.w2, acadoWorkspace.w1 );
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
acadoWorkspace.lb[20] = acadoVariables.lbValues[20] - acadoVariables.u[20];
acadoWorkspace.lb[21] = acadoVariables.lbValues[21] - acadoVariables.u[21];
acadoWorkspace.lb[22] = acadoVariables.lbValues[22] - acadoVariables.u[22];
acadoWorkspace.lb[23] = acadoVariables.lbValues[23] - acadoVariables.u[23];
acadoWorkspace.lb[24] = acadoVariables.lbValues[24] - acadoVariables.u[24];
acadoWorkspace.lb[25] = acadoVariables.lbValues[25] - acadoVariables.u[25];
acadoWorkspace.lb[26] = acadoVariables.lbValues[26] - acadoVariables.u[26];
acadoWorkspace.lb[27] = acadoVariables.lbValues[27] - acadoVariables.u[27];
acadoWorkspace.lb[28] = acadoVariables.lbValues[28] - acadoVariables.u[28];
acadoWorkspace.lb[29] = acadoVariables.lbValues[29] - acadoVariables.u[29];
acadoWorkspace.lb[30] = acadoVariables.lbValues[30] - acadoVariables.u[30];
acadoWorkspace.lb[31] = acadoVariables.lbValues[31] - acadoVariables.u[31];
acadoWorkspace.lb[32] = acadoVariables.lbValues[32] - acadoVariables.u[32];
acadoWorkspace.lb[33] = acadoVariables.lbValues[33] - acadoVariables.u[33];
acadoWorkspace.lb[34] = acadoVariables.lbValues[34] - acadoVariables.u[34];
acadoWorkspace.lb[35] = acadoVariables.lbValues[35] - acadoVariables.u[35];
acadoWorkspace.lb[36] = acadoVariables.lbValues[36] - acadoVariables.u[36];
acadoWorkspace.lb[37] = acadoVariables.lbValues[37] - acadoVariables.u[37];
acadoWorkspace.lb[38] = acadoVariables.lbValues[38] - acadoVariables.u[38];
acadoWorkspace.lb[39] = acadoVariables.lbValues[39] - acadoVariables.u[39];
acadoWorkspace.lb[40] = acadoVariables.lbValues[40] - acadoVariables.u[40];
acadoWorkspace.lb[41] = acadoVariables.lbValues[41] - acadoVariables.u[41];
acadoWorkspace.lb[42] = acadoVariables.lbValues[42] - acadoVariables.u[42];
acadoWorkspace.lb[43] = acadoVariables.lbValues[43] - acadoVariables.u[43];
acadoWorkspace.lb[44] = acadoVariables.lbValues[44] - acadoVariables.u[44];
acadoWorkspace.lb[45] = acadoVariables.lbValues[45] - acadoVariables.u[45];
acadoWorkspace.lb[46] = acadoVariables.lbValues[46] - acadoVariables.u[46];
acadoWorkspace.lb[47] = acadoVariables.lbValues[47] - acadoVariables.u[47];
acadoWorkspace.lb[48] = acadoVariables.lbValues[48] - acadoVariables.u[48];
acadoWorkspace.lb[49] = acadoVariables.lbValues[49] - acadoVariables.u[49];
acadoWorkspace.lb[50] = acadoVariables.lbValues[50] - acadoVariables.u[50];
acadoWorkspace.lb[51] = acadoVariables.lbValues[51] - acadoVariables.u[51];
acadoWorkspace.lb[52] = acadoVariables.lbValues[52] - acadoVariables.u[52];
acadoWorkspace.lb[53] = acadoVariables.lbValues[53] - acadoVariables.u[53];
acadoWorkspace.lb[54] = acadoVariables.lbValues[54] - acadoVariables.u[54];
acadoWorkspace.lb[55] = acadoVariables.lbValues[55] - acadoVariables.u[55];
acadoWorkspace.lb[56] = acadoVariables.lbValues[56] - acadoVariables.u[56];
acadoWorkspace.lb[57] = acadoVariables.lbValues[57] - acadoVariables.u[57];
acadoWorkspace.lb[58] = acadoVariables.lbValues[58] - acadoVariables.u[58];
acadoWorkspace.lb[59] = acadoVariables.lbValues[59] - acadoVariables.u[59];
acadoWorkspace.lb[60] = acadoVariables.lbValues[60] - acadoVariables.u[60];
acadoWorkspace.lb[61] = acadoVariables.lbValues[61] - acadoVariables.u[61];
acadoWorkspace.lb[62] = acadoVariables.lbValues[62] - acadoVariables.u[62];
acadoWorkspace.lb[63] = acadoVariables.lbValues[63] - acadoVariables.u[63];
acadoWorkspace.lb[64] = acadoVariables.lbValues[64] - acadoVariables.u[64];
acadoWorkspace.lb[65] = acadoVariables.lbValues[65] - acadoVariables.u[65];
acadoWorkspace.lb[66] = acadoVariables.lbValues[66] - acadoVariables.u[66];
acadoWorkspace.lb[67] = acadoVariables.lbValues[67] - acadoVariables.u[67];
acadoWorkspace.lb[68] = acadoVariables.lbValues[68] - acadoVariables.u[68];
acadoWorkspace.lb[69] = acadoVariables.lbValues[69] - acadoVariables.u[69];
acadoWorkspace.lb[70] = acadoVariables.lbValues[70] - acadoVariables.u[70];
acadoWorkspace.lb[71] = acadoVariables.lbValues[71] - acadoVariables.u[71];
acadoWorkspace.lb[72] = acadoVariables.lbValues[72] - acadoVariables.u[72];
acadoWorkspace.lb[73] = acadoVariables.lbValues[73] - acadoVariables.u[73];
acadoWorkspace.lb[74] = acadoVariables.lbValues[74] - acadoVariables.u[74];
acadoWorkspace.lb[75] = acadoVariables.lbValues[75] - acadoVariables.u[75];
acadoWorkspace.lb[76] = acadoVariables.lbValues[76] - acadoVariables.u[76];
acadoWorkspace.lb[77] = acadoVariables.lbValues[77] - acadoVariables.u[77];
acadoWorkspace.lb[78] = acadoVariables.lbValues[78] - acadoVariables.u[78];
acadoWorkspace.lb[79] = acadoVariables.lbValues[79] - acadoVariables.u[79];
acadoWorkspace.lb[80] = acadoVariables.lbValues[80] - acadoVariables.u[80];
acadoWorkspace.lb[81] = acadoVariables.lbValues[81] - acadoVariables.u[81];
acadoWorkspace.lb[82] = acadoVariables.lbValues[82] - acadoVariables.u[82];
acadoWorkspace.lb[83] = acadoVariables.lbValues[83] - acadoVariables.u[83];
acadoWorkspace.lb[84] = acadoVariables.lbValues[84] - acadoVariables.u[84];
acadoWorkspace.lb[85] = acadoVariables.lbValues[85] - acadoVariables.u[85];
acadoWorkspace.lb[86] = acadoVariables.lbValues[86] - acadoVariables.u[86];
acadoWorkspace.lb[87] = acadoVariables.lbValues[87] - acadoVariables.u[87];
acadoWorkspace.lb[88] = acadoVariables.lbValues[88] - acadoVariables.u[88];
acadoWorkspace.lb[89] = acadoVariables.lbValues[89] - acadoVariables.u[89];
acadoWorkspace.lb[90] = acadoVariables.lbValues[90] - acadoVariables.u[90];
acadoWorkspace.lb[91] = acadoVariables.lbValues[91] - acadoVariables.u[91];
acadoWorkspace.lb[92] = acadoVariables.lbValues[92] - acadoVariables.u[92];
acadoWorkspace.lb[93] = acadoVariables.lbValues[93] - acadoVariables.u[93];
acadoWorkspace.lb[94] = acadoVariables.lbValues[94] - acadoVariables.u[94];
acadoWorkspace.lb[95] = acadoVariables.lbValues[95] - acadoVariables.u[95];
acadoWorkspace.lb[96] = acadoVariables.lbValues[96] - acadoVariables.u[96];
acadoWorkspace.lb[97] = acadoVariables.lbValues[97] - acadoVariables.u[97];
acadoWorkspace.lb[98] = acadoVariables.lbValues[98] - acadoVariables.u[98];
acadoWorkspace.lb[99] = acadoVariables.lbValues[99] - acadoVariables.u[99];
acadoWorkspace.lb[100] = acadoVariables.lbValues[100] - acadoVariables.u[100];
acadoWorkspace.lb[101] = acadoVariables.lbValues[101] - acadoVariables.u[101];
acadoWorkspace.lb[102] = acadoVariables.lbValues[102] - acadoVariables.u[102];
acadoWorkspace.lb[103] = acadoVariables.lbValues[103] - acadoVariables.u[103];
acadoWorkspace.lb[104] = acadoVariables.lbValues[104] - acadoVariables.u[104];
acadoWorkspace.lb[105] = acadoVariables.lbValues[105] - acadoVariables.u[105];
acadoWorkspace.lb[106] = acadoVariables.lbValues[106] - acadoVariables.u[106];
acadoWorkspace.lb[107] = acadoVariables.lbValues[107] - acadoVariables.u[107];
acadoWorkspace.lb[108] = acadoVariables.lbValues[108] - acadoVariables.u[108];
acadoWorkspace.lb[109] = acadoVariables.lbValues[109] - acadoVariables.u[109];
acadoWorkspace.lb[110] = acadoVariables.lbValues[110] - acadoVariables.u[110];
acadoWorkspace.lb[111] = acadoVariables.lbValues[111] - acadoVariables.u[111];
acadoWorkspace.lb[112] = acadoVariables.lbValues[112] - acadoVariables.u[112];
acadoWorkspace.lb[113] = acadoVariables.lbValues[113] - acadoVariables.u[113];
acadoWorkspace.lb[114] = acadoVariables.lbValues[114] - acadoVariables.u[114];
acadoWorkspace.lb[115] = acadoVariables.lbValues[115] - acadoVariables.u[115];
acadoWorkspace.lb[116] = acadoVariables.lbValues[116] - acadoVariables.u[116];
acadoWorkspace.lb[117] = acadoVariables.lbValues[117] - acadoVariables.u[117];
acadoWorkspace.lb[118] = acadoVariables.lbValues[118] - acadoVariables.u[118];
acadoWorkspace.lb[119] = acadoVariables.lbValues[119] - acadoVariables.u[119];
acadoWorkspace.lb[120] = acadoVariables.lbValues[120] - acadoVariables.u[120];
acadoWorkspace.lb[121] = acadoVariables.lbValues[121] - acadoVariables.u[121];
acadoWorkspace.lb[122] = acadoVariables.lbValues[122] - acadoVariables.u[122];
acadoWorkspace.lb[123] = acadoVariables.lbValues[123] - acadoVariables.u[123];
acadoWorkspace.lb[124] = acadoVariables.lbValues[124] - acadoVariables.u[124];
acadoWorkspace.lb[125] = acadoVariables.lbValues[125] - acadoVariables.u[125];
acadoWorkspace.lb[126] = acadoVariables.lbValues[126] - acadoVariables.u[126];
acadoWorkspace.lb[127] = acadoVariables.lbValues[127] - acadoVariables.u[127];
acadoWorkspace.lb[128] = acadoVariables.lbValues[128] - acadoVariables.u[128];
acadoWorkspace.lb[129] = acadoVariables.lbValues[129] - acadoVariables.u[129];
acadoWorkspace.lb[130] = acadoVariables.lbValues[130] - acadoVariables.u[130];
acadoWorkspace.lb[131] = acadoVariables.lbValues[131] - acadoVariables.u[131];
acadoWorkspace.lb[132] = acadoVariables.lbValues[132] - acadoVariables.u[132];
acadoWorkspace.lb[133] = acadoVariables.lbValues[133] - acadoVariables.u[133];
acadoWorkspace.lb[134] = acadoVariables.lbValues[134] - acadoVariables.u[134];
acadoWorkspace.lb[135] = acadoVariables.lbValues[135] - acadoVariables.u[135];
acadoWorkspace.lb[136] = acadoVariables.lbValues[136] - acadoVariables.u[136];
acadoWorkspace.lb[137] = acadoVariables.lbValues[137] - acadoVariables.u[137];
acadoWorkspace.lb[138] = acadoVariables.lbValues[138] - acadoVariables.u[138];
acadoWorkspace.lb[139] = acadoVariables.lbValues[139] - acadoVariables.u[139];
acadoWorkspace.lb[140] = acadoVariables.lbValues[140] - acadoVariables.u[140];
acadoWorkspace.lb[141] = acadoVariables.lbValues[141] - acadoVariables.u[141];
acadoWorkspace.lb[142] = acadoVariables.lbValues[142] - acadoVariables.u[142];
acadoWorkspace.lb[143] = acadoVariables.lbValues[143] - acadoVariables.u[143];
acadoWorkspace.lb[144] = acadoVariables.lbValues[144] - acadoVariables.u[144];
acadoWorkspace.lb[145] = acadoVariables.lbValues[145] - acadoVariables.u[145];
acadoWorkspace.lb[146] = acadoVariables.lbValues[146] - acadoVariables.u[146];
acadoWorkspace.lb[147] = acadoVariables.lbValues[147] - acadoVariables.u[147];
acadoWorkspace.lb[148] = acadoVariables.lbValues[148] - acadoVariables.u[148];
acadoWorkspace.lb[149] = acadoVariables.lbValues[149] - acadoVariables.u[149];
acadoWorkspace.lb[150] = acadoVariables.lbValues[150] - acadoVariables.u[150];
acadoWorkspace.lb[151] = acadoVariables.lbValues[151] - acadoVariables.u[151];
acadoWorkspace.lb[152] = acadoVariables.lbValues[152] - acadoVariables.u[152];
acadoWorkspace.lb[153] = acadoVariables.lbValues[153] - acadoVariables.u[153];
acadoWorkspace.lb[154] = acadoVariables.lbValues[154] - acadoVariables.u[154];
acadoWorkspace.lb[155] = acadoVariables.lbValues[155] - acadoVariables.u[155];
acadoWorkspace.lb[156] = acadoVariables.lbValues[156] - acadoVariables.u[156];
acadoWorkspace.lb[157] = acadoVariables.lbValues[157] - acadoVariables.u[157];
acadoWorkspace.lb[158] = acadoVariables.lbValues[158] - acadoVariables.u[158];
acadoWorkspace.lb[159] = acadoVariables.lbValues[159] - acadoVariables.u[159];
acadoWorkspace.lb[160] = acadoVariables.lbValues[160] - acadoVariables.u[160];
acadoWorkspace.lb[161] = acadoVariables.lbValues[161] - acadoVariables.u[161];
acadoWorkspace.lb[162] = acadoVariables.lbValues[162] - acadoVariables.u[162];
acadoWorkspace.lb[163] = acadoVariables.lbValues[163] - acadoVariables.u[163];
acadoWorkspace.lb[164] = acadoVariables.lbValues[164] - acadoVariables.u[164];
acadoWorkspace.lb[165] = acadoVariables.lbValues[165] - acadoVariables.u[165];
acadoWorkspace.lb[166] = acadoVariables.lbValues[166] - acadoVariables.u[166];
acadoWorkspace.lb[167] = acadoVariables.lbValues[167] - acadoVariables.u[167];
acadoWorkspace.lb[168] = acadoVariables.lbValues[168] - acadoVariables.u[168];
acadoWorkspace.lb[169] = acadoVariables.lbValues[169] - acadoVariables.u[169];
acadoWorkspace.lb[170] = acadoVariables.lbValues[170] - acadoVariables.u[170];
acadoWorkspace.lb[171] = acadoVariables.lbValues[171] - acadoVariables.u[171];
acadoWorkspace.lb[172] = acadoVariables.lbValues[172] - acadoVariables.u[172];
acadoWorkspace.lb[173] = acadoVariables.lbValues[173] - acadoVariables.u[173];
acadoWorkspace.lb[174] = acadoVariables.lbValues[174] - acadoVariables.u[174];
acadoWorkspace.lb[175] = acadoVariables.lbValues[175] - acadoVariables.u[175];
acadoWorkspace.lb[176] = acadoVariables.lbValues[176] - acadoVariables.u[176];
acadoWorkspace.lb[177] = acadoVariables.lbValues[177] - acadoVariables.u[177];
acadoWorkspace.lb[178] = acadoVariables.lbValues[178] - acadoVariables.u[178];
acadoWorkspace.lb[179] = acadoVariables.lbValues[179] - acadoVariables.u[179];
acadoWorkspace.lb[180] = acadoVariables.lbValues[180] - acadoVariables.u[180];
acadoWorkspace.lb[181] = acadoVariables.lbValues[181] - acadoVariables.u[181];
acadoWorkspace.lb[182] = acadoVariables.lbValues[182] - acadoVariables.u[182];
acadoWorkspace.lb[183] = acadoVariables.lbValues[183] - acadoVariables.u[183];
acadoWorkspace.lb[184] = acadoVariables.lbValues[184] - acadoVariables.u[184];
acadoWorkspace.lb[185] = acadoVariables.lbValues[185] - acadoVariables.u[185];
acadoWorkspace.lb[186] = acadoVariables.lbValues[186] - acadoVariables.u[186];
acadoWorkspace.lb[187] = acadoVariables.lbValues[187] - acadoVariables.u[187];
acadoWorkspace.lb[188] = acadoVariables.lbValues[188] - acadoVariables.u[188];
acadoWorkspace.lb[189] = acadoVariables.lbValues[189] - acadoVariables.u[189];
acadoWorkspace.lb[190] = acadoVariables.lbValues[190] - acadoVariables.u[190];
acadoWorkspace.lb[191] = acadoVariables.lbValues[191] - acadoVariables.u[191];
acadoWorkspace.lb[192] = acadoVariables.lbValues[192] - acadoVariables.u[192];
acadoWorkspace.lb[193] = acadoVariables.lbValues[193] - acadoVariables.u[193];
acadoWorkspace.lb[194] = acadoVariables.lbValues[194] - acadoVariables.u[194];
acadoWorkspace.lb[195] = acadoVariables.lbValues[195] - acadoVariables.u[195];
acadoWorkspace.lb[196] = acadoVariables.lbValues[196] - acadoVariables.u[196];
acadoWorkspace.lb[197] = acadoVariables.lbValues[197] - acadoVariables.u[197];
acadoWorkspace.lb[198] = acadoVariables.lbValues[198] - acadoVariables.u[198];
acadoWorkspace.lb[199] = acadoVariables.lbValues[199] - acadoVariables.u[199];
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
acadoWorkspace.ub[20] = acadoVariables.ubValues[20] - acadoVariables.u[20];
acadoWorkspace.ub[21] = acadoVariables.ubValues[21] - acadoVariables.u[21];
acadoWorkspace.ub[22] = acadoVariables.ubValues[22] - acadoVariables.u[22];
acadoWorkspace.ub[23] = acadoVariables.ubValues[23] - acadoVariables.u[23];
acadoWorkspace.ub[24] = acadoVariables.ubValues[24] - acadoVariables.u[24];
acadoWorkspace.ub[25] = acadoVariables.ubValues[25] - acadoVariables.u[25];
acadoWorkspace.ub[26] = acadoVariables.ubValues[26] - acadoVariables.u[26];
acadoWorkspace.ub[27] = acadoVariables.ubValues[27] - acadoVariables.u[27];
acadoWorkspace.ub[28] = acadoVariables.ubValues[28] - acadoVariables.u[28];
acadoWorkspace.ub[29] = acadoVariables.ubValues[29] - acadoVariables.u[29];
acadoWorkspace.ub[30] = acadoVariables.ubValues[30] - acadoVariables.u[30];
acadoWorkspace.ub[31] = acadoVariables.ubValues[31] - acadoVariables.u[31];
acadoWorkspace.ub[32] = acadoVariables.ubValues[32] - acadoVariables.u[32];
acadoWorkspace.ub[33] = acadoVariables.ubValues[33] - acadoVariables.u[33];
acadoWorkspace.ub[34] = acadoVariables.ubValues[34] - acadoVariables.u[34];
acadoWorkspace.ub[35] = acadoVariables.ubValues[35] - acadoVariables.u[35];
acadoWorkspace.ub[36] = acadoVariables.ubValues[36] - acadoVariables.u[36];
acadoWorkspace.ub[37] = acadoVariables.ubValues[37] - acadoVariables.u[37];
acadoWorkspace.ub[38] = acadoVariables.ubValues[38] - acadoVariables.u[38];
acadoWorkspace.ub[39] = acadoVariables.ubValues[39] - acadoVariables.u[39];
acadoWorkspace.ub[40] = acadoVariables.ubValues[40] - acadoVariables.u[40];
acadoWorkspace.ub[41] = acadoVariables.ubValues[41] - acadoVariables.u[41];
acadoWorkspace.ub[42] = acadoVariables.ubValues[42] - acadoVariables.u[42];
acadoWorkspace.ub[43] = acadoVariables.ubValues[43] - acadoVariables.u[43];
acadoWorkspace.ub[44] = acadoVariables.ubValues[44] - acadoVariables.u[44];
acadoWorkspace.ub[45] = acadoVariables.ubValues[45] - acadoVariables.u[45];
acadoWorkspace.ub[46] = acadoVariables.ubValues[46] - acadoVariables.u[46];
acadoWorkspace.ub[47] = acadoVariables.ubValues[47] - acadoVariables.u[47];
acadoWorkspace.ub[48] = acadoVariables.ubValues[48] - acadoVariables.u[48];
acadoWorkspace.ub[49] = acadoVariables.ubValues[49] - acadoVariables.u[49];
acadoWorkspace.ub[50] = acadoVariables.ubValues[50] - acadoVariables.u[50];
acadoWorkspace.ub[51] = acadoVariables.ubValues[51] - acadoVariables.u[51];
acadoWorkspace.ub[52] = acadoVariables.ubValues[52] - acadoVariables.u[52];
acadoWorkspace.ub[53] = acadoVariables.ubValues[53] - acadoVariables.u[53];
acadoWorkspace.ub[54] = acadoVariables.ubValues[54] - acadoVariables.u[54];
acadoWorkspace.ub[55] = acadoVariables.ubValues[55] - acadoVariables.u[55];
acadoWorkspace.ub[56] = acadoVariables.ubValues[56] - acadoVariables.u[56];
acadoWorkspace.ub[57] = acadoVariables.ubValues[57] - acadoVariables.u[57];
acadoWorkspace.ub[58] = acadoVariables.ubValues[58] - acadoVariables.u[58];
acadoWorkspace.ub[59] = acadoVariables.ubValues[59] - acadoVariables.u[59];
acadoWorkspace.ub[60] = acadoVariables.ubValues[60] - acadoVariables.u[60];
acadoWorkspace.ub[61] = acadoVariables.ubValues[61] - acadoVariables.u[61];
acadoWorkspace.ub[62] = acadoVariables.ubValues[62] - acadoVariables.u[62];
acadoWorkspace.ub[63] = acadoVariables.ubValues[63] - acadoVariables.u[63];
acadoWorkspace.ub[64] = acadoVariables.ubValues[64] - acadoVariables.u[64];
acadoWorkspace.ub[65] = acadoVariables.ubValues[65] - acadoVariables.u[65];
acadoWorkspace.ub[66] = acadoVariables.ubValues[66] - acadoVariables.u[66];
acadoWorkspace.ub[67] = acadoVariables.ubValues[67] - acadoVariables.u[67];
acadoWorkspace.ub[68] = acadoVariables.ubValues[68] - acadoVariables.u[68];
acadoWorkspace.ub[69] = acadoVariables.ubValues[69] - acadoVariables.u[69];
acadoWorkspace.ub[70] = acadoVariables.ubValues[70] - acadoVariables.u[70];
acadoWorkspace.ub[71] = acadoVariables.ubValues[71] - acadoVariables.u[71];
acadoWorkspace.ub[72] = acadoVariables.ubValues[72] - acadoVariables.u[72];
acadoWorkspace.ub[73] = acadoVariables.ubValues[73] - acadoVariables.u[73];
acadoWorkspace.ub[74] = acadoVariables.ubValues[74] - acadoVariables.u[74];
acadoWorkspace.ub[75] = acadoVariables.ubValues[75] - acadoVariables.u[75];
acadoWorkspace.ub[76] = acadoVariables.ubValues[76] - acadoVariables.u[76];
acadoWorkspace.ub[77] = acadoVariables.ubValues[77] - acadoVariables.u[77];
acadoWorkspace.ub[78] = acadoVariables.ubValues[78] - acadoVariables.u[78];
acadoWorkspace.ub[79] = acadoVariables.ubValues[79] - acadoVariables.u[79];
acadoWorkspace.ub[80] = acadoVariables.ubValues[80] - acadoVariables.u[80];
acadoWorkspace.ub[81] = acadoVariables.ubValues[81] - acadoVariables.u[81];
acadoWorkspace.ub[82] = acadoVariables.ubValues[82] - acadoVariables.u[82];
acadoWorkspace.ub[83] = acadoVariables.ubValues[83] - acadoVariables.u[83];
acadoWorkspace.ub[84] = acadoVariables.ubValues[84] - acadoVariables.u[84];
acadoWorkspace.ub[85] = acadoVariables.ubValues[85] - acadoVariables.u[85];
acadoWorkspace.ub[86] = acadoVariables.ubValues[86] - acadoVariables.u[86];
acadoWorkspace.ub[87] = acadoVariables.ubValues[87] - acadoVariables.u[87];
acadoWorkspace.ub[88] = acadoVariables.ubValues[88] - acadoVariables.u[88];
acadoWorkspace.ub[89] = acadoVariables.ubValues[89] - acadoVariables.u[89];
acadoWorkspace.ub[90] = acadoVariables.ubValues[90] - acadoVariables.u[90];
acadoWorkspace.ub[91] = acadoVariables.ubValues[91] - acadoVariables.u[91];
acadoWorkspace.ub[92] = acadoVariables.ubValues[92] - acadoVariables.u[92];
acadoWorkspace.ub[93] = acadoVariables.ubValues[93] - acadoVariables.u[93];
acadoWorkspace.ub[94] = acadoVariables.ubValues[94] - acadoVariables.u[94];
acadoWorkspace.ub[95] = acadoVariables.ubValues[95] - acadoVariables.u[95];
acadoWorkspace.ub[96] = acadoVariables.ubValues[96] - acadoVariables.u[96];
acadoWorkspace.ub[97] = acadoVariables.ubValues[97] - acadoVariables.u[97];
acadoWorkspace.ub[98] = acadoVariables.ubValues[98] - acadoVariables.u[98];
acadoWorkspace.ub[99] = acadoVariables.ubValues[99] - acadoVariables.u[99];
acadoWorkspace.ub[100] = acadoVariables.ubValues[100] - acadoVariables.u[100];
acadoWorkspace.ub[101] = acadoVariables.ubValues[101] - acadoVariables.u[101];
acadoWorkspace.ub[102] = acadoVariables.ubValues[102] - acadoVariables.u[102];
acadoWorkspace.ub[103] = acadoVariables.ubValues[103] - acadoVariables.u[103];
acadoWorkspace.ub[104] = acadoVariables.ubValues[104] - acadoVariables.u[104];
acadoWorkspace.ub[105] = acadoVariables.ubValues[105] - acadoVariables.u[105];
acadoWorkspace.ub[106] = acadoVariables.ubValues[106] - acadoVariables.u[106];
acadoWorkspace.ub[107] = acadoVariables.ubValues[107] - acadoVariables.u[107];
acadoWorkspace.ub[108] = acadoVariables.ubValues[108] - acadoVariables.u[108];
acadoWorkspace.ub[109] = acadoVariables.ubValues[109] - acadoVariables.u[109];
acadoWorkspace.ub[110] = acadoVariables.ubValues[110] - acadoVariables.u[110];
acadoWorkspace.ub[111] = acadoVariables.ubValues[111] - acadoVariables.u[111];
acadoWorkspace.ub[112] = acadoVariables.ubValues[112] - acadoVariables.u[112];
acadoWorkspace.ub[113] = acadoVariables.ubValues[113] - acadoVariables.u[113];
acadoWorkspace.ub[114] = acadoVariables.ubValues[114] - acadoVariables.u[114];
acadoWorkspace.ub[115] = acadoVariables.ubValues[115] - acadoVariables.u[115];
acadoWorkspace.ub[116] = acadoVariables.ubValues[116] - acadoVariables.u[116];
acadoWorkspace.ub[117] = acadoVariables.ubValues[117] - acadoVariables.u[117];
acadoWorkspace.ub[118] = acadoVariables.ubValues[118] - acadoVariables.u[118];
acadoWorkspace.ub[119] = acadoVariables.ubValues[119] - acadoVariables.u[119];
acadoWorkspace.ub[120] = acadoVariables.ubValues[120] - acadoVariables.u[120];
acadoWorkspace.ub[121] = acadoVariables.ubValues[121] - acadoVariables.u[121];
acadoWorkspace.ub[122] = acadoVariables.ubValues[122] - acadoVariables.u[122];
acadoWorkspace.ub[123] = acadoVariables.ubValues[123] - acadoVariables.u[123];
acadoWorkspace.ub[124] = acadoVariables.ubValues[124] - acadoVariables.u[124];
acadoWorkspace.ub[125] = acadoVariables.ubValues[125] - acadoVariables.u[125];
acadoWorkspace.ub[126] = acadoVariables.ubValues[126] - acadoVariables.u[126];
acadoWorkspace.ub[127] = acadoVariables.ubValues[127] - acadoVariables.u[127];
acadoWorkspace.ub[128] = acadoVariables.ubValues[128] - acadoVariables.u[128];
acadoWorkspace.ub[129] = acadoVariables.ubValues[129] - acadoVariables.u[129];
acadoWorkspace.ub[130] = acadoVariables.ubValues[130] - acadoVariables.u[130];
acadoWorkspace.ub[131] = acadoVariables.ubValues[131] - acadoVariables.u[131];
acadoWorkspace.ub[132] = acadoVariables.ubValues[132] - acadoVariables.u[132];
acadoWorkspace.ub[133] = acadoVariables.ubValues[133] - acadoVariables.u[133];
acadoWorkspace.ub[134] = acadoVariables.ubValues[134] - acadoVariables.u[134];
acadoWorkspace.ub[135] = acadoVariables.ubValues[135] - acadoVariables.u[135];
acadoWorkspace.ub[136] = acadoVariables.ubValues[136] - acadoVariables.u[136];
acadoWorkspace.ub[137] = acadoVariables.ubValues[137] - acadoVariables.u[137];
acadoWorkspace.ub[138] = acadoVariables.ubValues[138] - acadoVariables.u[138];
acadoWorkspace.ub[139] = acadoVariables.ubValues[139] - acadoVariables.u[139];
acadoWorkspace.ub[140] = acadoVariables.ubValues[140] - acadoVariables.u[140];
acadoWorkspace.ub[141] = acadoVariables.ubValues[141] - acadoVariables.u[141];
acadoWorkspace.ub[142] = acadoVariables.ubValues[142] - acadoVariables.u[142];
acadoWorkspace.ub[143] = acadoVariables.ubValues[143] - acadoVariables.u[143];
acadoWorkspace.ub[144] = acadoVariables.ubValues[144] - acadoVariables.u[144];
acadoWorkspace.ub[145] = acadoVariables.ubValues[145] - acadoVariables.u[145];
acadoWorkspace.ub[146] = acadoVariables.ubValues[146] - acadoVariables.u[146];
acadoWorkspace.ub[147] = acadoVariables.ubValues[147] - acadoVariables.u[147];
acadoWorkspace.ub[148] = acadoVariables.ubValues[148] - acadoVariables.u[148];
acadoWorkspace.ub[149] = acadoVariables.ubValues[149] - acadoVariables.u[149];
acadoWorkspace.ub[150] = acadoVariables.ubValues[150] - acadoVariables.u[150];
acadoWorkspace.ub[151] = acadoVariables.ubValues[151] - acadoVariables.u[151];
acadoWorkspace.ub[152] = acadoVariables.ubValues[152] - acadoVariables.u[152];
acadoWorkspace.ub[153] = acadoVariables.ubValues[153] - acadoVariables.u[153];
acadoWorkspace.ub[154] = acadoVariables.ubValues[154] - acadoVariables.u[154];
acadoWorkspace.ub[155] = acadoVariables.ubValues[155] - acadoVariables.u[155];
acadoWorkspace.ub[156] = acadoVariables.ubValues[156] - acadoVariables.u[156];
acadoWorkspace.ub[157] = acadoVariables.ubValues[157] - acadoVariables.u[157];
acadoWorkspace.ub[158] = acadoVariables.ubValues[158] - acadoVariables.u[158];
acadoWorkspace.ub[159] = acadoVariables.ubValues[159] - acadoVariables.u[159];
acadoWorkspace.ub[160] = acadoVariables.ubValues[160] - acadoVariables.u[160];
acadoWorkspace.ub[161] = acadoVariables.ubValues[161] - acadoVariables.u[161];
acadoWorkspace.ub[162] = acadoVariables.ubValues[162] - acadoVariables.u[162];
acadoWorkspace.ub[163] = acadoVariables.ubValues[163] - acadoVariables.u[163];
acadoWorkspace.ub[164] = acadoVariables.ubValues[164] - acadoVariables.u[164];
acadoWorkspace.ub[165] = acadoVariables.ubValues[165] - acadoVariables.u[165];
acadoWorkspace.ub[166] = acadoVariables.ubValues[166] - acadoVariables.u[166];
acadoWorkspace.ub[167] = acadoVariables.ubValues[167] - acadoVariables.u[167];
acadoWorkspace.ub[168] = acadoVariables.ubValues[168] - acadoVariables.u[168];
acadoWorkspace.ub[169] = acadoVariables.ubValues[169] - acadoVariables.u[169];
acadoWorkspace.ub[170] = acadoVariables.ubValues[170] - acadoVariables.u[170];
acadoWorkspace.ub[171] = acadoVariables.ubValues[171] - acadoVariables.u[171];
acadoWorkspace.ub[172] = acadoVariables.ubValues[172] - acadoVariables.u[172];
acadoWorkspace.ub[173] = acadoVariables.ubValues[173] - acadoVariables.u[173];
acadoWorkspace.ub[174] = acadoVariables.ubValues[174] - acadoVariables.u[174];
acadoWorkspace.ub[175] = acadoVariables.ubValues[175] - acadoVariables.u[175];
acadoWorkspace.ub[176] = acadoVariables.ubValues[176] - acadoVariables.u[176];
acadoWorkspace.ub[177] = acadoVariables.ubValues[177] - acadoVariables.u[177];
acadoWorkspace.ub[178] = acadoVariables.ubValues[178] - acadoVariables.u[178];
acadoWorkspace.ub[179] = acadoVariables.ubValues[179] - acadoVariables.u[179];
acadoWorkspace.ub[180] = acadoVariables.ubValues[180] - acadoVariables.u[180];
acadoWorkspace.ub[181] = acadoVariables.ubValues[181] - acadoVariables.u[181];
acadoWorkspace.ub[182] = acadoVariables.ubValues[182] - acadoVariables.u[182];
acadoWorkspace.ub[183] = acadoVariables.ubValues[183] - acadoVariables.u[183];
acadoWorkspace.ub[184] = acadoVariables.ubValues[184] - acadoVariables.u[184];
acadoWorkspace.ub[185] = acadoVariables.ubValues[185] - acadoVariables.u[185];
acadoWorkspace.ub[186] = acadoVariables.ubValues[186] - acadoVariables.u[186];
acadoWorkspace.ub[187] = acadoVariables.ubValues[187] - acadoVariables.u[187];
acadoWorkspace.ub[188] = acadoVariables.ubValues[188] - acadoVariables.u[188];
acadoWorkspace.ub[189] = acadoVariables.ubValues[189] - acadoVariables.u[189];
acadoWorkspace.ub[190] = acadoVariables.ubValues[190] - acadoVariables.u[190];
acadoWorkspace.ub[191] = acadoVariables.ubValues[191] - acadoVariables.u[191];
acadoWorkspace.ub[192] = acadoVariables.ubValues[192] - acadoVariables.u[192];
acadoWorkspace.ub[193] = acadoVariables.ubValues[193] - acadoVariables.u[193];
acadoWorkspace.ub[194] = acadoVariables.ubValues[194] - acadoVariables.u[194];
acadoWorkspace.ub[195] = acadoVariables.ubValues[195] - acadoVariables.u[195];
acadoWorkspace.ub[196] = acadoVariables.ubValues[196] - acadoVariables.u[196];
acadoWorkspace.ub[197] = acadoVariables.ubValues[197] - acadoVariables.u[197];
acadoWorkspace.ub[198] = acadoVariables.ubValues[198] - acadoVariables.u[198];
acadoWorkspace.ub[199] = acadoVariables.ubValues[199] - acadoVariables.u[199];

}

void acado_expand(  )
{
int lRun1;
for (lRun1 = 0; lRun1 < 200; ++lRun1)
acadoVariables.u[lRun1] += acadoWorkspace.x[lRun1];

acadoWorkspace.sbar[0] = acadoWorkspace.Dx0[0];
acadoWorkspace.sbar[1] = acadoWorkspace.Dx0[1];
acadoWorkspace.sbar[2] = acadoWorkspace.Dx0[2];
acadoWorkspace.sbar[3] = acadoWorkspace.Dx0[3];
for (lRun1 = 0; lRun1 < 800; ++lRun1)
acadoWorkspace.sbar[lRun1 + 4] = acadoWorkspace.d[lRun1];

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
acado_expansionStep( &(acadoWorkspace.evGx[ 320 ]), &(acadoWorkspace.evGu[ 80 ]), &(acadoWorkspace.x[ 20 ]), &(acadoWorkspace.sbar[ 80 ]), &(acadoWorkspace.sbar[ 84 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 336 ]), &(acadoWorkspace.evGu[ 84 ]), &(acadoWorkspace.x[ 21 ]), &(acadoWorkspace.sbar[ 84 ]), &(acadoWorkspace.sbar[ 88 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 352 ]), &(acadoWorkspace.evGu[ 88 ]), &(acadoWorkspace.x[ 22 ]), &(acadoWorkspace.sbar[ 88 ]), &(acadoWorkspace.sbar[ 92 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 368 ]), &(acadoWorkspace.evGu[ 92 ]), &(acadoWorkspace.x[ 23 ]), &(acadoWorkspace.sbar[ 92 ]), &(acadoWorkspace.sbar[ 96 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 384 ]), &(acadoWorkspace.evGu[ 96 ]), &(acadoWorkspace.x[ 24 ]), &(acadoWorkspace.sbar[ 96 ]), &(acadoWorkspace.sbar[ 100 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 400 ]), &(acadoWorkspace.evGu[ 100 ]), &(acadoWorkspace.x[ 25 ]), &(acadoWorkspace.sbar[ 100 ]), &(acadoWorkspace.sbar[ 104 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 416 ]), &(acadoWorkspace.evGu[ 104 ]), &(acadoWorkspace.x[ 26 ]), &(acadoWorkspace.sbar[ 104 ]), &(acadoWorkspace.sbar[ 108 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 432 ]), &(acadoWorkspace.evGu[ 108 ]), &(acadoWorkspace.x[ 27 ]), &(acadoWorkspace.sbar[ 108 ]), &(acadoWorkspace.sbar[ 112 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 448 ]), &(acadoWorkspace.evGu[ 112 ]), &(acadoWorkspace.x[ 28 ]), &(acadoWorkspace.sbar[ 112 ]), &(acadoWorkspace.sbar[ 116 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 464 ]), &(acadoWorkspace.evGu[ 116 ]), &(acadoWorkspace.x[ 29 ]), &(acadoWorkspace.sbar[ 116 ]), &(acadoWorkspace.sbar[ 120 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 480 ]), &(acadoWorkspace.evGu[ 120 ]), &(acadoWorkspace.x[ 30 ]), &(acadoWorkspace.sbar[ 120 ]), &(acadoWorkspace.sbar[ 124 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 496 ]), &(acadoWorkspace.evGu[ 124 ]), &(acadoWorkspace.x[ 31 ]), &(acadoWorkspace.sbar[ 124 ]), &(acadoWorkspace.sbar[ 128 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 512 ]), &(acadoWorkspace.evGu[ 128 ]), &(acadoWorkspace.x[ 32 ]), &(acadoWorkspace.sbar[ 128 ]), &(acadoWorkspace.sbar[ 132 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 528 ]), &(acadoWorkspace.evGu[ 132 ]), &(acadoWorkspace.x[ 33 ]), &(acadoWorkspace.sbar[ 132 ]), &(acadoWorkspace.sbar[ 136 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 544 ]), &(acadoWorkspace.evGu[ 136 ]), &(acadoWorkspace.x[ 34 ]), &(acadoWorkspace.sbar[ 136 ]), &(acadoWorkspace.sbar[ 140 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 560 ]), &(acadoWorkspace.evGu[ 140 ]), &(acadoWorkspace.x[ 35 ]), &(acadoWorkspace.sbar[ 140 ]), &(acadoWorkspace.sbar[ 144 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.evGu[ 144 ]), &(acadoWorkspace.x[ 36 ]), &(acadoWorkspace.sbar[ 144 ]), &(acadoWorkspace.sbar[ 148 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 592 ]), &(acadoWorkspace.evGu[ 148 ]), &(acadoWorkspace.x[ 37 ]), &(acadoWorkspace.sbar[ 148 ]), &(acadoWorkspace.sbar[ 152 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 608 ]), &(acadoWorkspace.evGu[ 152 ]), &(acadoWorkspace.x[ 38 ]), &(acadoWorkspace.sbar[ 152 ]), &(acadoWorkspace.sbar[ 156 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 624 ]), &(acadoWorkspace.evGu[ 156 ]), &(acadoWorkspace.x[ 39 ]), &(acadoWorkspace.sbar[ 156 ]), &(acadoWorkspace.sbar[ 160 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 640 ]), &(acadoWorkspace.evGu[ 160 ]), &(acadoWorkspace.x[ 40 ]), &(acadoWorkspace.sbar[ 160 ]), &(acadoWorkspace.sbar[ 164 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 656 ]), &(acadoWorkspace.evGu[ 164 ]), &(acadoWorkspace.x[ 41 ]), &(acadoWorkspace.sbar[ 164 ]), &(acadoWorkspace.sbar[ 168 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 672 ]), &(acadoWorkspace.evGu[ 168 ]), &(acadoWorkspace.x[ 42 ]), &(acadoWorkspace.sbar[ 168 ]), &(acadoWorkspace.sbar[ 172 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 688 ]), &(acadoWorkspace.evGu[ 172 ]), &(acadoWorkspace.x[ 43 ]), &(acadoWorkspace.sbar[ 172 ]), &(acadoWorkspace.sbar[ 176 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 704 ]), &(acadoWorkspace.evGu[ 176 ]), &(acadoWorkspace.x[ 44 ]), &(acadoWorkspace.sbar[ 176 ]), &(acadoWorkspace.sbar[ 180 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 720 ]), &(acadoWorkspace.evGu[ 180 ]), &(acadoWorkspace.x[ 45 ]), &(acadoWorkspace.sbar[ 180 ]), &(acadoWorkspace.sbar[ 184 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 736 ]), &(acadoWorkspace.evGu[ 184 ]), &(acadoWorkspace.x[ 46 ]), &(acadoWorkspace.sbar[ 184 ]), &(acadoWorkspace.sbar[ 188 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 752 ]), &(acadoWorkspace.evGu[ 188 ]), &(acadoWorkspace.x[ 47 ]), &(acadoWorkspace.sbar[ 188 ]), &(acadoWorkspace.sbar[ 192 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 768 ]), &(acadoWorkspace.evGu[ 192 ]), &(acadoWorkspace.x[ 48 ]), &(acadoWorkspace.sbar[ 192 ]), &(acadoWorkspace.sbar[ 196 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 784 ]), &(acadoWorkspace.evGu[ 196 ]), &(acadoWorkspace.x[ 49 ]), &(acadoWorkspace.sbar[ 196 ]), &(acadoWorkspace.sbar[ 200 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 800 ]), &(acadoWorkspace.evGu[ 200 ]), &(acadoWorkspace.x[ 50 ]), &(acadoWorkspace.sbar[ 200 ]), &(acadoWorkspace.sbar[ 204 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 816 ]), &(acadoWorkspace.evGu[ 204 ]), &(acadoWorkspace.x[ 51 ]), &(acadoWorkspace.sbar[ 204 ]), &(acadoWorkspace.sbar[ 208 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 832 ]), &(acadoWorkspace.evGu[ 208 ]), &(acadoWorkspace.x[ 52 ]), &(acadoWorkspace.sbar[ 208 ]), &(acadoWorkspace.sbar[ 212 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 848 ]), &(acadoWorkspace.evGu[ 212 ]), &(acadoWorkspace.x[ 53 ]), &(acadoWorkspace.sbar[ 212 ]), &(acadoWorkspace.sbar[ 216 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 864 ]), &(acadoWorkspace.evGu[ 216 ]), &(acadoWorkspace.x[ 54 ]), &(acadoWorkspace.sbar[ 216 ]), &(acadoWorkspace.sbar[ 220 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 880 ]), &(acadoWorkspace.evGu[ 220 ]), &(acadoWorkspace.x[ 55 ]), &(acadoWorkspace.sbar[ 220 ]), &(acadoWorkspace.sbar[ 224 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 896 ]), &(acadoWorkspace.evGu[ 224 ]), &(acadoWorkspace.x[ 56 ]), &(acadoWorkspace.sbar[ 224 ]), &(acadoWorkspace.sbar[ 228 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 912 ]), &(acadoWorkspace.evGu[ 228 ]), &(acadoWorkspace.x[ 57 ]), &(acadoWorkspace.sbar[ 228 ]), &(acadoWorkspace.sbar[ 232 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 928 ]), &(acadoWorkspace.evGu[ 232 ]), &(acadoWorkspace.x[ 58 ]), &(acadoWorkspace.sbar[ 232 ]), &(acadoWorkspace.sbar[ 236 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 944 ]), &(acadoWorkspace.evGu[ 236 ]), &(acadoWorkspace.x[ 59 ]), &(acadoWorkspace.sbar[ 236 ]), &(acadoWorkspace.sbar[ 240 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 960 ]), &(acadoWorkspace.evGu[ 240 ]), &(acadoWorkspace.x[ 60 ]), &(acadoWorkspace.sbar[ 240 ]), &(acadoWorkspace.sbar[ 244 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 976 ]), &(acadoWorkspace.evGu[ 244 ]), &(acadoWorkspace.x[ 61 ]), &(acadoWorkspace.sbar[ 244 ]), &(acadoWorkspace.sbar[ 248 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 992 ]), &(acadoWorkspace.evGu[ 248 ]), &(acadoWorkspace.x[ 62 ]), &(acadoWorkspace.sbar[ 248 ]), &(acadoWorkspace.sbar[ 252 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1008 ]), &(acadoWorkspace.evGu[ 252 ]), &(acadoWorkspace.x[ 63 ]), &(acadoWorkspace.sbar[ 252 ]), &(acadoWorkspace.sbar[ 256 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1024 ]), &(acadoWorkspace.evGu[ 256 ]), &(acadoWorkspace.x[ 64 ]), &(acadoWorkspace.sbar[ 256 ]), &(acadoWorkspace.sbar[ 260 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1040 ]), &(acadoWorkspace.evGu[ 260 ]), &(acadoWorkspace.x[ 65 ]), &(acadoWorkspace.sbar[ 260 ]), &(acadoWorkspace.sbar[ 264 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1056 ]), &(acadoWorkspace.evGu[ 264 ]), &(acadoWorkspace.x[ 66 ]), &(acadoWorkspace.sbar[ 264 ]), &(acadoWorkspace.sbar[ 268 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1072 ]), &(acadoWorkspace.evGu[ 268 ]), &(acadoWorkspace.x[ 67 ]), &(acadoWorkspace.sbar[ 268 ]), &(acadoWorkspace.sbar[ 272 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1088 ]), &(acadoWorkspace.evGu[ 272 ]), &(acadoWorkspace.x[ 68 ]), &(acadoWorkspace.sbar[ 272 ]), &(acadoWorkspace.sbar[ 276 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1104 ]), &(acadoWorkspace.evGu[ 276 ]), &(acadoWorkspace.x[ 69 ]), &(acadoWorkspace.sbar[ 276 ]), &(acadoWorkspace.sbar[ 280 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1120 ]), &(acadoWorkspace.evGu[ 280 ]), &(acadoWorkspace.x[ 70 ]), &(acadoWorkspace.sbar[ 280 ]), &(acadoWorkspace.sbar[ 284 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1136 ]), &(acadoWorkspace.evGu[ 284 ]), &(acadoWorkspace.x[ 71 ]), &(acadoWorkspace.sbar[ 284 ]), &(acadoWorkspace.sbar[ 288 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1152 ]), &(acadoWorkspace.evGu[ 288 ]), &(acadoWorkspace.x[ 72 ]), &(acadoWorkspace.sbar[ 288 ]), &(acadoWorkspace.sbar[ 292 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1168 ]), &(acadoWorkspace.evGu[ 292 ]), &(acadoWorkspace.x[ 73 ]), &(acadoWorkspace.sbar[ 292 ]), &(acadoWorkspace.sbar[ 296 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1184 ]), &(acadoWorkspace.evGu[ 296 ]), &(acadoWorkspace.x[ 74 ]), &(acadoWorkspace.sbar[ 296 ]), &(acadoWorkspace.sbar[ 300 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1200 ]), &(acadoWorkspace.evGu[ 300 ]), &(acadoWorkspace.x[ 75 ]), &(acadoWorkspace.sbar[ 300 ]), &(acadoWorkspace.sbar[ 304 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1216 ]), &(acadoWorkspace.evGu[ 304 ]), &(acadoWorkspace.x[ 76 ]), &(acadoWorkspace.sbar[ 304 ]), &(acadoWorkspace.sbar[ 308 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1232 ]), &(acadoWorkspace.evGu[ 308 ]), &(acadoWorkspace.x[ 77 ]), &(acadoWorkspace.sbar[ 308 ]), &(acadoWorkspace.sbar[ 312 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1248 ]), &(acadoWorkspace.evGu[ 312 ]), &(acadoWorkspace.x[ 78 ]), &(acadoWorkspace.sbar[ 312 ]), &(acadoWorkspace.sbar[ 316 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1264 ]), &(acadoWorkspace.evGu[ 316 ]), &(acadoWorkspace.x[ 79 ]), &(acadoWorkspace.sbar[ 316 ]), &(acadoWorkspace.sbar[ 320 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1280 ]), &(acadoWorkspace.evGu[ 320 ]), &(acadoWorkspace.x[ 80 ]), &(acadoWorkspace.sbar[ 320 ]), &(acadoWorkspace.sbar[ 324 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1296 ]), &(acadoWorkspace.evGu[ 324 ]), &(acadoWorkspace.x[ 81 ]), &(acadoWorkspace.sbar[ 324 ]), &(acadoWorkspace.sbar[ 328 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1312 ]), &(acadoWorkspace.evGu[ 328 ]), &(acadoWorkspace.x[ 82 ]), &(acadoWorkspace.sbar[ 328 ]), &(acadoWorkspace.sbar[ 332 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1328 ]), &(acadoWorkspace.evGu[ 332 ]), &(acadoWorkspace.x[ 83 ]), &(acadoWorkspace.sbar[ 332 ]), &(acadoWorkspace.sbar[ 336 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1344 ]), &(acadoWorkspace.evGu[ 336 ]), &(acadoWorkspace.x[ 84 ]), &(acadoWorkspace.sbar[ 336 ]), &(acadoWorkspace.sbar[ 340 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1360 ]), &(acadoWorkspace.evGu[ 340 ]), &(acadoWorkspace.x[ 85 ]), &(acadoWorkspace.sbar[ 340 ]), &(acadoWorkspace.sbar[ 344 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1376 ]), &(acadoWorkspace.evGu[ 344 ]), &(acadoWorkspace.x[ 86 ]), &(acadoWorkspace.sbar[ 344 ]), &(acadoWorkspace.sbar[ 348 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1392 ]), &(acadoWorkspace.evGu[ 348 ]), &(acadoWorkspace.x[ 87 ]), &(acadoWorkspace.sbar[ 348 ]), &(acadoWorkspace.sbar[ 352 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1408 ]), &(acadoWorkspace.evGu[ 352 ]), &(acadoWorkspace.x[ 88 ]), &(acadoWorkspace.sbar[ 352 ]), &(acadoWorkspace.sbar[ 356 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1424 ]), &(acadoWorkspace.evGu[ 356 ]), &(acadoWorkspace.x[ 89 ]), &(acadoWorkspace.sbar[ 356 ]), &(acadoWorkspace.sbar[ 360 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1440 ]), &(acadoWorkspace.evGu[ 360 ]), &(acadoWorkspace.x[ 90 ]), &(acadoWorkspace.sbar[ 360 ]), &(acadoWorkspace.sbar[ 364 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1456 ]), &(acadoWorkspace.evGu[ 364 ]), &(acadoWorkspace.x[ 91 ]), &(acadoWorkspace.sbar[ 364 ]), &(acadoWorkspace.sbar[ 368 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1472 ]), &(acadoWorkspace.evGu[ 368 ]), &(acadoWorkspace.x[ 92 ]), &(acadoWorkspace.sbar[ 368 ]), &(acadoWorkspace.sbar[ 372 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1488 ]), &(acadoWorkspace.evGu[ 372 ]), &(acadoWorkspace.x[ 93 ]), &(acadoWorkspace.sbar[ 372 ]), &(acadoWorkspace.sbar[ 376 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1504 ]), &(acadoWorkspace.evGu[ 376 ]), &(acadoWorkspace.x[ 94 ]), &(acadoWorkspace.sbar[ 376 ]), &(acadoWorkspace.sbar[ 380 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1520 ]), &(acadoWorkspace.evGu[ 380 ]), &(acadoWorkspace.x[ 95 ]), &(acadoWorkspace.sbar[ 380 ]), &(acadoWorkspace.sbar[ 384 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1536 ]), &(acadoWorkspace.evGu[ 384 ]), &(acadoWorkspace.x[ 96 ]), &(acadoWorkspace.sbar[ 384 ]), &(acadoWorkspace.sbar[ 388 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1552 ]), &(acadoWorkspace.evGu[ 388 ]), &(acadoWorkspace.x[ 97 ]), &(acadoWorkspace.sbar[ 388 ]), &(acadoWorkspace.sbar[ 392 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1568 ]), &(acadoWorkspace.evGu[ 392 ]), &(acadoWorkspace.x[ 98 ]), &(acadoWorkspace.sbar[ 392 ]), &(acadoWorkspace.sbar[ 396 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1584 ]), &(acadoWorkspace.evGu[ 396 ]), &(acadoWorkspace.x[ 99 ]), &(acadoWorkspace.sbar[ 396 ]), &(acadoWorkspace.sbar[ 400 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1600 ]), &(acadoWorkspace.evGu[ 400 ]), &(acadoWorkspace.x[ 100 ]), &(acadoWorkspace.sbar[ 400 ]), &(acadoWorkspace.sbar[ 404 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1616 ]), &(acadoWorkspace.evGu[ 404 ]), &(acadoWorkspace.x[ 101 ]), &(acadoWorkspace.sbar[ 404 ]), &(acadoWorkspace.sbar[ 408 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1632 ]), &(acadoWorkspace.evGu[ 408 ]), &(acadoWorkspace.x[ 102 ]), &(acadoWorkspace.sbar[ 408 ]), &(acadoWorkspace.sbar[ 412 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1648 ]), &(acadoWorkspace.evGu[ 412 ]), &(acadoWorkspace.x[ 103 ]), &(acadoWorkspace.sbar[ 412 ]), &(acadoWorkspace.sbar[ 416 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1664 ]), &(acadoWorkspace.evGu[ 416 ]), &(acadoWorkspace.x[ 104 ]), &(acadoWorkspace.sbar[ 416 ]), &(acadoWorkspace.sbar[ 420 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1680 ]), &(acadoWorkspace.evGu[ 420 ]), &(acadoWorkspace.x[ 105 ]), &(acadoWorkspace.sbar[ 420 ]), &(acadoWorkspace.sbar[ 424 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1696 ]), &(acadoWorkspace.evGu[ 424 ]), &(acadoWorkspace.x[ 106 ]), &(acadoWorkspace.sbar[ 424 ]), &(acadoWorkspace.sbar[ 428 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1712 ]), &(acadoWorkspace.evGu[ 428 ]), &(acadoWorkspace.x[ 107 ]), &(acadoWorkspace.sbar[ 428 ]), &(acadoWorkspace.sbar[ 432 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1728 ]), &(acadoWorkspace.evGu[ 432 ]), &(acadoWorkspace.x[ 108 ]), &(acadoWorkspace.sbar[ 432 ]), &(acadoWorkspace.sbar[ 436 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1744 ]), &(acadoWorkspace.evGu[ 436 ]), &(acadoWorkspace.x[ 109 ]), &(acadoWorkspace.sbar[ 436 ]), &(acadoWorkspace.sbar[ 440 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1760 ]), &(acadoWorkspace.evGu[ 440 ]), &(acadoWorkspace.x[ 110 ]), &(acadoWorkspace.sbar[ 440 ]), &(acadoWorkspace.sbar[ 444 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1776 ]), &(acadoWorkspace.evGu[ 444 ]), &(acadoWorkspace.x[ 111 ]), &(acadoWorkspace.sbar[ 444 ]), &(acadoWorkspace.sbar[ 448 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1792 ]), &(acadoWorkspace.evGu[ 448 ]), &(acadoWorkspace.x[ 112 ]), &(acadoWorkspace.sbar[ 448 ]), &(acadoWorkspace.sbar[ 452 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1808 ]), &(acadoWorkspace.evGu[ 452 ]), &(acadoWorkspace.x[ 113 ]), &(acadoWorkspace.sbar[ 452 ]), &(acadoWorkspace.sbar[ 456 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1824 ]), &(acadoWorkspace.evGu[ 456 ]), &(acadoWorkspace.x[ 114 ]), &(acadoWorkspace.sbar[ 456 ]), &(acadoWorkspace.sbar[ 460 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1840 ]), &(acadoWorkspace.evGu[ 460 ]), &(acadoWorkspace.x[ 115 ]), &(acadoWorkspace.sbar[ 460 ]), &(acadoWorkspace.sbar[ 464 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1856 ]), &(acadoWorkspace.evGu[ 464 ]), &(acadoWorkspace.x[ 116 ]), &(acadoWorkspace.sbar[ 464 ]), &(acadoWorkspace.sbar[ 468 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1872 ]), &(acadoWorkspace.evGu[ 468 ]), &(acadoWorkspace.x[ 117 ]), &(acadoWorkspace.sbar[ 468 ]), &(acadoWorkspace.sbar[ 472 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1888 ]), &(acadoWorkspace.evGu[ 472 ]), &(acadoWorkspace.x[ 118 ]), &(acadoWorkspace.sbar[ 472 ]), &(acadoWorkspace.sbar[ 476 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1904 ]), &(acadoWorkspace.evGu[ 476 ]), &(acadoWorkspace.x[ 119 ]), &(acadoWorkspace.sbar[ 476 ]), &(acadoWorkspace.sbar[ 480 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1920 ]), &(acadoWorkspace.evGu[ 480 ]), &(acadoWorkspace.x[ 120 ]), &(acadoWorkspace.sbar[ 480 ]), &(acadoWorkspace.sbar[ 484 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1936 ]), &(acadoWorkspace.evGu[ 484 ]), &(acadoWorkspace.x[ 121 ]), &(acadoWorkspace.sbar[ 484 ]), &(acadoWorkspace.sbar[ 488 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1952 ]), &(acadoWorkspace.evGu[ 488 ]), &(acadoWorkspace.x[ 122 ]), &(acadoWorkspace.sbar[ 488 ]), &(acadoWorkspace.sbar[ 492 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1968 ]), &(acadoWorkspace.evGu[ 492 ]), &(acadoWorkspace.x[ 123 ]), &(acadoWorkspace.sbar[ 492 ]), &(acadoWorkspace.sbar[ 496 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 1984 ]), &(acadoWorkspace.evGu[ 496 ]), &(acadoWorkspace.x[ 124 ]), &(acadoWorkspace.sbar[ 496 ]), &(acadoWorkspace.sbar[ 500 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2000 ]), &(acadoWorkspace.evGu[ 500 ]), &(acadoWorkspace.x[ 125 ]), &(acadoWorkspace.sbar[ 500 ]), &(acadoWorkspace.sbar[ 504 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2016 ]), &(acadoWorkspace.evGu[ 504 ]), &(acadoWorkspace.x[ 126 ]), &(acadoWorkspace.sbar[ 504 ]), &(acadoWorkspace.sbar[ 508 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2032 ]), &(acadoWorkspace.evGu[ 508 ]), &(acadoWorkspace.x[ 127 ]), &(acadoWorkspace.sbar[ 508 ]), &(acadoWorkspace.sbar[ 512 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2048 ]), &(acadoWorkspace.evGu[ 512 ]), &(acadoWorkspace.x[ 128 ]), &(acadoWorkspace.sbar[ 512 ]), &(acadoWorkspace.sbar[ 516 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2064 ]), &(acadoWorkspace.evGu[ 516 ]), &(acadoWorkspace.x[ 129 ]), &(acadoWorkspace.sbar[ 516 ]), &(acadoWorkspace.sbar[ 520 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2080 ]), &(acadoWorkspace.evGu[ 520 ]), &(acadoWorkspace.x[ 130 ]), &(acadoWorkspace.sbar[ 520 ]), &(acadoWorkspace.sbar[ 524 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2096 ]), &(acadoWorkspace.evGu[ 524 ]), &(acadoWorkspace.x[ 131 ]), &(acadoWorkspace.sbar[ 524 ]), &(acadoWorkspace.sbar[ 528 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2112 ]), &(acadoWorkspace.evGu[ 528 ]), &(acadoWorkspace.x[ 132 ]), &(acadoWorkspace.sbar[ 528 ]), &(acadoWorkspace.sbar[ 532 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2128 ]), &(acadoWorkspace.evGu[ 532 ]), &(acadoWorkspace.x[ 133 ]), &(acadoWorkspace.sbar[ 532 ]), &(acadoWorkspace.sbar[ 536 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2144 ]), &(acadoWorkspace.evGu[ 536 ]), &(acadoWorkspace.x[ 134 ]), &(acadoWorkspace.sbar[ 536 ]), &(acadoWorkspace.sbar[ 540 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2160 ]), &(acadoWorkspace.evGu[ 540 ]), &(acadoWorkspace.x[ 135 ]), &(acadoWorkspace.sbar[ 540 ]), &(acadoWorkspace.sbar[ 544 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2176 ]), &(acadoWorkspace.evGu[ 544 ]), &(acadoWorkspace.x[ 136 ]), &(acadoWorkspace.sbar[ 544 ]), &(acadoWorkspace.sbar[ 548 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2192 ]), &(acadoWorkspace.evGu[ 548 ]), &(acadoWorkspace.x[ 137 ]), &(acadoWorkspace.sbar[ 548 ]), &(acadoWorkspace.sbar[ 552 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2208 ]), &(acadoWorkspace.evGu[ 552 ]), &(acadoWorkspace.x[ 138 ]), &(acadoWorkspace.sbar[ 552 ]), &(acadoWorkspace.sbar[ 556 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2224 ]), &(acadoWorkspace.evGu[ 556 ]), &(acadoWorkspace.x[ 139 ]), &(acadoWorkspace.sbar[ 556 ]), &(acadoWorkspace.sbar[ 560 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2240 ]), &(acadoWorkspace.evGu[ 560 ]), &(acadoWorkspace.x[ 140 ]), &(acadoWorkspace.sbar[ 560 ]), &(acadoWorkspace.sbar[ 564 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2256 ]), &(acadoWorkspace.evGu[ 564 ]), &(acadoWorkspace.x[ 141 ]), &(acadoWorkspace.sbar[ 564 ]), &(acadoWorkspace.sbar[ 568 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2272 ]), &(acadoWorkspace.evGu[ 568 ]), &(acadoWorkspace.x[ 142 ]), &(acadoWorkspace.sbar[ 568 ]), &(acadoWorkspace.sbar[ 572 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2288 ]), &(acadoWorkspace.evGu[ 572 ]), &(acadoWorkspace.x[ 143 ]), &(acadoWorkspace.sbar[ 572 ]), &(acadoWorkspace.sbar[ 576 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2304 ]), &(acadoWorkspace.evGu[ 576 ]), &(acadoWorkspace.x[ 144 ]), &(acadoWorkspace.sbar[ 576 ]), &(acadoWorkspace.sbar[ 580 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2320 ]), &(acadoWorkspace.evGu[ 580 ]), &(acadoWorkspace.x[ 145 ]), &(acadoWorkspace.sbar[ 580 ]), &(acadoWorkspace.sbar[ 584 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2336 ]), &(acadoWorkspace.evGu[ 584 ]), &(acadoWorkspace.x[ 146 ]), &(acadoWorkspace.sbar[ 584 ]), &(acadoWorkspace.sbar[ 588 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2352 ]), &(acadoWorkspace.evGu[ 588 ]), &(acadoWorkspace.x[ 147 ]), &(acadoWorkspace.sbar[ 588 ]), &(acadoWorkspace.sbar[ 592 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2368 ]), &(acadoWorkspace.evGu[ 592 ]), &(acadoWorkspace.x[ 148 ]), &(acadoWorkspace.sbar[ 592 ]), &(acadoWorkspace.sbar[ 596 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2384 ]), &(acadoWorkspace.evGu[ 596 ]), &(acadoWorkspace.x[ 149 ]), &(acadoWorkspace.sbar[ 596 ]), &(acadoWorkspace.sbar[ 600 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2400 ]), &(acadoWorkspace.evGu[ 600 ]), &(acadoWorkspace.x[ 150 ]), &(acadoWorkspace.sbar[ 600 ]), &(acadoWorkspace.sbar[ 604 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2416 ]), &(acadoWorkspace.evGu[ 604 ]), &(acadoWorkspace.x[ 151 ]), &(acadoWorkspace.sbar[ 604 ]), &(acadoWorkspace.sbar[ 608 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2432 ]), &(acadoWorkspace.evGu[ 608 ]), &(acadoWorkspace.x[ 152 ]), &(acadoWorkspace.sbar[ 608 ]), &(acadoWorkspace.sbar[ 612 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2448 ]), &(acadoWorkspace.evGu[ 612 ]), &(acadoWorkspace.x[ 153 ]), &(acadoWorkspace.sbar[ 612 ]), &(acadoWorkspace.sbar[ 616 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2464 ]), &(acadoWorkspace.evGu[ 616 ]), &(acadoWorkspace.x[ 154 ]), &(acadoWorkspace.sbar[ 616 ]), &(acadoWorkspace.sbar[ 620 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2480 ]), &(acadoWorkspace.evGu[ 620 ]), &(acadoWorkspace.x[ 155 ]), &(acadoWorkspace.sbar[ 620 ]), &(acadoWorkspace.sbar[ 624 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2496 ]), &(acadoWorkspace.evGu[ 624 ]), &(acadoWorkspace.x[ 156 ]), &(acadoWorkspace.sbar[ 624 ]), &(acadoWorkspace.sbar[ 628 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2512 ]), &(acadoWorkspace.evGu[ 628 ]), &(acadoWorkspace.x[ 157 ]), &(acadoWorkspace.sbar[ 628 ]), &(acadoWorkspace.sbar[ 632 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2528 ]), &(acadoWorkspace.evGu[ 632 ]), &(acadoWorkspace.x[ 158 ]), &(acadoWorkspace.sbar[ 632 ]), &(acadoWorkspace.sbar[ 636 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2544 ]), &(acadoWorkspace.evGu[ 636 ]), &(acadoWorkspace.x[ 159 ]), &(acadoWorkspace.sbar[ 636 ]), &(acadoWorkspace.sbar[ 640 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2560 ]), &(acadoWorkspace.evGu[ 640 ]), &(acadoWorkspace.x[ 160 ]), &(acadoWorkspace.sbar[ 640 ]), &(acadoWorkspace.sbar[ 644 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2576 ]), &(acadoWorkspace.evGu[ 644 ]), &(acadoWorkspace.x[ 161 ]), &(acadoWorkspace.sbar[ 644 ]), &(acadoWorkspace.sbar[ 648 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2592 ]), &(acadoWorkspace.evGu[ 648 ]), &(acadoWorkspace.x[ 162 ]), &(acadoWorkspace.sbar[ 648 ]), &(acadoWorkspace.sbar[ 652 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2608 ]), &(acadoWorkspace.evGu[ 652 ]), &(acadoWorkspace.x[ 163 ]), &(acadoWorkspace.sbar[ 652 ]), &(acadoWorkspace.sbar[ 656 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2624 ]), &(acadoWorkspace.evGu[ 656 ]), &(acadoWorkspace.x[ 164 ]), &(acadoWorkspace.sbar[ 656 ]), &(acadoWorkspace.sbar[ 660 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2640 ]), &(acadoWorkspace.evGu[ 660 ]), &(acadoWorkspace.x[ 165 ]), &(acadoWorkspace.sbar[ 660 ]), &(acadoWorkspace.sbar[ 664 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2656 ]), &(acadoWorkspace.evGu[ 664 ]), &(acadoWorkspace.x[ 166 ]), &(acadoWorkspace.sbar[ 664 ]), &(acadoWorkspace.sbar[ 668 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2672 ]), &(acadoWorkspace.evGu[ 668 ]), &(acadoWorkspace.x[ 167 ]), &(acadoWorkspace.sbar[ 668 ]), &(acadoWorkspace.sbar[ 672 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2688 ]), &(acadoWorkspace.evGu[ 672 ]), &(acadoWorkspace.x[ 168 ]), &(acadoWorkspace.sbar[ 672 ]), &(acadoWorkspace.sbar[ 676 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2704 ]), &(acadoWorkspace.evGu[ 676 ]), &(acadoWorkspace.x[ 169 ]), &(acadoWorkspace.sbar[ 676 ]), &(acadoWorkspace.sbar[ 680 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2720 ]), &(acadoWorkspace.evGu[ 680 ]), &(acadoWorkspace.x[ 170 ]), &(acadoWorkspace.sbar[ 680 ]), &(acadoWorkspace.sbar[ 684 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2736 ]), &(acadoWorkspace.evGu[ 684 ]), &(acadoWorkspace.x[ 171 ]), &(acadoWorkspace.sbar[ 684 ]), &(acadoWorkspace.sbar[ 688 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2752 ]), &(acadoWorkspace.evGu[ 688 ]), &(acadoWorkspace.x[ 172 ]), &(acadoWorkspace.sbar[ 688 ]), &(acadoWorkspace.sbar[ 692 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2768 ]), &(acadoWorkspace.evGu[ 692 ]), &(acadoWorkspace.x[ 173 ]), &(acadoWorkspace.sbar[ 692 ]), &(acadoWorkspace.sbar[ 696 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2784 ]), &(acadoWorkspace.evGu[ 696 ]), &(acadoWorkspace.x[ 174 ]), &(acadoWorkspace.sbar[ 696 ]), &(acadoWorkspace.sbar[ 700 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2800 ]), &(acadoWorkspace.evGu[ 700 ]), &(acadoWorkspace.x[ 175 ]), &(acadoWorkspace.sbar[ 700 ]), &(acadoWorkspace.sbar[ 704 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2816 ]), &(acadoWorkspace.evGu[ 704 ]), &(acadoWorkspace.x[ 176 ]), &(acadoWorkspace.sbar[ 704 ]), &(acadoWorkspace.sbar[ 708 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2832 ]), &(acadoWorkspace.evGu[ 708 ]), &(acadoWorkspace.x[ 177 ]), &(acadoWorkspace.sbar[ 708 ]), &(acadoWorkspace.sbar[ 712 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2848 ]), &(acadoWorkspace.evGu[ 712 ]), &(acadoWorkspace.x[ 178 ]), &(acadoWorkspace.sbar[ 712 ]), &(acadoWorkspace.sbar[ 716 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2864 ]), &(acadoWorkspace.evGu[ 716 ]), &(acadoWorkspace.x[ 179 ]), &(acadoWorkspace.sbar[ 716 ]), &(acadoWorkspace.sbar[ 720 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2880 ]), &(acadoWorkspace.evGu[ 720 ]), &(acadoWorkspace.x[ 180 ]), &(acadoWorkspace.sbar[ 720 ]), &(acadoWorkspace.sbar[ 724 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2896 ]), &(acadoWorkspace.evGu[ 724 ]), &(acadoWorkspace.x[ 181 ]), &(acadoWorkspace.sbar[ 724 ]), &(acadoWorkspace.sbar[ 728 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2912 ]), &(acadoWorkspace.evGu[ 728 ]), &(acadoWorkspace.x[ 182 ]), &(acadoWorkspace.sbar[ 728 ]), &(acadoWorkspace.sbar[ 732 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2928 ]), &(acadoWorkspace.evGu[ 732 ]), &(acadoWorkspace.x[ 183 ]), &(acadoWorkspace.sbar[ 732 ]), &(acadoWorkspace.sbar[ 736 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2944 ]), &(acadoWorkspace.evGu[ 736 ]), &(acadoWorkspace.x[ 184 ]), &(acadoWorkspace.sbar[ 736 ]), &(acadoWorkspace.sbar[ 740 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2960 ]), &(acadoWorkspace.evGu[ 740 ]), &(acadoWorkspace.x[ 185 ]), &(acadoWorkspace.sbar[ 740 ]), &(acadoWorkspace.sbar[ 744 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2976 ]), &(acadoWorkspace.evGu[ 744 ]), &(acadoWorkspace.x[ 186 ]), &(acadoWorkspace.sbar[ 744 ]), &(acadoWorkspace.sbar[ 748 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 2992 ]), &(acadoWorkspace.evGu[ 748 ]), &(acadoWorkspace.x[ 187 ]), &(acadoWorkspace.sbar[ 748 ]), &(acadoWorkspace.sbar[ 752 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 3008 ]), &(acadoWorkspace.evGu[ 752 ]), &(acadoWorkspace.x[ 188 ]), &(acadoWorkspace.sbar[ 752 ]), &(acadoWorkspace.sbar[ 756 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 3024 ]), &(acadoWorkspace.evGu[ 756 ]), &(acadoWorkspace.x[ 189 ]), &(acadoWorkspace.sbar[ 756 ]), &(acadoWorkspace.sbar[ 760 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 3040 ]), &(acadoWorkspace.evGu[ 760 ]), &(acadoWorkspace.x[ 190 ]), &(acadoWorkspace.sbar[ 760 ]), &(acadoWorkspace.sbar[ 764 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 3056 ]), &(acadoWorkspace.evGu[ 764 ]), &(acadoWorkspace.x[ 191 ]), &(acadoWorkspace.sbar[ 764 ]), &(acadoWorkspace.sbar[ 768 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 3072 ]), &(acadoWorkspace.evGu[ 768 ]), &(acadoWorkspace.x[ 192 ]), &(acadoWorkspace.sbar[ 768 ]), &(acadoWorkspace.sbar[ 772 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 3088 ]), &(acadoWorkspace.evGu[ 772 ]), &(acadoWorkspace.x[ 193 ]), &(acadoWorkspace.sbar[ 772 ]), &(acadoWorkspace.sbar[ 776 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 3104 ]), &(acadoWorkspace.evGu[ 776 ]), &(acadoWorkspace.x[ 194 ]), &(acadoWorkspace.sbar[ 776 ]), &(acadoWorkspace.sbar[ 780 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 3120 ]), &(acadoWorkspace.evGu[ 780 ]), &(acadoWorkspace.x[ 195 ]), &(acadoWorkspace.sbar[ 780 ]), &(acadoWorkspace.sbar[ 784 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 3136 ]), &(acadoWorkspace.evGu[ 784 ]), &(acadoWorkspace.x[ 196 ]), &(acadoWorkspace.sbar[ 784 ]), &(acadoWorkspace.sbar[ 788 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 3152 ]), &(acadoWorkspace.evGu[ 788 ]), &(acadoWorkspace.x[ 197 ]), &(acadoWorkspace.sbar[ 788 ]), &(acadoWorkspace.sbar[ 792 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 3168 ]), &(acadoWorkspace.evGu[ 792 ]), &(acadoWorkspace.x[ 198 ]), &(acadoWorkspace.sbar[ 792 ]), &(acadoWorkspace.sbar[ 796 ]) );
acado_expansionStep( &(acadoWorkspace.evGx[ 3184 ]), &(acadoWorkspace.evGu[ 796 ]), &(acadoWorkspace.x[ 199 ]), &(acadoWorkspace.sbar[ 796 ]), &(acadoWorkspace.sbar[ 800 ]) );
for (lRun1 = 0; lRun1 < 804; ++lRun1)
acadoVariables.x[lRun1] += acadoWorkspace.sbar[lRun1];

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
acadoVariables.lbValues[20] = -2.0000000000000000e+00;
acadoVariables.lbValues[21] = -2.0000000000000000e+00;
acadoVariables.lbValues[22] = -2.0000000000000000e+00;
acadoVariables.lbValues[23] = -2.0000000000000000e+00;
acadoVariables.lbValues[24] = -2.0000000000000000e+00;
acadoVariables.lbValues[25] = -2.0000000000000000e+00;
acadoVariables.lbValues[26] = -2.0000000000000000e+00;
acadoVariables.lbValues[27] = -2.0000000000000000e+00;
acadoVariables.lbValues[28] = -2.0000000000000000e+00;
acadoVariables.lbValues[29] = -2.0000000000000000e+00;
acadoVariables.lbValues[30] = -2.0000000000000000e+00;
acadoVariables.lbValues[31] = -2.0000000000000000e+00;
acadoVariables.lbValues[32] = -2.0000000000000000e+00;
acadoVariables.lbValues[33] = -2.0000000000000000e+00;
acadoVariables.lbValues[34] = -2.0000000000000000e+00;
acadoVariables.lbValues[35] = -2.0000000000000000e+00;
acadoVariables.lbValues[36] = -2.0000000000000000e+00;
acadoVariables.lbValues[37] = -2.0000000000000000e+00;
acadoVariables.lbValues[38] = -2.0000000000000000e+00;
acadoVariables.lbValues[39] = -2.0000000000000000e+00;
acadoVariables.lbValues[40] = -2.0000000000000000e+00;
acadoVariables.lbValues[41] = -2.0000000000000000e+00;
acadoVariables.lbValues[42] = -2.0000000000000000e+00;
acadoVariables.lbValues[43] = -2.0000000000000000e+00;
acadoVariables.lbValues[44] = -2.0000000000000000e+00;
acadoVariables.lbValues[45] = -2.0000000000000000e+00;
acadoVariables.lbValues[46] = -2.0000000000000000e+00;
acadoVariables.lbValues[47] = -2.0000000000000000e+00;
acadoVariables.lbValues[48] = -2.0000000000000000e+00;
acadoVariables.lbValues[49] = -2.0000000000000000e+00;
acadoVariables.lbValues[50] = -2.0000000000000000e+00;
acadoVariables.lbValues[51] = -2.0000000000000000e+00;
acadoVariables.lbValues[52] = -2.0000000000000000e+00;
acadoVariables.lbValues[53] = -2.0000000000000000e+00;
acadoVariables.lbValues[54] = -2.0000000000000000e+00;
acadoVariables.lbValues[55] = -2.0000000000000000e+00;
acadoVariables.lbValues[56] = -2.0000000000000000e+00;
acadoVariables.lbValues[57] = -2.0000000000000000e+00;
acadoVariables.lbValues[58] = -2.0000000000000000e+00;
acadoVariables.lbValues[59] = -2.0000000000000000e+00;
acadoVariables.lbValues[60] = -2.0000000000000000e+00;
acadoVariables.lbValues[61] = -2.0000000000000000e+00;
acadoVariables.lbValues[62] = -2.0000000000000000e+00;
acadoVariables.lbValues[63] = -2.0000000000000000e+00;
acadoVariables.lbValues[64] = -2.0000000000000000e+00;
acadoVariables.lbValues[65] = -2.0000000000000000e+00;
acadoVariables.lbValues[66] = -2.0000000000000000e+00;
acadoVariables.lbValues[67] = -2.0000000000000000e+00;
acadoVariables.lbValues[68] = -2.0000000000000000e+00;
acadoVariables.lbValues[69] = -2.0000000000000000e+00;
acadoVariables.lbValues[70] = -2.0000000000000000e+00;
acadoVariables.lbValues[71] = -2.0000000000000000e+00;
acadoVariables.lbValues[72] = -2.0000000000000000e+00;
acadoVariables.lbValues[73] = -2.0000000000000000e+00;
acadoVariables.lbValues[74] = -2.0000000000000000e+00;
acadoVariables.lbValues[75] = -2.0000000000000000e+00;
acadoVariables.lbValues[76] = -2.0000000000000000e+00;
acadoVariables.lbValues[77] = -2.0000000000000000e+00;
acadoVariables.lbValues[78] = -2.0000000000000000e+00;
acadoVariables.lbValues[79] = -2.0000000000000000e+00;
acadoVariables.lbValues[80] = -2.0000000000000000e+00;
acadoVariables.lbValues[81] = -2.0000000000000000e+00;
acadoVariables.lbValues[82] = -2.0000000000000000e+00;
acadoVariables.lbValues[83] = -2.0000000000000000e+00;
acadoVariables.lbValues[84] = -2.0000000000000000e+00;
acadoVariables.lbValues[85] = -2.0000000000000000e+00;
acadoVariables.lbValues[86] = -2.0000000000000000e+00;
acadoVariables.lbValues[87] = -2.0000000000000000e+00;
acadoVariables.lbValues[88] = -2.0000000000000000e+00;
acadoVariables.lbValues[89] = -2.0000000000000000e+00;
acadoVariables.lbValues[90] = -2.0000000000000000e+00;
acadoVariables.lbValues[91] = -2.0000000000000000e+00;
acadoVariables.lbValues[92] = -2.0000000000000000e+00;
acadoVariables.lbValues[93] = -2.0000000000000000e+00;
acadoVariables.lbValues[94] = -2.0000000000000000e+00;
acadoVariables.lbValues[95] = -2.0000000000000000e+00;
acadoVariables.lbValues[96] = -2.0000000000000000e+00;
acadoVariables.lbValues[97] = -2.0000000000000000e+00;
acadoVariables.lbValues[98] = -2.0000000000000000e+00;
acadoVariables.lbValues[99] = -2.0000000000000000e+00;
acadoVariables.lbValues[100] = -2.0000000000000000e+00;
acadoVariables.lbValues[101] = -2.0000000000000000e+00;
acadoVariables.lbValues[102] = -2.0000000000000000e+00;
acadoVariables.lbValues[103] = -2.0000000000000000e+00;
acadoVariables.lbValues[104] = -2.0000000000000000e+00;
acadoVariables.lbValues[105] = -2.0000000000000000e+00;
acadoVariables.lbValues[106] = -2.0000000000000000e+00;
acadoVariables.lbValues[107] = -2.0000000000000000e+00;
acadoVariables.lbValues[108] = -2.0000000000000000e+00;
acadoVariables.lbValues[109] = -2.0000000000000000e+00;
acadoVariables.lbValues[110] = -2.0000000000000000e+00;
acadoVariables.lbValues[111] = -2.0000000000000000e+00;
acadoVariables.lbValues[112] = -2.0000000000000000e+00;
acadoVariables.lbValues[113] = -2.0000000000000000e+00;
acadoVariables.lbValues[114] = -2.0000000000000000e+00;
acadoVariables.lbValues[115] = -2.0000000000000000e+00;
acadoVariables.lbValues[116] = -2.0000000000000000e+00;
acadoVariables.lbValues[117] = -2.0000000000000000e+00;
acadoVariables.lbValues[118] = -2.0000000000000000e+00;
acadoVariables.lbValues[119] = -2.0000000000000000e+00;
acadoVariables.lbValues[120] = -2.0000000000000000e+00;
acadoVariables.lbValues[121] = -2.0000000000000000e+00;
acadoVariables.lbValues[122] = -2.0000000000000000e+00;
acadoVariables.lbValues[123] = -2.0000000000000000e+00;
acadoVariables.lbValues[124] = -2.0000000000000000e+00;
acadoVariables.lbValues[125] = -2.0000000000000000e+00;
acadoVariables.lbValues[126] = -2.0000000000000000e+00;
acadoVariables.lbValues[127] = -2.0000000000000000e+00;
acadoVariables.lbValues[128] = -2.0000000000000000e+00;
acadoVariables.lbValues[129] = -2.0000000000000000e+00;
acadoVariables.lbValues[130] = -2.0000000000000000e+00;
acadoVariables.lbValues[131] = -2.0000000000000000e+00;
acadoVariables.lbValues[132] = -2.0000000000000000e+00;
acadoVariables.lbValues[133] = -2.0000000000000000e+00;
acadoVariables.lbValues[134] = -2.0000000000000000e+00;
acadoVariables.lbValues[135] = -2.0000000000000000e+00;
acadoVariables.lbValues[136] = -2.0000000000000000e+00;
acadoVariables.lbValues[137] = -2.0000000000000000e+00;
acadoVariables.lbValues[138] = -2.0000000000000000e+00;
acadoVariables.lbValues[139] = -2.0000000000000000e+00;
acadoVariables.lbValues[140] = -2.0000000000000000e+00;
acadoVariables.lbValues[141] = -2.0000000000000000e+00;
acadoVariables.lbValues[142] = -2.0000000000000000e+00;
acadoVariables.lbValues[143] = -2.0000000000000000e+00;
acadoVariables.lbValues[144] = -2.0000000000000000e+00;
acadoVariables.lbValues[145] = -2.0000000000000000e+00;
acadoVariables.lbValues[146] = -2.0000000000000000e+00;
acadoVariables.lbValues[147] = -2.0000000000000000e+00;
acadoVariables.lbValues[148] = -2.0000000000000000e+00;
acadoVariables.lbValues[149] = -2.0000000000000000e+00;
acadoVariables.lbValues[150] = -2.0000000000000000e+00;
acadoVariables.lbValues[151] = -2.0000000000000000e+00;
acadoVariables.lbValues[152] = -2.0000000000000000e+00;
acadoVariables.lbValues[153] = -2.0000000000000000e+00;
acadoVariables.lbValues[154] = -2.0000000000000000e+00;
acadoVariables.lbValues[155] = -2.0000000000000000e+00;
acadoVariables.lbValues[156] = -2.0000000000000000e+00;
acadoVariables.lbValues[157] = -2.0000000000000000e+00;
acadoVariables.lbValues[158] = -2.0000000000000000e+00;
acadoVariables.lbValues[159] = -2.0000000000000000e+00;
acadoVariables.lbValues[160] = -2.0000000000000000e+00;
acadoVariables.lbValues[161] = -2.0000000000000000e+00;
acadoVariables.lbValues[162] = -2.0000000000000000e+00;
acadoVariables.lbValues[163] = -2.0000000000000000e+00;
acadoVariables.lbValues[164] = -2.0000000000000000e+00;
acadoVariables.lbValues[165] = -2.0000000000000000e+00;
acadoVariables.lbValues[166] = -2.0000000000000000e+00;
acadoVariables.lbValues[167] = -2.0000000000000000e+00;
acadoVariables.lbValues[168] = -2.0000000000000000e+00;
acadoVariables.lbValues[169] = -2.0000000000000000e+00;
acadoVariables.lbValues[170] = -2.0000000000000000e+00;
acadoVariables.lbValues[171] = -2.0000000000000000e+00;
acadoVariables.lbValues[172] = -2.0000000000000000e+00;
acadoVariables.lbValues[173] = -2.0000000000000000e+00;
acadoVariables.lbValues[174] = -2.0000000000000000e+00;
acadoVariables.lbValues[175] = -2.0000000000000000e+00;
acadoVariables.lbValues[176] = -2.0000000000000000e+00;
acadoVariables.lbValues[177] = -2.0000000000000000e+00;
acadoVariables.lbValues[178] = -2.0000000000000000e+00;
acadoVariables.lbValues[179] = -2.0000000000000000e+00;
acadoVariables.lbValues[180] = -2.0000000000000000e+00;
acadoVariables.lbValues[181] = -2.0000000000000000e+00;
acadoVariables.lbValues[182] = -2.0000000000000000e+00;
acadoVariables.lbValues[183] = -2.0000000000000000e+00;
acadoVariables.lbValues[184] = -2.0000000000000000e+00;
acadoVariables.lbValues[185] = -2.0000000000000000e+00;
acadoVariables.lbValues[186] = -2.0000000000000000e+00;
acadoVariables.lbValues[187] = -2.0000000000000000e+00;
acadoVariables.lbValues[188] = -2.0000000000000000e+00;
acadoVariables.lbValues[189] = -2.0000000000000000e+00;
acadoVariables.lbValues[190] = -2.0000000000000000e+00;
acadoVariables.lbValues[191] = -2.0000000000000000e+00;
acadoVariables.lbValues[192] = -2.0000000000000000e+00;
acadoVariables.lbValues[193] = -2.0000000000000000e+00;
acadoVariables.lbValues[194] = -2.0000000000000000e+00;
acadoVariables.lbValues[195] = -2.0000000000000000e+00;
acadoVariables.lbValues[196] = -2.0000000000000000e+00;
acadoVariables.lbValues[197] = -2.0000000000000000e+00;
acadoVariables.lbValues[198] = -2.0000000000000000e+00;
acadoVariables.lbValues[199] = -2.0000000000000000e+00;
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
acadoVariables.ubValues[20] = 2.0000000000000000e+00;
acadoVariables.ubValues[21] = 2.0000000000000000e+00;
acadoVariables.ubValues[22] = 2.0000000000000000e+00;
acadoVariables.ubValues[23] = 2.0000000000000000e+00;
acadoVariables.ubValues[24] = 2.0000000000000000e+00;
acadoVariables.ubValues[25] = 2.0000000000000000e+00;
acadoVariables.ubValues[26] = 2.0000000000000000e+00;
acadoVariables.ubValues[27] = 2.0000000000000000e+00;
acadoVariables.ubValues[28] = 2.0000000000000000e+00;
acadoVariables.ubValues[29] = 2.0000000000000000e+00;
acadoVariables.ubValues[30] = 2.0000000000000000e+00;
acadoVariables.ubValues[31] = 2.0000000000000000e+00;
acadoVariables.ubValues[32] = 2.0000000000000000e+00;
acadoVariables.ubValues[33] = 2.0000000000000000e+00;
acadoVariables.ubValues[34] = 2.0000000000000000e+00;
acadoVariables.ubValues[35] = 2.0000000000000000e+00;
acadoVariables.ubValues[36] = 2.0000000000000000e+00;
acadoVariables.ubValues[37] = 2.0000000000000000e+00;
acadoVariables.ubValues[38] = 2.0000000000000000e+00;
acadoVariables.ubValues[39] = 2.0000000000000000e+00;
acadoVariables.ubValues[40] = 2.0000000000000000e+00;
acadoVariables.ubValues[41] = 2.0000000000000000e+00;
acadoVariables.ubValues[42] = 2.0000000000000000e+00;
acadoVariables.ubValues[43] = 2.0000000000000000e+00;
acadoVariables.ubValues[44] = 2.0000000000000000e+00;
acadoVariables.ubValues[45] = 2.0000000000000000e+00;
acadoVariables.ubValues[46] = 2.0000000000000000e+00;
acadoVariables.ubValues[47] = 2.0000000000000000e+00;
acadoVariables.ubValues[48] = 2.0000000000000000e+00;
acadoVariables.ubValues[49] = 2.0000000000000000e+00;
acadoVariables.ubValues[50] = 2.0000000000000000e+00;
acadoVariables.ubValues[51] = 2.0000000000000000e+00;
acadoVariables.ubValues[52] = 2.0000000000000000e+00;
acadoVariables.ubValues[53] = 2.0000000000000000e+00;
acadoVariables.ubValues[54] = 2.0000000000000000e+00;
acadoVariables.ubValues[55] = 2.0000000000000000e+00;
acadoVariables.ubValues[56] = 2.0000000000000000e+00;
acadoVariables.ubValues[57] = 2.0000000000000000e+00;
acadoVariables.ubValues[58] = 2.0000000000000000e+00;
acadoVariables.ubValues[59] = 2.0000000000000000e+00;
acadoVariables.ubValues[60] = 2.0000000000000000e+00;
acadoVariables.ubValues[61] = 2.0000000000000000e+00;
acadoVariables.ubValues[62] = 2.0000000000000000e+00;
acadoVariables.ubValues[63] = 2.0000000000000000e+00;
acadoVariables.ubValues[64] = 2.0000000000000000e+00;
acadoVariables.ubValues[65] = 2.0000000000000000e+00;
acadoVariables.ubValues[66] = 2.0000000000000000e+00;
acadoVariables.ubValues[67] = 2.0000000000000000e+00;
acadoVariables.ubValues[68] = 2.0000000000000000e+00;
acadoVariables.ubValues[69] = 2.0000000000000000e+00;
acadoVariables.ubValues[70] = 2.0000000000000000e+00;
acadoVariables.ubValues[71] = 2.0000000000000000e+00;
acadoVariables.ubValues[72] = 2.0000000000000000e+00;
acadoVariables.ubValues[73] = 2.0000000000000000e+00;
acadoVariables.ubValues[74] = 2.0000000000000000e+00;
acadoVariables.ubValues[75] = 2.0000000000000000e+00;
acadoVariables.ubValues[76] = 2.0000000000000000e+00;
acadoVariables.ubValues[77] = 2.0000000000000000e+00;
acadoVariables.ubValues[78] = 2.0000000000000000e+00;
acadoVariables.ubValues[79] = 2.0000000000000000e+00;
acadoVariables.ubValues[80] = 2.0000000000000000e+00;
acadoVariables.ubValues[81] = 2.0000000000000000e+00;
acadoVariables.ubValues[82] = 2.0000000000000000e+00;
acadoVariables.ubValues[83] = 2.0000000000000000e+00;
acadoVariables.ubValues[84] = 2.0000000000000000e+00;
acadoVariables.ubValues[85] = 2.0000000000000000e+00;
acadoVariables.ubValues[86] = 2.0000000000000000e+00;
acadoVariables.ubValues[87] = 2.0000000000000000e+00;
acadoVariables.ubValues[88] = 2.0000000000000000e+00;
acadoVariables.ubValues[89] = 2.0000000000000000e+00;
acadoVariables.ubValues[90] = 2.0000000000000000e+00;
acadoVariables.ubValues[91] = 2.0000000000000000e+00;
acadoVariables.ubValues[92] = 2.0000000000000000e+00;
acadoVariables.ubValues[93] = 2.0000000000000000e+00;
acadoVariables.ubValues[94] = 2.0000000000000000e+00;
acadoVariables.ubValues[95] = 2.0000000000000000e+00;
acadoVariables.ubValues[96] = 2.0000000000000000e+00;
acadoVariables.ubValues[97] = 2.0000000000000000e+00;
acadoVariables.ubValues[98] = 2.0000000000000000e+00;
acadoVariables.ubValues[99] = 2.0000000000000000e+00;
acadoVariables.ubValues[100] = 2.0000000000000000e+00;
acadoVariables.ubValues[101] = 2.0000000000000000e+00;
acadoVariables.ubValues[102] = 2.0000000000000000e+00;
acadoVariables.ubValues[103] = 2.0000000000000000e+00;
acadoVariables.ubValues[104] = 2.0000000000000000e+00;
acadoVariables.ubValues[105] = 2.0000000000000000e+00;
acadoVariables.ubValues[106] = 2.0000000000000000e+00;
acadoVariables.ubValues[107] = 2.0000000000000000e+00;
acadoVariables.ubValues[108] = 2.0000000000000000e+00;
acadoVariables.ubValues[109] = 2.0000000000000000e+00;
acadoVariables.ubValues[110] = 2.0000000000000000e+00;
acadoVariables.ubValues[111] = 2.0000000000000000e+00;
acadoVariables.ubValues[112] = 2.0000000000000000e+00;
acadoVariables.ubValues[113] = 2.0000000000000000e+00;
acadoVariables.ubValues[114] = 2.0000000000000000e+00;
acadoVariables.ubValues[115] = 2.0000000000000000e+00;
acadoVariables.ubValues[116] = 2.0000000000000000e+00;
acadoVariables.ubValues[117] = 2.0000000000000000e+00;
acadoVariables.ubValues[118] = 2.0000000000000000e+00;
acadoVariables.ubValues[119] = 2.0000000000000000e+00;
acadoVariables.ubValues[120] = 2.0000000000000000e+00;
acadoVariables.ubValues[121] = 2.0000000000000000e+00;
acadoVariables.ubValues[122] = 2.0000000000000000e+00;
acadoVariables.ubValues[123] = 2.0000000000000000e+00;
acadoVariables.ubValues[124] = 2.0000000000000000e+00;
acadoVariables.ubValues[125] = 2.0000000000000000e+00;
acadoVariables.ubValues[126] = 2.0000000000000000e+00;
acadoVariables.ubValues[127] = 2.0000000000000000e+00;
acadoVariables.ubValues[128] = 2.0000000000000000e+00;
acadoVariables.ubValues[129] = 2.0000000000000000e+00;
acadoVariables.ubValues[130] = 2.0000000000000000e+00;
acadoVariables.ubValues[131] = 2.0000000000000000e+00;
acadoVariables.ubValues[132] = 2.0000000000000000e+00;
acadoVariables.ubValues[133] = 2.0000000000000000e+00;
acadoVariables.ubValues[134] = 2.0000000000000000e+00;
acadoVariables.ubValues[135] = 2.0000000000000000e+00;
acadoVariables.ubValues[136] = 2.0000000000000000e+00;
acadoVariables.ubValues[137] = 2.0000000000000000e+00;
acadoVariables.ubValues[138] = 2.0000000000000000e+00;
acadoVariables.ubValues[139] = 2.0000000000000000e+00;
acadoVariables.ubValues[140] = 2.0000000000000000e+00;
acadoVariables.ubValues[141] = 2.0000000000000000e+00;
acadoVariables.ubValues[142] = 2.0000000000000000e+00;
acadoVariables.ubValues[143] = 2.0000000000000000e+00;
acadoVariables.ubValues[144] = 2.0000000000000000e+00;
acadoVariables.ubValues[145] = 2.0000000000000000e+00;
acadoVariables.ubValues[146] = 2.0000000000000000e+00;
acadoVariables.ubValues[147] = 2.0000000000000000e+00;
acadoVariables.ubValues[148] = 2.0000000000000000e+00;
acadoVariables.ubValues[149] = 2.0000000000000000e+00;
acadoVariables.ubValues[150] = 2.0000000000000000e+00;
acadoVariables.ubValues[151] = 2.0000000000000000e+00;
acadoVariables.ubValues[152] = 2.0000000000000000e+00;
acadoVariables.ubValues[153] = 2.0000000000000000e+00;
acadoVariables.ubValues[154] = 2.0000000000000000e+00;
acadoVariables.ubValues[155] = 2.0000000000000000e+00;
acadoVariables.ubValues[156] = 2.0000000000000000e+00;
acadoVariables.ubValues[157] = 2.0000000000000000e+00;
acadoVariables.ubValues[158] = 2.0000000000000000e+00;
acadoVariables.ubValues[159] = 2.0000000000000000e+00;
acadoVariables.ubValues[160] = 2.0000000000000000e+00;
acadoVariables.ubValues[161] = 2.0000000000000000e+00;
acadoVariables.ubValues[162] = 2.0000000000000000e+00;
acadoVariables.ubValues[163] = 2.0000000000000000e+00;
acadoVariables.ubValues[164] = 2.0000000000000000e+00;
acadoVariables.ubValues[165] = 2.0000000000000000e+00;
acadoVariables.ubValues[166] = 2.0000000000000000e+00;
acadoVariables.ubValues[167] = 2.0000000000000000e+00;
acadoVariables.ubValues[168] = 2.0000000000000000e+00;
acadoVariables.ubValues[169] = 2.0000000000000000e+00;
acadoVariables.ubValues[170] = 2.0000000000000000e+00;
acadoVariables.ubValues[171] = 2.0000000000000000e+00;
acadoVariables.ubValues[172] = 2.0000000000000000e+00;
acadoVariables.ubValues[173] = 2.0000000000000000e+00;
acadoVariables.ubValues[174] = 2.0000000000000000e+00;
acadoVariables.ubValues[175] = 2.0000000000000000e+00;
acadoVariables.ubValues[176] = 2.0000000000000000e+00;
acadoVariables.ubValues[177] = 2.0000000000000000e+00;
acadoVariables.ubValues[178] = 2.0000000000000000e+00;
acadoVariables.ubValues[179] = 2.0000000000000000e+00;
acadoVariables.ubValues[180] = 2.0000000000000000e+00;
acadoVariables.ubValues[181] = 2.0000000000000000e+00;
acadoVariables.ubValues[182] = 2.0000000000000000e+00;
acadoVariables.ubValues[183] = 2.0000000000000000e+00;
acadoVariables.ubValues[184] = 2.0000000000000000e+00;
acadoVariables.ubValues[185] = 2.0000000000000000e+00;
acadoVariables.ubValues[186] = 2.0000000000000000e+00;
acadoVariables.ubValues[187] = 2.0000000000000000e+00;
acadoVariables.ubValues[188] = 2.0000000000000000e+00;
acadoVariables.ubValues[189] = 2.0000000000000000e+00;
acadoVariables.ubValues[190] = 2.0000000000000000e+00;
acadoVariables.ubValues[191] = 2.0000000000000000e+00;
acadoVariables.ubValues[192] = 2.0000000000000000e+00;
acadoVariables.ubValues[193] = 2.0000000000000000e+00;
acadoVariables.ubValues[194] = 2.0000000000000000e+00;
acadoVariables.ubValues[195] = 2.0000000000000000e+00;
acadoVariables.ubValues[196] = 2.0000000000000000e+00;
acadoVariables.ubValues[197] = 2.0000000000000000e+00;
acadoVariables.ubValues[198] = 2.0000000000000000e+00;
acadoVariables.ubValues[199] = 2.0000000000000000e+00;
return ret;
}

void acado_initializeNodesByForwardSimulation(  )
{
int index;
for (index = 0; index < 200; ++index)
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
for (index = 0; index < 200; ++index)
{
acadoVariables.x[index * 4] = acadoVariables.x[index * 4 + 4];
acadoVariables.x[index * 4 + 1] = acadoVariables.x[index * 4 + 5];
acadoVariables.x[index * 4 + 2] = acadoVariables.x[index * 4 + 6];
acadoVariables.x[index * 4 + 3] = acadoVariables.x[index * 4 + 7];
}

if (strategy == 1 && xEnd != 0)
{
acadoVariables.x[800] = xEnd[0];
acadoVariables.x[801] = xEnd[1];
acadoVariables.x[802] = xEnd[2];
acadoVariables.x[803] = xEnd[3];
}
else if (strategy == 2) 
{
state[0] = acadoVariables.x[800];
state[1] = acadoVariables.x[801];
state[2] = acadoVariables.x[802];
state[3] = acadoVariables.x[803];
if (uEnd != 0)
{
state[24] = uEnd[0];
}
else
{
state[24] = acadoVariables.u[199];
}

acado_integrate(state, 1);

acadoVariables.x[800] = state[0];
acadoVariables.x[801] = state[1];
acadoVariables.x[802] = state[2];
acadoVariables.x[803] = state[3];
}
}

void acado_shiftControls( real_t* const uEnd )
{
int index;
for (index = 0; index < 199; ++index)
{
acadoVariables.u[index] = acadoVariables.u[index + 1];
}

if (uEnd != 0)
{
acadoVariables.u[199] = uEnd[0];
}
}

real_t acado_getKKT(  )
{
real_t kkt;

int index;
real_t prd;

kkt = + acadoWorkspace.g[0]*acadoWorkspace.x[0] + acadoWorkspace.g[1]*acadoWorkspace.x[1] + acadoWorkspace.g[2]*acadoWorkspace.x[2] + acadoWorkspace.g[3]*acadoWorkspace.x[3] + acadoWorkspace.g[4]*acadoWorkspace.x[4] + acadoWorkspace.g[5]*acadoWorkspace.x[5] + acadoWorkspace.g[6]*acadoWorkspace.x[6] + acadoWorkspace.g[7]*acadoWorkspace.x[7] + acadoWorkspace.g[8]*acadoWorkspace.x[8] + acadoWorkspace.g[9]*acadoWorkspace.x[9] + acadoWorkspace.g[10]*acadoWorkspace.x[10] + acadoWorkspace.g[11]*acadoWorkspace.x[11] + acadoWorkspace.g[12]*acadoWorkspace.x[12] + acadoWorkspace.g[13]*acadoWorkspace.x[13] + acadoWorkspace.g[14]*acadoWorkspace.x[14] + acadoWorkspace.g[15]*acadoWorkspace.x[15] + acadoWorkspace.g[16]*acadoWorkspace.x[16] + acadoWorkspace.g[17]*acadoWorkspace.x[17] + acadoWorkspace.g[18]*acadoWorkspace.x[18] + acadoWorkspace.g[19]*acadoWorkspace.x[19] + acadoWorkspace.g[20]*acadoWorkspace.x[20] + acadoWorkspace.g[21]*acadoWorkspace.x[21] + acadoWorkspace.g[22]*acadoWorkspace.x[22] + acadoWorkspace.g[23]*acadoWorkspace.x[23] + acadoWorkspace.g[24]*acadoWorkspace.x[24] + acadoWorkspace.g[25]*acadoWorkspace.x[25] + acadoWorkspace.g[26]*acadoWorkspace.x[26] + acadoWorkspace.g[27]*acadoWorkspace.x[27] + acadoWorkspace.g[28]*acadoWorkspace.x[28] + acadoWorkspace.g[29]*acadoWorkspace.x[29] + acadoWorkspace.g[30]*acadoWorkspace.x[30] + acadoWorkspace.g[31]*acadoWorkspace.x[31] + acadoWorkspace.g[32]*acadoWorkspace.x[32] + acadoWorkspace.g[33]*acadoWorkspace.x[33] + acadoWorkspace.g[34]*acadoWorkspace.x[34] + acadoWorkspace.g[35]*acadoWorkspace.x[35] + acadoWorkspace.g[36]*acadoWorkspace.x[36] + acadoWorkspace.g[37]*acadoWorkspace.x[37] + acadoWorkspace.g[38]*acadoWorkspace.x[38] + acadoWorkspace.g[39]*acadoWorkspace.x[39] + acadoWorkspace.g[40]*acadoWorkspace.x[40] + acadoWorkspace.g[41]*acadoWorkspace.x[41] + acadoWorkspace.g[42]*acadoWorkspace.x[42] + acadoWorkspace.g[43]*acadoWorkspace.x[43] + acadoWorkspace.g[44]*acadoWorkspace.x[44] + acadoWorkspace.g[45]*acadoWorkspace.x[45] + acadoWorkspace.g[46]*acadoWorkspace.x[46] + acadoWorkspace.g[47]*acadoWorkspace.x[47] + acadoWorkspace.g[48]*acadoWorkspace.x[48] + acadoWorkspace.g[49]*acadoWorkspace.x[49] + acadoWorkspace.g[50]*acadoWorkspace.x[50] + acadoWorkspace.g[51]*acadoWorkspace.x[51] + acadoWorkspace.g[52]*acadoWorkspace.x[52] + acadoWorkspace.g[53]*acadoWorkspace.x[53] + acadoWorkspace.g[54]*acadoWorkspace.x[54] + acadoWorkspace.g[55]*acadoWorkspace.x[55] + acadoWorkspace.g[56]*acadoWorkspace.x[56] + acadoWorkspace.g[57]*acadoWorkspace.x[57] + acadoWorkspace.g[58]*acadoWorkspace.x[58] + acadoWorkspace.g[59]*acadoWorkspace.x[59] + acadoWorkspace.g[60]*acadoWorkspace.x[60] + acadoWorkspace.g[61]*acadoWorkspace.x[61] + acadoWorkspace.g[62]*acadoWorkspace.x[62] + acadoWorkspace.g[63]*acadoWorkspace.x[63] + acadoWorkspace.g[64]*acadoWorkspace.x[64] + acadoWorkspace.g[65]*acadoWorkspace.x[65] + acadoWorkspace.g[66]*acadoWorkspace.x[66] + acadoWorkspace.g[67]*acadoWorkspace.x[67] + acadoWorkspace.g[68]*acadoWorkspace.x[68] + acadoWorkspace.g[69]*acadoWorkspace.x[69] + acadoWorkspace.g[70]*acadoWorkspace.x[70] + acadoWorkspace.g[71]*acadoWorkspace.x[71] + acadoWorkspace.g[72]*acadoWorkspace.x[72] + acadoWorkspace.g[73]*acadoWorkspace.x[73] + acadoWorkspace.g[74]*acadoWorkspace.x[74] + acadoWorkspace.g[75]*acadoWorkspace.x[75] + acadoWorkspace.g[76]*acadoWorkspace.x[76] + acadoWorkspace.g[77]*acadoWorkspace.x[77] + acadoWorkspace.g[78]*acadoWorkspace.x[78] + acadoWorkspace.g[79]*acadoWorkspace.x[79] + acadoWorkspace.g[80]*acadoWorkspace.x[80] + acadoWorkspace.g[81]*acadoWorkspace.x[81] + acadoWorkspace.g[82]*acadoWorkspace.x[82] + acadoWorkspace.g[83]*acadoWorkspace.x[83] + acadoWorkspace.g[84]*acadoWorkspace.x[84] + acadoWorkspace.g[85]*acadoWorkspace.x[85] + acadoWorkspace.g[86]*acadoWorkspace.x[86] + acadoWorkspace.g[87]*acadoWorkspace.x[87] + acadoWorkspace.g[88]*acadoWorkspace.x[88] + acadoWorkspace.g[89]*acadoWorkspace.x[89] + acadoWorkspace.g[90]*acadoWorkspace.x[90] + acadoWorkspace.g[91]*acadoWorkspace.x[91] + acadoWorkspace.g[92]*acadoWorkspace.x[92] + acadoWorkspace.g[93]*acadoWorkspace.x[93] + acadoWorkspace.g[94]*acadoWorkspace.x[94] + acadoWorkspace.g[95]*acadoWorkspace.x[95] + acadoWorkspace.g[96]*acadoWorkspace.x[96] + acadoWorkspace.g[97]*acadoWorkspace.x[97] + acadoWorkspace.g[98]*acadoWorkspace.x[98] + acadoWorkspace.g[99]*acadoWorkspace.x[99] + acadoWorkspace.g[100]*acadoWorkspace.x[100] + acadoWorkspace.g[101]*acadoWorkspace.x[101] + acadoWorkspace.g[102]*acadoWorkspace.x[102] + acadoWorkspace.g[103]*acadoWorkspace.x[103] + acadoWorkspace.g[104]*acadoWorkspace.x[104] + acadoWorkspace.g[105]*acadoWorkspace.x[105] + acadoWorkspace.g[106]*acadoWorkspace.x[106] + acadoWorkspace.g[107]*acadoWorkspace.x[107] + acadoWorkspace.g[108]*acadoWorkspace.x[108] + acadoWorkspace.g[109]*acadoWorkspace.x[109] + acadoWorkspace.g[110]*acadoWorkspace.x[110] + acadoWorkspace.g[111]*acadoWorkspace.x[111] + acadoWorkspace.g[112]*acadoWorkspace.x[112] + acadoWorkspace.g[113]*acadoWorkspace.x[113] + acadoWorkspace.g[114]*acadoWorkspace.x[114] + acadoWorkspace.g[115]*acadoWorkspace.x[115] + acadoWorkspace.g[116]*acadoWorkspace.x[116] + acadoWorkspace.g[117]*acadoWorkspace.x[117] + acadoWorkspace.g[118]*acadoWorkspace.x[118] + acadoWorkspace.g[119]*acadoWorkspace.x[119] + acadoWorkspace.g[120]*acadoWorkspace.x[120] + acadoWorkspace.g[121]*acadoWorkspace.x[121] + acadoWorkspace.g[122]*acadoWorkspace.x[122] + acadoWorkspace.g[123]*acadoWorkspace.x[123] + acadoWorkspace.g[124]*acadoWorkspace.x[124] + acadoWorkspace.g[125]*acadoWorkspace.x[125] + acadoWorkspace.g[126]*acadoWorkspace.x[126] + acadoWorkspace.g[127]*acadoWorkspace.x[127] + acadoWorkspace.g[128]*acadoWorkspace.x[128] + acadoWorkspace.g[129]*acadoWorkspace.x[129] + acadoWorkspace.g[130]*acadoWorkspace.x[130] + acadoWorkspace.g[131]*acadoWorkspace.x[131] + acadoWorkspace.g[132]*acadoWorkspace.x[132] + acadoWorkspace.g[133]*acadoWorkspace.x[133] + acadoWorkspace.g[134]*acadoWorkspace.x[134] + acadoWorkspace.g[135]*acadoWorkspace.x[135] + acadoWorkspace.g[136]*acadoWorkspace.x[136] + acadoWorkspace.g[137]*acadoWorkspace.x[137] + acadoWorkspace.g[138]*acadoWorkspace.x[138] + acadoWorkspace.g[139]*acadoWorkspace.x[139] + acadoWorkspace.g[140]*acadoWorkspace.x[140] + acadoWorkspace.g[141]*acadoWorkspace.x[141] + acadoWorkspace.g[142]*acadoWorkspace.x[142] + acadoWorkspace.g[143]*acadoWorkspace.x[143] + acadoWorkspace.g[144]*acadoWorkspace.x[144] + acadoWorkspace.g[145]*acadoWorkspace.x[145] + acadoWorkspace.g[146]*acadoWorkspace.x[146] + acadoWorkspace.g[147]*acadoWorkspace.x[147] + acadoWorkspace.g[148]*acadoWorkspace.x[148] + acadoWorkspace.g[149]*acadoWorkspace.x[149] + acadoWorkspace.g[150]*acadoWorkspace.x[150] + acadoWorkspace.g[151]*acadoWorkspace.x[151] + acadoWorkspace.g[152]*acadoWorkspace.x[152] + acadoWorkspace.g[153]*acadoWorkspace.x[153] + acadoWorkspace.g[154]*acadoWorkspace.x[154] + acadoWorkspace.g[155]*acadoWorkspace.x[155] + acadoWorkspace.g[156]*acadoWorkspace.x[156] + acadoWorkspace.g[157]*acadoWorkspace.x[157] + acadoWorkspace.g[158]*acadoWorkspace.x[158] + acadoWorkspace.g[159]*acadoWorkspace.x[159] + acadoWorkspace.g[160]*acadoWorkspace.x[160] + acadoWorkspace.g[161]*acadoWorkspace.x[161] + acadoWorkspace.g[162]*acadoWorkspace.x[162] + acadoWorkspace.g[163]*acadoWorkspace.x[163] + acadoWorkspace.g[164]*acadoWorkspace.x[164] + acadoWorkspace.g[165]*acadoWorkspace.x[165] + acadoWorkspace.g[166]*acadoWorkspace.x[166] + acadoWorkspace.g[167]*acadoWorkspace.x[167] + acadoWorkspace.g[168]*acadoWorkspace.x[168] + acadoWorkspace.g[169]*acadoWorkspace.x[169] + acadoWorkspace.g[170]*acadoWorkspace.x[170] + acadoWorkspace.g[171]*acadoWorkspace.x[171] + acadoWorkspace.g[172]*acadoWorkspace.x[172] + acadoWorkspace.g[173]*acadoWorkspace.x[173] + acadoWorkspace.g[174]*acadoWorkspace.x[174] + acadoWorkspace.g[175]*acadoWorkspace.x[175] + acadoWorkspace.g[176]*acadoWorkspace.x[176] + acadoWorkspace.g[177]*acadoWorkspace.x[177] + acadoWorkspace.g[178]*acadoWorkspace.x[178] + acadoWorkspace.g[179]*acadoWorkspace.x[179] + acadoWorkspace.g[180]*acadoWorkspace.x[180] + acadoWorkspace.g[181]*acadoWorkspace.x[181] + acadoWorkspace.g[182]*acadoWorkspace.x[182] + acadoWorkspace.g[183]*acadoWorkspace.x[183] + acadoWorkspace.g[184]*acadoWorkspace.x[184] + acadoWorkspace.g[185]*acadoWorkspace.x[185] + acadoWorkspace.g[186]*acadoWorkspace.x[186] + acadoWorkspace.g[187]*acadoWorkspace.x[187] + acadoWorkspace.g[188]*acadoWorkspace.x[188] + acadoWorkspace.g[189]*acadoWorkspace.x[189] + acadoWorkspace.g[190]*acadoWorkspace.x[190] + acadoWorkspace.g[191]*acadoWorkspace.x[191] + acadoWorkspace.g[192]*acadoWorkspace.x[192] + acadoWorkspace.g[193]*acadoWorkspace.x[193] + acadoWorkspace.g[194]*acadoWorkspace.x[194] + acadoWorkspace.g[195]*acadoWorkspace.x[195] + acadoWorkspace.g[196]*acadoWorkspace.x[196] + acadoWorkspace.g[197]*acadoWorkspace.x[197] + acadoWorkspace.g[198]*acadoWorkspace.x[198] + acadoWorkspace.g[199]*acadoWorkspace.x[199];
kkt = fabs( kkt );
for (index = 0; index < 200; ++index)
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

for (lRun1 = 0; lRun1 < 200; ++lRun1)
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
acadoWorkspace.objValueIn[0] = acadoVariables.x[800];
acadoWorkspace.objValueIn[1] = acadoVariables.x[801];
acadoWorkspace.objValueIn[2] = acadoVariables.x[802];
acadoWorkspace.objValueIn[3] = acadoVariables.x[803];
acado_evaluateLSQEndTerm( acadoWorkspace.objValueIn, acadoWorkspace.objValueOut );
acadoWorkspace.DyN[0] = acadoWorkspace.objValueOut[0] - acadoVariables.yN[0];
acadoWorkspace.DyN[1] = acadoWorkspace.objValueOut[1] - acadoVariables.yN[1];
acadoWorkspace.DyN[2] = acadoWorkspace.objValueOut[2] - acadoVariables.yN[2];
acadoWorkspace.DyN[3] = acadoWorkspace.objValueOut[3] - acadoVariables.yN[3];
objVal = 0.0000000000000000e+00;
for (lRun1 = 0; lRun1 < 200; ++lRun1)
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

