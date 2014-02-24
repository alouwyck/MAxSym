/****************************************************************************************************************************************************/
/*

    C source file for mex function to solve finite-difference equations for 
    2D axi-symmetric non-uniform unconfined ground water flow using ADI
    See mexADIunconf.m for syntax discription

    Copyright (c) 2011, Ghent University (Andy Louwyck), Belgium
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
          notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
          notice, this list of conditions and the following disclaimer in the
          documentation and/or other materials provided with the distribution.
        * Neither the name of the Ghent University nor the
          names of its contributors may be used to endorse or promote products
          derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL GHENT UNIVERSITY BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
/****************************************************************************************************************************************************/
/* header */
#include <mex.h>
#include <matrix.h>
#include <stdlib.h>
#include <math.h>

/****************************************************************************************************************************************************/
/* solve equations using ADI for unconfined system */
void ADIunconf(double *s, double *niter, // output 
               mwSize nz, mwSize nr, mwSize nper, mwSize *nt, double *dt, 
               double *qrc, double *qzc, double *qssc, double *qsyc, double D1, double *hskz,
               mwSize nq, mwSize *idq, double *q, mwSize ns0, mwSize *ids0, double *s0, signed char *ibound,
               double delta, mwSize mni)
{
    
    // declarations
    mwSize i, j, p, t, ni; // i = layer index; j = row index; p = stress period index; t = time step index; ni = iteration index
    mwSize ij, ij2, ijt, tt; // ij = linear index for (i,j); ij2 = linear index for active; ijt = linear index for s(i,j,t) = ij+tt
    mwSize nznr, nz2, k, m, n; // nznr = nz*nr; nz2 = nz+2; k = index for v and e; m = nz+2 + 2*nr-1; n = (nr+1)*(nz+2)-1
    unsigned char rowwise, ok; // flag indicating if iteration is row-wise (rowwise=1) or column-wise (rowwise=0); flag to continue iterations
    double *rhs, *qp; // equation's right hand side; stress period discharges
    double *e, *v, d; // ADI variables
    double difs, tmp; // maximum absolute drawdown difference; temporary variable

    // set variables
    nz2 = nz+2;
    m = nz2 + 2*nr-1;
    n = (nr+1)*nz2 - 1;
    nznr = nz*nr;
    
    // allocate arrays
    rhs = (double *) mxMalloc(nznr*sizeof(double));
    qp = (double *) mxMalloc(nznr*sizeof(double));
    e = (double *) mxMalloc(nznr*sizeof(double));
    v = (double *) mxMalloc(nznr*sizeof(double));
    
    // stress period loop
    tt = 0;
    p = 0;
    ok = 1;
    while (p<nper && ok>0) {
        
        // add s0 to initial drawdown
        if (ns0 > 0) { 
            k = 0;
            ij = nznr * *(nt+p);
            while (k < ns0 && p >= *(ids0+k)) {
                if (p == *(ids0+k)) *(s+ij+*(ids0+ns0+k)) += *(s0+k);
                k++;
            }
        }
        
        // determine discharge matrix
        for (ij=0; ij<nznr; ij++) *(qp+ij) = 0.0;
        if (nq > 0) { 
            k = 0;
            while (k < nq && p >= *(idq+k)) {
                if (p == *(idq+k)) *(qp+*(idq+nq+k)) = *(q+k);
                k++;
            }
        }
        
        // time step loop
        t = *(nt+p);
        while (t<*(nt+p+1) && ok>0){

            // calculate rhs and copy s(t-1) to s(t)
            ij = 0;
            ij2 = nz2;
            for (j=0; j<nr; j++) {
                ij2++;
                for (i=0; i<nz; i++) {
                    ijt = ij+tt;
                    if (*(ibound+ij2) > 0) {
                        *(qssc+ij) /= *(dt+t);
                        if (i>0) {
                            *(rhs+ij) = -*(qssc+ij) * *(s+ijt) + *(qp+ij);
                        } else {
                            *(qsyc+j) /= *(dt+t);
                            *(rhs+ij) = -(*(qssc+ij)+*(qsyc+j)) * *(s+ijt) + *(qp+ij);
                        }
                    }
                    if (*(ibound+ij2)!=0) *(s+nznr+ijt) = *(s+ijt);
                    ij++;
                    ij2++;
                }
                ij2++;
            } 
            tt += nznr;

            // iteration loop
            ni = 0; // iteration index
            rowwise = 1;
            difs = delta + 1.0;

            while (difs>delta && ni<mni && ok>0){

                // row-wise
                if (rowwise > 0) { 

                    // forward substitution 
                    k = 0;
                    for (i=0; i<nz; i++) {
                        ij = i;
                        ij2 = nz2+i+1;

                        for (j=0; j<nr; j++) {
                            ijt = ij+tt;

                            if (*(ibound+ij2) > 0) {
                                if (i>0) {
                                    d = -*(qssc+ij); // E
                                    *(v+k) = *(rhs+ij); // Rhs
                                } else {
                                    d = -*(qsyc+j) - *(qssc+ij) * (1.0 + *(s+ijt) / D1); // E
                                    *(v+k) = *(rhs+ij) - *(s+ijt) * *(qssc+ij) * *(s+ijt-nznr) / D1; // Rhs
                                }
                                if (nz>1 && *(ibound+ij2-1)!=0) {
                                    tmp = *(qzc+ij); // B
                                    if (i==1 && *(hskz+j)>0) tmp = 1.0 / (1.0 / tmp + *(s+ijt-1) / *(hskz+j)); // B
                                    *(v+k) -= tmp * *(s+ijt-1); // Rhs
                                    d -= tmp; // E
                                }
                                if (nz>1 && *(ibound+ij2+1)!=0) {
                                    tmp = *(qzc+ij+1); // H
                                    if (i==0 && *(hskz+j)>0) tmp = 1.0 / (1.0 / tmp + *(s+ijt) / *(hskz+j)); // H
                                    *(v+k) -= tmp * *(s+ijt+1); // Rhs
                                    d -= tmp; // E
                                }
                                if (*(ibound+ij2-nz2)!=0) {
                                    tmp = *(qrc+ij); // D
                                    if (i==0) tmp *= 1.0 + (*(s+ijt)+*(s+ijt-nz)) / 2.0 / D1; // D
                                    d -= tmp * (1.0 + *(e+k-1)); // d = E - D*e(k-1)
                                    *(v+k) -= tmp * *(v+k-1); // v(k) = Rhs - D*v(k-1))
                                }
                                if (*(ibound+ij2+nz2)!=0) {
                                    tmp = *(qrc+ij+nz); // F
                                    if (i==0) tmp *= 1.0 + (*(s+ijt)+*(s+ijt+nz)) / 2.0 / D1; // F
                                    d -= tmp; // E
                                    *(e+k) =  tmp / d; // e(k) = F/d;
                                } else *(e+k) = 0.0;
                                *(v+k) /= d; // v/d

                            } else { // not active
                                *(e+k) = 0.0;
                                *(v+k) = 0.0;
                            }    

                            k++;
                            ij += nz;
                            ij2 += nz2;
                        }
                    } // end forward loops

                    // backward substitution 
                    difs = 0.0;
                    k = nznr-1;
                    ij = k;
                    for (i=nz-1; i>=0; i--) {
                        ijt = ij+tt;
                        ij2 = ij+m;

                        for (j=nr-1; j>=0; j--) {

                            if (*(ibound+ij2) > 0 && ok>0) {
                                if (j<nr-1) *(v+k) -= *(e+k)**(v+k+1);
                                if (fabs(*(v+k)-*(s+ijt))>difs) difs = fabs(*(v+k)-*(s+ijt));
                                *(s+ijt) = *(v+k);
                                if (i==0 && D1 <= -*(s+ijt)) {
                                    ok = 0;
                                }
                            }

                            ijt -= nz;
                            ij2 -= nz2;
                            k--;
                        }
                        ij--;
                    } // end backward loop

                // column-wise    
                } else {

                    // forward substitution
                    ij = 0;
                    ijt = tt;
                    ij2 = nz2;
                    for (j=0; j<nr; j++) {
                        ij2++;

                        for (i=0; i<nz; i++) {

                            if (*(ibound+ij2) > 0) {
                                if (i>0) {
                                    d = -*(qssc+ij); // E
                                    *(v+ij) = *(rhs+ij); // Rhs
                                } else {
                                    d = -*(qsyc+j) - *(qssc+ij) * (1.0 + *(s+ijt) / D1); // E
                                    *(v+ij) = *(rhs+ij) - *(s+ijt) * *(qssc+ij) * *(s+ijt-nznr) / D1; // Rhs
                                }
                                if (*(ibound+ij2-nz2)!=0) {
                                    tmp = *(qrc+ij); // B
                                    if (i==0) tmp *= 1.0 + (*(s+ijt) + *(s+ijt-nz)) / 2.0 / D1; // B
                                    d -= tmp; // E
                                    *(v+ij) -= tmp * *(s+ijt-nz); // Rhs
                                }
                                if (*(ibound+ij2+nz2)!=0) {
                                    tmp = *(qrc+ij+nz); // H
                                    if (i==0) tmp *= 1.0 + (*(s+ijt) + *(s+ijt+nz)) / 2.0 / D1; // H
                                    d -= tmp; // E
                                    *(v+ij) -= tmp * *(s+ijt+nz); // Rhs
                                }
                                if (*(ibound+ij2-1)!=0) {
                                    tmp = *(qzc+ij); // D
                                    if (i==1 && *(hskz+j)>0) tmp = 1.0 / (1.0 / tmp + *(s+ijt-1) / *(hskz+j)); // D
                                    d -= tmp * (1.0 + *(e+ij-1)); // d = E - D*e(ij-1)
                                    *(v+ij) -= tmp * *(v+ij-1); // v(ij) = (Rhs - D*v(ij-1))
                                }
                                if (*(ibound+ij2+1)!=0) {
                                    tmp = *(qzc+ij+1); // F
                                    if (i==0 && *(hskz+j)>0) tmp = 1.0 / (1.0 / tmp + *(s+ijt) / *(hskz+j)); // F
                                    d -= tmp; // E
                                    *(e+ij) = tmp / d; // e(ij) = F/d
                                } else *(e+ij) = 0.0;
                                *(v+ij) /= d; // v/d

                            } else { // not active
                                *(e+ij) = 0.0;
                                *(v+ij) = 0.0;
                            }    

                            ij++;
                            ijt++;
                            ij2++;                    
                        }
                        ij2++;
                    } // end forward loops

                    // backward substitution 
                    difs = 0.0;
                    ij = nznr-1;
                    ijt = tt+ij;
                    ij2 = n;
                    for (j=nr-1; j>=0; j--) {
                        ij2--;

                        for (i=nz-1; i>=0; i--) {

                            if (*(ibound+ij2) > 0 && ok>0) {
                                if (i<nz-1) *(v+ij) -= *(e+ij)**(v+ij+1);
                                if (fabs(*(v+ij)-*(s+ijt))>difs) difs = fabs(*(v+ij)-*(s+ijt));
                                *(s+ijt) = *(v+ij);
                                if (i==0 && D1 <= -*(s+ijt)) {
                                    ok = 0;
                                }
                            }

                            ij--;
                            ijt--;
                            ij2--;
                        }
                        ij2--;
                    } // end backward loop

                } // end row-wise branch

                // update counters
                ni++;
                if (nz>1) rowwise = 1 - rowwise;

            } // end of iteration loop

            // copy number of iterations
            if (ok>0) *(niter+t) = (double) ni;

            // augment time step index t
            ij = 0;
            for (j=0; j<nr; j++) {
                for (i=0; i<nz; i++) {
                    *(qssc+ij) *= *(dt+t);
                    if (i==0) *(qsyc+j) *= *(dt+t);
                    ij++;
                }
            } 
            if (ok>0) t++;
        
        } // end of time step loop
 
        // augment stress period index p
        if (ok>0) p++;
        
    } // end of stress period loop
    
    // free memory
    mxFree(rhs);
    mxFree(qp);
    mxFree(e);
    mxFree(v);
    
}


/****************************************************************************************************************************************************/
/* Gateway function */
void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{

    // Declarations
    mwSize i, n, ndim[3];
    mwSize nz, nr, nper, *nt; 
    mwSize nq, ns0, *idq, *ids0;
    mwSize mni;
    signed char *ibound;
    double *dt, *qrc, *qzc, *qrzc, *qssc, *qsyc;
    double *q, *s0;
    double D1, *hskz;
    double delta;
    double *s, *niter;

    // Input arguments
    nz = (mwSize) *mxGetPr(prhs[0]); // number of layers
    nr = (mwSize) *mxGetPr(prhs[1]); // number of rings
    nper = (mwSize) *mxGetPr(prhs[2]); // number of stress periods
    nt = mxGetPr(prhs[3]); // time step indices indicating the end of a each stress period (nper+1 values)
    dt = mxGetPr(prhs[4]); // time steps (nper values)
    qrc = mxGetPr(prhs[5]); // radial conductance (nz x nr matrix)
    qzc = mxGetPr(prhs[6]); // vertical conductance (nz x nr matrix)
    qssc = mxGetPr(prhs[7]); // storage rate constant related to specific elastic storage (nz x nr matrix)
    qsyc = mxGetPr(prhs[8]); // storage rate constant related to specific yield (nz x nr matrix)
    D1 = *mxGetPr(prhs[9]); // thickness of top layer
    hskz = mxGetPr(prhs[10]); // correction factor for vertical flow between top layers (nr values) - if hskz < 0, no correction is made
    nq = (mwSize) *mxGetPr(prhs[11]); // number of non-zero discharges
    idq = mxGetPr(prhs[12]); // stress period and cell indices of discharges (nq x 2 matrix)
    q = mxGetPr(prhs[13]); // discharges (nq values)
    ns0 = (mwSize) *mxGetPr(prhs[14]); // number of non-zero initial head changes
    ids0 = mxGetPr(prhs[15]); // stress period and cell indices of initial head changes (ns0 x 2 matrix)
    s0 = mxGetPr(prhs[16]); // initial head changes (ns0 values)
    ibound = mxGetPr(prhs[17]); // ibound matrix with cell characteristics: 1 is active, 0 is inactive, and -1 is constant (nz x nr matrix)
    delta = *mxGetPr(prhs[18]); // criterion for convergence
    mni = (mwSize) *mxGetPr(prhs[19]); // maximum number of iterations
    
    // Output arguments
    ndim[0] = nz;
    ndim[1] = nr;
    ndim[2] = *(nt+nper)+1;
    plhs[0] = mxCreateNumericArray(3,&ndim,mxDOUBLE_CLASS,mxREAL);
    s = mxGetPr(plhs[0]); // drawdown
    plhs[1] = mxCreateDoubleMatrix(*(nt+nper),1,mxREAL);
    niter = mxGetPr(plhs[1]); // number of iterations

    // Call solver algorithm
    ADIunconf(s,niter,nz,nr,nper,nt,dt,qrc,qzc,qssc,qsyc,D1,hskz,nq,idq,q,ns0,ids0,s0,ibound,delta,mni);
    
}
