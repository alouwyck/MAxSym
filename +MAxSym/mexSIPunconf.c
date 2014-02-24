/****************************************************************************************************************************************************/
/*

    C source file for mex function to solve finite-difference equations for 
    2D axi-symmetric non-uniform unconfined ground water flow using SIP
    See mexSIPunconf.m for syntax discription

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
/* solve equations using SIP for unconfined system */
void SIPunconf(double *s, double *niter, // output 
               mwSize nz, mwSize nr, mwSize nper, mwSize *nt, double *dt, 
               double *qrc, double *qzc, double *qssc, double *qsyc, double D1, double *hskz,
               mwSize nq, mwSize *idq, double *q, mwSize ns0, mwSize *ids0, double *s0, signed char *ibound,
               double delta, mwSize mni,mwSize nparm, double wseed, double accl)
{
    
    // declarations
    mwSize i, j, p, t, ni; // i = layer index; j = row index; p = stress period index; t = time step index; ni = iteration index
    mwSize ij, ij2, ijt, tt; // ij = linear index for (i,j); ij2 = linear index for active; ijt = linear index for s(i,j,t) = ij+tt
    mwSize nznr, nz2, k, m, n; // nznr = nz*nr; nz2 = nz+2; k = index for v and e; m = nz+2 + 2*nr-1; n = (nr+1)*(nz+2)-1
    int iw; // omega index
    unsigned char rowwise, ok; // flag indicating if iteration is row-wise (rowwise=1) or column-wise (rowwise=0); flag to continue iterations
    double *rhs, *qp; // equation's right hand side; stress period discharges
    double *e, *v, d, *f, C, G, b, c, *omega, w; // SIP variables
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
    omega = (double *) mxMalloc(nparm*sizeof(double));
    
    // calculate omega values
    if (nparm > 1) { 
        for (iw=0; iw<nparm; iw++) *(omega+iw) = 1.0 - pow(wseed,((double)iw)/((double)(nparm-1)));
    } else *(omega) = 1.0 - wseed;
    
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
            iw = 0; // omega index

            while (difs>delta && ni<mni && ok>0){

                // omega
                w = *(omega+iw);

                // row-wise
                if (rowwise > 0) { 

                    // allocate f
                    if (nz > 1) f = (double *) mxMalloc((nznr-nr)*sizeof(double));

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
                                if (*(ibound+ij2-nz2)!=0) {
                                    c = *(qrc+ij); // D
                                    if (i==0) c *= 1.0 + (*(s+ijt) + *(s+ijt-nz)) / 2.0 / D1; // D
                                    *(v+k) -= c * *(s+ijt-nz); // Rhs 
                                    d -= c; // E
                                } else c = 0.0; // D
                                if (*(ibound+ij2+nz2)!=0) {
                                    *(e+k) = *(qrc+ij+nz); // F
                                    if (i==0) *(e+k) *= 1.0 + (*(s+ijt) + *(s+ijt+nz)) / 2.0 / D1; // F
                                    *(v+k) -= *(e+k) * *(s+ijt+nz); // Rhs
                                    d -= *(e+k); // E
                                } else *(e+k) = 0.0; // F
                                if (nz>1 && *(ibound+ij2-1)!=0) {
                                    b = *(qzc+ij); // B
                                    if (i==1 && *(hskz+j)>0) b = 1.0 / (1.0 / b + *(s+ijt-1) / *(hskz+j)); // B
                                    *(v+k) -= b * *(s+ijt-1); // Rhs
                                    d -= b; // E
                                } else if (nz>1) b = 0.0; // B
                                if (nz>1 && *(ibound+ij2+1)!=0) {
                                    tmp = *(qzc+ij+1); // H
                                    if (i==0 && *(hskz+j)>0) tmp = 1.0 / (1.0 / tmp + *(s+ijt) / *(hskz+j)); // H
                                    d -= tmp; // E
                                    *(v+k) -= tmp * *(s+ijt+1); // Rhs
                                    if (i<nz-1) *(f+k) = tmp;
                                } else if (nz>1 && i<nz-1) *(f+k) = 0.0;
                                *(v+k) = accl * (*(v+k) - d * *(s+ijt)); // Res

                                if (nz > 1 && i > 0 && j < nr-1) {
                                    b /= 1.0 + w * *(e+k-nr);
                                    C = *(e+k-nr) * b * w;
                                    d += C;
                                    *(e+k) -= C;
                                } 
                                if (nz > 1 && i < nz-1 && j > 0) {
                                    c /= 1.0 + w * *(f+k-1);
                                    G = *(f+k-1) * c * w;
                                    d += G; 
                                    *(f+k) -= G;
                                }
                                if (k > 0) {
                                    d -= c * *(e+k-1);
                                    *(v+k) -= c * *(v+k-1);
                                }
                                if (nz > 1 && i > 0) {
                                    d -= b * *(f+k-nr);
                                    *(v+k) -= b * *(v+k-nr);
                                }
                                *(v+k) /= d;
                                *(e+k) /= d;
                                if (nz > 1 && i < nz-1) *(f+k) /= d;

                            } else { // not active
                                *(e+k) = 0.0;
                                *(v+k) = 0.0;
                                if (nz > 1 && i < nz-1) *(f+k) = 0.0;
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
                                if (j < nr-1) *(v+k) -= *(e+k) * *(v+k+1);
                                if (i < nz-1) *(v+k) -= *(f+k) * *(v+k+nr);
                                *(s+ijt) += *(v+k);
                                if (fabs(*(v+k)) > difs) difs = fabs(*(v+k));
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

                    // allocate f
                    f = (double *) mxMalloc((nznr-nz)*sizeof(double));

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
                                if (*(ibound+ij2-1)!=0) {
                                    c = *(qzc+ij); // D
                                    if (i==1 && *(hskz+j)>0) c = 1.0 / (1.0 / c + *(s+ijt-1) / *(hskz+j)); // D
                                    *(v+ij) -= c * *(s+ijt-1); // Rhs 
                                    d -= c; // E
                                } else c = 0.0; // D
                                if (*(ibound+ij2+1)!=0) {
                                    *(e+ij) = *(qzc+ij+1); // F
                                    if (i==0 && *(hskz+j)>0) *(e+ij) = 1.0 / (1.0 / *(e+ij) + *(s+ijt) / *(hskz+j)); // F
                                    *(v+ij) -= *(e+ij) * *(s+ijt+1); // Rhs
                                    d -= *(e+ij); // E
                                } else *(e+ij) = 0.0; // F
                                if (*(ibound+ij2-nz2)!=0) {
                                    b = *(qrc+ij); // B
                                    if (i==0) b *= 1.0 + (*(s+ijt) + *(s+ijt-nz)) / 2.0 / D1; // B
                                    *(v+ij) -= b * *(s+ijt-nz); // Rhs
                                    d -= b; // E
                                } else b = 0.0; // B
                                if (*(ibound+ij2+nz2)!=0) {
                                    tmp = *(qrc+ij+nz); // H
                                    if (i==0) tmp *= 1.0 + (*(s+ijt) + *(s+ijt+nz)) / 2.0 / D1; // H
                                    *(v+ij) -= tmp * *(s+ijt+nz); // Rhs
                                    d -= tmp; // E
                                    if (j < nr-1) *(f+ij) = tmp;
                                } else if (j < nr-1) *(f+ij) = 0.0;

                                *(v+ij) = accl * (*(v+ij) - d * *(s+ijt)); // Res
                                
                                if (j > 0 && i < nz-1) {
                                    b /= 1.0 + w * *(e+ij-nz);
                                    C = *(e+ij-nz) * b * w;
                                    d += C;
                                    *(e+ij) -= C;
                                } 
                                if (j < nr-1 && i > 0) { 
                                    c /= 1.0 + w * *(f+ij-1);
                                    G = *(f+ij-1) * c * w;
                                    d += G; 
                                    *(f+ij) -= G;
                                }
                                if (ij > 0) {
                                    d -= c * *(e+ij-1);
                                    *(v+ij) -= c * *(v+ij-1);
                                }
                                if (j > 0) {
                                    d -= b * *(f+ij-nz);
                                    *(v+ij) -= b * *(v+ij-nz);
                                }
                                *(v+ij) /= d;
                                *(e+ij) /= d;
                                if (j < nr-1) *(f+ij) /= d;

                            } else { // not active
                                *(e+ij) = 0.0;
                                *(v+ij) = 0.0;
                                if (j < nr-1) *(f+ij) = 0.0;
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
                                if (j<nr-1) *(v+ij) -= *(f+ij)**(v+ij+nz);
                                *(s+ijt) += *(v+ij);
                                if (fabs(*(v+ij))>difs) difs = fabs(*(v+ij));
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

                // update counters and free memory
                ni++;
                if (nz>1) rowwise = 1 - rowwise;
                iw++;
                if (iw == nparm) iw = 0;
                if (nz>1) mxFree(f);

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
    mxFree(omega);
    
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
    mwSize mni, nparm;
    signed char *ibound;
    double *dt, *qrc, *qzc, *qrzc, *qssc, *qsyc;
    double D1, *hskz;
    double *q, *s0;
    double delta, wseed, accl;
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
    nparm = (mwSize) *mxGetPr(prhs[20]); // number of iteration parameters to be used
    wseed = *mxGetPr(prhs[21]); // seed for calculaton of iteration parameter omega
    accl = *mxGetPr(prhs[22]); // acceleration parameter
    
    // Output arguments
    ndim[0] = nz;
    ndim[1] = nr;
    ndim[2] = *(nt+nper)+1;
    plhs[0] = mxCreateNumericArray(3,&ndim,mxDOUBLE_CLASS,mxREAL);
    s = mxGetPr(plhs[0]); // drawdown
    plhs[1] = mxCreateDoubleMatrix(*(nt+nper),1,mxREAL);
    niter = mxGetPr(plhs[1]); // number of iterations

    // Call solver algorithm
    SIPunconf(s,niter,nz,nr,nper,nt,dt,qrc,qzc,qssc,qsyc,D1,hskz,nq,idq,q,ns0,ids0,s0,ibound,delta,mni,nparm,wseed,accl);
    
}
