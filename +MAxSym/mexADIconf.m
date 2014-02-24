% Mex ADI routine to solve finite-difference equations for 2D axi-symmetric non-uniform confined ground water flow
%
% syntax: [s,niter] = mexADIconf(nz,nr,nper,nt,dt,qrc,qzc,qrzc,qssc,nq,idq,q,ns0,ids0,s0,ibound,delta,mni)
%  nz is number of layers (double)
%  nr is number of rings (double)
%  nper is number of stress periods (double)
%  nt is int32 vector with zero based time step indices indicating the end of each stress period, with length(nt) = nper+1
%   nt(1) = 0 and indicates the initial drawdown condition at the start of the simulation
%   the number of time steps for each stress period is diff(nt)
%   the total number of time steps is nt(end)
%  dt is a double vector with time steps [T], with length(dt) = nt(end)
%  qrc is nz x nr double matrix with radial conductances [L²/T] where qrc(:,1) = 0 corresponds to the inner model boundary
%  qzc is nz x nr double matrix with vertical conductances [L²/T] where qzc(1:,) = 0 corresponds to the upper model boundary
%  qrzc is nz x nr double matrix with total conductances [L²/T]
%   total conductance of a ring is the sum of upper and lower vertical conductances and inner and outer conductances
%  qssc is nz x nr double matrix with storage decrease terms [L²]
%  nq is total number of nonzero discharges (double)
%  idq is nq x 2 int32 matrix with zero based stress period and grid indices of discharged rings
%   or idq = [iper(:), igrid(:)] where iper is stress period index 
%   and igrid is linear grid index equal to nr*iring + ilay with iring and ilay the ring and layer indices
%  q is nq double vector with discharges [L³/T] corresponding to idq
%   a positive discharge value means that water is withdrawn
%  ns0 is total number initial head changes (double)
%  ids0 is ns0 x 2 int32 matrix with zero based stress period and grid indices of rings with nonzero initial head change
%   or ids0 = [iper(:), igrid(:)] where iper is stress period index 
%   and igrid is linear grid index equal to nr*iring + ilay with iring and ilay the ring and layer indices
%  s0 is ns0 double vector with initial head changes [L] corresponding to ids0
%   a positive initial head change value means that there is an initial rising of head
%  ibound is nz x nr int8 matrix indicating the head property of each ring:
%   ibound > 0 indicates an active, variable-head ring
%   ibound = 0 indicates an inactive, impervious ring
%   ibound < 0 indicates an active, constant-head ring
%  delta is criterion for convergence [L] (double)
%  mni is maximum number of iterations for one time step iteration (double)
%  s is nz x nr x (nt(end)+1) drawdown [L] matrix (double)
%   a positive drawdown value means that head has risen since the start of the simulation
%   note that s(:,:,1) = s0{1}
%   and s(:,:,nt(i)+1) = s(:,:,nt(i)) + s0{i}
%  niter is nt(end) double vector with number of iterations performed during each time step
% 
% remark: a steady state simulation is performed by setting all qssc terms equal to zero and nt equal to [0 1]
%  an arbitrary non zero value is assigned to dt in this case
%
% attention!
% - indices are zero based (see nt, idq and ids0)
% - empty matrices are not allowed as input and must be replaced by zero (e.g. qzc if nz = 1)
% - some arguments need to be converted to int32 (see nt, idq and ids0) or int8 (see ibound)
%   other arguments are doubles (Matlab's default numeric data type)
%
%
% Copyright (c) 2011, Ghent University (Andy Louwyck), Belgium
% All rights reserved.
 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of the Ghent University nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL GHENT UNIVERSITY BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
