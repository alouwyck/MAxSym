function [s,niter] = solveconfined(nz,nr,nper,ndt,dt,qrc,qzc,qssc,q,s0,inactive,constant,delta,mni,nparm,wseed,accl,mx)
% solve finite-difference equations for 2D axi-symmetric non-uniform confined ground water flow using ADI or SIP
%
% syntax: [s,niter] = solveconfined(nz,nr,nper,ndt,dt,qrc,qzc,qssc,q,s0,inactive,constant,delta,mni,nparm,wseed,accl,mx)
%  nz is number of layers
%  nr is number of rings
%  nper is number of stress periods
%  ndt is vector with number of time steps for each stress period, with length(ndt) = nper
%   if ndt is empty or smaller than 1, a steady state simulation is performed
%  dt is vector with time steps [T], with length(dt) = sum(ndt)
%   if the simulation is steady state, dt is ignored
%  qrc is nz x (nr-1) matrix with radial conductances [L²/T]
%  qzc is (nz-1) x nr matrix with vertical conductances [L²/T]
%   if the model has one layer only (nz = 1), qzc is ignored
%  qssc is nz x nr matrix with storage change terms [L²]
%   if the simulation is steady state, the qssc terms are ignored
%  q is nper cell where q{i} is nz x nr matrix with ring discharges [L³/T] for stress period i
%   a positive discharge value means that water is withdrawn
%   q is empty if all discharges are zero throughout the simulation 
%   q{i} is empty if all discharges are zero for stress period i 
%  s0 is nper cell where s0{i} is nz x nr matrix with initial head changes [L] at the start of stress period i 
%   a positive initial head change value means that there is an initial rising of head
%   s0 is empty if all initial head changes are zero throughout the simulation 
%   s0{i} is empty if all initial head changes are zero for stress period i 
%  inactive is nz x nr logical matrix indicating inactive rings
%   if inactive is empty, there are no inactive rings
%  constant is nz x nr logical matrix indicating constant-head rings
%   if constant is empty, there are no constant-head rings
%  delta is criterion for convergence [L]
%  mni is maximum number of iterations for one time step iteration
%  nparm is number of SIP iteration parameters
%   if nparm is empty or smaller than 1, ADI is used
%  wseed is seed for calculaton of SIP iteration parameters
%   if wseed is empty or smaller than zero, wseed is derived from the conductance matrices
%  accl is SIP accelerator, which is a value between 0 and 1
%   if accl is empty, accl is equal to 1
%  mx is logical indicating whether the mex routine must be called (true) or the matlab routine (false)
%  s is nz x nr x (sum(ndt)+1) drawdown [L] matrix
%   a positive drawdown value means that head has risen since the start of the simulation
%   note that s(:,:,1) = s0{1}
%   if the simulation is steady state, s is nz x nr matrix
%  niter is sum(ndt) vector with number of iterations performed during each time step
%   if the simulation is steady state, niter is scalar
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



% steady state?
if isscalar(ndt) && ndt < 1
    ndt = [];
end

% modify flow constants
[qrc,qzc,qrzc,qssc] = flowconstants(nz,nr,ndt,qrc,qzc,qssc);

% calculate ibound
ibound = getibound(nz,nr,inactive,constant);

% calculate nt
[nt,dt] = getnt(ndt,dt);

% convert q and s0 to list if mex
if mx
    [nq,idq,q] = numcell2list(q);
    [ns0,ids0,s0] = numcell2list(s0);
end

% ADI
if isempty(nparm) || nparm < 1
    
    % call mex file
    if mx
        [s,niter] = MAxSym.mexADIconf(nz,nr,nper,nt-1,dt,qrc,qzc,qrzc,qssc,nq,idq,q,ns0,ids0,s0,ibound,delta,mni);
    % call matlab routine
    else
        [s,niter] = MAxSym.ADIconf(nz,nr,nper,nt,dt,qrc,qzc,qrzc,qssc,q,s0,ibound,delta,mni);
    end

% SIP    
else
    
    % wseed and accl
    if isempty(wseed) || wseed < 0
        wseed = getwseed(nr,nz,qrc,qzc);
    end
    if isempty(accl)
        accl = 1;
    end
    
    % call mex file
    if mx
        [s,niter] = MAxSym.mexSIPconf(nz,nr,nper,nt-1,dt,qrc,qzc,qrzc,qssc,nq,idq,q,ns0,ids0,s0,ibound,delta,mni,nparm,wseed,accl);
    % call matlab routine
    else
        [s,niter] = MAxSym.SIPconf(nz,nr,nper,nt,dt,qrc,qzc,qrzc,qssc,q,s0,ibound,delta,mni,nparm,wseed,accl);
    end
    
end

% steady state?
if isempty(ndt)
    s(:,:,1) = [];
end

% inactive cells
if ~isempty(inactive)
    inactive = repmat(inactive,[1 1 size(s,3)]);
    s(inactive) = NaN;
end

end


function [qrc,qzc,qrzc,qssc] = flowconstants(nz,nr,ndt,qrc,qzc,qssc)
% modify or calculate flow constants
% nr is number of rings
% nz is number of layers
% ndt is nper vector with number of time steps for each stress period
% qrc (input) is nz x (nr-1) radial conductance matrix
% qrc (output) is nz x nr radial conductance matrix where qrc(:,1) = 0
% qzc (input) is (nz-1) x nr vertical conductance matrix
% qzc (output) is nz x nr vertical conductance matrix where qzc(1,:) = 0
% if nz = 1, qzc (output) = 0
% qssc is nz x nr matrix with storage decrease terms
% if steady state simulation, qssc (output) = zeros(nz,nr)
% qrzc is nz x nr matrix with total conductance:
% qrzc(i,j) = qrc(i,j-1) + qrc(i,j) + qzc(i-1,j) + qzc(i,j)
    z = zeros(nz,1); 
    qrc = [z, qrc]; 
    qrzc = [qrc, z]; 
    qrzc(isnan(qrzc)) = 0;
    qrzc = qrzc(:,1:end-1) + qrzc(:,2:end);  
    if nz > 1 
        z = zeros(1,nr); 
        qzc = [z; qzc]; 
        z = [qzc; z]; 
        z(isnan(z)) = 0;
        qrzc = qrzc + z(1:end-1,:) + z(2:end,:); 
    else
        qzc = 0;
    end 
    if isempty(ndt) % steady state 
        qssc = zeros(nz,nr); 
    end 
end


function ibound = getibound(nz,nr,inactive,constant)
% get nz x nr int8 ibound matrix
% ibound > 0: ring is active
% ibound = 0: ring is inactive
% ibound < 0: ring has constant head
% nr is number of rings
% nz is number of layers
% inactive is nz x nr logical matrix indicating inactive rings
% constant is nz x nr logical matrix indicating constant head rings
    if isempty(inactive) 
        inactive = false(nz,nr); 
    end 
    if isempty(constant) 
        constant = false(nz,nr); 
    end 
    inactive = [true(nz+2,1), [true(1,nr); inactive; true(1,nr)], true(nz+2,1)]; 
    constant = [false(nz+2,1), [false(1,nr); constant; false(1,nr)], false(nz+2,1)]; 
    ibound = int8(ones(nz+2,nr+2)); 
    ibound(constant) = int8(-1); 
    ibound(inactive) = int8(0); 
end


function [nt,dt] = getnt(ndt,dt)
% get time step variable nt
% ndt is nper vector with number of time steps for each stress period
% dt is vector with time steps
% nt is int32 vector with time step indices indicating the end of a each stress period
% nt(1) = 1 and indicates the initial drawdown condition at the start of the simulation
    if ~isempty(ndt) % transient
        nt = cumsum([1; ndt(:)]); 
        dt = dt(:); 
    else % steady state 
        nt = [1 2]; 
        dt = 1;
    end 
    nt = int32(nt);
end


function [n,id,v] = numcell2list(c)
% convert numeric q or s0 cell to list
% c is q or s0 cell
% n is number of list records
% id is n x 2 int32 matrix [iper, igrid]
%  where iper is stress period index
%  and igrid is linear grid index minus 1 (because it's to be used in C routine)
% v is column vector of length n with corresponding q or s0 values
% if n is zero, id and v are zero too
    v = [];
    id = int32([]);
    for i = 1:length(c)
        if ~isempty(c{i}) && any(c{i}(:))
            b = c{i} ~= 0;
            id = [id; int32([(i-1)*ones(sum(b(:)),1), find(b(:))-1])];
            v = [v; reshape(c{i}(b),[],1)];
        end
    end
    n = length(v);
    if ~n
        id = int32(0);
        v = 0;
    end
end


function wseed = getwseed(nr,nz,qrc,qzc)
% calculate wseed
% nr is number of rings
% nz is number of layers
% qrc is nz x nr radial conductance matrix
% qzc is nz x nr vertical conductance matrix
    if nz > 2
        zz = zeros(1,nr);
        zr = zeros(nz,1);
        iz = Inf(1,nr);
        ir = Inf(nz,1);
        qz = [qzc; zz];
        qr = [ir, qrc(:,2:end), ir];
        r1 = max(qz(1:end-1,:),qz(2:end,:))./...
             min(qr(:,1:end-1),qr(:,2:end));
        qz = [iz; qzc(2:end,:); iz];
        qr = [qrc, zr];
        r2 = max(qr(:,1:end-1),qr(:,2:end))./...
             min(qz(1:end-1,:),qz(2:end,:));
        b = ~isnan(r1) & ~isnan(r2) & ~isinf(r1) & ~isinf(r2) & r1 > 0 & r2 > 0;
        if any(b(:))
            wseed = min(pi^2./(1+r1(b))/2/nr^2,pi^2./(1+r2(b))/2/nz^2);
            wseed = mean(wseed);
        else
            wseed = 1e-7;
        end
    elseif nz == 2
        ir = Inf(nz,1);
        qr = [ir, qrc(:,2:end), ir];
        r = qzc./min(qr(:,1:end-1),qr(:,2:end));
        b = ~isnan(r) & ~isinf(r) & r > 0;
        if any(b(:))
            wseed = pi^2./(1+r(b))/2/nr^2;
            wseed = mean(wseed);
        else
            wseed = 1e-7;            
        end
    else
        wseed = 1e-7;
    end
end