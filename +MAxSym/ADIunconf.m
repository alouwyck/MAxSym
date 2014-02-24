function [s,niter] = ADIunconf(nz,nr,nper,nt,dt,qrc,qzc,qssc,qsyc,D1,hskz,q,s0,ibound,delta,mni)
% Matlab ADI routine to solve finite-difference equations for 2D axi-symmetric non-uniform unconfined ground water flow
%
% syntax: [s,niter] = ADIunconf(nz,nr,nper,nt,dt,qrc,qzc,qssc,qsyc,D1,hskz,q,s0,ibound,delta,mni)
%  nz is number of layers
%  nr is number of rings
%  nper is number of stress periods
%  nt is vector with time step indices indicating the end of a each stress period, with length(nt) = nper+1
%   nt(1) = 1 and indicates the initial drawdown condition at the start of the simulation
%   the number of time steps for each stress period is diff(nt)
%   the total number of time steps is nt(end)-1
%  dt is vector with time steps [T], with length(dt) = nt(end)-1
%  qrc is nz x nr matrix with radial conductances [L²/T] where qrc(:,1) = 0 corresponds to the inner model boundary
%  qzc is nz x nr matrix with vertical conductances [L²/T] where qzc(1:,) = 0 corresponds to the upper model boundary
%  qssc is nz x nr matrix with storage decrease terms [L²]
%  qsyc is 1 x nr vector with specific yield terms [L²]
%  D1 is thickness [L] of upper phreatic layer
%  hskz is 1 x nr vector with vertical flow correction terms for the top layer
%   hskz(j) = 2 * hs(j) * kz(1,j) where hs(j) is the horizontal surface area [L²] of ring j 
%   and kz(1,j) the vertical conductivity [L/T] of ring j in the top layer
%   if hskz(j) <= 0, vertical flow in ring j between the top layer and the underlying layer is not corrected
%  q is nper cell where q{i} is nz x nr matrix with ring discharges [L³/T] for stress period i
%   a positive discharge value means that water is withdrawn
%   q is empty if all discharges are zero throughout the simulation 
%   q{i} is empty if all discharges are zero for stress period i 
%  s0 is nper cell where s0{i} is nz x nr matrix with initial head changes [L] at the start of stress period i 
%   a positive initial head change value means that there is an initial rising of head
%   s0 is empty if all initial head changes are zero throughout the simulation 
%   s0{i} is empty if all initial head changes are zero for stress period i 
%  ibound is nz x nr matrix indicating the head property of each ring:
%   ibound > 0 indicates an active, variable-head ring
%   ibound = 0 indicates an inactive, impervious ring
%   ibound < 0 indicates an active, constant-head ring
%  delta is criterion for convergence [L]
%  mni is maximum number of iterations for one time step iteration
%  s is nz x nr x nt(end) drawdown [L] matrix
%   a positive drawdown value means that head has risen since the start of the simulation
%   note that s(:,:,1) = s0{1}
%   and s(:,:,nt(i)+1) = s(:,:,nt(i)) + s0{i}
%  niter is nt(end)-1 vector with number of iterations performed during each time step
% 
% remark: a steady state simulation is performed by setting all qssc and qsyc terms equal to zero and nt equal to [1 2]
%  an arbitrary non zero value is assigned to dt in this case 
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



    % initialize s and niter
    s = zeros(nz,nr,nt(end));
    niter = zeros(nt(end)-1,1);

    % total number of rings
    nznr = nz*nr;
    
    % stress period loop
    ok = true;
    p = 1;
    while ok &&  p <= nper

        % initialization
        if ~(isempty(s0) || isempty(s0{p}))
            s(:,:,nt(p)) = s(:,:,nt(p)) + s0{p};
        end

        % time step loop
        t = nt(p) + 1;
        while ok && t <= nt(p+1)

            % calculate rhs and copy s(t-1) to s(t)
            qssc = qssc/dt(t-1);
            qsyc = qsyc/dt(t-1);
            rhs = -qssc .* s(:,:,t-1);
            rhs(1,:) = rhs(1,:) - qsyc .* s(1,:,t-1);            
            if ~(isempty(q) || isempty(q{p}))
                rhs = rhs + q{p};
            end
            s(:,:,t) = s(:,:,t-1);

            % iteration loop
            ni = 0; % iteration index
            rowwise = true;
            difs = delta + 1; % drawdown difference

            while (difs>delta && ni<mni && ok)

                % row-wise
                if rowwise

                    % allocate v and e
                    [v,e] = deal(zeros(nznr,1));

                    % forward substitution 
                    k = 1;
                    for i = 1:nz
                        for j = 1:nr

                            if ibound(i+1,j+1) > 0

                                if i > 1
                                    E = -qssc(i,j);
                                    R = rhs(i,j);
                                else
                                    E = -qsyc(j) - qssc(i,j) * (1.0 + s(i,j,t)/D1);
                                    R = rhs(i,j) - s(i,j,t)*qssc(i,j)*s(i,j,t-1)/D1;
                                end
                                if ibound(i+1,j)
                                    D = qrc(i,j);
                                    if i == 1
                                        D = D * (1.0 + (s(i,j,t) + s(i,j-1,t))/2/D1);
                                    end
                                    E = E - D;
                                else
                                    D = 0.0;
                                end
                                if ibound(i+1,j+2)
                                    F = qrc(i,j+1);
                                    if i == 1
                                        F = F * (1.0 + (s(i,j,t) + s(i,j+1,t))/2/D1);
                                    end
                                    E = E - F;
                                else
                                    F = 0.0;
                                end
                                if nz>1 && ibound(i,j+1)
                                    B = qzc(i,j);
                                    if i == 2 && hskz(j) > 0
                                        B = 1.0 / (1.0/B + s(i-1,j,t)/hskz(j));
                                    end
                                    E = E - B;
                                    B = B * s(i-1,j,t);
                                else
                                    B = 0.0;
                                end
                                if nz>1 && ibound(i+2,j+1)
                                    H = qzc(i+1,j);
                                    if i==1 && hskz(j) > 0
                                        H = 1.0 / (1.0/H + s(i,j,t)/hskz(j));
                                    end
                                    E = E - H;
                                    H = H * s(i+1,j,t);
                                else
                                    H = 0.0;
                                end
                                R = R - B - H;

                                if k>1
                                    d = E - D*e(k-1);
                                else
                                    d = E;
                                end
                                e(k) = F/d;
                                if k>1
                                    v(k) = (R - D*v(k-1)) / d;
                                else
                                    v(k) = R/d;
                                end

                            end 

                            k = k+1;

                        end
                    end % end forward loops

                    % backward substitution 
                    difs = 0.0;
                    k = nznr;
                    for i = nz:-1:1
                        for j = nr:-1:1

                            if ibound(i+1,j+1) > 0
                                if j < nr
                                    v(k) = v(k) - e(k)*v(k+1);
                                end
                                if abs(v(k)-s(i,j,t)) > difs 
                                    difs = abs(v(k)-s(i,j,t));
                                end
                                s(i,j,t) = v(k);
                                if i == 1 && D1 <= -s(i,j,t)
                                    disp(['Watertable is below top layer! (time step ' int2str(t-1) ')'])
                                    ok = false;
                                    s(:,:,t:end) = NaN;
                                    break
                                end
                            end

                            k = k-1;

                        end
                    end % end backward loop

                % column-wise    
                else

                    % allocate v and e
                    [v,e] = deal(zeros(nznr,1));

                    % forward substitution 
                    k = 1;
                    for j = 1:nr
                        for i = 1:nz

                            if ibound(i+1,j+1) > 0

                                if i > 1
                                    E = -qssc(i,j);
                                    R = rhs(i,j);
                                else
                                    E = -qsyc(j) - qssc(i,j) * (1.0 + s(i,j,t)/D1);
                                    R = rhs(i,j) - s(i,j,t)*qssc(i,j)*s(i,j,t-1)/D1;
                                end
                                if ibound(i,j+1)
                                    D = qzc(i,j);
                                    if i == 2 && hskz(j) > 0
                                        D = 1.0 / (1.0/D + s(i-1,j,t)/hskz(j));
                                    end
                                    E = E - D;
                                else
                                    D = 0.0;
                                end
                                if ibound(i+2,j+1)
                                    F = qzc(i+1,j);
                                    if i == 1 && hskz(j) > 0
                                        F = 1.0 / (1.0/F + s(i,j,t)/hskz(j));
                                    end
                                    E = E - F;
                                else
                                    F = 0.0;
                                end
                                if ibound(i+1,j)
                                    B = qrc(i,j);
                                    if i == 1
                                        B = B * (1.0 + (s(i,j,t) + s(i,j-1,t))/2/D1);
                                    end
                                    E = E - B;
                                    B = B * s(i,j-1,t);
                                else
                                    B = 0.0;
                                end
                                if ibound(i+1,j+2)
                                    H = qrc(i,j+1);
                                    if i == 1
                                        H = H * (1.0 + (s(i,j,t) + s(i,j+1,t))/2/D1);
                                    end
                                    E = E - H;
                                    H = H * s(i,j+1,t);
                                else
                                    H = 0.0;
                                end
                                R = R - B - H;

                                if k>1
                                    d = E - D*e(k-1);
                                else
                                    d = E;
                                end
                                e(k) = F/d;
                                if k>1
                                    v(k) = (R - D*v(k-1)) / d;
                                else
                                    v(k) = R/d;
                                end

                            end    

                            k = k+1;

                        end
                    end % end forward loops

                    % backward substitution 
                    difs = 0.0;
                    k = nznr;
                    for j = nr:-1:1
                        for i = nz:-1:1

                            if ibound(i+1,j+1) > 0

                                if i < nz
                                    v(k) = v(k) - e(k)*v(k+1);
                                end
                                if abs(v(k)-s(i,j,t)) > difs 
                                    difs = abs(v(k)-s(i,j,t));
                                end
                                s(i,j,t) = v(k);
                                if i == 1 && D1 <= -s(i,j,t)
                                    disp(['Watertable is below top layer! (time step ' int2str(t-1) ')'])
                                    ok = false;
                                    s(:,:,t:end) = NaN;
                                    break
                                end
                                
                            end

                            k = k-1;

                        end
                    end % end backward loop

                end % end row-wise branch

                % update counters
                ni = ni+1;
                if (nz>1)
                    rowwise = ~rowwise;
                end

            end % end of iteration loop

            % copy number of iterations and augment time step counter
            if ok
                niter(t-1,1) = ni;
            end
            qssc = qssc * dt(t-1);
            qsyc = qsyc * dt(t-1);
            t = t + 1;

        end % end of time step loop
        
        p = p + 1;
        
    end % end of stress period loop

        
end % end of function
