function [s,niter] = SIPconf(nz,nr,nper,nt,dt,qrc,qzc,qrzc,qssc,q,s0,ibound,delta,mni,nparm,wseed,accl)
% Matlab SIP routine to solve finite-difference equations for 2D axi-symmetric non-uniform confined ground water flow
%
% syntax: [s,niter] = SIPconf(nz,nr,nper,nt,dt,qrc,qzc,qrzc,qssc,q,s0,ibound,delta,mni,nparm,wseed,accl)
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
%  qrzc is nz x nr matrix with total conductances [L²/T]
%   total conductance of a ring is the sum of upper and lower vertical conductances and inner and outer conductances
%  qssc is nz x nr matrix with storage decrease terms [L²]
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
%  nparm is number of iteration parameters (5 is sufficient in most cases)
%  wseed is seed for calculaton of iteration parameters
%  accl is accelerator (value between 0 and 1, in most cases 1 or close to 1)
%  s is nz x nr x nt(end) drawdown [L] matrix
%   a positive drawdown value means that head has risen since the start of the simulation
%   note that s(:,:,1) = s0{1}
%   and s(:,:,nt(i)+1) = s(:,:,nt(i)) + s0{i}
%  niter is nt(end)-1 vector with number of iterations performed during each time step
% 
% remark: a steady state simulation is performed by setting all qssc terms equal to zero and nt equal to [1 2]
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



    % initialize s
    s = zeros(nz,nr,nt(end));

    % total number of rings
    nznr = nz*nr;

    % calculate gamma values
    if nparm > 1
        g = 1.0 - wseed.^((0:nparm-1)/(nparm-1));
    else
        g = 1.0 - wseed;
    end
    

    % stress period loop
    for p = 1:nper

        % initialization
        if ~(isempty(s0) || isempty(s0{p}))
            s(:,:,nt(p)) = s(:,:,nt(p)) + s0{p};
        end
    
        % time step loop
        for t = nt(p)+1:nt(p+1)

            % calculate rhs and copy s(t-1) to s(t)
            qssc = qssc/dt(t-1);
            rhs = -qssc .* s(:,:,t-1);
            if ~(isempty(q) || isempty(q{p}))
                rhs = rhs + q{p};
            end 
            s(:,:,t) = s(:,:,t-1);

            % iteration loop
            ni = 0; % iteration index
            rowwise = true;
            difs = delta + 1; % drawdown difference
            iw = 1; % gamma index

            while (difs>delta && ni<mni)

                % gamma
                w = g(iw);

                % row-wise
                if rowwise

                    % allocate v, e, and f
                    [v,e] = deal(zeros(nznr,1));
                    f = zeros(nznr-nr,1);

                    % forward substitution 
                    k = 1;
                    for i = 1:nz
                        for j = 1:nr

                            if ibound(i+1,j+1) > 0

                                E = -qrzc(i,j) - qssc(i,j);
                                if ibound(i+1,j)
                                    D = qrc(i,j);
                                    sD = s(i,j-1,t);
                                else
                                    if ibound(i+1,j) < 0
                                        E = E + qrc(i,j);
                                    end
                                    D = 0.0;
                                    sD = 0.0;
                                end
                                if ibound(i+1,j+2)
                                    F = qrc(i,j+1);
                                    sF = s(i,j+1,t);
                                else
                                    if ibound(i+1,j+2) < 0 && j < nr
                                        E = E + qrc(i,j+1);
                                    end
                                    F = 0.0;
                                    sF = 0.0;
                                end
                                if nz > 1 && ibound(i,j+1)
                                    B = qzc(i,j);
                                    sB = s(i-1,j,t);
                                else
                                    if ibound(i,j+1) < 0 && nz > 1
                                        E = E + qzc(i,j);
                                    end
                                    B = 0.0;
                                    sB = 0.0;
                                end
                                if nz > 1 && ibound(i+2,j+1)
                                    H = qzc(i+1,j);
                                    sH = s(i+1,j,t);
                                else
                                    if ibound(i+2,j+1) < 0 && nz > 1 && i < nz
                                        E = E + qzc(i+1,j);
                                    end
                                    H = 0.0;
                                    sH = 0.0;
                                end

                                v(k) = accl * (rhs(i,j) -E*s(i,j,t) - D*sD - F*sF - B*sB - H*sH);

                                d = E;
                                if i > 1
                                    if j < nr
                                        b = B / (1.0 + w*e(k-nr));
                                        C = e(k-nr)*b*w;
                                        d = d + C;
                                        e(k) = F - C;
                                    else
                                        b = B;
                                        e(k) = F;
                                    end
                                else
                                    e(k) = F;
                                end
                                if i < nz
                                    if j > 1
                                        c = D / (1.0 + w*f(k-1));
                                        G = f(k-1)*c*w;
                                        d = d + G; 
                                        f(k) = H - G;
                                    else
                                        f(k) = H;
                                        c = 0.0;
                                    end
                                else
                                    c = D;
                                end
                                if k > 1
                                    d = d - c*e(k-1);
                                    v(k) = v(k) - c*v(k-1);
                                end
                                if i > 1
                                    d = d - b*f(k-nr);
                                    v(k) = v(k) - b*v(k-nr);
                                end
                                v(k) = v(k)/d;
                                e(k) = e(k)/d;
                                if i < nz
                                    f(k) = f(k)/d;
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
                                if i < nz
                                    v(k) = v(k) - f(k)*v(k+nr);
                                end
                                s(i,j,t) = v(k) + s(i,j,t);
                                if abs(v(k)) > difs 
                                    difs = abs(v(k));
                                end

                            end

                            k = k-1;

                        end
                    end % end backward loop

                % column-wise    
                else

                    % allocate v, e, and f
                    [v,e] = deal(zeros(nznr,1));
                    f = zeros(nznr-nz,1);

                    % forward substitution
                    k = 1;
                    for j = 1:nr
                        for i = 1:nz

                            if ibound(i+1,j+1) > 0

                                E = -qrzc(i,j) - qssc(i,j);
                                if ibound(i,j+1)
                                    D = qzc(i,j);
                                    sD = s(i-1,j,t);
                                else
                                    if ibound(i,j+1) < 0
                                        E = E + qzc(i,j);
                                    end
                                    D = 0.0;
                                    sD = 0.0;
                                end
                                if ibound(i+2,j+1)
                                    F = qzc(i+1,j);
                                    sF = s(i+1,j,t);
                                else
                                    if ibound(i+2,j+1) < 0 && i < nz
                                        E = E + qzc(i+1,j);
                                    end
                                    F = 0.0;
                                    sF = 0.0;
                                end
                                if ibound(i+1,j)
                                    B = qrc(i,j);
                                    sB = s(i,j-1,t);
                                else
                                    if ibound(i+1,j) < 0
                                        E = E + qrc(i,j);
                                    end
                                    B = 0.0;
                                    sB = 0.0;
                                end
                                if ibound(i+1,j+2)
                                    H = qrc(i,j+1);
                                    sH = s(i,j+1,t);
                                else
                                    if ibound(i+1,j+2) < 0 && j < nr
                                        E = E + qrc(i,j+1);
                                    end
                                    H = 0.0;
                                    sH = 0.0;
                                end

                                v(k) = accl * (rhs(i,j) -E*s(i,j,t) - D*sD - F*sF - B*sB - H*sH);

                                d = E;
                                if j > 1
                                    if i < nz
                                        b = B / (1.0 + w*e(k-nz));
                                        C = e(k-nz)*b*w;
                                        d = d + C;
                                        e(k) = F - C;
                                    else
                                        b = B;
                                        e(k) = F;
                                    end
                                else
                                    e(k) = F;
                                end
                                if j < nr
                                    if i > 1
                                        c = D / (1.0 + w*f(k-1));
                                        G = f(k-1)*c*w;
                                        d = d + G; 
                                        f(k) = H - G;
                                    else
                                        f(k) = H;
                                        c = 0.0;
                                    end
                                else
                                    c = D;
                                end
                                if k > 1
                                    d = d - c*e(k-1);
                                    v(k) = v(k) - c*v(k-1);
                                end
                                if j > 1
                                    d = d - b*f(k-nz);
                                    v(k) = v(k) - b*v(k-nz);
                                end
                                v(k) = v(k)/d;
                                e(k) = e(k)/d;
                                if j < nr
                                    f(k) = f(k)/d;
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
                                if j < nr
                                    v(k) = v(k) - f(k)*v(k+nz);
                                end
                                s(i,j,t) = v(k) + s(i,j,t);
                                if abs(v(k)) > difs 
                                    difs = abs(v(k));
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
                iw = iw+1;
                if iw == nparm+1 
                    iw = 1;
                end

            end % end of iteration loop

            % copy number of iterations
            niter(t-1,1) = ni;
            qssc = qssc * dt(t-1);

        end % end of time step loop

    end % end of stress period loop


end % function end
