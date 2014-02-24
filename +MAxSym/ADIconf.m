function [s,niter] = ADIconf(nz,nr,nper,nt,dt,qrc,qzc,qrzc,qssc,q,s0,ibound,delta,mni)
% Matlab ADI routine to solve finite-difference equations for 2D axi-symmetric non-uniform confined ground water flow
%
% syntax: [s,niter] = ADIconf(nz,nr,nper,nt,dt,qrc,qzc,qrzc,qssc,q,s0,ibound,delta,mni)
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
%  qssc is nz x nr matrix with storage change terms [L²]
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

            while (difs>delta && ni<mni)

                % row-wise
                if rowwise

                    % allocate v and e
                    [v,e] = deal(zeros(nznr,1));

                    % forward substitution 
                    k = 1;
                    for i = 1:nz
                        for j = 1:nr

                            if ibound(i+1,j+1) > 0

                                E = -qrzc(i,j) - qssc(i,j);
                                if ibound(i+1,j)
                                    D = qrc(i,j);
                                else
                                    if ibound(i+1,j) < 0
                                        E = E + qrc(i,j);
                                    end
                                    D = 0.0;
                                end
                                if ibound(i+1,j+2)
                                    F = qrc(i,j+1);
                                else
                                    if ibound(i+1,j+2) < 0 && j < nr
                                        E = E + qrc(i,j+1);
                                    end
                                    F = 0.0;
                                end
                                if nz > 1 && ibound(i,j+1)
                                    B = qzc(i,j)*s(i-1,j,t);
                                else
                                    if ibound(i,j+1) < 0 && nz > 1
                                        E = E + qzc(i,j);
                                    end
                                    B = 0.0;
                                end
                                if nz > 1 && ibound(i+2,j+1)
                                    H = qzc(i+1,j)*s(i+1,j,t);
                                else
                                    if ibound(i+2,j+1) < 0 && nz > 1 && i < nz
                                        E = E + qzc(i+1,j);
                                    end
                                    H = 0.0;
                                end
                                R = rhs(i,j) - B - H;
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

                                E = -qrzc(i,j) - qssc(i,j);
                                if ibound(i,j+1)
                                    D = qzc(i,j);
                                else
                                    if ibound(i,j+1) < 0
                                        E = E + qzc(i,j);
                                    end
                                    D = 0.0;
                                end
                                if ibound(i+2,j+1)
                                    F = qzc(i+1,j);
                                else
                                    if ibound(i+2,j+1) < 0 && i < nz
                                        E = E + qzc(i+1,j);
                                    end
                                    F = 0.0;
                                end
                                if ibound(i+1,j)
                                    B = qrc(i,j)*s(i,j-1,t);
                                else
                                    if ibound(i+1,j) < 0
                                        E = E + qrc(i,j);
                                    end
                                    B = 0.0;
                                end
                                if ibound(i+1,j+2)
                                    H = qrc(i,j+1)*s(i,j+1,t);
                                else
                                    if ibound(i+1,j+2) < 0 && j < nr
                                        E = E + qrc(i,j+1);
                                    end
                                    H = 0.0;
                                end
                                R = rhs(i,j) - B - H;
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

            % copy number of iterations
            niter(t-1,1) = ni;
            qssc = qssc * dt(t-1);

        end % end of time step loop

    end % end of stress period loop

    
end % function end
