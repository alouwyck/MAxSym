classdef Model < MAxSym.Check
% class to build and run axi-symmetric model MAxSym
%
% 1.  create MAxSym.Model object: m = MAxSym.Model;
% 2.  define stress periods and time steps: m.<a href="matlab:help MAxSym.Model.settime">settime</a>(dt,ndt);
% 3.  define model grid: m.<a href="matlab:help MAxSym.Model.setgrid">setgrid</a>(rb,D,confined);
% 4.  define hydraulic parameters:
%     - radial conductivity: m.par.<a href="matlab:help MAxSym.Parameters.kr">kr</a> = ...;
%     - radial resistance: m.par.<a href="matlab:help MAxSym.Parameters.cr">cr</a> = ...;
%     - vertical conductivity: m.par.<a href="matlab:help MAxSym.Parameters.kz">kz</a> = ...;
%     - vertical resistance: m.par.<a href="matlab:help MAxSym.Parameters.cz">cz</a> = ...;
%     - specific elastic storage coefficient: m.par.<a href="matlab:help MAxSym.Parameters.ss">ss</a> = ...;
%     - specific yield: m.par.<a href="matlab:help MAxSym.Parameters.sy">sy</a> = ...;
% 5.  define inactive and/or constant head rings if any:
%     - inactive rings: m.par.<a href="matlab:help MAxSym.Parameters.inactive">inactive</a> = ...;
%     - constant head rings: m.par.<a href="matlab:help MAxSym.Parameters.constant">constant</a> = ...;
% 6.  define stresses for each stress period k: 
%     - discharge: m.stress(k).<a href="matlab:help MAxSym.Stresses.q">q</a> = ...;
%     - initial head change: m.stress(k).<a href="matlab:help MAxSym.Stresses.s0">s0</a> = ...;
% 7.  define solver: m.<a href="matlab:help MAxSym.Model.setsolver">setsolver</a>(delta,mni,nparm);
% 8.  run model: m.<a href="matlab:help MAxSym.Model.run">run</a>;
% 9.  check flow and volumetric budget terms: 
%     - radial flow terms: m.<a href="matlab:help MAxSym.Model.qr">qr</a>
%     - vertical flow terms: m.<a href="matlab:help MAxSym.Model.qz">qz</a>
%     - storage change terms: m.<a href="matlab:help MAxSym.Model.qs">qs</a>
%     - volumetric budget terms: m.<a href="matlab:help MAxSym.Model.bud">bud</a>
%     - total volumetric budget: m.<a href="matlab:help MAxSym.Model.totbud">totbud</a>
% 10. get drawdown matrix: m.<a href="matlab:help MAxSym.Model.s">s</a>
% 11. apply linear interpolation to get drawdown at arbitrary distances and times:
%     - linear interpolation in dimension of log(r): s = m.<a href="matlab:help MAxSym.Model.interp1r">interp1r</a>(ilay,r,it,rw);
%     - linear interpolation in dimension of log(t): s = m.<a href="matlab:help MAxSym.Model.interp1t">interp1t</a>(ilay,ir,t);
%     - bilinear interpolation in dimensions of log(r) and log(t): s = m.<a href="matlab:help MAxSym.Model.interp2">interp2</a>(ilay,r,t,rw);
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



    % input
    properties (SetAccess = protected)
        time = MAxSym.TimeSteps.empty; % time steps
        grid = MAxSym.Grid.empty; % model grid
        par = MAxSym.Parameters.empty; % hydraulic parameters
        stress = MAxSym.Stresses.empty; % stresses
        solver = MAxSym.Solver.empty; % solver
    end
    
    % output
    properties (SetAccess = protected)
        s % drawdown [L] - nz x nr x nt array, where nt is the number of simulation times
        niter % number of iterations for each time step - ndt vector, where ndt is the number of time steps
    end
    
    % output
    properties (Dependent)
        qr  % radial flow [L³/T] - nz x (nr+1) x nt array, where nt is the number of simulation times
        qz  % vertical flow [L³/T] - (nz+1) x nr x nt array, where nt is the number of simulation times, empty if nz = 1
        qs  % storage change [L³/T] - nz x nr x (nt-1) array, where nt is the number of simulation times, empty if steady state
        bud % volumetric budget [L³/T] - nz x nr x (nt-1) array, where nt is the number of simulation times
        totbud % total volumetric budget [L³/T] - (nt-1) x 2 array: column 1 refers to variable head rings, column 2 to constant head rings
    end
    
    % flow constants
    properties (Access = protected)
        qrc  % radial conductances [L²/T] - nz x nr-1 matrix
        qzc  % vertical conductances [L²/T] - nz-1 x nr matrix - empty if nz = 1
        qssc % storage change terms [L²] - nz x nr matrix - empty if steady state
        qsyc % specific yield terms [L²] - nr vector - empty if steady state or confined
        hskz % correction terms for vertical flow between top layers [L³/T] - nr vector - empty if nz = 1 or confined
    end
    
    methods
        
        function settime(obj,varargin)
        % define time steps
        %
        % if transient simulation:
        % syntax: obj.settime(dt,ndt)
        %  dt is time step vector [T]
        %   time steps must be infinite and strictly positive
        %  ndt is vector with number of time steps for each stress period
        %   note that length(dt) must be equal to sum(ndt)
        %   ndt is only required if more than 1 stress period
        %
        % if steady state simulation:
        % syntax: obj.settime(0)
        %
        % settime creates a TimeSteps object, which is saved in obj.time
        % if obj.time is already set, its contents is replaced by the new TimeSteps object
        % if a grid is defined, settime creates a Parameters object, which is saved in obj.par
        % also, a Stresses object is defined for each stress period, and the resulting Stresses vector is stored in obj.stress
        % if obj.par is already set, it is reset if the simulation state (steady or unsteady) was modified
        % if obj.stress is already set, it is reset if the number of stress periods was changed
        
            if ~isempty(obj.time)
                steady = obj.time.steady;
            else
                steady = obj.equal(varargin{1}(1),0);
            end
            obj.time = MAxSym.TimeSteps(varargin{:});
            if ~isempty(obj.grid)
                if isempty(obj.par) || (steady && ~obj.time.steady) || (~steady && obj.time.steady)
                    obj.setpar;
                end
                if isempty(obj.stress) || length(obj.stress) ~= obj.time.nper
                    obj.setstress;
                end
            end
        end
        
        function setgrid(obj,varargin)
        % define model grid
        %
        % syntax: obj.setgrid(rb,D,confined)
        %  rb is vector with radii of ring boundaries [L]
        %   rb must be strictly monotonically increasing
        %   rb must contain finite, strictly positive elements
        %   rb must contain two elements at least
        %  D is vector with layer thicknesses [L]
        %   D must contain finite, strictly positive elements
        %   D must contain one element at least
        %   D(1) is top layer thickness, D(end) is bottom layer thickness
        %  confined is logical indicating whether the aquifer system is confined (true) or not (false)
        %
        % setgrid creates a Grid object, which is saved in obj.grid
        % if obj.grid is already set, its contents is replaced by the new Grid object
        % if time steps are defined, setgrid creates a Parameters object, which is saved in obj.par
        % also, a Stresses object is defined for each stress period, and the resulting Stresses vector is stored in obj.stress
        % if obj.par or obj.stress are already set, they are reset if different grid dimensions [nz nr] were defined
        
            obj.grid = MAxSym.Grid(varargin{:});
            if ~isempty(obj.time)
                if isempty(obj.par) || obj.grid.nr ~= obj.par.nr || obj.grid.nz ~= obj.par.nz || ...
                   xor(obj.grid.confined,obj.par.confined)
                    obj.setpar;
                end
                if isempty(obj.stress) || obj.stress(1).nr ~= obj.grid.nr || obj.stress(1).nz ~= obj.grid.nz
                    obj.setstress;
                end
            end
        end
        
        function setsolver(obj,varargin)
        % define solver
        %
        % syntax: obj.setsolver(delta,mni,nparm)
        %  delta is criterion for convergence [L]
        %  mni is maximum number of iterations for one time step
        %  nparm (optional) is number of SIP iteration parameters (5 is sufficient in most cases)
        %   if nparm is given and greater than 0, SIP is used, otherwise, ADI is used
        %
        % setsolver creates a Solver object, which is stored in obj.solver
            obj.solver = MAxSym.Solver(varargin{:});
        end
        
        function run(obj)
        % run model
        %
        % syntax: obj.run
            
            % check input
            obj.checkinput;
            
            % set flow constants
            obj.setflowconstants;
            
            % solve equations
            obj.solve;
                        
        end
        
        function s = interp1r(obj,ilay,r,it,rw)
        % interpolate drawdown in dimension of log(r)
        %
        % if simulation is transient:
        % syntax: s = obj.interp1r(ilay,r,it,rw)
        %  ilay is vector with indices of considered layers
        %   if ilay is empty, all layers are considered or ilay = 1:obj.grid.nz
        %  r is vector with radial distances at which drawdown is interpolated
        %   radial distances must be between zero and obj.grid.r(end)
        %  it is vector with indices of considered simulation times
        %   if it is empty, all simulation times are considered or it = 1:length(obj.time.nt)
        %  rw (optional) is a vector of length(ilay) with well radius for each layer
        %   if rw is empty or not given, the well radius is assumed to coincide with the first nodal circle
        %   if rw is a scalar, all layers are assumed to have the same well radius rw
        %   for distances r smaller than or equal to rw, drawdown s is not interpolated,
        %    but equal to drawdown in the first ring or: if r <= rw(n), then s = obj.s(ilay(n),1,:)
        %  s is the interpolated drawdown array
        %   size(s) = [length(ilay), length(r), length(it)]
        %
        % if simulation is steady state
        % syntax: s = obj.interp1r(ilay,r,rw)
        %  ilay is vector with indices of considered layers
        %   if ilay is empty, all layers are considered or ilay = 1:obj.grid.nz
        %  r is vector with radial distances at which drawdown is interpolated
        %   radial distances must be between zero and obj.grid.r(end)
        %  rw (optional) is a vector of length(ilay) with well radius for each layer
        %   if rw is empty or not given, the well radius is assumed to coincide with the first nodal circle
        %   if rw is a scalar, all layers are assumed to have the same well radius rw
        %   for distances r smaller than or equal to rw, drawdown s is not interpolated,
        %    but equal to drawdown in the first ring or: if r <= rw(n), then s = obj.s(ilay(n),1)
        %  s is the interpolated drawdown matrix
        %   size(s) = [length(ilay), length(r)]
        
            % check number of input arguments
            if obj.time.steady
                if ~any(nargin == 3:4)
                    error('Method ''interp1r'' requires 3 or 4 input arguments (Model object included) if model is steady state!')
                end
            elseif ~any(nargin == 4:5)
                error('Method ''interp1r'' requires 4 or 5 input arguments (Model object included) if model is transient!')                
            end
            
            % discretization parameters
            nz = obj.grid.nz;
            nt = size(obj.s,3);
        
            % input ilay
            if ~(obj.vector(ilay,true,true) && obj.integer(ilay) && obj.between(ilay,[1 nz],[true true]))
                error(['''ilay'' must be integer vector with elements between 1 and  ' int2str(nz) '!'])
            end
            if isempty(ilay)
                ilay = 1:nz;
            end
            
            % input r
            if ~(obj.vector(r,true,false) && isnumeric(r) && obj.between(r,[0 obj.grid.r(end)],[true true]))
                error('''r'' must be numeric vector with elements between 0 and radius of outer nodal circle!')
            end
            
            % input it
            if ~obj.time.steady
                if ~(obj.vector(it,true,true) && obj.integer(it) && obj.between(it,[1 nt],[true true]))
                    error(['''it'' must be integer vector with elements between 1 and  ' int2str(nt) '!'])
                end
                if isempty(it)
                    it = 1:nt;
                end
            else
                if nargin > 3
                    rw = it;
                end
                it = 1;
            end
            
            % input rw
            if (obj.time.steady && nargin == 4) || nargin == 5
                if ~(obj.vector(rw,true,true) && isnumeric(rw) && obj.length(rw,length(ilay),true,true) && ...
                     obj.between(rw,[0 obj.grid.r(end)],[true true]))
                    error(['''rw'' must be numeric vector of length ' int2str(length(ilay)) ' with elements between 0 and radius of outer nodal circle!'])
                end
                rw = obj.repmat(rw,size(ilay));
            else
                rw = obj.grid.r(1) * ones(size(ilay));
            end
            
            % allocate s
            s = NaN(length(r),length(ilay),length(it));
            
            % interpolate
            if any(r > min(rw))
                s = interp1(log10(obj.grid.r),permute(obj.s(ilay,:,it),[2 1 3]),log10(r));
            end
            
            % permute
            s = permute(s,[2 1 3]);

            % well radius
            for n = 1:length(ilay)
                b = r <= rw(n);
                if any(b)
                    s(n,b,:) = repmat(obj.s(ilay(n),1,it),[1 sum(b) 1]);
                end
            end
            
        end
        
        function s = interp1t(obj,ilay,ir,t)
        % interpolate drawdown in dimension of log(t)
        %
        % syntax: s = obj.interp1t(ilay,ir,t)
        %  ilay is vector with indices of considered layers
        %   if ilay is empty, all layers are considered or ilay = 1:obj.grid.nz
        %  ir is vector with indices of considered rings
        %   if ir is empty, all rings are considered or ir = 1:obj.grid.nr
        %  t is vector with times at which drawdown is interpolated
        %   times must be between zero and obj.time.t(end)
        %   for times t smaller than obj.time.t(2), drawdown s is not interpolated,
        %    but equal to initial drawdown or: if t < obj.time.t(2), then s = obj.stress(1).s0
        %  s is the interpolated drawdown array
        %   size(s) = [length(ilay), length(ir), length(t)]
        %
        %  remark: method interp1t can only be used if the model is transient
            
            % check if not steady state
            if obj.time.steady
                error('Applying linear interpolation in time dimension is useless if model is steady state!')
            end
            
            % discretization parameters
            nz = obj.grid.nz;
            nr = obj.grid.nr;
        
            % input ilay
            if ~(obj.vector(ilay,true,true) && obj.integer(ilay) && obj.between(ilay,[1 nz],[true true]))
                error(['''ilay'' must be integer vector with elements between 1 and  ' int2str(nz) '!'])
            end
            if isempty(ilay)
                ilay = 1:nz;
            end
            
            % input ir
            if ~(obj.vector(ir,true,true) && obj.integer(ir) && obj.between(ir,[1 nr],[true true]))
                error(['''ir'' must be integer vector with elements between 1 and  ' int2str(nr) '!'])
            end
            if isempty(ir)
                ir = 1:nr;
            end
            
            % input t
            if ~(obj.vector(t,true,false) && isnumeric(t) && obj.between(t,[0 obj.time.t(end)],[true true]))
                error('''t'' must be numeric vector with elements between 0 and last simulation time!')
            end
            
            % allocate s
            s = NaN(length(t),length(ilay),length(ir));            
            t0 = obj.time.t(2);
            
            % instantaneous head changes
            [b,ti,si] = obj.times4s0(t);
            if any(b)
                s(b,:,:) = interp1(log10(ti),permute(si(ilay,ir,:),[3 1 2]),log10(t(b)));
            end
            
            % interpolate
            j = t < t0;
            b = ~(b|j);
            if any(b)
                s(b,:,:) = interp1(log10(obj.time.t(2:end)),permute(obj.s(ilay,ir,2:end),[3 1 2]),log10(t(b)));
            end
            
            % permute
            s = permute(s,[2 3 1]);
            
            % initial times t < t(2)
            if any(j)
               s(:,:,j) = repmat(obj.s(ilay,ir,1),[1 1 sum(j)]);
            end
            
        end
        
        function s = interp2(obj,ilay,r,t,rw)
        % interpolate drawdown in dimension of log(r) and dimension of log(t)
        %
        % syntax: s = obj.interp2(ilay,r,t,rw)
        %  ilay is vector with indices of considered layers
        %   if ilay is empty, all layers are considered or ilay = 1:obj.grid.nz
        %  r is vector with radial distances at which drawdown is interpolated
        %   radial distances must be between zero and obj.grid.r(end)
        %  t is vector with times at which drawdown is interpolated
        %   times must be between zero and obj.time.t(end)
        %   for times t smaller than obj.time.t(2), drawdown s is not interpolated,
        %    but equal to initial drawdown or: if t < obj.time.t(2), then s = obj.stress(1).s0
        %  rw (optional) is a vector of length(ilay) with well radius for each layer
        %   if rw is empty or not given, the well radius is assumed to coincide with the first nodal circle
        %   if rw is a scalar, all layers are assumed to have the same well radius rw
        %   for distances r smaller than or equal to rw, drawdown s is not interpolated,
        %    but equal to drawdown in the first ring or: if r <= rw(n), then s = obj.s(ilay(n),1,:)
        %  s is the interpolated drawdown array
        %   size(s) = [length(ilay), length(r), length(t)]
        %
        %  remark: method interp2 can only be used if the model is transient
            
            % discretization parameters
            nz = obj.grid.nz;
        
            % input ilay
            if ~(obj.vector(ilay,true,true) && obj.integer(ilay) && obj.between(ilay,[1 nz],[true true]))
                error(['''ilay'' must be integer vector with elements between 1 and  ' int2str(nz) '!'])
            end
            if isempty(ilay)
                ilay = 1:nz;
            end
            
            % input r
            if ~(obj.vector(r,true,false) && isnumeric(r) && obj.between(r,[0 obj.grid.r(end)],[true true]))
                error('''r'' must be numeric vector with elements between 0 and radius of outer nodal circle!')
            end
            
            % input t
            if ~(obj.vector(t,true,false) && isnumeric(t) && obj.between(t,[0 obj.time.t(end)],[true true]))
                error('''t'' must be numeric vector with elements between 0 and last simulation time!')
            end
            
            % input rw
            if nargin > 4
                if ~(obj.vector(rw,true,true) && isnumeric(rw) && obj.length(rw,length(ilay),true,true) && ...
                     obj.between(rw,[0 obj.grid.r(end)],[true true]))
                    error(['''rw'' must be numeric vector of length ' int2str(length(ilay)) ' with elements between 0 and radius of outer nodal circle!'])
                end
                rw = obj.repmat(rw,size(ilay));
            else
                rw = obj.grid.r(1) * ones(size(ilay));
            end
            
            % initialize variables
            s = NaN(length(t),length(r),length(ilay));            
            j = t < obj.time.t(2);
            [bt,ti,si] = obj.times4s0(t);
            [rax,tax] = deal(log10(obj.grid.r),log10(obj.time.t(2:end)));
            [r,t] = deal(log10(r(:)'),log10(t(:)));
            rw = log10(rw);

            % instantaneous head changes
            if any(bt)
                ti = log10(ti);
                for n = 1:length(ilay)
                    s(bt,:,n) = interp2(rax,ti,permute(si(ilay(n),:,:),[3 2 1]),r,t(bt));
                end
            end
            
            % initial times t < t(2)
            if any(j)
                tmp = interp1(rax,permute(obj.s(ilay,:,1),[2 1 3]),r);
                s(j,:,:) = permute(repmat(tmp,[1 1 sum(j)]),[3 1 2]);
            end
            
            % interpolate
            bt = ~(j|bt);
            for n = 1:length(ilay)
                b = r > rw(n);
                % points outside well
                if any(b) && any(bt)
                    s(bt,b,n) = interp2(rax,tax,permute(obj.s(ilay(n),:,2:end),[3 2 1]),r(b),t(bt));
                end
                % points inside well
                b = ~b;
                if any(b)
                    if any(~j)
                        s(~j,b,n) = repmat(interp1(tax,squeeze(obj.s(ilay(n),1,2:end)),t(~j)),[1 sum(b)]);
                    end
                    if any(j)
                        s(j,b,n) = s(ilay(n),1,1);
                    end
                end
            end
            
            % permute
            s = permute(s,[3 2 1]);

        end
        
        function qr = get.qr(obj)
            if isempty(obj.s)
                qr = [];
            else
                nt = size(obj.s,3);
                qr = repmat(obj.qrc,[1 1 nt]);
                if ~obj.grid.confined
                    qr(1,:,:) = qr(1,:,:) .* (1 + (obj.s(1,1:end-1,:)+obj.s(1,2:end,:))/2/obj.grid.D(1));
                end
                qr = qr .* -diff(obj.s,1,2);
                z = zeros(obj.par.nz,1,nt);
                qr = cat(2,z,qr,z);
                qr(isnan(qr)) = 0;
            end
        end
        
        function qz = get.qz(obj)
            if isempty(obj.s) || obj.par.nz == 1
                qz = [];
            else
                nt = size(obj.s,3);
                qz = repmat(obj.qzc,[1 1 nt]);
                if ~obj.grid.confined
                    b = obj.hskz > 0;
                    if any(b)
                        qz(1,b,:) = 1 ./ (1./qz(1,b,:) + obj.s(1,b,:)./repmat(obj.hskz(b),[1 1 nt]));
                    end
                end
                qz = qz .* diff(obj.s,1,1);
                z = zeros(1,obj.par.nr,nt);
                qz = cat(1,z,qz,z);
                qz(isnan(qz)) = 0;
            end
        end
        
        function qs = get.qs(obj)
            
            % no storage
            if isempty(obj.s) || obj.time.steady
                qs = [];
            
            % transient simulation
            else
                
                nt = size(obj.s,3) - 1;
                [nz,nr] = deal(obj.par.nz,obj.par.nr);
                qs = repmat(obj.qssc,[1 1 nt]);
                if ~obj.grid.confined
                     qs(1,:,:) = repmat(obj.qsyc,[1 1 nt]) + qs(1,:,:) .* (1 + obj.s(1,:,2:end)/obj.grid.D(1));
                end
                qs = qs .* diff(obj.s,1,3) ./ repmat(reshape(obj.time.dt,[1 1 nt]),nz,nr);
                
                % take into account initial head changes
                nt = cumsum(obj.time.ndt(:));
                for k = 1:length(nt)-1
                    s0 = obj.stress(k+1).s0;
                    if ~isempty(s0) && any(s0(:))
                        dsdt = obj.repmat(s0,[nz,nr])/obj.time.dt(nt(k));
                        qs(:,:,nt(k)) = qs(:,:,nt(k)) - obj.qssc .* dsdt;
                        if ~obj.grid.confined
                            qs(1,:,nt(k)) = qs(1,:,nt(k)) - obj.qsyc .* dsdt(1,:);
                        end
                    end
                end
                qs(isnan(qs)) = 0;
                
            end
            
        end
        
        function bud = get.bud(obj)
            
            % flow components
            [nz,nr] = deal(obj.par.nz,obj.par.nr);
            qr = obj.qr;
            bud = qr(:,2:end,:)-qr(:,1:end-1,:);
            if nz > 1
                qz = obj.qz;
                bud = bud - qz(2:end,:,:)+qz(1:end-1,:,:);
            end
            
            % variable head rings
            b = true(nz,nr);
            if ~isempty(obj.par.inactive)
                b = b & ~obj.repmat(obj.par.inactive,[nz,nr]);
            end
            if ~isempty(obj.par.constant)
                b = b & ~obj.repmat(obj.par.constant,[nz,nr]);
            end
            
            % steady state: sources and sinks
            if obj.time.steady
                if ~isempty(obj.stress.q)
                    q = full(obj.repmat(obj.stress.q,[nz,nr]));
                    bud(b) = bud(b) + q(b);
                end
                
            % transient: storage decrease and sources/sinks
            else
                bud = bud(:,:,2:end) + obj.qs;
                ndt = obj.time.ndt;
                nt = cumsum([0; ndt(:)]);
                b = ~b;
                for n = 1:obj.time.nper
                    if ~isempty(obj.stress(n).q)
                        q = full(obj.repmat(obj.stress(n).q,[nz,nr]));
                        q(b) = 0;
                        bud(:,:,nt(n)+1:nt(n+1)) = bud(:,:,nt(n)+1:nt(n+1)) + repmat(q,[1 1 ndt(n)]);
                    end
                end
            end
        end
        
        function totbud = get.totbud(obj)
            b = obj.repmat(obj.par.constant,[obj.grid.nz,obj.grid.nr]);
            bud = obj.bud;
            for n = 1:size(bud,3)
                tmp = bud(:,:,n);
                totbud(n,:) = [sum(tmp(~b)),sum(tmp(b))];
            end
        end
        
    end
    
    methods (Access = protected)
        
        function setpar(obj)
        % set obj.par
            obj.par = MAxSym.Parameters(obj.grid.nz,obj.grid.nr,obj.grid.confined,obj.time.steady);
        end
        
        function setstress(obj)
        % set obj.stress    
            for n = 1:obj.time.nper
                obj.stress(n) = MAxSym.Stresses(obj.grid.nz,obj.grid.nr);
            end
        end
        
        function checkinput(obj)
        % check model input    
            if isempty(obj.time)
                error('Time steps are not defined!')
            end
            if isempty(obj.grid)
                error('Grid is not defined!')
            end
            if ~obj.par.isdefined
                error('Some parameters are missing!')
            end
            for n = 1:length(obj.stress)
                b(n) = obj.stress(n).isdefined;
            end
            if ~any(b)
                error('No sources or sinks defined!')
            elseif obj.par.steady && (isempty(obj.par.constant) || ~any(obj.par.constant(:)))
                q = obj.stress.q;
                if isempty(q) || abs(sum(q(:)))>sum(q(q>0))*1.001
                    error('Steady state model is ill-posed!')
                end
            end
            if isempty(obj.solver)
                error('Solver is not defined!')
            end
        end
        
        function setflowconstants(obj)
        % set flow constants required for solver routine:
        %  obj.qrc, obj.qzc, obj.qssc, obj.qsyc, obj.hskz
        
            % get grid properties
            nr = obj.par.nr;
            nz = obj.par.nz;
            steady = obj.par.steady;
            confined = obj.par.confined;
            rb = obj.repmat(obj.grid.rb,[nz,nr+1]);
            D  = obj.repmat(obj.grid.D,[nz,nr]);
            hs = obj.repmat(obj.grid.hs,[nz,nr]);
            vol = obj.repmat(obj.grid.vol,[nz,nr]);
            ia = obj.repmat(obj.par.inactive,[nz,nr]);
            cst = obj.repmat(obj.par.constant,[nz,nr]);
            if isempty(ia)
                ia = false(nz,nr);
            end
            if isempty(cst)
                cst = false(nz,nr);
            end
            
            % set qrc
            kr = obj.repmat(obj.par.kr,[nz,nr]);
            cr = full(obj.repmat(obj.par.cr,[nz,nr-1]));
            if ~isempty(kr)
                dr = log(rb(:,2:end)./rb(:,1:end-1));
                obj.qrc = 4*pi * D(:,1:end-1) ./ (dr(:,1:end-1)./kr(:,1:end-1) + dr(:,2:end)./kr(:,2:end));
                if ~isempty(cr)
                    b = cr > 0 & ~isnan(cr);
                    dr = 2*pi * rb(:,2:end-1) .* D(:,1:end-1);
                    obj.qrc(b) = dr(b) ./ cr(b);
                end
            elseif ~isempty(cr)
                obj.qrc = 2*pi * rb(:,2:end-1) .* D(:,1:end-1) ./ obj.cr;
            else
                error('Radial conductivity is not defined!')
            end
            b = ~(ia(:,1:end-1) | ia(:,2:end));
            if any(isnan(obj.qrc(b)) | obj.qrc(b) < 0 | isinf(obj.qrc(b)))
                error('Radial conductivity is not defined properly!')
            end
            obj.qrc(~b) = 0;
            
            % set qzc
            if nz > 1
                kz = obj.repmat(obj.par.kz,[nz,nr]);
                cz = full(obj.repmat(obj.par.cz,[nz-1,nr]));
                if ~isempty(kz)
                    obj.qzc = 2*hs(1:end-1,:) ./ (D(1:end-1,:)./kz(1:end-1,:) + D(2:end,:)./kz(2:end,:));
                    if ~isempty(cz)
                        b = cz > 0 & ~isnan(cz);
                        tmp = hs(1:end-1,:);
                        obj.qzc(b) = tmp(b) ./ cz(b);
                    end
                elseif ~isempty(cz)
                    obj.qzc = hs(1:end-1,:) ./ cz;                
                else
                    error('Vertical conductivity is not defined!')
                end
                b = ~(ia(1:end-1,:) | ia(2:end,:));
                if any(isnan(obj.qzc(b)) | obj.qzc(b) < 0 | isinf(obj.qzc(b)))
                    error('Vertical conductivity is not defined properly!')
                end
                obj.qzc(~b) = 0;
            else
                obj.qzc = [];
            end
            
            % set qssc
            if ~steady
                ss  = obj.repmat(obj.par.ss,[nz,nr]);        
                if ~isempty(ss)
                    obj.qssc = vol .* ss;
                else
                    error('Specific elastic storage is not defined!')
                end
                b = ~(ia | cst);
                if any(isnan(obj.qssc(b)) | obj.qssc(b) < 0 | isinf(obj.qssc(b)))
                    error('Specific elastic storage is not defined properly!')
                end
                obj.qssc(~b) = 0;
            else
                obj.qssc = [];
            end
            
            % set qsyc
            if ~(steady || confined)
                sy  = obj.repmat(obj.par.sy,[1,nr]);
                if ~isempty(sy)
                    obj.qsyc = hs(1,:) .* sy;
                else
                    error('Specific yield is not defined!')
                end
                b = ~(ia(1,:) | cst(1,:));
                if any(isnan(obj.qsyc(b)) | obj.qsyc(b) < 0 | isinf(obj.qsyc(b)))
                    error('Specific yield is not defined properly!')
                end
                obj.qsyc(~b) = 0;
            else
                obj.qsyc = [];
            end
            
            % set hskz
            if nz > 1 && ~confined 
                if ~isempty(kz)
                    obj.hskz = 2 * hs(1,:) .* kz(1,:);
                    if ~isempty(cz)
                        b = cz(1,:) > 0 & ~isnan(cz(1,:));
                        obj.hskz(b) = -1;
                    end
                else
                    obj.hskz = -ones(1,nr);
                end
            else
                obj.hskz = [];
            end
            
        end
        
        function solve(obj)
        % calls appropriate solver routine:
        %  solveconfined if confined system
        %  solveunconfined if unconfined system
            [nz,nr] = deal(obj.par.nz,obj.par.nr);
            if obj.par.confined
                [obj.s,obj.niter] = MAxSym.solveconfined(nz,nr,obj.time.nper,obj.time.ndt,obj.time.dt,...
                                                         obj.qrc,obj.qzc,obj.qssc,...
                                                         obj.repcellelements({obj.stress.q},[nz,nr]),obj.repcellelements({obj.stress.s0},[nz,nr]),...
                                                         obj.repmat(obj.par.inactive,[nz,nr]),obj.repmat(obj.par.constant,[nz,nr]),...
                                                         obj.solver.delta,obj.solver.mni,obj.solver.nparm,obj.solver.wseed,obj.solver.accl,obj.solver.mex);            
            else
                [obj.s,obj.niter] = MAxSym.solveunconfined(nz,nr,obj.time.nper,obj.time.ndt,obj.time.dt,...
                                                           obj.qrc,obj.qzc,obj.qssc,obj.qsyc,obj.grid.D(1),obj.hskz,...
                                                           obj.repcellelements({obj.stress.q},[nz,nr]),obj.repcellelements({obj.stress.s0},[nz,nr]),...
                                                           obj.repmat(obj.par.inactive,[nz,nr]),obj.repmat(obj.par.constant,[nz,nr]),...
                                                           obj.solver.delta,obj.solver.mni,obj.solver.nparm,obj.solver.wseed,obj.solver.accl,obj.solver.mex);            
            end
        end
                
        function [b,ti,si] = times4s0(obj,t)
        % selects times for which interpolated drawdown must be corrected for initial drawdown s0
        %
        % syntax: [b,ti,si] = times4s0(obj,t)
        % t is vector with times at which drawdown is interpolated
        % b indicates which times t are inside time steps ending with an instantaneous head change
        % ti is vector with relevant simulation times 
        % si is array with relevant drawdowns corrected for instantaneous head changes
        %  size(si) = [m.grid.nz, m.grid.nr, length(ti)]
        %
        % times4s0 is called by methods interp1t and interp2
            b = false(size(t));
            [ti,si] = deal([]);
            it = cumsum([1;obj.time.ndt(:)]);
            for n = 2:length(obj.stress)
                s0 = obj.repmat(obj.stress(n).s0,[obj.grid.nz,obj.grid.nr]);
                if any(abs(s0(:))>0)
                    k = [it(n)-1, it(n)];
                    j = t > obj.time.t(k(1)) & t < obj.time.t(k(2));
                    if any(j)
                        b = b|j;
                        ti(end+1:end+2,1) = obj.time.t(k);
                        si = cat(3,si,obj.s(:,:,k(1)),obj.s(:,:,k(2))-s0);
                    end
                end
            end
        end
        
    end
    
end