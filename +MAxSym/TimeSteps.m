classdef TimeSteps < MAxSym.Check    
% class to define stress periods and time steps in axi-symmetric model MAxSym  
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



    properties (Dependent)
        steady % steady state simulation? true or false
        t      % simulation times [T] - zero if steady state - vector of length nt
        dt     % time steps [T] - empty if steady state - vector of length (nt-1)
        ndt    % number of time steps in each stress period - empty if steady state - vector of length nper 
        tper   % stress period times [T] - zero if steady state - vector of length (nper+1)
        dtper  % stress period lengths [T] - empty if steady state - vector of length nper
        nper   % number of stress periods - one if steady state
    end
    
    properties (Access = protected)
        tback % backing field for t
        ndtback % backing field for ndt
    end
    
    methods
        
        function obj = TimeSteps(dt,ndt)
        % create TimeSteps object
        % 
        % syntax (1): obj = TimeSteps(0)
        % creates object for steady state simulation
        %
        % syntax (2): obj = TimeSteps(dt,ndt)
        % creates object for transient simulation
        %  dt is time step vector [T]
        %   time steps must be infinite and strictly positive
        %  ndt (optional) is vector with number of time steps for each stress period
        %   note that length(dt) must be equal to sum(ndt)
        %   if ndt is not given, only one stress period is considered
        %  obj is TimeSteps object
        
            % input
            if nargin
                % check dt
                if obj.vector(dt,true,false) && isnumeric(dt)
                    % steady state
                    if isscalar(dt) && obj.equal(dt,0) 
                        obj.tback = 0;
                    % transient    
                    elseif obj.between(dt,[0 Inf],false(1,2))
                        obj.tback = [0; cumsum(dt(:))];
                        % ndt given
                        if nargin > 1 && ~isempty(ndt)
                            if obj.integer(ndt) && obj.between(ndt,[0 Inf],false(1,2)) && obj.length(dt,sum(ndt))
                                obj.ndtback = ndt(:);
                            else
                                error('Input ndt must be vector with positive integers and sum(ndt) must be equal to length(dt)!')
                            end
                        % ndt not given: 1 stress period    
                        else
                            obj.ndtback = length(dt);
                        end
                    else
                        error('Time steps must be finite and strictly positive!')
                    end
                else
                    error('Input dt must be non-empty numeric vector!')
                end
            end
            
        end
        
        function tf = get.steady(obj)
            tf = getsteady(obj);
        end
        
        function t = get.t(obj)
            t = obj.gett;
        end
        
        function dt = get.dt(obj)
            dt = obj.getdt;
        end
        
        function ndt = get.ndt(obj)
            ndt = obj.getndt;
        end
        
        function tper = get.tper(obj)
            tper = obj.gettper;
        end
        
        function dtper = get.dtper(obj)
            dtper = obj.getdtper;
        end
        
        function nper = get.nper(obj)
            nper = obj.getnper;
        end
        
    end
    
    methods (Access = protected)
        
        function tf = getsteady(obj)
            tf = length(obj.tback) == 1;
        end
        
        function t = gett(obj)
            t = obj.tback;
        end
        
        function dt = getdt(obj)
            dt = diff(obj.t);
        end
        
        function ndt = getndt(obj)
            ndt = obj.ndtback;
        end
        
        function tper = gettper(obj)
            tper = obj.t(cumsum([1; obj.ndt]));
        end
        
        function dtper = getdtper(obj)
            dtper = diff(obj.tper);
        end
        
        function nper = getnper(obj)
            if obj.steady
                nper = 1;
            else
                nper = length(obj.ndt);
            end
        end
        
    end
    
end 