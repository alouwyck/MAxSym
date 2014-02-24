classdef Solver < MAxSym.Check
% class to define solver for simulation of axi-symmetric model MAxSym  
%
% remark:
%  ADI is used if property nparm is empty or smaller than 1
%  otherwise SIP is used
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
        delta % criterion for convergence [L] (e.g. 1e-5)
        mni % maximum number of iterations for one time step
        nparm % number of SIP iteration parameters (e.g. 5) - if empty or <= 0, ADI will be used
        wseed = -1; % seed for calculaton of SIP iteration parameters (e.g. 1e-5) - optional
        accl = 1; % SIP accelerator (value between 0 and 1) - optional, default is 1
        mex = true; % call mex file to solve equations? true (default) or false
    end
    
    properties (Access = protected)
        deltaback % backing field for delta
        mniback % backing field for mni
        nparmback % backing field for nparm
        wseedback = -1; % backing field for wseed
        acclback = 1; % backing field for accl
        mexback = true; % backing field for mex
    end
    
    methods
        
        function obj = Solver(delta,mni,nparm)
        % create Solver object
        %
        % syntax: obj = Solver(delta,mni,nparm)
        %  delta is criterion for convergence [L]
        %  mni is maximum number of iterations for one time step
        %  nparm (optional) is number of SIP iteration parameters (5 is sufficient in most cases)
        %   if nparm is given and greater than 0, SIP is used, otherwise, ADI is used
            
            if any(nargin == [2 3])
                obj.delta = delta;
                obj.mni = mni;
                if nargin == 3
                    obj.nparm = nparm;
                end
            elseif nargin
                error('Wrong number of input arguments!')
            end
        end
        
        function delta = get.delta(obj)
            delta = obj.getdelta;
        end
        
        function set.delta(obj,delta)
            obj.setdelta(delta);
        end
        
        function mni = get.mni(obj)
            mni = obj.getmni;
        end
        
        function set.mni(obj,mni)
            obj.setmni(mni);
        end
        
        function nparm = get.nparm(obj)
            nparm = obj.getnparm;
        end
        
        function set.nparm(obj,nparm)
            obj.setnparm(nparm);
        end
        
        function wseed = get.wseed(obj)
            wseed = obj.getwseed;
        end
        
        function set.wseed(obj,wseed)
            obj.setwseed(wseed);
        end
        
        function accl = get.accl(obj)
            accl = obj.getaccl;
        end
        
        function set.accl(obj,accl)
            obj.setaccl(accl);
        end
        
        function mx = get.mex(obj)
            mx = obj.getmex;
        end
        
        function set.mex(obj,mx)
            obj.setmex(mx);
        end
        
    end
    
    methods (Access = protected)
        
        function delta = getdelta(obj)
            delta = obj.deltaback;
        end
        
        function setdelta(obj,delta)
            if obj.scalar(delta,false) && isnumeric(delta) && obj.between(delta,[0 Inf],false(1,2))
                obj.deltaback = delta;
            else
                error('delta must be a strictly positive numeric value!')
            end
        end
        
        function mni = getmni(obj)
            mni = obj.mniback;
        end
        
        function setmni(obj,mni)
            if obj.scalar(mni,false) && obj.integer(mni) && obj.between(mni,[0 Inf],false(1,2))
                obj.mniback = mni;
            else
                error('mni must be a strictly positive integer value!')
            end
        end
        
        function nparm = getnparm(obj)
            nparm = obj.nparmback;
        end
        
        function setnparm(obj,nparm)
            if obj.scalar(nparm,true) && obj.integer(nparm)
                obj.nparmback = nparm;
            else
                error('nparm must be an integer value!')
            end
        end
        
        function wseed = getwseed(obj)
            wseed = obj.wseedback;
        end
        
        function setwseed(obj,wseed)
            if obj.scalar(wseed,true) && isnumeric(wseed)
                obj.wseedback = wseed;
            else
                error('wseed must be a numeric value!')
            end
        end
        
        function accl = getaccl(obj)
            accl = obj.acclback;
        end
        
        function setaccl(obj,accl)
            if obj.scalar(accl,true) && isnumeric(accl) && obj.between(accl,[0 1],[false true])
                obj.acclback = accl;
            else
                error('accl must be a numeric value between 0 and 1 (0 not included)!')
            end
        end
        
        function mx = getmex(obj)
            mx = obj.mexback;
        end
        
        function setmex(obj,mx)
            if obj.scalar(mx,false) && obj.logical(mx)
                obj.mexback = logical(mx);
            else
                error('mex must be true or false!')
            end
        end
                
    end
    
end