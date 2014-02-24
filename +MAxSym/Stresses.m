classdef Stresses < MAxSym.Check
% class to define stresses for one stress period in axi-symmetric model MAxSym  
%
% remark:
%  singleton dimensions are allowed for stress matrices
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

    
    
    properties (SetAccess = protected)
        nz % total number of layers
        nr % total number of rings
    end
    
    properties (Dependent) 
        q  % discharges [L³/T], q > 0 means extraction of water - nz x nr matrix - optional, q = 0 for all rings by default 
        s0 % instantaneous head change [L], s0 > 0 means head rise - nz x nr matrix - optional, s0 = 0 for all rings by default
    end

    properties (Access = protected)
        qback = 0;  % backing field for q
        s0back = 0; % backing field for s0
    end
    
    methods
        
        function obj = Stresses(nz,nr)
        % create Stresses object
        %
        % syntax: obj = Stresses(nz,nr)
        %  nz is number of layers
        %  nr is number of rings

            if nargin == 2
                % nr
                if obj.scalar(nr,false) && obj.integer(nr) && obj.between(nr,[0 Inf],false(1,2))
                    obj.nr = nr;
                else
                    error('Input nr must be a strictly positive integer value!')
                end
                % nz
                if obj.scalar(nz,false) && obj.integer(nz) && obj.between(nz,[0 Inf],false(1,2))
                    obj.nz = nz;
                else
                    error('Input nz must be a strictly positive integer value!')
                end
            elseif nargin
                error('Wrong number of input arguments!')
            end
        end
        
        function set.q(obj,q)
            obj.setq(q);
        end
        
        function q = get.q(obj)
            q = obj.getq;
        end
        
        function set.s0(obj,s0)
            obj.sets0(s0);
        end
        
        function s0 = get.s0(obj)
            s0 = obj.gets0;
        end
        
        function tf = isdefined(obj)
        % check if any stress is defined
        % 
        % syntax: tf = isdefined(obj)
        %  tf is true or false
            tf = ~(isempty(obj.q) & isempty(obj.s0));
            if tf
                tf = any(obj.q(:)) | any(obj.s0(:));
            end
        end
        
    end
    
    methods (Access = protected)
        
        function setq(obj,q)
            if obj.size(q,[obj.nz,obj.nr],true,true) && isnumeric(q) && ~any(isnan(q(:))) && obj.between(q,[-Inf Inf],[false false])
                if isempty(q)
                    obj.qback = 0;
                else
                    obj.qback = q;
                end
            else
                error(['''q'' must be numeric matrix of size ' mat2str([obj.nz obj.nr]) ' with finite non NaN elements!'])
            end
        end
        
        function q = getq(obj)
            q = obj.qback;
        end
        
        function sets0(obj,s0)
            if obj.size(s0,[obj.nz,obj.nr],true,true) && isnumeric(s0) && ~any(isnan(s0(:))) && obj.between(s0,[-Inf Inf],[false false])
                if isempty(s0)
                    obj.s0back = 0;
                else
                    obj.s0back = s0;
                end
            else
                error(['''s0'' must be numeric matrix of size ' mat2str([obj.nz obj.nr]) ' with finite non NaN elements!'])
            end
        end
        
        function s0 = gets0(obj)
            s0 = obj.s0back;
        end
        
    end
    
end