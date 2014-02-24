classdef Grid < MAxSym.Check
% class to define grid in axi-symmetric model MAxSym  
%
% remark:
%  - rings are counted from inner to outer grid boundary
%  - layers are counted from top to bottom
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
        confined % aquifer system is confined? true or false
        rb  % radii of ring boundaries [L] - 1 x (nr+1) vector
        r   % radii of nodal circles [L] - 1 x nr vector
        dr  % ring widths [L] - 1 x nr vector
        nr  % number of rings
        zb  % levels of layer boundaries [L] - (nz+1) x 1 vector
        z   % levels of layer middles [L] - nz x 1 vector
        D   % layer thicknesses [L] - nz x 1 vector
        nz  % number of layers
        hs  % horizontal ring surface areas [L²] - 1 x nr vector
        vol % ring volumes [L³] - nz x nr matrix
    end
    
    properties (Access = protected)
        confinedback % backing field for confined
        rbback % backing field for rb
        zbback % backing field for zb
    end
    
    methods
        
        function obj = Grid(rb,D,confined)
        % create Grid object
        %
        % syntax: obj = Grid(rb,D,confined)
        %  rb is vector with radii of ring boundaries [L]
        %   rb must be strictly monotonically increasing
        %   rb must contain finite, strictly positive elements
        %   rb must contain two elements at least
        %  D is vector with layer thickesses [L]
        %   D must contain finite, strictly positive elements
        %   D must contain one element at least
        %   D(1) is top layer thickness, D(end) is bottom layer thickness
        %  confined is logical indicating whether the aquifer system is confined (true) or not (false)
        %  obj is Grid object
        
            if nargin == 3
                % rb
                if obj.vector(rb,false,false) && isnumeric(rb) && obj.between(rb,[0 Inf],false(1,2))
                    if all(diff(rb) >= obj.eps)
                        obj.rbback = rb(:)';
                    else
                        error('Input rb must be strictly monotonically increasing')
                    end
                else
                    error('Input rb must be non-empty numeric vector with finite and strictly positive elements (two at least)!')
                end
                % zb
                if obj.vector(D,true,false) && isnumeric(D) && obj.between(D,[0 Inf],false(1,2)) 
                    obj.zbback = flipud(cumsum([0; flipud(D(:))]));
                else
                    error('Input D must be non-empty numeric vector with finite and strictly positive elements!')
                end
                % confined
                if obj.scalar(confined,false) && obj.logical(confined)
                    obj.confinedback = logical(confined);
                else
                    error('Input ''confined'' must be logical scalar!')
                end
            elseif nargin
                error('Wrong number of input arguments!')
            end
        end
        
        function tf = get.confined(obj)
            tf = obj.getconfined;
        end
        
        function rb = get.rb(obj)
            rb = obj.getrb;
        end
        
        function r = get.r(obj)
            r = obj.getr;
        end
        
        function dr = get.dr(obj)
            dr = obj.getdr;
        end
        
        function nr = get.nr(obj)
            nr = obj.getnr;
        end
        
        function zb = get.zb(obj)
            zb = obj.getzb;
        end
        
        function D = get.D(obj)
            D = obj.getD;
        end
        
        function z = get.z(obj)
            z = obj.getz;
        end
        
        function nz = get.nz(obj)
            nz = obj.getnz;
        end
        
        function hs = get.hs(obj)
            hs = obj.geths;
        end
        
        function vol = get.vol(obj)
            vol = obj.getvol;
        end
        
    end
    
    methods (Access = protected)
        
        function tf = getconfined(obj)
            tf = obj.confinedback;
        end
        
        function rb = getrb(obj)
            rb = obj.rbback;
        end
        
        function r = getr(obj)
            r = obj.rb;
            r = sqrt(r(1:end-1).*r(2:end));
        end
        
        function dr = getdr(obj)
            dr = diff(obj.rb);
        end
        
        function nr = getnr(obj)
            nr = length(obj.rb) - 1;
        end
        
        function zb = getzb(obj)
            zb = obj.zbback;
        end
        
        function z = getz(obj)
            z = obj.zb;
            z = (z(1:end-1) + z(2:end)) / 2;
        end
        
        function D = getD(obj)
            D = -diff(obj.zb);
        end
        
        function nz = getnz(obj)
            nz = length(obj.zb) - 1;
        end
        
        function hs = geths(obj)
            hs = obj.rb;
            hs = pi * (hs(2:end).^2 - hs(1:end-1).^2);
        end
        
        function vol = getvol(obj)
            vol = obj.D * obj.hs;
        end
        
    end
    
end