classdef Parameters < MAxSym.Check
% class to define parameters for axi-symmetric model MAxSym  
%
% remark:
%  singleton dimensions are allowed for parameter matrices
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
        confined % aquifer system is confined? true or false
        steady % simulation is steady state? true or false
        nz % total number of layers
        nr % total number of rings
    end
       
    properties (Dependent)
        inactive % inactive rings - nz x nr logical matrix - optional, all rings are active by default
        constant % constant head rings - nz x nr logical matrix - optional, all rings have variable head by default
        kr % radial component of hydraulic conductivity [L/T] - nz x nr matrix, NaN values are ignored - required if not all cr are defined
        cr % radial resistance between rings [T] - nz x (nr-1) matrix, NaN values and values <= 0 are ignored - required if kr is not set
        kz % vertical component of hydraulic conductivity [L/T] - nz x nr matrix, NaN values are ignored - required if nz > 1 and not all cz are defined
        cz % vertical resistance between layers [T] - (nz-1) x nr matrix, NaN values and values <= 0 are ignored - required if nz > 1 and kz is not set
        ss % specific elastic storage coefficient [1/L] - nz x nr matrix, NaN values are ignored - required if transient simulation
        sy % specific yield [-] - nr vector, NaN values are ignored - required if transient simulation and aquifer system is not confined
    end
    
    properties (Access = protected)
        inactiveback = false; % backing field for inactive
        constantback = false; % backing field for constant
        krback % backing field for kr
        crback % backing field for cr
        kzback % backing field for kz
        czback % backing field for cz
        ssback % backing field for ss
        syback % backing field for sy
    end
    
    methods
        
        function obj = Parameters(nz,nr,confined,steady)
        % create Parameters object
        %
        % syntax: obj = Parameters(nz,nr,confined,steady)
        %  nz is number of layers
        %  nr is number of rings
        %  confined is logical indicating whether the aquifer system is confined (true) or unconfined (false)
        %  steady is logical indicating whether the simulation is steady state (true) or transient (false)
        
            if nargin == 4
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
                % confined
                if obj.scalar(confined,false) && obj.logical(confined)
                    obj.confined = confined;
                else
                    error('Input ''confined'' must be true or false!')
                end
                % steady
                if obj.scalar(steady,false) && obj.logical(steady)
                    obj.steady = steady;
                else
                    error('Input ''steady'' must be true or false!')
                end
            elseif nargin
                error('Wrong number of input arguments!')
            end
        end
        
        function set.inactive(obj,inactive)
            obj.setinactive(inactive);
        end
        
        function inactive = get.inactive(obj)
            inactive = obj.getinactive;
        end
        
        function set.constant(obj,constant)
            obj.setconstant(constant);
        end
        
        function constant = get.constant(obj)
            constant = obj.getconstant;
        end
        
        function set.kr(obj,kr)
            obj.setkr(kr);
        end
        
        function kr = get.kr(obj)
            kr = obj.getkr;
        end
        
        function set.cr(obj,cr)
            obj.setcr(cr);
        end
        
        function cr = get.cr(obj)
            cr = obj.getcr;
        end
        
        function set.kz(obj,kz)
            obj.setkz(kz);
        end
        
        function kz = get.kz(obj)
            kz = obj.getkz;
        end
        
        function set.cz(obj,cz)
            obj.setcz(cz);
        end
        
        function cz = get.cz(obj)
            cz = obj.getcz;
        end
        
        function set.ss(obj,ss)
            obj.setss(ss);
        end
        
        function ss = get.ss(obj)
            ss = obj.getss;
        end
        
        function set.sy(obj,sy)
            obj.setsy(sy);
        end
        
        function sy = get.sy(obj)
            sy = obj.getsy;
        end
        
        function tf = isdefined(obj)
        % check if all required parameters are defined
        %
        % syntax: tf = isdefined(obj)
        %  tf is true or false
        
            tf = ~(isempty(obj.kr) & isempty(obj.cr));
            if tf && obj.nz > 1
                tf = ~(isempty(obj.kz) & isempty(obj.cz));
            end
            if tf && ~obj.steady
                tf = ~isempty(obj.ss);
                if tf && ~obj.confined
                    tf = ~isempty(obj.sy);
                end
            end
        end
        
    end
    
    methods (Access = protected)
        
        function setinactive(obj,inactive)
            if obj.size(inactive,[obj.nz,obj.nr],true,true) && obj.logical(inactive)
                if isempty(inactive)
                    obj.inactiveback = false;
                elseif ~all(inactive(:))
                    obj.inactiveback = inactive;
                else
                    error('setting all rings inactive is useless!')
                end
            else
                error(['''inactive'' must be logical matrix of size ' mat2str([obj.nz obj.nr]) '!'])
            end
        end
        
        function inactive = getinactive(obj)
            inactive = obj.inactiveback;
        end
        
        function setconstant(obj,constant)
            if obj.size(constant,[obj.nz,obj.nr],true,true) && obj.logical(constant)
                if isempty(constant)
                    obj.constantback = false;
                elseif ~all(constant(:))
                    obj.constantback = constant;
                else
                    error('setting all rings constant is useless!')
                end
            else
                error(['''constant'' must be logical matrix of size ' mat2str([obj.nz obj.nr]) '!'])
            end
        end
        
        function constant = getconstant(obj)
            constant = obj.constantback;
        end
                
        function setkr(obj,kr)
            if obj.size(kr,[obj.nz,obj.nr],true,true) && isnumeric(kr) && obj.between(kr(~isnan(kr)),[0 Inf],[true false])
                obj.krback = kr;
            else
                error(['''kr'' must be numeric matrix of size ' mat2str([obj.nz obj.nr]) ' with positive elements!'])
            end
        end
        
        function kr = getkr(obj)
            kr = obj.krback;
        end
        
        function setcr(obj,cr)
            if obj.size(cr,[obj.nz,obj.nr-1],true,true) && isnumeric(cr)
                obj.crback = cr;
            else
                error(['''cr'' must be numeric matrix of size ' mat2str([obj.nz obj.nr-1]) '!'])
            end
        end
        
        function cr = getcr(obj)
            cr = obj.crback;
        end
                
        function setkz(obj,kz)
            if obj.nz == 1 && ~isempty(kz)
                error('One layer only! kz is not needed!')
            elseif obj.size(kz,[obj.nz,obj.nr],true,true) && isnumeric(kz) && obj.between(kz(~isnan(kz)),[0 Inf],[true false])
                obj.kzback = kz;
            else
                error(['''kz'' must be numeric matrix of size ' mat2str([obj.nz obj.nr]) ' with positive elements!'])
            end
        end
        
        function kz = getkz(obj)
            kz = obj.kzback;
        end
        
        function setcz(obj,cz)
            if obj.nz == 1 && ~isempty(cz)
                error('One layer only! cz is not needed!')
            elseif obj.size(cz,[obj.nz-1,obj.nr],true,true) && isnumeric(cz)
                obj.czback = cz;
            else
                error(['''cz'' must be numeric matrix of size ' mat2str([obj.nz-1 obj.nr]) '!'])
            end
        end
        
        function cz = getcz(obj)
            cz = obj.czback;
        end
                
        function setss(obj,ss)
            if obj.steady && ~isempty(ss)
                error('Steady state simulation! ss is not needed!')
            elseif obj.size(ss,[obj.nz,obj.nr],true,true) && isnumeric(ss) && obj.between(ss(~isnan(ss)),[0 Inf],[true false])
                obj.ssback = ss;
            else
                error(['''ss'' must be numeric matrix of size ' mat2str([obj.nz obj.nr]) ' with positive elements!'])
            end
        end
        
        function ss = getss(obj)
            ss = obj.ssback;
        end
        
        function setsy(obj,sy)
            if obj.steady && ~isempty(sy)
                error('Steady state simulation! sy is not needed!')
            elseif obj.confined && ~isempty(sy)
                error('Confined aquifer system! sy is not needed!')
            elseif obj.length(sy,obj.nr,true,true) && isnumeric(sy) && obj.between(sy(~isnan(sy)),[0 Inf],[true false])
                obj.syback = sy(:)';
            else
                error(['''sy'' must be numeric vector of length ' int2str(obj.nr) ' with positive elements!'])
            end
        end
        
        function sy = getsy(obj)
            sy = obj.syback;
        end
        
    end
    
end
    