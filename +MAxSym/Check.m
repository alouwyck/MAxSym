classdef Check < handle
% static class to perform variable checking
% all methods are protected, so methods are accessible only to subclasses of Check
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



    properties (GetAccess = protected, Constant)
        eps = 1e-15; % maximum absolute difference used to test equality of floating point numbers
    end
    
    methods (Static, Access = protected)
        
        function ok = scalar(x,allowempty)
        % check if array is scalar
        %
        % syntax: ok = scalar(x,allowempty)
        %  x is array to check
        %  allowempty (optional) is logical indicating whether value may be empty (true) or not (false, default) 
            if nargin < 2
                allowempty = false;
            end
            ok = isempty(x) & allowempty;
            if ~(ok || isempty(x))
                ok = isscalar(x);
            end
        end
        
        function ok = vector(v,allowscalar,allowempty)
        % check if array is vector
        %
        % syntax: ok = vector(v,allowscalar,allowempty)
        %  v is array to check
        %  allowscalar (optional) is logical indicating whether vector may be scalar (true) or not (false, default)
        %  allowempty (optional) is logical indicating whether vector may be empty (true) or not (false, default)
            if nargin < 2
                allowscalar = false;
            end
            if nargin < 3
                allowempty = false;
            end
            ok = isempty(v) & allowempty;
            if ~(ok || isempty(v))
                ok = isscalar(v) & allowscalar;
                if ~(ok || isscalar(v))
                    ok = isvector(v);
                end
            end
        end
        
        function ok = length(v,len,allowscalar,allowempty)
        % check vector length
        %
        % syntax: ok = length(v,len,allowscalar,allowempty)
        %  v is vector to check
        %  len is required vector length
        %  allowscalar (optional) is logical indicating whether vector may be scalar (true) or not (false, default)
        %  allowempty (optional) is logical indicating whether vector may be empty (true) or not (false, default)
            if nargin < 3
                allowscalar = false;
            end
            if nargin < 4
                allowempty = false;
            end
            ok = isempty(v) & (allowempty | ~len);
            if ~(ok || isempty(v))
                ok = isscalar(v) & (allowscalar | len==1);
                if ~(ok || isscalar(v))
                    ok = isvector(v) & length(v) == len;
                end
            end
        end
        
        function ok = size(a,siz,allowsingleton,allowempty)
        % check array size
        %
        % syntax: ok = size(a,siz,allowsingleton,allowempty)
        %  a is array to check
        %  siz is required array size
        %  allowsingleton (optional) is logical indicating whether array may have singleton dimensions (true) or not (false, default)
        %  allowempty (optional) is logical indicating whether array may be empty (true) or not (false, default)
            if nargin < 3
                allowsingleton = false;
            end
            if nargin < 4
                allowempty = false;
            end
            ok = isempty(a) & (allowempty | any(siz==0));
            if ~(ok || isempty(a))
                s = size(a);
                m = length(s);
                n = length(siz);
                ok = m == n | (allowsingleton & m < n);
                if ok
                    s = [s ones(1,n-m)];
                    ok = all(s == siz | (allowsingleton & s == 1));
                end
            end
        end
        
        function ok = integer(x)
        % check if array has integer elements
        %
        % syntax: ok = integer(x)
            ok = isempty(x);
            if ~ok
                ok = isnumeric(x);
                if ok
                    ok = all(MAxSym.Check.equal(rem(x(:),1),0));
                end
            end
        end
        
        function ok = logical(tf)
        % check if array is logical
        % numeric array is also allowed
        %
        % syntax: ok = logical(tf)
            ok = isempty(tf);
            if ~ok
                ok = islogical(tf) | isnumeric(tf);
            end
        end
        
        function ok = between(x,xlim,incl)
        % check if array elements are inside given interval
        %
        % syntax: ok = between(x,xlim,incl)
        %  xlim is vector [xmin xmax] indicating the interval boundaries
        %  incl is logical vector [includexmin includexmax] 
        %   indicating whether interval boundaries must be included or not
            ok = isempty(x);
            if ~ok
                b = x(:) > xlim(1) & x(:) < xlim(2);
                ok = all(b);
                if ~ok && any(incl)
                    if incl(1)
                        b(~b) = MAxSym.Check.equal(x(~b),xlim(1));
                        ok = all(b);
                    end
                    if ~ok && incl(2)
                        b(~b) = MAxSym.Check.equal(x(~b),xlim(2));
                        ok = all(b);
                    end
                end
            end
        end
        
        function m = repmat(m,siz)
        % replicate input array to get full size array
        % 
        % syntax: f = repmat(a,siz)
        %  a is array which may contain singleton dimensions
        %  siz is required full array size
        %  f is full array in which singleton dimensions are repeated according to given size
        %
        % remark: if a is empty, an empty array f will be returned 
            if ~isempty(m)
                s = size(m);
                s = [s ones(1,length(siz)-length(s))];
                b = s==1 & siz>1;
                s(b) = siz(b);
                s(~b) = 1;
                m = repmat(m,s);
            end
        end
        
        function c = repcellelements(c,siz)
        % get cell with full size array elements (where each element has the same size)
        % 
        % syntax: cout = repcellelements(cin,siz)
        %  cout{i} = repmat(cin{i},siz)
        %
        % see Check.repmat
            for n = 1:numel(c)
                c{n} = MAxSym.Check.repmat(c{n},siz);
            end
        end
        
        function tf = equal(x,y)
        % check if elements of x are equal to corresponding elements of y
        %
        % syntax: tf = equal(x,y)
        %  x and y are numeric arrays of same size
        %  one of both arrays may be scalar
        %
        % equal(x,y) is equivalent to x == y but is less strict
        % since equal(x,y) checks if the absolute difference between 
        % x and y is less than property eps = 1e-15
        % equal(x,y) returns true if this is the case
        %
        % e.g. equal(1e-16,1e-17) returns true
        %      1e-16 == 1e-17 returns false
            tf = abs(x-y) < MAxSym.Check.eps;
        end
        
        function tf = equalsize(x,y,allowsingleton,allowempty)
        % check if x and y have equal size
        %
        % syntax: tf = equalsize(x,y,allowsingleton,allowempty)
        %  x and y are arrays
        %  allowsingleton (optional) is logical indicating whether arrays may have singleton dimensions (true) or not (false, default)
        %   singleton dimensions are not considered if allowsingleton is true 
        %  allowempty (optional) is logical indicating whether arrays may be empty (true) or not (false, default)
        %   equalsize always returns true if x or y is empty and allowempty is true
            if nargin < 3
                allowsingleton = false;
            end
            if nargin < 4
                allowempty = false;
            end
            b = isempty(x) | isempty(y);
            tf = allowempty & b;
            if ~(tf || b)
                x = size(x);
                y = size(y);
                tf = length(x) == length(y);
                if tf
                    tf = all(x == y);
                    if ~tf && allowsingleton
                        b = x > 1 & y > 1;
                        tf = all(x(b) == y(b));
                    end
                elseif allowsingleton
                    if length(x) < length(y)
                        x = [x, ones(1,length(y)-length(x))];
                    else
                        y = [y, ones(1,length(x)-length(y))];
                    end
                    b = x > 1 & y > 1;
                    tf = all(x(b) == y(b));
                end
            end
        end
        
    end
    
end
