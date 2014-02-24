classdef Excel < handle
    
    properties (SetAccess = protected)
        filename
        sheets
        model
    end
    
    properties (Access = protected)
        steady
        confined
        s
        qr
        qz
        qs
        bud
        totbud
        niter
    end
    
    methods
        
        function obj = Excel(varargin)
            obj.setfilename(varargin{:});
            if ~isempty(obj.filename)
                obj.setmodel;
                obj.settime;
                obj.setgrid;
                obj.setpar;
                obj.setstress;
                obj.setsolver;
                obj.model.run;
                obj.setoutput;
                obj.setobs;
            end
        end
        
    end
    
    methods (Access = protected)
        
        function setfilename(obj,filename)
            if nargin == 1
                [filename, pathname] = uigetfile({'*.xls;*.xlsx;*.xlsm;*.xlsb;*.ods','MicroSoft Excel Spreadsheet (*.xls, *.xlsx, *.xlsm, *.xlsb, *.ods)';...
                                                  '*.*','All Files (*.*)'},'Open Spreadsheet File');
                if isequal(filename,0)
                    disp('Cancelled...')
                    return
                else
                    filename = fullfile(pathname,filename);
                end
            else
                [pathname,filename,ext] = fileparts(filename);
                if isempty(pathname)
                    pathname = pwd;
                end
                filename = fullfile(pathname,[filename ext]);
            end
            [filetype,sheets] = xlsfinfo(filename);
            if isempty(filetype)
                error('Selected file is not a MicroSoft Excel Spreadsheet!')
            else
                obj.filename = filename;
                obj.sheets = sheets;
            end
        end
        
        function setmodel(obj)
            obj.model = MAxSym.Model;
            num = obj.readsheet('model',true);
            num(isnan(num)) = [];
            obj.steady = logical(num(1));
            obj.confined = logical(num(2));
            obj.s = num(3);
            obj.qr = num(4);
            obj.qz = num(5);
            obj.qs = num(6);
            obj.bud = num(7);
            obj.totbud = num(8);
            obj.niter = num(9);            
        end
        
        function settime(obj)
            if obj.steady
                obj.model.settime(0);
            else
                dt = obj.readsheet('dt',true);
                if isvector(dt)
                    obj.model.settime(dt);
                else
                    [~,ndt] = unique(dt(:,2),'last');
                    obj.model.settime(dt(:,1),[ndt(1);diff(ndt(:))]);
                end
            end
        end
        
        function setgrid(obj)
            rb = obj.readsheet('rb',true);
            D = obj.readsheet('D',true);
            obj.model.setgrid(rb(:)',D(:),obj.confined);
        end
        
        function setpar(obj)
            obj.setparmatrix('kr',false);
            obj.setparmatrix('cr',false);
            if obj.model.grid.nz > 1
                obj.setparmatrix('kz',false);
                obj.setparmatrix('cz',false);
            end
            if ~obj.model.time.steady 
                obj.setparmatrix('ss',true);
                if ~obj.model.grid.confined
                    obj.setparmatrix('sy',true);
                end
            end
            obj.setparmatrix('inactive',false,true);
            obj.setparmatrix('constant',false,true);
        end
        
        function setparmatrix(obj,par,required,isbool)
            m = obj.readsheet(par,required);
            if ~isempty(m)
                if nargin == 4 && isbool
                    m = logical(m);
                end
                obj.model.par.(par) = m;
            end
        end
        
        function setstress(obj)
            obj.setstressmatrices('q');
            obj.setstressmatrices('s0');
        end
        
        function setstressmatrices(obj,stress)
            m = obj.readsheet(stress,false);
            if ~isempty(m)
                n = find(isnan(m(:,1)));
                im1 = n + 1;
                im2 = n - 1;
                im1(ismember(im1,n)) = [];
                im2(ismember(im2,n)) = [];
                im1 = [1; im1];
                im2 = [im2; size(m,1)];
                for n = 1:length(im1)
                    i = m(im1(n),1);
                    s = m(im1(n)+1:im2(n),:);
                    b = ~all(isnan(s),1);
                    obj.model.stress(i).(stress) = s(:,b);
                end
            end
        end
        
        function setsolver(obj)
            num = obj.readsheet('solver',true);
            if length(num) < 2
                error('Not all solver parameters are given!')
            elseif length(num) == 2 || isnan(num(3))
                obj.model.setsolver(num(1),num(2));
            else
                obj.model.setsolver(num(1),num(2),num(3));                
            end
            if length(num) > 3 && ~isnan(num(4))
                obj.model.solver.wseed = num(4);
            end
            if length(num) > 4 && ~isnan(num(5))
                obj.model.solver.wseed = num(5);
            end
        end
        
        function setoutput(obj)
            if obj.s
                obj.array2sheet('s','s');
            end
            if obj.qr
                obj.array2sheet('qr','qr');
            end
            if obj.qz && obj.model.grid.nz > 1
                obj.array2sheet('qz','qz');
            end
            if obj.qs && obj.steady
                obj.array2sheet('qs','qs');
            end
            if obj.bud
                obj.array2sheet('bud','bud');
            end
            if obj.totbud
                obj.array2sheet('totbud','totbud');
            end
            if obj.niter
                obj.array2sheet('niter','niter');
            end
        end
        
        function setobs(obj)

            if obj.steady
                % s vs r
                c = obj.readsheet('s_vs_r',false,true);
                [irow,icol] = find(strcmpi(c,'ilay'));
                for n = 1:length(irow)
                    [i,j] = deal(irow(n),icol(n));
                    ilay = c{i,j+1};
                    r = cell2mat(c((i+3):end,j));
                    r(isnan(r)) = [];
                    s = obj.model.interp1r(ilay,r);
                    c(i+2+(1:length(s)),j+1) = num2cell(s(:));
                end
                xlswrite(obj.filename,c,'s_vs_r');
                
            else
                
                % s vs t
                c = obj.readsheet('s_vs_t',false,true);
                [irow,icol] = find(strcmpi(c,'ilay'));
                for n = 1:length(irow)
                    [i,j] = deal(irow(n),icol(n));
                    ilay = c{i,j+1};
                    r = c{i+1,j+1};
                    t = cell2mat(c((i+3):end,j));
                    t(isnan(t)) = [];
                    s = obj.model.interp2(ilay,r,t);
                    c(i+2+(1:length(s)),j+1) = num2cell(s(:));
                end
                xlswrite(obj.filename,c,'s_vs_t');
                
                % s vs r
                c = obj.readsheet('s_vs_r',false,true);
                [irow,icol] = find(strcmpi(c,'ilay'));
                for n = 1:length(irow)
                    [i,j] = deal(irow(n),icol(n));
                    ilay = c{i,j+1};
                    t = c{i+1,j+1};
                    r = cell2mat(c((i+3):end,j));
                    r(isnan(r)) = [];
                    s = obj.model.interp2(ilay,r,t);
                    c(i+2+(1:length(s)),j+1) = num2cell(s(:));
                end
                xlswrite(obj.filename,c,'s_vs_r');

            end
            
        end
        
        function num = readsheet(obj,sheetname,required,getraw)
            b = strcmpi(obj.sheets,sheetname);
            if any(b)
                if nargin == 3 || ~getraw
                    num = xlsread(obj.filename,obj.sheets{b});
                else
                    [~,~,num] = xlsread(obj.filename,obj.sheets{b});
                end
            elseif required
                error(['Sheet ''' sheetname ''' is not found!'])
            else
                num = [];
            end
        end
        
        function array2sheet(obj,prop,sheetname)
            a = obj.model.(prop);
            a(end+1,:,:) = NaN;
            a = reshape(permute(a,[1 3 2]),[],size(a,2));
            xlswrite(obj.filename,a(1:end-1,:),sheetname);
        end
        
    end
    
end