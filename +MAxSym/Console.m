classdef Console < dynamicprops & MAxSym.Model

    
    
    properties (Access = protected)
        hwindow__
        hinput__
        houtput__
        height__
        print__
    end
    
    methods
        
        function obj = Console
            obj.hwindow__ = figure('resizefcn',@(~,~) obj.resize,'handlevisib','off',...
                                   'numbertitle','off','name','MAxSym Console','menu','no','toolbar','no');
            p = get(obj.hwindow__,'position');
            obj.height__ = 25;
            lst = {'parent',obj.hwindow__,'style','edit','fontname','lucida console','fontsize',10,'horizontalalign','left',...
                   'backgroundcolor','k','foregroundcolor','w','string',{}};
            obj.hinput__ = uicontrol(lst{:},'position',[0, p(4)-obj.height__, p(3), obj.height__],'callback',@(~,~) obj.evaluate);
            obj.houtput__ = uicontrol(lst{:},'max',2,'position',[0, 0, p(3), p(4)-obj.height__],'enable','inactive');
            obj.print__ = true;
        end
        
        function clc(obj)
            set(obj.houtput__,'string','')
            obj.print__ = false;
        end
        
        function addprop(obj,propname)
            addprop@dynamicprops(obj,propname);
        end
        
        function script(obj,fname)
            if nargin == 1
                [fname,pname,id] = uigetfile('*.*','Open script...');
                if id
                    fname = strcat(pname,filesep,fname);
                else
                    disp('Cancelled...')
                end
            end
            fid = fopen(fname,'r');
            c = textscan(fid,'%s','delimiter',char(10));
            fclose(fid);
            for cmd = c{1}(:)'
                cmd = strtrim(cmd{1});
                i = strfind(cmd,'%');
                if ~isempty(i)
                    cmd(i(1):end) = '';
                end
                if ~isempty(cmd)
                    obj.eval(cmd);
                end
            end
            obj.print__ = false;
        end
        
        function setforegroundcolor(obj,col)
            if nargin == 1
                col = uisetcolor;
            end
            if isscalar(col) && col == 0
                disp('Cancelled...')
            else
                set([obj.hinput__,obj.houtput__],'foregroundcolor',col)
            end
        end
        
        function setbackgroundcolor(obj,col)
            if nargin == 1
                col = uisetcolor;
            end
            if isscalar(col) && col == 0
                disp('Cancelled...')
            else
                set([obj.hinput__,obj.houtput__],'backgroundcolor',col)
            end
        end
        
        function setfont(obj,fontsct)
            if nargin == 1
                fontsct = uisetfont(obj.hinput__);
            end
            set([obj.hinput__,obj.houtput__],fontsct)
        end
        
    end
    
    methods (Access = protected)
        
        function resize(obj)
            p = get(obj.hwindow__,'position');
            set(obj.hinput__,'position',[0 p(4)-obj.height__ p(3) obj.height__])
            set(obj.houtput__,'position',[0 0 p(3) p(4)-obj.height__])
        end
    
        function evaluate(m)
            cmd = char(get(m.hinput__,'string'));
            set(m.hinput__,'string',{});
            m.eval(cmd);
        end
        
        function eval(m,cmd)
            try
                out = evalc(cmd);
                if ~isempty(out)
                    out = cellstr(out);
                else
                    out = {};
                end
            catch err
                out = cellstr(err.message);
            end
            out = [out(:);{''}];
            cmd = ['>> ' cmd];
            str = cellstr(get(m.houtput__,'string'));
            str = [{cmd};out;str];
            if m.print__
                set(m.houtput__,'string',str)
            else
                m.print__ = true;
            end
        end
        
    end
    
end