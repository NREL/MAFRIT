classdef st1a < Components
  % based on the block diagram for the ST1A model given in IEEE Std
  % 421.5-2005.
    properties
        tr;
        vimax;
        vimin;
        tc;
        tb;
        ka;
        ta;
        vrmax
        vrmin
        kc;
        kf;
        tf;
        tc1;
        tb1;
        vamax
        vamin
        ilr;
        klr;
        extStates
        extOutputLocation;
        stWrite;
    end
    
    methods
        function ext=st1a()
        end
        
        function ext=AssignParams(ext)
            ext.tr=ext.compParams(1);
            ext.vimax=ext.compParams(2);
            ext.vimin=ext.compParams(3);
            ext.tc=ext.compParams(4);
            ext.tb=ext.compParams(5);
            ext.ka=ext.compParams(6);
            ext.ta=ext.compParams(7);
            ext.vrmax=ext.compParams(8);
            ext.vrmin=ext.compParams(9);
            ext.kc=ext.compParams(10);
            ext.kf=ext.compParams(11);
            ext.tf=ext.compParams(12);
            ext.tc1=ext.compParams(13);
            ext.tb1=ext.compParams(14);
            ext.vamax=ext.compParams(15);
            ext.vamin=ext.compParams(16);
            ext.ilr=ext.compParams(17);
            ext.klr=ext.compParams(18);
            ext.stWrite=fopen(ext.extOutputLocation,'w');
            if ext.stWrite<0
                sprintf('error opening file')
            else
                fprintf(ext.stWrite,'%8s %8s %8s \r\n',['State 1' 'State 2' 'State 3']);
            end  
        end
        %% Writing exciter states    
        function ext=WriteData(ext)
                fprintf(ext.stWrite,'%8.4f %8.4f %8.4f \r\n',[ext.extStates(1) ext.extStates(2) ext.extStates(3)]);  
        end
        %%
        function [ext,vref]=Init(ext,efd,vt)
            u=efd;
            x=u/ext.ka;
            if (ext.tb~=0 && ext.tc~=0)
                y1=(ext.tb-ext.tc)/ext.tb*x;
                z=(x-y1)*(ext.tb/ext.tc);
            else
                y1=0;
                z=x;
            end
            y=vt;
            vref=z+y;
            ext.extStates=[y;y1;u]; % 3 X 1 vector of initial conditions
        end
        %%
        function [ext,efd]=newStates(ext,tspan,vref,ifd,vpss,vt1)
            options=odeset('RelTol',1e-9);
            [T,X]=ode23tb(@(t,x) ext_dyn(t,x,vref,vpss,vt1,ext.tr,ext.vimax,ext.vimin,ext.tc,ext.tb,ext.ka,ext.ta,ext.kf,ext.tf,ext.tc1,ext.tb1,ext.vamax,ext.vamin,ext.ilr,ext.klr),tspan,ext.extStates,options);
            
            function dy=ext_dyn(t1,y,vref,vpss,vt,tr,vimax,vimin,tc,tb,ka,ta,kf,tf,tc1,tb1,vamax,vamin,ilr,klr)
                 dy=zeros(3,1);     
                % Transducer Dynamics
                
                if (tr~=0)
                    dy(1)=-y(1)/tr+vt/tr;
                else
                    dy(1)=0;
                end
                % Wind-Up Limiter
                
                if (vref-y(1)>=vimax)
                    z=vimax;
                elseif (vref-y(1)<=vimin)
                    z=vimin;
                else
                    z=vref-y(1)+vpss;
                end
                
                % Loop Shaping Transfer Function Dynamics to Reduce Transient Gain
                
                if (tb~=0 && tc~=0)
                    dy(2)=-y(2)/tb+(1-tc/tb)*z/tb;
                    xx=(tc/tb)*z+y(2);
                else
                    dy(2)=0;
                    xx=z;
                end
                
                % Regulator + Static Exciter Dynamics
                
                dy(3)=-y(3)/ta+xx*ka/ta;
                
                if (y(3)>=vamax && dy(3)>0)
                    dy(3)=0;
                elseif (y(3)<=vamin && dy(3)<0)
                    dy(3)=0;
                end
            end
            
            ext.extStates=X(end,:);
            
            % Calculation of field voltage Efd to be used by the generator model as input
            
            efd=ext.extStates(3);
            if (efd> ext.vrmax*vt1-ext.kc*ifd)
                efd=ext.vrmax*vt1-ext.kc*ifd;
            elseif (efd< ext.vrmin*vt1)
                efd=ext.vrmin*vt1;
            end
        end
        
    end
    
end

