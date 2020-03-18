classdef repc_a < Components
    % REPC_A represents the plant controller module for Type 3 and 4 Wind
    % Turbines and Central PV Stations
    % Based on the model description given here: https://www.wecc.org/Reliability/WECC-Second-Generation-Wind-Turbine-Models-012314.pdf        
    properties
        mvab;
        refflag;
        vcompflag;
        freqflag;
        tfltr;
        vbus;
        frombus;
        tobus;
        ckt;
        rc;
        xc;
        kc;
        dbd;
        emax;
        emin;
        kp;
        ki;
        qmax;
        qmin;
        vfrz;
        tft;
        tfv;
        fdbd1;
        fdbd2;
        ddn;
        dup;
        tp;
        femax;
        femin;
        kpg;
        kig;
        pmax;
        pmin;
        tlag;
        vref;
        plantpref;
        freqref;
        qref;
        repcStates;
        repcOutputLocation;
        stWrite;
    end
    
    methods
        function repc=repc_a()
        end
        
        function repc=AssignParams(repc)
            repc.mvab=repc.compParams(1);
            repc.tfltr=repc.compParams(2);
            repc.kp=repc.compParams(3);
            repc.ki=repc.compParams(4);
            repc.tft=repc.compParams(5);
            repc.tfv=repc.compParams(6);
            repc.refflag=repc.compParams(7);
            repc.vfrz=repc.compParams(8);
            repc.rc=repc.compParams(9);
            repc.xc=repc.compParams(10);
            repc.kc=repc.compParams(11);
            repc.vcompflag=repc.compParams(12);
            repc.emax=repc.compParams(13);
            repc.emin=repc.compParams(14);
            repc.dbd=repc.compParams(15);
            repc.qmax=repc.compParams(16);
            repc.qmin=repc.compParams(17);
            repc.kpg=repc.compParams(18);
            repc.kig=repc.compParams(19);
            repc.tp=repc.compParams(20);
            repc.fdbd1=repc.compParams(21);
            repc.fdbd2=repc.compParams(22);
            repc.femax=repc.compParams(23);
            repc.femin=repc.compParams(24);
            repc.pmax=repc.compParams(25);
            repc.pmin=repc.compParams(26);
            repc.tlag=repc.compParams(27);
            repc.ddn=repc.compParams(28);
            repc.dup=repc.compParams(29);
            repc.freqflag=repc.compParams(30);
            repc.stWrite=fopen(repc.repcOutputLocation,'w');
            if repc.stWrite<0
                sprintf('error opening file')
            else
                fprintf(repc.stWrite,'%8s %8s %8s %8s %8s %8s %8s \r\n',['State 1' 'State 2' 'State 3' 'State 4' 'State 5' 'State 6' 'State 7']);
            end  
        end
        %% Writing exciter states    
        function repc=WriteData(repc)
                fprintf(repc.stWrite,'%8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f \r\n',[repc.repcStates(1) repc.repcStates(2) repc.repcStates(3) repc.repcStates(4) repc.repcStates(5) repc.repcStates(6) repc.repcStates(7)]);  
        end
        %%
        function repc=Init(repc,pref,qext,ibranch,vreg,qbranch,pbranch,freq)
            if repc.refflag==1
                if repc.vcompflag==1
                    x1=abs(vreg-(repc.rc+1i*repc.xc)*ibranch);
                else
                    x1=vreg+qbranch*repc.kc;
                end
                x2=x1;
                x4=0;
                repc.vref=x2;
                repc.qref=0;
            elseif repc.refflag==0
                x2=0;
                x4=qbranch;
                repc.vref=0;
                repc.qref=x4;
            end
            x9=qext;
            x91=x9;
            qext1=x9*(1-repc.tft/repc.tfv);
            
            if repc.freqflag==1
                z1=pbranch;
                z5=pref;
                z41=z5;
                repc.plantpref=z1;
                repc.freqref=freq;
            else
                z1=0;
                z5=0;
                z41=0;
                repc.plantpref=0;
                repc.freqref=0;
            end
            
            repc.repcStates=[x2;x4;x91;qext1;z1;z41;z5]; % 7 X 1 vector of initial conditions
        end
        %%
        function [repc,pref,qext]=newStates(repc,tspan,prefrepc,ibranch,vreg,qbranch,pbranch,freq,dPlantPref)
            options=odeset('RelTol',1e-9);
            [T,X]=ode23t(@(t,x) repc_dyn(t,x,ibranch,vreg,qbranch,pbranch,freq,dPlantPref,repc.vref,repc.plantpref,repc.freqref,repc.qref,repc.refflag,repc.vcompflag,repc.tfltr,repc.rc,repc.xc,repc.kc,repc.dbd,repc.emax,repc.emin,repc.kp,repc.ki,repc.qmax,repc.qmin,repc.vfrz,repc.tft,repc.tfv,repc.fdbd1,repc.fdbd2,repc.ddn,repc.dup,repc.tp,repc.femax,repc.femin,repc.kpg,repc.kig,repc.pmax,repc.pmin,repc.tlag,repc.freqflag),tspan,repc.repcStates,options);
            
            function dy=repc_dyn(t1,y,ibranch,vreg,qbranch,pbranch,freq,dPlantPref,vref,plantpref,freqref,qref,refflag,vcompflag,tfltr,rc,xc,kc,dbd,emax,emin,kp,ki,qmax,qmin,vfrz,tft,tfv,fdbd1,fdbd2,ddn,dup,tp,femax,femin,kpg,kig,pmax,pmin,tlag,freqflag)
                
                dy=zeros(7,1);    
               
                if refflag==1
                    if vcompflag==1
                        x1=abs(vreg-(rc+1i*xc)*ibranch);
                    else
                        x1=vreg+qbranch*kc;
                    end
                    dy(1)=(x1-y(1))/tfltr;
                    x6=-y(1)+vref;
                    dy(2)=0;
                elseif refflag==0
                    dy(1)=0;
                    dy(2)=(qbranch-y(2))/tfltr;
                    x6=-y(2)+qref;
                end
                
                if x6>=-dbd && x6<=dbd
                    x7=0;
                else
                    x7=x6;
                end
                
                x8=x7;
                if x7>=emax
                    x8=emax;
                elseif x7<=emin
                    x8=emin;
                end
                
                x9=kp*x8+y(3);
                dy(3)=ki*x8;
                if x9>=qmax
                    dy(3)=0;
                    x9=qmax;
                elseif x9<=qmin
                    dy(3)=0;
                    x9=qmin;
                end
                
                if vreg<vfrz
                    dy(3)=0;
                end
                
                dy(4)=(x9*(1-tft/tfv)-y(4))/tfv; 
                if freqflag==1
                    dy(5)=(pbranch-y(5))/tp;
                    z6=-freq+freqref;
                    if (z6<=0 && z6>=fdbd1) || (z6> 0 && z6<=fdbd2)
                        z7=0;
                    else
                        z7=z6;
                    end
                    
                    z8=z7*ddn;
                    z9=z7*dup;
                    if z8>=0
                        z8=0;
                    end
                    
                    if z9<=0
                        z9=0;
                    end
                    z10=z8+z9;
                    
                    z2=plantpref-y(5)+z10+dPlantPref;
                    if z2>=femax
                        z3=femax;
                    elseif z2<=femin
                        z3=femin;
                    else
                        z3=z2;
                    end
                    
                    dy(6)=kig*z3;
                    z4=z3*kpg+y(6);
                    if z4>=pmax
                        dy(6)=0;
                        z4=pmax;
                    elseif z4<=pmin
                        dy(6)=0;
                        z4=pmin;
                    end
                    
                    dy(7)=(z4-y(7))/tlag;
                else
                    dy(6)=0;
                    dy(7)=0;
                end
            end
            
            repc.repcStates=X(end,:);   
            %%
                if repc.refflag==1
                    x6_1=-repc.repcStates(1)+repc.vref;
                elseif repc.refflag==0
                    x6_1=-repc.repcStates(2)+repc.qref;
                end
                
                if x6_1>=-repc.dbd && x6_1<=repc.dbd
                    x7_1=0;
                else
                    x7_1=x6_1;
                end
                
                x8_1=x7_1;
                if x7_1>=repc.emax
                    x8_1=repc.emax;
                elseif x7_1<=repc.emin
                    x8_1=repc.emin;
                end
                
                x9_1=repc.kp*x8_1+repc.repcStates(3);
                qext=x9_1*repc.tft/repc.tfv+repc.repcStates(4);               
                if repc.freqflag==0
                    pref=prefrepc; % governor control is disabled
                elseif repc.freqflag==1
                    pref=repc.repcStates(7);
                end
        end
        
    end
    
end

