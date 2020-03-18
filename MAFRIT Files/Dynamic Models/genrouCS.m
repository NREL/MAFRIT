classdef genrouCS < Components
    %genrouCS Initializes all generators of type genrouCS and updates generator states during simulation
    % instead of passing voltages to the power flow at each integration step, it passes current injections.
    % based on the generator model defined in the book "Power System
    % Dynamics and Stability" by Peter W. Sauer and M.A. Pai
    properties
        td01;
        td011;
        tq01;
        tq011;
        h;
        d;
        xd;
        xq;
        xd1;
        xq1;
        xd11;
        xq11;
        x1;
        s1_0;
        s1_2;
        ra;
        rcomp;
        xcomp;
        w01;
        a;
        b;
        sbgen;
        genStates;
        genOutputLocation;
        stWrite;
    end
    
    methods
        %% Constructor
        function gen=genrouCS()
            gen.w01=120*pi;
        end
        %% Assigning generator parameters to variables
        function gen=AssignParams(gen)
            gen.sbgen=gen.compParams(1);
            gen.td01=gen.compParams(2);
            gen.td011=gen.compParams(3);
            gen.tq01=gen.compParams(4);
            gen.tq011=gen.compParams(5);
            gen.h=gen.compParams(6);
            gen.d=gen.compParams(7);
            gen.xd=gen.compParams(8);
            gen.xq=gen.compParams(9);
            gen.xd1=gen.compParams(10);
            gen.xq1=gen.compParams(11);
            gen.xd11=gen.compParams(12);
            gen.xq11=gen.compParams(13);
            gen.x1=gen.compParams(14);
            gen.s1_0=gen.compParams(15);
            gen.s1_2=gen.compParams(16);
            gen.ra=gen.compParams(17);
            gen.rcomp=gen.compParams(18);
            gen.xcomp=gen.compParams(19);
            gen.a=1.2-(1.0-1.2)/(-1+sqrt(1.0*gen.s1_0/(1.2*gen.s1_2))); % saturation parameter A calculated from the generaror parameters specified in the parameter text file
            gen.b=(1.0*gen.s1_0)/(1.0-gen.a)^2;                         % saturation parameter A calculated from the generaror parameters specified in the parameter text file
                                                                        % quadratic saturation function given in the PSS/E manual is used
            gen.stWrite=fopen(gen.genOutputLocation,'w');
            if gen.stWrite<0
                sprintf('error opening file')
            else
                fprintf(gen.stWrite,'%8s %8s %8s %8s %8s %8s \r\n',['State 1' 'State 2' 'State 3' 'State 4' 'State 5' 'State 6']);
            end           
        end
       %% Writing generator states    
        function gen=WriteData(gen)
            fprintf(gen.stWrite,'%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \r\n',[gen.genStates(1) gen.genStates(2) gen.genStates(3) gen.genStates(4) gen.genStates(5) gen.genStates(6)]);     
        end
       %% Initilizing generator states. Since the current version of MATLAB does not have the symbolic toolbox, fzero had to be used within a parameterized function to obtain the initial rotor angle. 
       % These equations are based on the power system dynamics book by Sauer and Pai
       function [gen,pmech,ifd]=Init(gen,eD,eQ,iD,iQ)
            y=SubstSymbolic(gen.xd,gen.xq,gen.a,gen.b,gen.xd11,gen.xq11,gen.ra,gen.x1);
            function y=SubstSymbolic(xd,xq,a,b,xd11,xq11,ra,x1)
                y=fzero(@delta,0);
                function y=delta(x)
                    psiq11=(-eD-iD*ra+iQ*xq11)*sin(x)+(eQ+iD*xq11+iQ*ra)*cos(x);
                    psid11=(eQ+iQ*ra+iD*xd11)*sin(x)+(eD-iQ*xq11+iD*ra)*cos(x);
                    psi11=sqrt(psiq11^2+psid11^2);
                    se=(b*(psi11-a)^2)/psi11;
                    s1q=psiq11*se*(xq-x1)/(xd-x1);
                    y=s1q-(eD*sin(x)-eQ*cos(x)+ra*(iD*sin(x)-iQ*cos(x))-xq*(iD*cos(x)+iQ*sin(x)));
                end
            end
            if y<0%atan2(eQ,eD)
                y=pi+y; % this is because a generator's rotor angle cannot be less than the angle of its terminal voltage
            end
            y-atan2(eQ,eD);
            % calculate state and algebraic variables after obtaining delta
            ed=eD*sin(y)-eQ*cos(y);
            eq=eD*cos(y)+eQ*sin(y);
            id=iD*sin(y)-iQ*cos(y);
            iq=iD*cos(y)+iQ*sin(y);
            psiq11=(-eD-iD*gen.ra+iQ*gen.xq11)*sin(y)+(eQ+iD*gen.xq11+iQ*gen.ra)*cos(y);
            psid11=(eQ+iQ*gen.ra+iD*gen.xd11)*sin(y)+(eD-iQ*gen.xq11+iD*gen.ra)*cos(y);
            psi11=sqrt(psiq11^2+psid11^2);
            se=(gen.b*(psi11-gen.a)^2)/psi11;
            sfd=psid11*se;
            s1q=psiq11*se*(gen.xq-gen.x1)/(gen.xd-gen.x1);
            ed1=iq*(gen.xq-gen.xq1)+s1q;
            eq1=eq+id*gen.xd1+gen.ra*iq;
            psikd=eq1-id*(gen.xd1-gen.x1);
            psi2q=-ed1-iq*(gen.xq1-gen.x1);
            psi1q=-ed1;
            efd=eq1+id*(gen.xd-gen.xd1)+sfd; % 8
            psid=-gen.xd11*id+psid11; % d axis flux linkage
            psiq=-gen.xq11*iq+psiq11; % q axis flux linkage
            pelec=psid*iq-psiq*id;
            pmech=pelec;
            ifd=efd;
            gen.genStates=[eq1;psikd;psi1q;psi2q;120*pi;y]; % 6 X 1 vector of initial conditions
        end
        %% function is called at every simulation iteration to calculate new generator states.
        function [gen,id2,iq2,ifd] = newStates(gen,tspan,pmech1,id1,iq1,ed,eq,efd1)
            options=odeset('RelTol',1e-9);
            [T,X]=ode23t(@(t,x) gen_dyn(t,x,gen.td01,gen.td011,gen.tq01,gen.tq011,gen.h,gen.d,gen.xd,gen.xq,gen.xd1,gen.xq1,gen.xd11,gen.xq11,gen.x1,gen.a,gen.b,pmech1,id1,iq1,efd1,gen.w01),tspan,gen.genStates,options);
            function dx = gen_dyn(t,x,td01,td011,tq01,tq011,h,d,xd,xq,xd1,xq1,xd11,xq11,x1,a,b,pmech,id,iq,efd,w0)
                dx=zeros(6,1);
                psid11=((xd1-xd11)/(xd1-x1))*x(2)+((xd11-x1)/(xd1-x1))*x(1); % d axis subtransient flux: 'psi-d-double-dash'
                psiq11=((xq1-xq11)/(xq1-x1))*x(4)+((xq11-x1)/(xq1-x1))*x(3); % q axis subtransient flux: 'psi-q-double-dash'
                psi11=sqrt(psid11^2+psiq11^2); % total subtransient flux
                se=b*(psi11-a)^2/psi11; % saturation due to subtransient flux
                sfd=psid11*se; % saturation in the field circuit
                s1q=psiq11*se*(xq-x1)/(xd-x1); % saturation in the 1st damper winding circuit
                psid=-xd11*id+psid11; % d axis flux linkage
                psiq=-xq11*iq+psiq11; % q axis flux linkage
                te=psid*iq-psiq*id;   % electric torque generated by the machine
                dx(1)=((xd-xd1)*(xd1-xd11)/(td01*(xd1-x1)^2))*x(2)-(x(1)/td01)*(1+(xd-xd1)*(xd1-xd11)/(xd1-x1)^2)-id*((xd-xd1)/td01)*(1-(xd1-xd11)/(xd1-x1))+efd/td01-sfd/td01; %field circuit state equation
                dx(2)=-x(2)/td011+x(1)/td011-id*(xd1-x1)/td011; % d-axis damper winding circuit state equation
                dx(3)=((xq-xq1)*(xq1-xq11)/(tq01*(xq1-x1)^2))*x(4)-(x(3)/tq01)*(1+(xq-xq1)*(xq1-xq11)/(xq1-x1)^2)-iq*((xq-xq1)/tq01)*(1-(xq1-xq11)/(xq1-x1))-s1q/tq01; % 1st q-axis damper winding circuit state equation
                dx(4)=-x(4)/tq011+x(3)/tq011-iq*(xq1-x1)/tq011; % 2nd q-axis damper winding circuit state equation
                dx(5)=((pmech-d*(x(5)/w0-1))/(x(5)/w0)-te)*(w0/(2*h)); % rotor speed state equation
                dx(6)=x(5)-w0; % rotor angle state equation
            end
            
            gen.genStates=X(end,:); % Since MATLAB solves the differential equations a number of times, we pick up the last state only for next iteration and for calculating the new d-q voltages
            
            % Voltage update based on new values of 'gen.genStates'
            psid110=((gen.xd1-gen.xd11)/(gen.xd1-gen.x1))*gen.genStates(2)+((gen.xd11-gen.x1)/(gen.xd1-gen.x1))*gen.genStates(1); % last step - d axis subtransient flux: 'psi-d-double-dash'
            psiq110=((gen.xq1-gen.xq11)/(gen.xq1-gen.x1))*gen.genStates(4)+((gen.xq11-gen.x1)/(gen.xq1-gen.x1))*gen.genStates(3); % last step - q axis subtransient flux: 'psi-q-double-dash'
            psi110=sqrt(psid110^2+psiq110^2); % last step - total subtransient flux
            se0=gen.b*(psi110-gen.a)^2/psi110; % last step - saturation due to subtransient flux
            sfd0=psid110*se0; % last step - saturation in the field circuit
            ifd=gen.genStates(1)+(gen.xd-gen.xd1)*(id1+(gen.genStates(1)-gen.genStates(2)-id1*(gen.xd1-gen.x1))*((gen.xd1-gen.xd11)/(gen.xd1-gen.x1)^2))+sfd0; % last step - represents field current that is used in the exciter equations to account for voltage drop in AC excitation systems due to field circuit loading
            idq=(1/(gen.xd11*gen.xq11+gen.ra^2))*[-gen.xq11 gen.ra;-gen.ra -gen.xd11]*([eq; -ed]-[psid110;psiq110]); % we should make wr=1 when neglecting the time constant 1/ws
            id2=idq(1);
            iq2=idq(2);
        end
    end
    
end

