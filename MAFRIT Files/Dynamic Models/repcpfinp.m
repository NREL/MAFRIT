function [vreg,ibranch,pbranch,qbranch,freq]=repcpfinp(fbus,tbus,ckt,V,Ybus,deltat)
%repcpfinp calculates the branch current and real and reactive power for repc_a module
%   When both fbus and tbus are non zero, it is assumed here that our
%   intention is to regulate the voltage of the tbus and fbus is the
%   monitored bus

persistent vangold;
persistent deltaw;

if fbus==0 || tbus==0
    % Regulate To Bus Voltage and pbranch, ibranch, qbranch, freq ar not
    % used. In this case refflag should be set to 1 and freqflag should be
    % set to 0 because any other flag setting won't make sense.
    if fbus==0
        vreg=abs(V(tbus));
    elseif tbus==0
        vreg=abs(V(fbus));
    end
    ibranch=0;
    pbranch=0;
    qbranch=0;
    freq=0;
else
    vreg=abs(V(fbus));
    ibranch=full((V(fbus)-V(tbus))*(-Ybus(fbus,tbus)));
    qbranch=full(imag(V(fbus)*conj(ibranch)));
    pbranch=full(real(V(fbus)*conj(ibranch)));
    if isempty(vangold)
        vangold=angle(V(fbus));
        deltaw=0;
    else
        deltaw=deltaw+((angle(V(fbus))-vangold)/(deltat*120*pi)-deltaw)*deltat/0.05; % forward euler is used to calculate the deviation in p.u. frequency 
                                                                                     % assuming a filter time constant of 0.05 seconds. 
                                                                                     
    end
    freq=1+deltaw;
end
end

