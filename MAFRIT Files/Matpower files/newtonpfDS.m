function [V, converged, i] = newtonpfDS(Ybus, Sbus, V0, ref, pv, pq, mpopt,busnumI)
%NEWTONPF  Solves the power flow using a full Newton's method.
%   [V, CONVERGED, I] = NEWTONPF(YBUS, SBUS, V0, REF, PV, PQ, MPOPT)
%   solves for bus voltages given the full system admittance matrix (for
%   all buses), the complex bus power injection vector (for all buses),
%   the initial vector of complex bus voltages, and column vectors with
%   the lists of bus indices for the swing bus, PV buses, and PQ buses,
%   respectively. The bus voltage vector contains the set point for
%   generator (including ref bus) buses, and the reference angle of the
%   swing bus, as well as an initial guess for remaining magnitudes and
%   angles. MPOPT is a MATPOWER options struct which can be used to 
%   set the termination tolerance, maximum number of iterations, and 
%   output options (see MPOPTION for details). Uses default options if
%   this parameter is not given. Returns the final complex voltages, a
%   flag which indicates whether it converged or not, and the number of
%   iterations performed.
%
%   See also RUNPF.

%   EDITS FOR DYNAMICS SIMULATION: During dynamics simulation, under the
%   assumption of a voltage source model of the generator, the voltage
%   magnitude and angle of the PV buses at the kth iteration will be
%   supplied by the genrou model. Hence, voltage angle update for the PV
%   buses need not be performed by the power flow. Therefore, the Jacobian
%   elements corresponding to the PV buses are removed in this file. 

%   MATPOWER
%   Copyright (c) 1996-2015 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   $Id: newtonpf.m 2644 2015-03-11 19:34:22Z ray $
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default arguments
if nargin < 7
    mpopt = mpoption;
end

%% options
tol     = mpopt.pf.tol;
max_it  = mpopt.pf.nr.max_it;
max_it=500;

%% initialize
converged = 0;
i = 0;
V = V0;
Va = angle(V);
Vm = abs(V);
%% set up indexing for updating V
npv = length(pv);
npq = length(pq);
j1 = 1;         j2 = npv;           %% j1:j2 - V angle of pv buses
j3 = j2 + 1;    j4 = j2 + npq;      %% j3:j4 - V angle of pq buses
j5 = j4 + 1;    j6 = j4 + npq;      %% j5:j6 - V mag of pq buses

%% evaluate F(x0)
if ~isempty(busnumI)            % This condition if true implies that current sources are present and net power injection at those buses should be based on current supplied by the wind generator model in the main program
    Sbus(busnumI)=Sbus(busnumI).*V(busnumI)./V0(busnumI);
    Sbus0=Sbus;
end
mis = V .* conj(Ybus * V) - Sbus;
    F = [   real(mis(pq));
            imag(mis(pq))   ];
%% check tolerance
normF = norm(F, inf);
if mpopt.verbose > 1
    fprintf('\n it    max P & Q mismatch (p.u.)');
    fprintf('\n----  ---------------------------');
    fprintf('\n%3d        %10.3e', i, normF);
end
if normF < tol
    converged = 1;
    if mpopt.verbose > 1
        fprintf('\nConverged!\n');
    end
end

%% do Newton iterations
while (~converged && i < max_it)
    %% update iteration counter
    i = i + 1;
    
    %% evaluate Jacobian
    [dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V);
    
    j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
    j12 = real(dSbus_dVm([pv; pq], pq));
    j21 = imag(dSbus_dVa(pq, [pv; pq]));
    j22 = imag(dSbus_dVm(pq, pq));
    
    j11(1:npv,:)=[];j11(:,1:npv)=[];
    j12(1:npv,:)=[];
    j21(:,1:npv)=[];
    
    J = [   j11 j12;
            j21 j22;    ];

    %% compute update step
    dx = -(J \ F);

    %% update voltage
        % Update to voltage source bus voltage angle is not needed as these
        % are updated by the dynamic models
%     if npv
%         Va(pv) = Va(pv) + dx(j1:j2);
%     end
    if npq
        Va(pq) = Va(pq) + dx(j3-j2:j4-j2);
        Vm(pq) = Vm(pq) + dx(j5-j2:j6-j2);
    end
    V = Vm .* exp(1j * Va);
    Vm = abs(V);            %% update Vm and Va again in case
    Va = angle(V);          %% we wrapped around with a negative Vm
    %% evalute F(x)
    if ~isempty(busnumI)  % This condition if true implies that current sources are present and net power injection at those buses should be based on current supplied by the wind generator model in the main program
        Sbus(busnumI)=Sbus0(busnumI).*V(busnumI)./V0(busnumI); % This expression ensures that current is kept constant during Newton Raphson iterations    end
    end
    mis = V .* conj(Ybus * V) - Sbus;
    F = [   real(mis(pq));
            imag(mis(pq))   ];

    %% check for convergence
    normF = norm(F, inf);
    if mpopt.verbose > 1
        fprintf('\n%3d        %10.3e', i, normF);
    end
    if normF < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\nNewton''s method power flow converged in %d iterations.\n', i);
        end
    end
end

if mpopt.verbose
    if ~converged
        fprintf('\nNewton''s method power flow did not converge in %d iterations.\n', i);
    end
end
