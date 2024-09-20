function [Ky,gopt] = sof_prempain(sys, ios, opt)

% Static output feedback design based on 
%
%   E. Prempain, “Static output feedback stabilisation with Hinf performance 
%   for a class of plants,” Syst. Control Lett., vol. 43, no. 3, pp. 159–166, 
%   Jul. 2001.
%
% Conditions:
%   - C2*B2 must be full row rank
%   - B2 and C2 must be parameter independet
%
% Use:
%    [Ksof,gopt] = sofsyn(sys, ios, [opt])
% 
% argument within [] are optional
%
% Inputs:
%   sys: augmented plant (LTI or LPV)
%   ios: - 1x2 matrix with number of measures ny and number of controls nu
%        ([ny nu])
%        - struct defined as
%           ios.perf = {'y(1)','y(2)','y(3)','y(4)'};
%           ios.dist = {'u(1)','u(2)','u(3)','u(4)'};
%           ios.meas = {'y(5)'};
%           ios.ctrl = {'u(5)','u(6)'};
%          or numerical arrays with the input/output indices 
%   opt: struct produced by sofsettings
%
% Outputs:
%   Ky: static output gain
%   gopt: optimal performance level

% fbianchi - 2024-08-02


% options
if (nargin < 3)
    opt = sofsettings();
end
etol = opt.etol;          % eigenvalue bound
penalty = opt.penalty;    % norm penalization
beta = opt.beta;          % lower eigenvalue bound to improve condition of Q1
minDecay = opt.minDecay;  % minimun -Re(poles)   
% solver options
solverOpt = sdpsettings('solver',opt.solver,...
                        'verbose',opt.verb,...
                        'removeequalities', opt.removeequalities); 

% convert to 3D sys
[sys,sysInfo] = standardizeSys(sys);
nv = sysInfo.nv;
% dimensions
[ns,nz,nw,ny,nu] = parsysdata(sys(:,:,1),ios,'dims');

% check if the conditions are satisfied
[opc,msg] = check_conditions(sys,ios);
if strcmp(opc,'no')
    error(msg);
else
    disp(msg);
end


% in case performance level is set
if isinf(opt.gamma)
    g = sdpvar(1);
else
    g = opt.gamma;
end

% transformation
C2 = parsysdata(sys(:,:,1),ios,'c2');
Ti = [pinv(C2), null(C2)];

% transformed matrices
B2 = parsysdata(sys(:,:,1),ios,'b2');
B2c = Ti\B2;
C2c = C2*Ti;

B2c1 = B2c(1:ny,:); B2c2 = B2c(ny+1:end,:);
N = (B2c2*B2c1'/(B2c1*B2c1'))';
Tn = [eye(ny), zeros(ny,ns-ny);
      N',      eye(ns-ny)];

% optimization variables
Q1 = sdpvar(ny);
Q2 = sdpvar(ns-ny);
Qd = blkdiag(Q1,Q2);
Y = sdpvar(nu,ny,'full');

lmis = [];
for ii = 1:nv
    
    [A,B1,~,C1,~,D11,D12] = parsysdata(sys(:,:,ii),ios);
    A = A + minDecay*eye(ns);
    
    % transformed system
    Ac  = Ti\A*Ti;
    B1c = Ti\B1;
    C1c = C1*Ti;

    lmis = [];
    brl = blkvar;
    M = (Tn\Ac*Tn*Qd + Tn\B2c*Y*C2c);
    brl(1,1) = M + M';
    brl(2,1) = C1c*Tn*Qd + D12*Y*C2c;
    brl(3,1) = (Tn\B1c)';
    brl(2,2) = -g*eye(nw);
    brl(3,2) = D11';
    brl(3,3) = -g*eye(nz);
    brl = sdpvar(brl);
  
    lmis = [lmis, brl <= -etol*eye(size(brl))];
    
end
lmis = [lmis, Q1 >= beta*eye(ny), Qd >= etol*eye(ns)];

% Objetive
if isinf(opt.gamma)
    obj = g + penalty*(trace(Qd) + norm(Y));
else
    obj = penalty*(trace(Qd) + norm(Y));
end
% solver
diagnostics = optimize(lmis,obj,solverOpt);
if (diagnostics.problem == 1)
    % Infeasible
    disp('Infeasible')
    Ky = [];
    gopt = inf;
else
    % Feasible
    Q1 = value(Q1);
    gopt = value(g);
    Ky = value(Y)/Q1;
    if (opt.verb)
        fprintf('---------------------------------------\n')
        fprintf('Eigenvalues of Q1: min %4.3f, max %4.3f\n',...
            min(eig(Q1)),max(eig(Q1)));
        for ii = 1:nv
            [A,B2,C2] = ssdata(sys(ios.meas,ios.ctrl,ii));
            eigCL = eig(A + B2*Ky*C2);
            fprintf('Closed-loop Eigenvalues: min %4.3f, max %4.3f\n',...
                min(real(eigCL)),max(real(eigCL)));
        end
        fprintf('Performance Hinf: %2.4f\n',gopt)
        fprintf('---------------------------------------\n\n')
    end
end


end


% #######################################################################
% Local functions
function [opc,msg] = check_conditions(sys,ios)

% to check the applicability conditions
[~,~,B2,~,C2,~,~,D21,D22] = parsysdata(sys(:,:,1),ios);
[~,~,~,ny,nu] = parsysdata(sys(:,:,1),ios,'dims');
nv = size(sys,3);

% checking if D21 and D22 are null
if (norm(D21) > 0) || (norm(D22) > 0)
    opc = 'no';
    msg = 'D21 and D22 must be zero matrices';
    return
end

% checking if B2 and C2 are paremeter dependent
for ii = 2:nv
    opc = 'yes';
    msg = 'B2 and C2 are parameter independent';

    [~,~,B2x,~,C2x,~,D12,D21] = parsysdata(sys(:,:,1),ios);
    if norm(C2 - C2x) > 0
        opc = 'no';
        msg = 'C2 is parameter dependent';
        return
    end
    if norm(B2 - B2x) > 0
        opc = 'no';
        msg = 'B2 is parameter dependent';
        return
    end
end

% checking the rank of B2*C2
if (rank(C2*B2) == ny) 
    opc = 'yes';
    msg = 'The system satisfies the condition';
else
    opc = 'no';
    msg = 'C2B2 is not full row rank';
end

end