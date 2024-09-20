function [Ky,gopt] = sof_mattei(sys, ios, opt)

% Static output feedback design based on 
%
% M. Mattei, ‘Robust multivariable PID control for linear parameter 
% varying systems’, Automatica, vol. 37, no. 12, pp. 1997–2003, 2001
%
% Conditions:
%   - B2, C2, D12, and D21 must be parameter independet
%   - B2 must be full row rank or
%   - C2 must be full column rank
%
% Similar to:
%  C.Crusius and A. Trofino, ‘Sufficient LMI conditions for output feedback 
%  control problems’, IEEE Transactions on Automatic Control, vol. 44, 
%  no. 5, pp. 1053–1057, May 1999, doi: 10.1109/9.763227.
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

% fbianchi - 2024-08-02 - rev


% options
if (nargin < 3)
    opt = sofsettings();
end
etol = opt.etol;          % eigenvalue bound
penalty = opt.penalty;    % norm penalization
minDecay = opt.minDecay;  % minimun -Re(poles)   
beta = opt.beta;          % to improve condition of Ky
if isempty(opt.case)
    mat = 'any';
else
    mat = opt.case;
end
% solver options
solverOpt = sdpsettings('solver',opt.solver,...
                        'verbose',opt.verb,...
                        'removeequalities', opt.removeequalities);                    
                    
% in case no solution
Ky = [];

% convert to 3D sys
[sys,sysInfo] = standardizeSys(sys);
nv = sysInfo.nv;
% dimensions
[ns,nz,nw,ny,nu] = parsysdata(sys(:,:,1),ios,'dims');
  
% check if the conditions are satisfied
[mat,msg] = check_conditions(sys,ios,mat);
if strcmp(mat,'no')
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
  
% optimization problem
switch mat
    case 'b2'   % case B2 full rank
        
        % transformation
        B2 = parsysdata(sys(:,:,1),ios,'b2');
        Ti = [pinv(B2); null(B2')'];

        % optimization variables
        X1 = sdpvar(nu);
        X2 = sdpvar(ns-nu);
        X  = blkdiag(X1,X2);
        Wp = sdpvar(nu,ny,'full');
        W  = [Wp; zeros(ns-nu,ny)];

        lmis = [];
        for ii = 1:nv
            
            [A,B1,~,C1,C2,D11,~,D21] = parsysdata(sys(:,:,ii),ios);
            
            % transformed system
            As  = Ti*(A + minDecay*eye(ns))/Ti; 
            C2s = C2/Ti; 
            B1s = Ti*B1; 
            C1s = C1/Ti;
            
            brl = blkvar;
            brl(1,1) = X*As+W*C2s+(X*As+W*C2s)';
            brl(1,2) = X*B1s+W*D21;
            brl(1,3) = C1s';
            brl(2,2) = -g*eye(nw);
            brl(2,3) = D11';
            brl(3,3) = -g*eye(nz);
            brl = sdpvar(brl);
            
            lmis = [lmis, brl <= -etol*eye(size(brl))];
            
        end
        lmis = [lmis, X >= etol*eye(ns), X1 >= beta*eye(ny)];
        
        % Objetive
        if isinf(opt.gamma)
            obj = g + penalty*(trace(X) + norm(Wp));
        else
            obj = penalty*(trace(X) + norm(Wp));
        end
        
        % solver
        diagnostics = optimize(lmis,obj,solverOpt);
        if (diagnostics.problem == 1)
            % Infeasible
            disp('Infeasible')
            gopt = inf;
        else
            % Feasible
            X = value(X);
            Ky = value(X1)\value(Wp);
            gopt = value(g);
            if (opt.verb)
                fprintf('---------------------------------------\n')
                fprintf('Eigenvalues of X1: min %4.3f, max %4.3f\n',...
                    min(eig(X)),max(eig(X)));
            end
            flag = 0;
            for ii = 1:nv
                [A,~,B2,~,C2] = parsysdata(sys(:,:,ii),ios);
                eigCL = eig(A+B2*Ky*C2);
                eigCLmx = max(real(eigCL));
                eigCLmn = min(real(eigCL));
                if (opt.verb)
                    fprintf('Closed-loop Eigenvalues: min %4.3f, max %4.3f\n',...
                        eigCLmn,eigCLmx);
                end
                if (eigCLmx > 0)
                    flag = 1;
                end
            end
            if (opt.verb)
                fprintf('Performance Hinf: %2.4f\n',gopt)
                fprintf('---------------------------------------\nºn')
            end
            if (flag == 1)
                warning('Numerical problem -> unreliable solution')
            end
        end

    
    case 'c2'   % case C2 full rank
        
        % transformation
        C2 = parsysdata(sys(:,:,1),ios,'c2');
        Ti = [pinv(C2), null(C2)];

        % optimization variables
        Y1 = sdpvar(ny);
        Y2 = sdpvar(ns-ny);
        Y  = blkdiag(Y1,Y2);
        Vp = sdpvar(nu,ny,'full');
        V  = [Vp, zeros(nu,ns-ny)];
    
        lmis = [];
        for ii = 1:nv
            
            [A,B1,B2,C1,~,D11,D12] = parsysdata(sys(:,:,ii),ios);
            
            % transformed system
            As  = Ti\(A + minDecay*eye(ns))*Ti; 
            B2s = Ti\B2; 
            B1s = Ti\B1; 
            C1s = C1*Ti;
            
            brl = blkvar;
            brl(1,1) = (As*Y + B2s*V) + (As*Y + B2s*V)';
            brl(1,2) = B1s;
            brl(1,3) = (C1s*Y + D12*V)';
            brl(2,2) = -g*eye(nw);
            brl(2,3) = D11';
            brl(3,3) = -g*eye(nz);
            brl = sdpvar(brl);

            lmis = [lmis, brl <= -etol*eye(size(brl))];
            
        end
        lmis = [lmis, Y >= etol*eye(ns), Y1 >= beta*eye(ny)];
        
        % objective
        if isinf(opt.gamma)
            obj = g + penalty*(trace(Y) + norm(Vp));
        else
            obj = penalty*(trace(Y) + norm(Vp));
        end
        % solver
        diagnostics = optimize(lmis,obj,solverOpt);
        if (diagnostics.problem == 1)
            % Infeasible
            disp('Infeasible')
            gopt = inf;
        else
            % Feasible
            Y = value(Y);
            Ky = value(Vp)/value(Y1);
            gopt = value(g);
            flag = 0;
            if (opt.verb)
                fprintf('---------------------------------------\n')
                fprintf('Eigenvalues of Y: min %4.3f, max %4.3f\n',...
                    min(eig(Y)),max(eig(Y)));
            end
            for ii = 1:nv
                [A,~,B2,~,C2] = parsysdata(sys(:,:,ii),ios);
                eigCL = eig(A+B2*Ky*C2);
                eigCLmx = max(real(eigCL));
                eigCLmn = min(real(eigCL));
                if (opt.verb)
                    fprintf('Closed-loop Eigenvalues: min %4.3f, max %4.3f\n',...
                        eigCLmn,eigCLmx);
                end
                if (eigCLmx > 0)
                    flag = 1;
                end
            end
            if (opt.verb)
                fprintf('Performance Hinf: %2.4f\n',gopt)
                fprintf('---------------------------------------\n\n')
            end
            if (flag == 1)
                warning('Numerical problem -> unreliable solution')
            end
        end
end
        


end


% ########################################################################
% Local functions

function [opc,msg] = check_conditions(sys,ios,opc)

% to check the applicability conditions

[~,~,B2,~,C2,~,D12,D21,D22] = parsysdata(sys(:,:,1),ios);
if (norm(D22) > 0)
    opc = 'no';
    msg = 'D22 must be zero';
    return
end
    
[~,~,~,ny,nu] = parsysdata(sys(:,:,1),ios,'dims');
nv = size(sys,3);

if (ny == rank(C2)) && (norm(D21) == 0)
    opcC2 = 'yes';
    msgC2 = 'C2 is full row rank, with D21=0';
    for ii = 2:nv
        C2x = parsysdata(sys(:,:,ii),ios,'c2');
        if norm(C2 - C2x) > 0
            opcC2 = 'no';
            msgC2 = 'C2 is parameter dependent';
            break
        end
    end
else
    opcC2 = 'no';
    msgC2 = 'C2 does not satisfy conditions';
end

if (nu == rank(B2)) && (norm(D12) == 0)
    opcB2 = 'yes';
    msgB2 = 'B2 is full column rank, with D12=0';
    for ii = 2:nv
        B2x = parsysdata(sys(:,:,ii),ios,'b2');
        if norm(B2 - B2x) > 0
            opcB2 = 'no';
            msgB2 = 'B2 is parameter dependent';
            break
        end
    end
else
    opcB2 = 'no';
    msgB2 = 'B2 does not satisfy conditions';
end

if strcmp(opc,'c2')
    if strcmp(opcC2,'yes')
        opc = 'c2'; msg = msgC2;
    else
        opc = 'no'; msg = msgC2;
    end
    
elseif strcmp(opc,'b2')
    if strcmp(opcB2,'yes')
        opc = 'b2'; msg = msgB2;
    else
        opc = 'no'; msg = msgB2;
    end
    
else
    if strcmp(opcC2,'yes')
        opc = 'c2'; msg = msgC2;
    elseif strcmp(opcB2,'yes')
        opc = 'b2'; msg = msgB2;
    else
        opc = 'no';
        msg = 'Neither B2 nor C2 satisfies conditions';
    end
end

end