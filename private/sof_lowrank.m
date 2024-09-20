function [Ky,gopt] = sof_lowrank(sys,ios,opt)

% Static output feedback design based on 
%
% L. El Ghaoui, F. Oustry, and M. AitRami, ‘A cone complementarity 
% linearization algorithm for static output-feedback and related problems’, 
% IEEE Trans. on Automatic Control, vol. 42, no. 8, p. 1171, 1997, 
% doi: 10.1109/9.618250.
%
% using the toolbox lmirank to impose the rank constraints
%
% Conditions:
%   - B2, C2, D12, and D21 must be parameter independet
%   - D22 must be zero
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


% fbianchi - 2024-09-18 - rev

% options
if (nargin < 3)
    opt = sofsettings();
end
etol = opt.etol;                % eigenvalue bound
rtol = opt.rtol;                % rank tolerance
gamma = opt.gamma;              % performance level
tol = opt.lowrank.tol;          % solution tolerance (bisection)
maxiter = opt.maxiter;          % max iteration (bisection)
minDecay = opt.minDecay;        % minimum decay rate
% solver options
solverOpt = sdpsettings('solver',opt.solver,...
                        'verbose', 0,...
                        'removeequalities', opt.removeequalities); 

% default solution
Ky = [];
gopt = 0;

% convert to 3D sys
[sys,sysInfo] = standardizeSys(sys);
nv = sysInfo.nv;
if (nv == 1)
    opt.type = 'rob';
elseif (nargin < 3)
    opt.type = 'gs';
end
% dimensions
[ns,nz,nw,ny,nu] = parsysdata(sys,ios,'dims');

% check if conditions are satisfied
[opc,msg] = check_conditions(sys,ios);
if strcmp(opc,'no')
    error(msg);
else
    disp(msg);
end


% =========================================================================
% Lmis set-up

% lyapunov functions
X = sdpvar(ns);
Y = sdpvar(ns);

% null spaces
[~,~,B2,~,C2,~,D12,D21] = parsysdata(sys(:,:,1),ios);
Ny = blkdiag(null([B2; D12]'),eye(nw));
Nx = blkdiag(null([C2, D21]),eye(nz));

% Bisection algorithm
g_lower = 0;
iter = 1;
if isinf(gamma)
    g = g_lower;
    g_upper = 1e4;
else
    g = gamma;
    g_upper = max(2.1*gamma,1e4);
end
while ((g_upper - g_lower) > tol) && (iter < maxiter)
    
    g_new = (g_lower + g_upper)/2;
    
    % stop condition
    if (iter > maxiter)
        fprintf('Max iteration number reached\n')
        break
    
    elseif (abs(g_new - g) <= tol)
        fprintf('Hinf Performance: %2.4f\n\n',g)
        break,
    else
        g = g_new;
    end
    
    lmis = [];
    for ii=1:nv
        
        [A,B1,~,C1,~,D11] = parsysdata(sys(:,:,ii),ios);
        A = A + minDecay*eye(ns);
        
        mat1 = blkvar;
        mat1(1,1) = A*Y + (A*Y)';
        mat1(1,2) = Y*C1';
        mat1(1,3) = B1;
        mat1(2,2) = -g*eye(nz);
        mat1(2,3) = D11;
        mat1(3,3) = -g*eye(nw);
        mat1 = sdpvar(mat1);
        mat1 = Ny'*mat1*Ny;
        if ~isnumeric(mat1)
            lmis = [lmis, mat1 <= -etol*eye(size(Ny,2))];
        else
            disp('LMI 1 is not necessary')
        end
        
        mat2 = blkvar;
        mat2(1,1) = X*A + (X*A)';
        mat2(1,2) = X*B1;
        mat2(1,3) = C1';
        mat2(2,2) = -g*eye(nw);
        mat2(2,3) = D11';
        mat2(3,3) = -g*eye(nz);
        mat2 = sdpvar(mat2);
        mat2 = Nx'*mat2*Nx;
        if ~isnumeric(mat2)
            lmis = [lmis, mat2 <= -etol*eye(size(Nx,2))];
        else
            disp('LMI 2 is not necessary')
        end
        
    end
    % lmis = [lmis, [X eye(ns);eye(ns) Y] >= -etol*eye(ns)]; % not needed,
    % it is added by the solver
    lmis = [lmis, rank([X eye(ns);eye(ns) Y]) <= ns];
     
    % solving
    warning('off', 'all')
    diagnostics = optimize(lmis, [], sdpsettings('lmirank.solver','sedumi','verbose',0));
    warning('on', 'all')
    P = value([X eye(ns);eye(ns) Y]);
    
    % checking solution
    if (diagnostics.problem ~= 1)
        rnc = rank(P, rtol);
        if (rnc == ns)
            g_upper = g;
            gopt = g;
            Xopt = value(X);
        else
            g_lower = g;
        end
    else
        if any(any(isnan(P)))
            rnc = inf;
        else
            rnc = rank(P, rtol);
        end
        g_lower = g;
    end
    fprintf('Iter: %3.0f, Rank: %2.0f (order=%2.0f), Norm-inf: %2.4f\n',...
        iter, rnc, ns, g);
    
    iter = iter + 1;
    
end

if (gopt == 0)
    Ky = [];
    fprintf('No solution was found for gamma < %4.2f\n',g_upper);
    return
end
    

% Controller computation
lmis = [];
md = sdpvar(1);
switch opt.type
    case 'gs'
        Ky = sdpvar(nu,ny,nv,'full');
        
        for ii = 1:nv

            [A,B1,B2,C1,C2,D11,D12,D21] = parsysdata(sys(:,:,ii),ios);
            
            phi = [2*md*Xopt + (Xopt*A) + (Xopt*A)', Xopt*B1, C1'; 
                   B1'*Xopt, -gopt*eye(nw), D11';
                   C1, D11, -gopt*eye(nz)];
            Px = [Xopt*B2; zeros(nw,nu); D12]';
            Q = [C2, D21, zeros(ny,nz)];

            lmis = [lmis, phi + Q'*Ky(:,:,ii)'*Px + (Q'*Ky(:,:,ii)'*Px)' <= 0];
        end
        
    case 'rob'
        Ky = sdpvar(nu,ny,'full');
        for ii = 1:nv
            
            [A,B1,B2,C1,C2,D11,D12,D21] = parsysdata(sys(:,:,ii),ios);
            
            phi = [2*md*Xopt + (Xopt*A) + (Xopt*A)', Xopt*B1, C1'; 
                   B1'*Xopt, -gopt*eye(nw), D11';
                   C1, D11, -gopt*eye(nz)];
            Px = [Xopt*B2; zeros(nw,nu); D12]';
            Q = [C2, D21, zeros(ny,nz)];

            lmis = [lmis, phi + Q'*Ky'*Px + (Q'*Ky'*Px)' <= 0];
        end
end
lmis = [lmis, md >= 0];
diagnostics = optimize(lmis, -md, solverOpt);

if (diagnostics.problem == 1)
    % Infeasible
    disp('The controller computation is infeasible')
    if strcmp(opt.type,'rob') && (nv > 1)
        disp('Try with option.type = ''gs''')
    end
    gopt = inf;
else
    Ky = value(Ky);
end

return


% ########################################################################
% Local functions

function [opc,msg] = check_conditions(sys,ios)

nv = size(sys,3);

if (nv == 1)
    opc = 'yes';
    msg = 'Conditions satisfied';
    
else
    
    [~,~,B2,~,C2,~,D12,D21,D22] = parsysdata(sys(:,:,1),ios);
    % check that b2,c2,d12,d21,d22 cte
    for ii = 2:nv
        [~,~,b2x,~,c2x,~,d12x,d21x,d22x] = parsysdata(sys(:,:,ii),ios);
        if (norm(B2 - b2x) > 0) || (norm(D12 - d12x) > 0)
            opc = 'no';
            msg = 'B2 and D12 must be parameter independent';

        elseif (norm(C2 - c2x) > 0) || (norm(D21 - d21x) > 0)
            opc = 'no';
            msg = 'C2 and D21 must be parameter independent';

        elseif (norm(D22-d22x) > 0)
            opc = 'no';
            msg = 'D22 must be zero';
            
        else
            opc = 'yes';
            msg = 'Conditions satisfied';

        end
    end
end

