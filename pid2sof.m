function syspid = pid2sof(sys,ios,str,Gi,tau)

% PID2SOF returns the augmented plant to compute a PID controller as a 
% static-output feedback problem
%
% Use:
%   syspid = pid2sof(sys,ios,str,tau)
%
% Input:
%   sys: plant model (LTI or LPv object)
%   ios: 1x2 matrix as [ny nu]
%        struct ios.perf: string arrar with element of z
%               ios.dist: string arrar with element of w
%               ios.meas: string arrar with element of y (error)
%               ios.ctrl: string arrar with element of u (control)
%   str: controller structure: pi, pi2dof, pid, d-pi, pid2dof
%
% The derivative term is approximated by s/(tau*s+1), by default 
% tau = 1/10*min(real(eig(sys)))

% fbianchi - 2024-07-29 - rev


% derivative term approximation
if (nargin < 5)
    pd = 10*min(real(eig(sys)));
else
    pd = 1/tau;
end

% old formats to new ones
if ispsys(sys)
    % psys systems
    [type,nv,ns,nu,ny] = psinfo(sys);
    if strcmp(type,'aff')
        sys = aff2pol(sys);
    end
    
elseif isnumeric(sys)
    % lmitool/robust (old) objects
    sys = mat2lti(sys);
end
% plant dimensions
[no,ni] = iosize(sys);
if isstruct(ios)
    if isfield(ios,'ctrl')
        nu = length(ios.ctrl);
    else
        error('IOS must have a field ctrl')
    end
    if isfield(ios,'meas')
        ny = length(ios.meas);
    else
        error('IOS must have a field meas')
    end
    if isfield(ios,'perf')
        nz = length(ios.perf);
    else
        error('IOS must have a field perf')
    end
    if isfield(ios,'dist')
        nw = length(ios.dist);
    else
        error('IOS must have a field dist')
    end
    
elseif isnumeric(ios) && (length(ios) == 2) 
    ny = ios(1);
    nz = no - ny;
    nu = ios(2);
    nw = ni - nu;
    ios = struct();
    ios.meas = sprintfc('y(%d)',1:ny);
    ios.ctrl = sprintfc('u(%d)',1:nu);
    ios.perf = sprintfc('z(%d)',1:nz);
    ios.dist = sprintfc('w(%d)',1:nw);
end
sys.u = [ios.dist, ios.ctrl];
sys.y = [ios.perf, ios.meas];

if (nargin < 4)
    Gi = eye(ny);
end
ns = size(Gi,1);

% building the Gi block
switch str
    case 'pi'
        
        Ac = zeros(ns);
        Bc = [Gi -Gi];
        Cc = [zeros(ny,ns); 
              eye(ns)];
        Dc = [eye(ny) -eye(ny);
              zeros(ns,2*ny)];
        
    case 'pi2dof'
        
        Ac = zeros(ny);
        Bc = [eye(ny), -eye(ny)];
        Cc = [zeros(ny);
              zeros(ny);
              eye(ny)];
        Dc = [eye(ny), zeros(ny);
              zeros(ny), -eye(ny);
              zeros(ny,2*ny)];

    case 'pid'
        
        Ac = blkdiag(zeros(ny),-pd*eye(ny));
        Bc = [eye(ny) -eye(ny);
              pd*eye(ny); -pd*eye(ny)];
        Cc = [zeros(ny,2*ny);
              eye(ny), zeros(ny);
              zeros(ny), pd*eye(ny)];
        Dc = [eye(ny), -eye(ny);
              zeros(ny,2*ny);
              pd*eye(ny), -pd*eye(ny)];
        
    case 'd-pi'
        
        Ac = blkdiag(zeros(ny),-pd*eye(ny));
        Bc = [eye(ny), -eye(ny);
              zeros(ny), -pd*eye(ny)];
        Cc = [zeros(ny,2*ny);
              eye(ny), zeros(ny);
              zeros(ny), pd*eye(ny)];
        Dc = [eye(ny), -eye(ny);
              zeros(ny,2*ny);
              zeros(ny), pd*eye(ny)];

    case 'pid2dof'
        
        Ac = blkdiag(zeros(nc),-pd*eye(nc));
        Bc = [eye(ny), -eye(nc); 
              pd*eye(ny), -pd*eye(ny)];
        Cc = [zeros(2*ny,2*ny);
              eye(ny), zeros(ny);
              zeros(ny), pd*eye(ny)];
        Dc = [eye(ny), zeros(ny);
              zeros(ny), -eye(ny);
              zeros(ny), zeros(ny);
              zeros(ny), pd*eye(ny)];
       
    otherwise
        error('Invalid controller structure')
end
Kstr = ss(Ac, Bc, Cc, Dc);
Kstr.u = [sprintfc('r(%d)',1:ny), ios.meas];
Kstr.y = sprintfc('v(%d)',1:size(Cc,1));

% augmented plant for synthesis
in = [ios.dist, sprintfc('r(%d)',1:ny), ios.ctrl];
out = [ios.perf, Kstr.y'];
syspid = connect(sys,Kstr,in,out);

