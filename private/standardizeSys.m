function [sys,sysInfo] = standardizeSys(sys)

% STANDARDIZESYS serves to standarize the systems in design functions
% the output is a lti object of nv x 1 dimension

% fbianchi - 2021-07-31

if isa(sys,'pass') || isa(sys,'ppss')
    % lpv objects
    sysInfo.sys = sys;
    sysInfo.typ = 'lpv';
    [ny,nu,ns] = size(sys);
    sysInfo.ns = ns;
    sysInfo.ny = ny;
    sysInfo.nu = nu;
    sys = ss(sys);
    sysInfo.nv = size(sys,3);
    
elseif (isa(sys,'ss') || isa(sys,'tf') || isa(sys,'zpk'))
    % ss/tf/zpk systems
    sysInfo.typ = 'lti';
    sysInfo.ns = order(sys);
    [ny,nu] = iosize(sys);
    sysInfo.ny = ny;
    sysInfo.nu = nu;
    nd = size(sys);
    if (length(nd) > 2)
        % convert into nv x 1 array
        sysInfo.nv = prod(nd(3:end));
        sys = reshape(sys,[sysInfo.nv, 1]);
    else
        sysInfo.nv = 1;
        sys = ss(minreal(sys));
    end
    
elseif ispsys(sys)
    % psys systems
    sysInfo.typ = 'psys';
    [type,nv,ns,nu,ny] = psinfo(sys);
    if strcmp(type,'aff')
        sys = aff2pol(sys);
    end
    for ii = 1:nv
        sys_aux(:,:,ii) = mat2lti(psinfo(sys,'sys',ii));
    end
    sys = sys_aux;
    sysInfo.ns = ns;
    sysInfo.ny = ny;
    sysInfo.nu = nu;
    sysInfo.nv = nv;
    
elseif isnumeric(sys)
    % lmitool/robust (old) objects
    sysInfo.typ = 'mat';
    sys = mat2lti(sys);
    sysInfo.ns = order(sys);
    [nu,ny] = iosize(sys);
    sysInfo.ny = ny;
    sysInfo.nu = nu;
    sysInfo.nv = 1;
else
    error('STANDARDIZESYS:InputError',...
        'SYS is not a valid system description')
end
