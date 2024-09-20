function varargout = parsysdata(sys,ios,flag)

% PARSYSDATA returns the system matrices partitioned according to
%
%	 dx/dt = A  * x + B1  * w + B2  * u
%    z     = C1 * x + D11 * w + D12 * u
%	 y     = C2 * x + D21 * w + D22 * u
%
% Use:
%   - [a,b1,b2,c1,c2,d11,d12,d21,d22] = PARSYSDATA(sys,ios,flag)
%   - data = PARSYSDATA(sys,ios,flag)
%   - [ns,nz,nw,ny,nu] = PARSYSDATA(sys,ios,'dims')
% 
% Input:
% 	- sys:  system (mat or lti format)
%   - ios:  vector [ny nu], ny is # rows of C2 and nu # of columns of B2
%           struct with fields: 
%                   ios.perf = cell array or numeric vector for signals in Z
%                   ios.dist = cell array or numeric vector for signals in W
%                   ios.meas = cell array or numeric vector for signals in Y
%                   ios.ctrl = cell array or numeric vector for signals in U
%
%      flag: (optional)
%            'dims': it returns the dimensions of the partition matrices
%            'a','b1',... returns only the correspondig matrix, e.g. 
%                  c1 = parsysdata(P,r,'c1')
%

% fbianchi - 04/09/2013 - from hinfpar
% fbianchi - 16/06/2014 - rev 1.1
% fbianchi - 2020-07-01 - rev 2.0


if (nargin < 2)
    error('PARSYSDATA:inputError',...
        'Use [a,b1,b2,c1,c2,d11,d12,d21,d22] = PARSYSDATA(sys,ios)');

elseif (nargin < 3)
    flag = 'none';
    
end

% I/O maps
%
if isnumeric(ios) && isvector(ios) && (length(ios) == 2)
    % ios numeric vector
    
    % system information
    if (isnumeric(sys) && sys(end,end) == -Inf)
        % sys an lmitool object
        [~,ni,no] = sinfo(sys);
        [a,b,c,d,~] = ltiss(sys);

    elseif isa(sys,'ss') || isa(sys,'tf') || isa(sys,'zpk')
        % lti control object
        [no,ni] = iosize(sys);
        [a,b,c,d] = ssdata(ss(sys));

    else
        error('PARSYSDATA:inputError',...
            'SYS must be an LTI model')

    end

    % ios is a numeric vector
    ny = ios(1);
    nu = ios(2);
    nz = no - ny;
    nw = ni - nu;

    % indicess corresponding to i/o signals
    y_ind = nz+1:no;
    u_ind = nw+1:ni;
    z_ind = 1:nz;
    w_ind = 1:nw;

    
elseif isstruct(ios) && isfield(ios,'dist')  && isfield(ios,'perf')
    % ios struct
    
    % system information
    if isa(sys,'ss') || isa(sys,'tf') || isa(sys,'zpk')
        % lti control object
        [no,ni] = iosize(sys);
        [a,b,c,d] = ssdata(ss(sys));

    else
        error('PARSYSDATA:inputError',...
            'When IOS is a struct, SYS must be a SS, TF or ZPK model')

    end
                    
    % finding signal indices
    % u -> control 
    if ~isfield(ios,'ctrl')
        if isfield(sys.InputGroup,'ctrl')
            ios.ctrl = sys.InputGroup.ctrl;
        else
            error('PARSYSDATA:inputError',...
                    'There is no information about the control input U')
        end
    end
    if isnumeric(ios.ctrl)
        if (max(ios.ctrl) > ni) && (min(ios.ctrl) <= 0)
            error('PARSYSDATA:inputError',...
                'Control indices exceed the number of pdG inputs')
        else
            u_ind = ios.ctrl;
        end
        
    elseif iscellstr(ios.ctrl)
        % list of names
        nu = length(ios.ctrl);
        u_ind(1,nu) = 0;
        i_names = sys.u;
        for ii = 1:nu
            aux = find(strcmp(i_names, ios.ctrl{ii}));
            if isempty(aux)
                error('PARSYSDATA:inputError',...
                    'At least one control name is not an input of SYS')
            else
                u_ind(ii) = aux;
            end
        end
        
    elseif ischar(ios.ctrl)
        % only 1 name
        i_names = sys.u;
        u_ind = find(strcmp(i_names, ios.ctrl));
        if isempty(u_ind)
            error('PARSYSDATA:inputError',...
                'At least one control name is not an input of SYS')
        end
        
    else
        error('PARSYSDATA:inputError',...
            'IOS.CTRL must be vector with signal indices or names')
    end    
    
    % y -> measure
    if ~isfield(ios,'meas')
        if isfield(sys.OutputGroup,'meas')
            ios.meas = sys.OutputGroup.meas;
        else
            error('There is no information about measured output Y')
        end
    end
    if isnumeric(ios.meas)
        if (max(ios.meas) > ni) && (min(ios.meas) <= 0)
            error('PARSYSDATA:inputError',...
                'The indices of Y exceed the number of SYS outputs')
        else
            y_ind = ios.meas;
        end
        
    elseif iscellstr(ios.meas)
        % list of names
        ny = length(ios.meas);
        y_ind(1,ny) = 0;
        i_names = sys.y;
        for ii = 1:ny
            aux = find(strcmp(i_names, ios.meas{ii}));
            if isempty(aux)
                error('PARSYSDATA:inputError',...
                    'At least one measured output name is not an output of SYS')
            else
                y_ind(ii) = aux;
            end
        end
        
    elseif ischar(ios.meas)
        % only 1 name
        i_names = sys.y;
        y_ind = find(strcmp(i_names, ios.meas));
        if isempty(y_ind)
            error('PARSYSDATA:inputError',...
                'At least one measured output name is not an output of SYS')
        end
        
    else
        error('PARSYSDATA:inputError',...
            'IOS.MEAS must be vector with signal indices or names')
    end    

    % w -> disturbance
    if isnumeric(ios.dist)
        if (max(ios.dist) > ni) && (min(ios.dist) <= 0)
            error('PARSYSDATA:inputError',...
                'The indices of W exceed the number of SYS inputs')
        else
            w_ind = ios.dist;
        end
        
    elseif iscellstr(ios.dist)
        % list of names
        nw = length(ios.dist);
        w_ind(1,nw) = 0;
        i_names = sys.u;
        for ii = 1:nw
            aux = find(strcmp(i_names, ios.dist{ii}),1,'first');
            if isempty(aux)
                error('PARSYSDATA:inputError',...
                    'At least one disturbance input name is not an input of SYS')
            else
                w_ind(ii) = aux;
            end
        end
        
    elseif ischar(ios.dist)
        % only 1 name
        i_names = sys.u;
        w_ind = find(strcmp(i_names, ios.dist),1,'first');
        if isempty(w_ind)
            error('PARSYSDATA:inputError',...
                'At least one disturbance input name is not an input of SYS')
        end
        
    else
        error('PARSYSDATA:inputError',...
            'IOS.DIST must be vector with signal indices or names')
    end    

    % z -> performance
    if isnumeric(ios.perf)
        if (max(ios.perf) > ni) && (min(ios.perf) <= 0)
            error('PARSYSDATA:inputError',...
                'The indices of Z exceed the number of SYS outputs')
        else
            z_ind = ios.perf;
        end
        
    elseif iscellstr(ios.perf)
        % list of names
        nz = length(ios.perf);
        z_ind(1,nz) = 0;
        i_names = sys.y;
        for ii = 1:nz
            aux = find(strcmp(i_names, ios.perf{ii}),1,'first');
            if isempty(aux)
                error('PARSYSDATA:inputError',...
                    'At least one performance output name is not an output of SYS')
            else
                z_ind(ii) = aux;
            end
        end
        
    elseif ischar(ios.dist)
        % only 1 name
        i_names = sys.y;
        z_ind = find(strcmp(i_names, ios.perf),1,'first');
        if isempty(z_ind)
            error('PARSYSDATA:inputError',...
                'At least one performance output name is not an output of SYS')
        end
        
    else
        error('PARSYSDATA:inputError',...
            'IOS.PERF must be vector with indices or signal names')
    end    
   
else
    varargout = [];
end


switch flag
    
    case 'dims'
        
        % returns on the dimensions
        
        varargout{1} = order(sys(:,:,1));      % ns
        varargout{2} = length(z_ind);   % nz
        varargout{3} = length(w_ind);   % nw
        varargout{4} = length(y_ind);   % ny
        varargout{5} = length(u_ind);   % nu

    case 'none'
        
        % returns all matrices
        
        b1  = b(:,w_ind); 
        b2  = b(:,u_ind);
        c1  = c(z_ind,:); 
        c2  = c(y_ind,:); 
        d11 = d(z_ind,w_ind); 
        d12 = d(z_ind,u_ind);
        d21 = d(y_ind,w_ind); 
        d22 = d(y_ind,u_ind);
        
        varargout = {a,b1,b2,c1,c2,d11,d12,d21,d22};
        
    case 'a'
        varargout = {a};
        
    case 'b1'
        varargout = {b(:,w_ind)};
        
    case 'b2'
        varargout = {b(:,u_ind)};
        
    case 'c1'
        varargout = {c(z_ind,:)};
        
    case 'c2'
        varargout = {c(y_ind,:)}; 
        
    case 'd11'
        varargout = {d(z_ind,w_ind)};
        
    case 'd12'
        varargout = {d(z_ind,u_ind)};
        
    case 'd21'
        varargout = {d(y_ind,w_ind)};
        
    case 'd22'
        varargout = {d(y_ind,u_ind)};
        
    otherwise
        error('Unknown option')
end
        
