function newOpts = sofsettings(varargin)

% SOFSETTINGS allows to set the synthesis settings
%
% Use:
%   newOpts = SOFSETTINGS(Name, Value, ...)
%
% where opts is an struct with the some of the following fields
%
%     opts.etol = 1e-6;              % min eigenvalue to decide positive definitiveness
%     opts.penalty = 1e-4;           % penalty on Lyapunov matrices
%     opts.minDecay = 1e-2;          % minimum decay rate  
%     opts.gamma = inf;              % performance level
%     opts.type = 'rob';             % controller type: rob/gs (robust/gain-sheduled)  
% 
%     opts.beta = 1e-1;              % lower eigenvalue bound to improve solutions
%     
%     % specific for mattei/crucius
%     opts.case = 'any';             % to force c2 or b2 design
%     
%     % specific for lowrank algorithm
%     opts.rtol = 1e-6;              % rank tolerance  
%     opts.maxiter = 25;
%     opts.lowrank.tol = 0.01;       % tolerance in bisection
%     
%     % display options
%     opts.verb = true;              % verbose mode
% 
%     % solver options
%     opts.solver = 'mosek';         % solver for sdp: sedumi, sdpt3, mosek
%     opts.solTol = 1e-6;            % solver tolerance
%     opts.maxiter = 100;            % max solver iterations
%     opts.dualize = false;          % check YALMIP help
%     opts.removeequalities = true;  % check YALMIP help
%
% calling the function without arguments returns the default settings

% fbianchi - 2024-09-18


if (nargin == 0)
   % no arguments => return default settings
   
    newOpts.etol = 1e-6;              % min eigenvalue to decide positive definitiveness
    newOpts.penalty = 1e-4;           % penalty on Lyapunov matrices
    newOpts.minDecay = 0;             % minimum decay rate  
    newOpts.gamma = inf;              % performance level
    newOpts.type = 'rob';             % controller type: rob/gs (robust/gain-sheduled)  

    newOpts.beta = 1e-1;              % lower eigenvalue bound to improve solutions
    
    % specific for mattei/crucius
    newOpts.case = 'any';             % to force c2 or b2 design
    
    % specific for lowrank algorithm
    newOpts.rtol = 1e-6;              % rank tolerance  
    newOpts.maxiter = 25;
    newOpts.lowrank.tol = 0.01;       % tolerance in bisection
    
    % display options
    newOpts.verb = true;                % verbose mode

    % solver options
    if exist('mosekdiag','file')      % solver for sdp: sedumi, sdpt3, mosek  
        newOpts.solver = 'mosek';       
    else
        newOpts.solver = 'sedumi';    
    end
    newOpts.solTol = 1e-6;            % solver tolerance
    newOpts.maxiter = 100;            % max solver iterations
    newOpts.dualize = false;          % check YALMIP help
    newOpts.removeequalities = true;  % check YALMIP help
    
else
    
    % new setting struct with default values
    newOpts = sofsettings;
    
    for ii = 1:2:length(varargin)
        
        if isfield(newOpts, varargin{ii})
            
            % check for valid names
            checkField(varargin{ii}, varargin{ii+1});
            
            % changing values
            newOpts.(varargin{ii}) = varargin{ii+1};
            
        else
            error('SOFSYNSETTINGS:inputError',...
                '%s is not an acceptable setting', varargin{ii})
            
        end
        
    end
        
    
end


% ----------------------------------------------------------------------

function checkField(field, val)

    % checking for valid name-value pairs

    bool = false;
    switch field
        case {'etol', 'rtol', 'penalty', 'maxiter', 'solTol',...
                'minDecay', 'gamma', 'beta','lowrank.tol'}
            
            if ~isnumeric(val)
                msg = 'a numeric value';
                bool = true;
            end
            
        case {'verb'}
            if ~(isnumeric(val) || islogical(val))
                msg = 'logical value';
                bool = true;
            end
            
        case {'type'}
            
            if ~ismember(val,{'rob', 'gs'})
                msg = 'rob or gs';
                bool = true;
            end
            
        case {'case'}
            
            if ~ismember(val,{'any', 'b2', 'c2'})
                msg = 'any, b2 or c2';
                bool = true;
            end
            
        case {'solver'}
            
            if ~ismember(val,{'sedumi', 'sdpt3', 'mosek'})
                msg = 'sedumi, sdpt3 or mosek';
                bool = true;
            end
            
        case {'removeequalities'}
            
            if ~ismember(val,{-1, 0, 1, 2})
                msg = '-1, 0, 1, or 2, check sdpsettings for help';
                bool = true;
            end
            
        case 'dualize'
            
            if ~(islogical(val) || (val == 1) || (val == 0))
                msg = 'must be true o false';
                bool = true;
            end
            
        otherwise
            error('Unexpected field')    
        
    end

    if bool
        error('SOFSYNSETTINGS:inputError','%s must be %s',...
            field, msg)
    end


