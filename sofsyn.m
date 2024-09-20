function [Ksof,gopt] = sofsyn(sys,ios, strmethod, opt)

% Static output feedback design for a Hinf performance objective using 
% several approaches.
%
% Use:
%    [Ksof,gopt] = sofsyn(sys,ios, [method], [opt])
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
%   method: - mattei (default) using the relaxation proposed in        
%               M. Mattei, ‘Robust multivariable PID control for linear 
%               parameter varying systems’, Automatica, vol. 37, no. 12, 
%               pp. 1997–2003, 2001
%             the following conditions must hold:
%               - B2, C2, D12, and D21 must be parameter independet
%               - B2 must be full row rank or C2 must be full column rank
%
%           - crucius: using the relaxation proposed in 
%               C.Crusius and A. Trofino, ‘Sufficient LMI conditions for 
%               output feedback control problems’, IEEE Trans. on Automatic 
%               Control, vol. 44, no. 5, pp. 1053–1057, 1999, doi: 10.1109/9.763227.
%             the same conditions applies.
%
%           - prempain using relaxation proposed in
%               E. Prempain, “Static output feedback stabilisation with Hinf 
%               performance  for a class of plants,” Syst. Control Lett., 
%               vol. 43, no. 3, pp. 159–166, Jul. 2001.
%             the following conditions must hold:
%               - C2*B2 must be full row rank
%               - B2 and C2 must be parameter independet
%
%           - lowrank using the approach based on rank constraints proposed in
%               L. El Ghaoui, F. Oustry, and M. AitRami, ‘A cone complementarity 
%               linearization algorithm for static output-feedback and related problems’, 
%               IEEE Trans. on Automatic Control, vol. 42, no. 8, p. 1171, 1997, 
%
%           - systune an implementaition solving the nonlinear programming
%               problem using systune
%
% Outputs:
%   Ky: static output gain
%   gopt: optimal performance level

% fbianchi - 2024-09-20 - rev 1



if (nargin < 3)
    strmethod = 'mattei';
elseif (nargin < 4)
    strmethod = 'mattei';
    opt = sofsetting();
end


% =========================================================================
% cases:
switch strmethod
    case 'mattei'
        [Ksof,gopt] = sof_mattei(sys, ios, opt);

    case 'crucius'
        [Ksof,gopt] = sof_crucius(sys, ios, opt);

    case 'prempain'
        [Ksof,gopt] = sof_prempain(sys, ios, opt);

    case 'lowrank'
        [Ksof,gopt] = sof_lowrank(sys, ios, opt);

    case 'systune'
        [Ksof,gopt] = sof_systune(sys, ios, opt);

    otherwise
        error('Method not available')
        
end

