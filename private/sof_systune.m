function [Ky,gopt] = sof_systune(sys,ios,opt)

% Static output feedback design based on systune.
%
% see systune for more details% 
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
% solver options
if (opt.verb)
    verb = 'iter';
else
    verb = 'final';
end
if (opt.minDecay <= 0)
    solverOpt = systuneOptions('Display',verb);
else
    solverOpt = systuneOptions('Display',verb,'MinDecay', opt.minDecay);
end                    

% convert to 3D sys
[sys,sysInfo] = standardizeSys(sys);
nv = sysInfo.nv;
% dimensions
[~,~,~,ny,nu] = parsysdata(sys(:,:,1),ios,'dims');
  
% gain to be found
if strcmp(opt.type,'rob')
    % robust case
    Ky = ltiblock.gain('Ky',nu,ny);
else
    % gain-scheduled case
    for ii = 1:nv
        name = sprintf('Ky_%d',ii);
        ky(:,:,ii) = ltiblock.gain(name,nu,ny);
    end
    Ky = genss(ky);
end
Ky.y = ios.ctrl;
Ky.u = ios.meas;

% closed-loop systems with tunable gains
sysGcl = connect(sys,Ky,ios.dist,ios.perf);

% objective
SGoal = TuningGoal.Gain(ios.dist,ios.perf,1);

% design
[~,gopt,~,info] = systune(sysGcl,SGoal,[],solverOpt);

% results
if strcmp(opt.type,'rob')
    % robust case
    aux = ss(info.Blocks.Ky);
else
    % gain-scheduled case
    f = fieldnames(info.Blocks);
    for ii = 1:nv
        aux(:,:,ii) = ss(info.Blocks.(f{ii}));
    end
end
Ky = aux.D;

