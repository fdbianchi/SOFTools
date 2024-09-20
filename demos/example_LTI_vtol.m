
% Example of optimal design of a static output feedback gain.
%
% Case LTI: longitudinal motion of a VTOL helicopter
%
% Example from:
%   E. Prempain, “Static output feedback stabilisation with Hinf performance 
%   for a class of plants,” Syst. Control Lett., vol. 43, no. 3, pp. 159–166, 
%   Jul. 2001. doi: 10.1016/S0167-6911(01)00087-1

% fbianchi - 2024-08-01 - rev

% cleaning
clearvars
close all
clc

% --------------------------------------------------------------------
% Model
A   = [-0.0366  0.0271  0.0188 -0.4555;
        0.0482 -1.0100  0.0024 -4.0208;
        0.1002  0.3681 -0.7070  1.4200;
        0       0       1       0];
B1  = eye(4);
B2  = [-0.4422  0.1761;
        3.5446 -7.5922;    
       -5.5200  4.4900;
        0       0];
C2  = [0 1 0 0];
C1  = blkdiag(eye(2), zeros(2));
D11 = zeros(4);
D12 = [zeros(2); eye(2)];
D21 = zeros(1,4);
D22 = zeros(1,2);
sys = ss(A,[B1 B2],[C1; C2], [D11 D12; D21 D22]);
sys.y = 'y';
sys.u = 'u';

% i/o map
ios.perf = {'y(1)','y(2)','y(3)','y(4)'};
ios.dist = {'u(1)','u(2)','u(3)','u(4)'};
ios.meas = {'y(5)'};
ios.ctrl = {'u(5)','u(6)'};


% -------------------------------------------------------------------------
% controller designs
%
% full order
[Kfull,~,gfull] = hinfsyn(sys,1,2,'method','lmi'); 
Gclfull = lft(sys,Kfull);
w = logspace(-2,2,100);
sv0 = sigma(Gclfull,w);
eigmx0 = max(real(eig(Gclfull)));

% using hinfstrcut
Ksof1 = ltiblock.gain('Ky',2,1);
sysCl = lft(sys,Ksof1);
[CL,gsof1,INFO] = hinfstruct(sysCl);
Ksof1 = ss(INFO.TunedBlocks.Ky); Ksof1 = Ksof1.d;
Gclsof1 = lft(sys,Ksof1);
sv1 = sigma(Gclsof1,w);
eigmx1 = max(real(eig(Gclsof1)));

% using paper's method
opt = sofsettings('minDecay',0.01,'beta',0.01);
[Ksof2,gsof2] = sofsyn(sys,ios,'prempain',opt);
Gclsof2 = lft(sys,Ksof2);
sv2 = sigma(Gclsof2,w);
eigmx2 = max(real(eig(Gclsof2)));

% using mattei method
[Ksof3,gsof3] = sofsyn(sys,ios,'mattei',opt);
Gclsof3 = lft(sys,Ksof3);
sv3 = sigma(Gclsof3,w);
eigmx3 = max(real(eig(Gclsof3)));

% using lmirank method
[Ksof4,gsof4] = sofsyn(sys,ios,'lowrank',opt);
Gclsof4 = lft(sys,Ksof4);
sv4 = sigma(Gclsof4,w);
eigmx4 = max(real(eig(Gclsof4)));


fprintf('\n')
fprintf('Performance comparison\n')
fprintf('Full order: %6.2f, max eig: %04.2f\n',gfull,eigmx0)
fprintf('Hinfstruct: %6.2f, max eig: %04.2f\n',gsof1,eigmx1)
fprintf('Prempain:   %6.2f, max eig: %04.2f\n',gsof2,eigmx2)
fprintf('Mattei:     %6.2f, max eig: %04.2f\n',gsof3,eigmx3)
fprintf('Low rank:   %6.2f, max eig: %04.2f\n',gsof4,eigmx4)
fprintf('\n')

% frequency response
figure
semilogx(w,20*log10(sv0(1,:)),w,20*log10(sv1(1,:)),w,20*log10(sv2(1,:)),...
         w,20*log10(sv3(1,:)),w,20*log10(sv4(1,:)))
legend('Full','Hinfstruct','Prempain','Mattei','Lmirank')

% step response
figure
t = linspace(0,40,400);
yfull = step(Gclfull,t);
ysof1 = step(Gclsof1,t);
ysof2 = step(Gclsof2,t);
ysof3 = step(Gclsof3,t);
ysof4 = step(Gclsof4,t);

for h = 1:4
    subplot(4,1,h);
    plot(t,yfull(:,h,h),t,ysof1(:,h,h),t,ysof2(:,h,h),...
        t,ysof4(:,h,h))
    ylabel('y_1')
end
xlabel('time (sec)')
legend('Full','Hinfstruct','Prempain','Lmirank')
