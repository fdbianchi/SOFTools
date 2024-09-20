
% Example of optimal design of a static output feedback gain.
%
% Case LPV: quadruple-tank process
%
% Example from:
% F.Bianchi, R.Mantz, C.Christiansen, Multivariable PID control with 
% set-point weighting via BMI optimisation, Automatica, Volume 44, Issue 2,
% 2008, pp 472-478, doi:10.1016/j.automatica.2007.05.021.systems,”

% fbianchi - 2024-09-19 - rev


% cleaning
clearvars
clc
close all


% Operating points
Vset = [1.6 3.6];   
pv = pset.Box(Vset);

% plant parameters
A1 = 28;    A3 = A1;
A2 = 32;    A4 = A2;
a1 = 0.071; a3 = a1; 
a2 = 0.057; a4 = a2; 
kc = 0.50;
g  = 981;
k1 = 3.33;  k2 = 3.35;
g1 = 0.70;  g2 = 0.60;

% tank heights
for ii = 1:2
    
    v0 = Vset(ii);
    h10 = ((g1*k1*v0+(1-g2)*k2*v0)/a1/sqrt(2*g))^2;
    h20 = (((1-g1)*k1*v0+g2*k2*v0)/a2/sqrt(2*g))^2;
    h30 = ((1-g2)/a3/sqrt(2*g)*k2*v0)^2;
    h40 = ((1-g1)/a4/sqrt(2*g)*k1*v0)^2;
    
    % time constants
    T1 = A1/a1*sqrt(2*h10/g); T3 = A3/a3*sqrt(2*h30/g);
    T2 = A2/a2*sqrt(2*h20/g); T4 = A4/a4*sqrt(2*h40/g);
    
    % model matrices
    A(:,:,ii) = [-1/T1   0      A3/A1/T3    0;
                  0     -1/T2   0           A4/A2/T4;
                  0      0     -1/T3        0;
                  0      0      0          -1/T4];
end
B = [g1*k1/A1       0;
     0              g2*k2/A2;
     0              (1-g2)*k2/A3;
     (1-g1)*k1/A1   0];
C = [kc  0   0   0;0  kc   0   0]; 
D = zeros(2);
pdG = ppss(A,B,C,D,pv);
pdG.y = 'y';
pdG.u = 'v';

% augmented plant
sb1 = sumblk('v = u + d',2);
sb2 = sumblk('e = r - y',2);
Gi = ss(zeros(2),eye(2),eye(2),zeros(2));
Gi.u = 'e'; Gi.y = 'ei';
pdGau = connect(pdG,Gi,sb1,sb2,...
                {'d','r','u'},...
                {'u','ei','r','ei','y'});
Gau = ss(pdGau);
            
% weigths
As = 20/100; wb = 100;
wui = tf(As*[1/(wb*0.1) 1],[1/(wb*10) 1]);
Wui = tf(As,1);
Wu  = append(Wui,Wui);
We  = ss(eye(2)*0.4);
Wi1 = ss(0.95*eye(2));
Wi2 = ss(0.05*eye(2));

Win = append(Wi1, Wi2, eye(2));
Wout = append(Wu,We,eye(6));
pdGaw = Wout*pdGau*Win;
pdGaw.u = {'dw(1)','dw(2)','rw(1)','rw(2)','u(1)','u(2)'};
pdGaw.y = [{'uw(1)','uw(2)','ew(1)','ew(2)'},...
           sprintfc('uc(%d)',1:6)];

% i/o map
ios.perf = pdGaw.y(1:4);
ios.dist = pdGaw.u(1:4);
ios.meas = pdGaw.y(5:10);
ios.ctrl = pdGaw.u(5:6);

% full order

[pdKfull,~,gfull] = lpvsyn(pdGaw,[6 2]); 
Kfull = ss(pdKfull);

% using lowrank
opt = sofsettings('type','gs','minDecay',0.01);
[Ksof1,gsof1] = sofsyn(pdGaw,ios,'lowrank',opt);
Ksof1 = ss([],[],[],Ksof1);

% using systune
[Ksof2,gsof2] = sofsyn(pdGaw,ios,'systune',opt);
Ksof2 = ss([],[],[],Ksof2);

% -----------------------------------------------------------------------
% checking results

sb4 = sumblk('v = u',2);
sb3 = sumblk('e = r - y',2);
pdGau = connect(pdG,Gi,sb3,sb4,...
                {'r','u'},...
                {'y','r','ei','y'});
Gau = ss(pdGau);
Gclfull = lft(Gau,Kfull);
eigmx0 = max(real(eig(Gclfull)));
Gclsof1 = lft(Gau,Ksof1);
eigmx1 = max(real(eig(Gclsof1)));
Gclsof2 = lft(Gau,Ksof2);
eigmx2 = max(real(eig(Gclsof2)));

% closed loop eigenvalues
fprintf('\n')
fprintf('Performance comparison\n')
fprintf('Full order: %6.2f, max eig: %04.2f\n',gfull,eigmx0(:,:,1))
fprintf('Low rank:   %6.2f, max eig: %04.2f\n',gsof1,eigmx1(:,:,1))
fprintf('Systune:    %6.2f, max eig: %04.2f\n',gsof2,eigmx2(:,:,1))
fprintf('\n')

% step response
np = 1000;
yfull(np,2,2) = 0;
ysof1(np,2,2) = 0;
ysof2(np,2,2) = 0;

t = linspace(0,100,np);
u(:,1) = max(0,min(t-5,1));
u(:,2) = max(0,min(t-5,1));

for ii = 1:2
    yfull(:,:,ii) = lsim(Gclfull(:,:,ii),u,t);
    ysof1(:,:,ii) = lsim(Gclsof1(:,:,ii),u,t);
    ysof2(:,:,ii) = lsim(Gclsof2(:,:,ii),u,t);
end

% plots
clines = lines(3);
figure

h = 1;
subplot(2,1,h);
plot(t,u(:,h),'Color',0.5*[1 1 1]);
hold on
plot(t,yfull(:,h),'Color',clines(1,:));
plot(t,ysof1(:,h),'Color',clines(2,:));
plot(t,ysof2(:,h),'Color',clines(3,:));
ylabel('y_1')

h = 2;
subplot(2,1,h);
plot(t,u(:,h),'Color',0.5*[1 1 1]);
hold on
plot(t,yfull(:,h),'Color',clines(1,:));
plot(t,ysof1(:,h),'Color',clines(2,:));
plot(t,ysof2(:,h),'Color',clines(3,:));
ylabel('y_2')

xlabel('time (sec)')
legend('Ref','Full','Lowrank','Systune')


figure
sigma(Gclfull,Gclsof1,Gclsof2)
legend('Full','Lowrank','Systune')



return

% =========================================================================
% Diseño del controlador

tau=1/100; beta = 0.01;

tip = 'rob'; gam = 11.5/11.5;
tip = 'gs';  gam = 10;

% -------------------------------------------------------------------------
% funciona
As = 1/100; wb = 30; we  = eye(nc)*0.4; wem  = eye(nc)*0.8;
% -------------------------------------------------------------------------

% funciones de peso
As = 20/100; wb = 100;
wui = ltisys('tf',As*[1/(wb*0.1) 1],[1/(wb*10) 1]);
wui = ltisys('tf',As,1);
wu  = sdiag(wui,wui);
we  = eye(nc)*0.4;
wi1 = 0.95*eye(2);
wi2 = 0.05*eye(2);

% PID 2DOF: método propuesto
tic;
Ky_lpv = pid_2dof_lpv(Gplpv,[2 2],tau,gam,'pi','gs',wi1,wi2,we,wu,beta);
toc
if strcmp(tip,'rob')
    Ky_f = psinfo(Ky_lpv,'sys',1);
else
    alpha = (Vset(2)-v0)*Vset(1)/(Vset(2)-Vset(1))/v0;
    Ky_f = psinfo(Ky_lpv,'eval',[alpha 1-alpha]);
end
[xx,xx,xx,Ky] = ltiss(Ky_f);

% PID 1DOF: Mattei
As = 20/100; wb = 10000;
wuim = ltisys('tf',As*[1/(wb*0.1) 1],[1/(wb*10) 1]);
wuim = ltisys('tf',As,[1/(wb) 1]);
wum  = sdiag(wuim,wuim);
wem  = eye(nc)*0.2;
% wum  = 0*eye(nc)*0.2;
% Ky_mt = pid_mattei(Gplpv,[2 2],eye(2),gam,'pi',wem,wum,beta);
Ky_lpv0 = pid_1dof_lpv(Gplpv,[2 2],tau,gam,'pi','gs',we,wu,beta);
if strcmp(tip,'rob')
    Ky_f = psinfo(Ky_lpv,'sys',1);
else
    alpha = (Vset(2)-v0)*Vset(1)/(Vset(2)-Vset(1))/v0;
    Ky_f = psinfo(Ky_lpv0,'eval',[alpha 1-alpha]);
end
[xx,xx,xx,Ky_mt] = ltiss(Ky_f);
Ky0   = [zeros(2),zeros(2),Ky_mt(:,1:2),zeros(2),Ky_mt(:,5:6)];

% =========================================================================
% sintesis controlador GS_Hinf
m = ltisys('tf',1,[1 0]);
M = sdiag(m,m);
we  = eye(nc)*0.4;

inputs  = 'r(2);w(2)';
outputs = 'M;K';
Kin     = 'K:r;M';
Plpvin  = 'Plpv:w+K';
Min     = 'M:r-Plpv'
[Glpv,rn] = sconnect(inputs,outputs,Kin,Plpvin,Plpv,Min,M);

% Glpv = smult(sdiag(wi1,wi2,eye(2)),Glpv,sdiag(we,wu,eye(2)));
Glpv = smult(sdiag(wi1,wi2,eye(2)),Glpv,sdiag(we,wu,eye(4)));

[ggs,Kgs] = hinfgs(Glpv,rn); ggs
[ak,bk,ck,dk] = ltiss(psinfo(smult(sdiag(eye(2),M),Kgs),'eval',[alpha 1-alpha]));


return


% =========================================================================
% Diseño del controlador

tau=1/100; beta = 0.01;

tip = 'rob'; gam = 11.5/11.5;
tip = 'gs';  gam = 10;

% -------------------------------------------------------------------------
% funciona
As = 1/100; wb = 30; we  = eye(nc)*0.4; wem  = eye(nc)*0.8;
% -------------------------------------------------------------------------

% funciones de peso
As = 20/100; wb = 100;
wui = ltisys('tf',As*[1/(wb*0.1) 1],[1/(wb*10) 1]);
wui = ltisys('tf',As,1);
wu  = sdiag(wui,wui);
we  = eye(nc)*0.4;
wi1 = 0.95*eye(2);
wi2 = 0.05*eye(2);

% PID 2DOF: método propuesto
tic;
Ky_lpv = pid_2dof_lpv(Gplpv,[2 2],tau,gam,'pi','gs',wi1,wi2,we,wu,beta);
toc
if strcmp(tip,'rob')
    Ky_f = psinfo(Ky_lpv,'sys',1);
else
    alpha = (V(2)-v0)*V(1)/(V(2)-V(1))/v0;
    Ky_f = psinfo(Ky_lpv,'eval',[alpha 1-alpha]);
end
[xx,xx,xx,Ky] = ltiss(Ky_f);

% PID 1DOF: Mattei
As = 20/100; wb = 10000;
wuim = ltisys('tf',As*[1/(wb*0.1) 1],[1/(wb*10) 1]);
wuim = ltisys('tf',As,[1/(wb) 1]);
wum  = sdiag(wuim,wuim);
wem  = eye(nc)*0.2;
% wum  = 0*eye(nc)*0.2;
% Ky_mt = pid_mattei(Gplpv,[2 2],eye(2),gam,'pi',wem,wum,beta);
Ky_lpv0 = pid_1dof_lpv(Gplpv,[2 2],tau,gam,'pi','gs',we,wu,beta);
if strcmp(tip,'rob')
    Ky_f = psinfo(Ky_lpv,'sys',1);
else
    alpha = (V(2)-v0)*V(1)/(V(2)-V(1))/v0;
    Ky_f = psinfo(Ky_lpv0,'eval',[alpha 1-alpha]);
end
[xx,xx,xx,Ky_mt] = ltiss(Ky_f);
Ky0   = [zeros(2),zeros(2),Ky_mt(:,1:2),zeros(2),Ky_mt(:,5:6)];

% =========================================================================
% sintesis controlador GS_Hinf
m = ltisys('tf',1,[1 0]);
M = sdiag(m,m);
we  = eye(nc)*0.4;

inputs  = 'r(2);w(2)';
outputs = 'M;K';
Kin     = 'K:r;M';
Plpvin  = 'Plpv:w+K';
Min     = 'M:r-Plpv'
[Glpv,rn] = sconnect(inputs,outputs,Kin,Plpvin,Plpv,Min,M);

% Glpv = smult(sdiag(wi1,wi2,eye(2)),Glpv,sdiag(we,wu,eye(2)));
Glpv = smult(sdiag(wi1,wi2,eye(2)),Glpv,sdiag(we,wu,eye(4)));

[ggs,Kgs] = hinfgs(Glpv,rn); ggs
[ak,bk,ck,dk] = ltiss(psinfo(smult(sdiag(eye(2),M),Kgs),'eval',[alpha 1-alpha]));


return
% Simulación en varios puntos
v0 = 1.6;
sys = cs2lmi(tanques_nl(v0,v0));
[A,B,C,D] = ltiss(sys);
if strcmp(tip,'rob')
    Ky_f = psinfo(Ky_lpv,'sys',1);
    Ky_f0 = psinfo(Ky_lpv0,'sys',1);
else
    alpha = (V(2)-v0)*V(1)/(V(2)-V(1))/v0;
    Ky_f = psinfo(Ky_lpv,'eval',[alpha 1-alpha]);
    Ky_f0 = psinfo(Ky_lpv0,'eval',[alpha 1-alpha]);
end
[xx,xx,xx,Ky] = ltiss(Ky_f);
[ak,bk,ck,dk] = ltiss(psinfo(smult(sdiag(eye(2),M),Kgs),'eval',[alpha 1-alpha]));
[xx,xx,xx,Ky_mt] = ltiss(Ky_f0);
Ky0   = [zeros(2),zeros(2),Ky_mt(:,1:2),zeros(2),Ky_mt(:,5:6)];
% ~~~~~~~~~~~~~~~~
v0 = 2.6;
sys = cs2lmi(tanques_nl(v0,v0));
[A,B,C,D] = ltiss(sys);
if strcmp(tip,'rob')
    Ky_f = psinfo(Ky_lpv,'sys',1);
    Ky_f0 = psinfo(Ky_lpv0,'sys',1);
else
    alpha = (V(2)-v0)*V(1)/(V(2)-V(1))/v0;
    Ky_f = psinfo(Ky_lpv,'eval',[alpha 1-alpha]);
    Ky_f0 = psinfo(Ky_lpv0,'eval',[alpha 1-alpha]);
end
[xx,xx,xx,Ky] = ltiss(Ky_f);
[ak,bk,ck,dk] = ltiss(psinfo(smult(sdiag(eye(2),M),Kgs),'eval',[alpha 1-alpha]));
[xx,xx,xx,Ky_mt] = ltiss(Ky_f0);
Ky0   = [zeros(2),zeros(2),Ky_mt(:,1:2),zeros(2),Ky_mt(:,5:6)];
% ~~~~~~~~~~~~~~~~
v0 = 3.6;
sys = cs2lmi(tanques_nl(v0,v0));
[A,B,C,D] = ltiss(sys);
if strcmp(tip,'rob')
    Ky_f = psinfo(Ky_lpv,'sys',1);
    Ky_f0 = psinfo(Ky_lpv0,'sys',1);
else
    alpha = (V(2)-v0)*V(1)/(V(2)-V(1))/v0;
    Ky_f = psinfo(Ky_lpv,'eval',[alpha 1-alpha]);
    Ky_f0 = psinfo(Ky_lpv0,'eval',[alpha 1-alpha]);
end
[xx,xx,xx,Ky] = ltiss(Ky_f);
[ak,bk,ck,dk] = ltiss(psinfo(smult(sdiag(eye(2),M),Kgs),'eval',[alpha 1-alpha]));
[xx,xx,xx,Ky_mt] = ltiss(Ky_f0);
Ky0   = [zeros(2),zeros(2),Ky_mt(:,1:2),zeros(2),Ky_mt(:,5:6)];


return
% -------------------------------------------------------------------------


% % pesos de performance
% wu = diag([0.3,0.15]);
% wu = diag([0.2 0.2]);
% we = eye(nc)*0.4;
% Ky_lpv = pid_2dof_lpv(Gplpv,[2 2],tau,gam,'pi','gs',we,wu,beta);
% 
% % Ky=pid_1dof(Gp,[2 2],10000,'pi');
% if strcmp(tip,'rob')
%     Ky_f = psinfo(Ky_lpv,'sys',1);
% else
%     alpha = (V(2)-v0)*V(1)/(V(2)-V(1))/v0;
%     Ky_f = psinfo(Ky_lpv,'eval',[alpha 1-alpha]);
% end
% [xx,xx,xx,Ky] = ltiss(Ky_f);
% 
% % mattei
% we  = eye(nc)*0.2;
% Ky_mt=pid_mattei(Gplpv,[2 2],eye(2),'pi',we,wu,0.01);
% Ky0 = [zeros(2),zeros(2),-Ky_mt(:,1:2),zeros(2),-Ky_mt(:,3:4)];
% % Ky0 = [-Ky_mt(:,1:2),zeros(2),-Ky_mt(:,1:2),zeros(2),-Ky_mt(:,3:4)];
% 
% % Valores de Johansson
% % Ky0 = [diag([3.0 2.7]),zeros(2),diag([3.0 2.7]),zeros(2),diag([3.0/30 2.7/40])];

gam = 10;
% PID 2DOF: método propuesto
As = 0.3/100; wb = 10;
wui = ltisys('tf',As*[1/(wb*0.1) 1],[1/(wb*10) 1]);
wu  = sdiag(wui,wui);
we  = eye(nc)*0.4;
wi1 = 0.95*eye(2);
wi2 = 0.05*eye(2);
Ky_lpv = pid_2dof_lpv(Gplpv,[2 2],tau,gam,'pi','gs',wi1,wi2,we,wu,beta);
if strcmp(tip,'rob')
    Ky_f = psinfo(Ky_lpv,'sys',1);
else
    alpha = (V(2)-v0)*V(1)/(V(2)-V(1))/v0;
    Ky_f = psinfo(Ky_lpv,'eval',[alpha 1-alpha]);
end
[xx,xx,xx,Ky] = ltiss(Ky_f);

% PID 1DOF
Ky_mt=pid_mattei(Gplpv,[2 2],eye(2),gam,'pi',we,wu,beta);
Ky0 = [zeros(2),zeros(2),-Ky_mt(:,1:2),zeros(2),-Ky_mt(:,3:4)];


% =========================================================================
% sintesis controlador GS_Hinf
m = tf(1,[1 0]);
M = 1*append(m,m);
Wu= diag([0.5,0.5]);
Gii = [];
for ii=1:2
    P = lmi2cs(psinfo(Plpv,'sys',ii));
    % Interconección de la planta aumentada
    warning off
    systemnames  = 'P M Wu';
    inputvar     = '[r{2};w{2};u{2}]';
    outputvar    = '[M;Wu;M]';
    input_to_M   = '[r-P]';
    input_to_P   = '[w+u]';
    input_to_Wu  = '[u]';
    Gii = [Gii cs2lmi(sysic)];
    warning on
end
Glpv = psys(Gii);
we = eye(nc)*0.4;
Glpv = smult(Glpv,sdiag(we,eye(4)));

[ggs,Kgs] = hinfgs(Glpv,[2,2]); ggs
[ak,bk,ck,dk] = ltiss(psinfo(smult(cs2lmi(M),Kgs),'eval',[alpha 1-alpha]));
