
% Example of optimal design of a static output feedback gain.
%
% Case LPV: 
%
% Example from:
% M. Mattei, “Robust multivariable PID control for linear parameter varying
% systems,” Automatica, vol. 37, no. 12, pp. 1997–2003, Dec. 2001.
%

% fbianchi - 2015-03-25
% fbianchi - 2024-08-02 - rev


% cleaning
clearvars
close all
clc

% Integrator structure
Gi = [1 0 0;0 1 0];

% plant dimensions
ns = 4; 
nu = 2; ny = 3;
nw = 2; nz = 2;
ni = 2; n  = ns + ni;

% plant model
A(:,:,1) = [-4.0 -0.4 -4.4  4.9; 
         1.5 -4.1  3.5  0.6; 
         1.1 -4.0 -0.5 -2.8; 
         3.1  3.7 -2.5  3.8];
A(:,:,2) = [ 2.3 -0.3  1.8 -2.6; 
         3.5  2.3  3.8  0.7; 
        -3.2  4.8  0.1  2.6; 
         0.4  1.7 -0.3 -1.5];
A(:,:,3) = [-0.5 -2.0 -0.7  2.5;
         0.2 -3.5  3.3  1.2;
        -3.4  4.4 -3.3 -2.5;
         4.8  1.7 -0.6 -1.5];     
Bw(:,:,1) = [-1.0 -2.8; 
         4.8 -3.6; 
        -2.0  4.8; 
         0.7 -2.9];
Bw(:,:,2) = [-0.4 -0.1;
        -0.4  0.0;
        -0.2  0.0; 
         0.3 -0.4];
Bw(:,:,3) = [ 0.5 -0.3;
        -0.1 -0.3;
        -0.2 -0.2;
         0.4  0.2];
Bu(:,:,1) = [ 1.1  0.2; 
        -1.6  1.3; 
         1.4  4.5; 
         0.2 -2.8];
Bu(:,:,2) = [ 4.9  4.2;
         3.8 -3.7;
         1.3  4.7;
        -2.9  1.2];     
Bu(:,:,3) = [ 1.4 -5.0;
        -4.9  0.2;
         2.9  1.2;
         4.0 -2.3];     
Cz(:,:,1) = [ 0.10 -0.13 -0.200 -0.030; 
        -0.24 -0.05  0.019 -0.095];
Cz(:,:,2) = [0.15  0.00 0.20 -0.05;
        0.10 -0.15 0.15  0.0];    
Cz(:,:,3) = [ 0.05  0.15  0.05 -0.05;
        -0.10 -0.05 -0.10 -0.15];    
Cy = cat(3,[ 1.0  3.3 -1.1 -2.2; 
        -1.1 -4.5 -0.5 -2.1; 
         4.4  3.3 -4.9 -0.7],zeros(3,4,2));
Dzw(:,:,1) = [-0.175 0.205; 
        -0.145 0.040];
Dzw(:,:,2) = [ 0.1  0.0;
         0.2 -0.1];
Dzw(:,:,3) = [-0.05 0.0;
        -0.20 0.25];     
Dzu = zeros(2,2,3); 
Dyw = zeros(3,2,3);
Dyu(:,:,1) = [-1.7 -3.4; 
        -1.6 -3.4; 
         3.0 -3.4];
Dyu(:,:,2) = [ 0.4 -0.2;
         0.3  0.1;
        -0.1  0.5];    
Dyu(:,:,3) = [ 0.1 -0.5;
        -0.3  0.3;
         0.3 -0.3];
% parameter set
pv = pset.Box([-0.1 0.1;-0.1 0.1]);

% LPV model
pdG = pass(A,[Bw Bu],[Cz; Cy],[Dzw Dzu; Dyw Dyu], pv);
pdG.y = [sprintfc('z(%d)',1:2), sprintfc('y(%d)',1:3)];
pdG.u = [sprintfc('w(%d)',1:2), sprintfc('u(%d)',1:2)];

% augmented plant: plant + PI structure
Gi = ss(zeros(2),eye(2,3),[zeros(3,2); eye(2)],...
        [eye(3); zeros(2,3)]);
Gi.y = 'v';
Gi.u = 'y';

in = [sprintfc('w(%d)',1:2) sprintfc('u(%d)',1:2)];
out = [sprintfc('z(%d)',1:2) sprintfc('v(%d)',1:5)];
pdGau = connect(pdG,Gi,in,out);
Gau = ss(pdGau);

% Note: the PI design ignores Dyu in the desing of Ky and it is considers
% in the PI implementation
pdGdes = pdGau;
Daux = pdGdes.D;
Daux(3:7,3:4,:) = 0;
pdGdes.D = Daux;

% io map
ios.perf = pdGdes.y(1:2);
ios.meas = pdGdes.y(3:7);
ios.dist = pdGdes.u(1:2);
ios.ctrl = pdGdes.u(3:4);

% control designs
opt = sofsettings('minDecay',0.4,'beta',1e0);
[Ky1,gopt1] = sofsyn(pdGdes, ios, 'mattei', opt);

[Ky2,gopt2] = sofsyn(pdGdes, ios, 'crucius', opt);

[Ky3,gopt3] = sofsyn(pdGdes, ios, 'systune', opt);

% gain from paper
Ky4 = [15.36  13.39 -3.85  0.0057  0.081;
      -17.06 -15.01  4.50 -0.0020 -0.091];


Pset = pv.points;
for ii = 1:size(Pset,2)
    % Dyu(p)
    p = Pset(:,ii);
    dyui = Dyu(:,:,1) + Dyu(:,:,2)*p(1) + Dyu(:,:,3)*p(2);
    
    % PI for 1
    k1 = Ky1(:,1:3); k2 = Ky1(:,4:5);
    Dk = eye(3) + dyui*k1;
    Kpi = [k1/Dk (eye(2) - k1/Dk*dyui)*k2];
    Kpi1(:,:,ii) = Kpi*Gi;
    
    % PI for 2
    k1 = Ky2(:,1:3); k2 = Ky2(:,4:5);
    Dk = eye(3) + dyui*k1;
    Kpi = [k1/Dk (eye(2) - k1/Dk*dyui)*k2];
    Kpi2(:,:,ii) = Kpi*Gi;
    
    % PI for 3
    k1 = Ky3(:,1:3); k2 = Ky3(:,4:5);
    Dk = eye(3) + dyui*k1;
    Kpi = [k1/Dk (eye(2) - k1/Dk*dyui)*k2];
    Kpi3(:,:,ii) = Kpi*Gi;
    
    % PI for 4
    k1 = Ky4(:,1:3); k2 = Ky4(:,4:5);
    Dk = eye(3) + dyui*k1;
    Kpi = [k1/Dk (eye(2) - k1/Dk*dyui)*k2];
    Kpi4(:,:,ii) = Kpi*Gi;

end

fprintf('\nResults:\n')

G = ss(pdG);
Gcl1 = lft(G,Kpi1);
E1 = eig(Gcl1);
fprintf('Using Mattei:\n')
for ii = 1:size(Pset,2)
    fprintf(' min(real(eig_cl)) = %d,\tmax(real(eig_cl)) = %d\n',...
        min(real(E1(:,:,ii))),max(real(E1(:,:,ii))));
end
fprintf(' performance = %6.3f\n',gopt1);

Gcl2 = lft(G,Kpi2);
E2 = eig(Gcl2);
fprintf('Using Crucius:\n')
for ii = 1:size(Pset,2)
    fprintf(' min(real(eig_cl)) = %d,\tmax(real(eig_cl)) = %d\n',...
        min(real(E2(:,:,ii))),max(real(E2(:,:,ii))));
end
fprintf(' performance = %6.3f\n',gopt2);

Gcl3 = lft(G,Kpi3);
E3 = eig(Gcl3);
fprintf('Using Systune:\n')
for ii = 1:size(Pset,2)
    fprintf(' min(real(eig_cl)) = %d,\tmax(real(eig_cl)) = %d\n',...
        min(real(E3(:,:,ii))),max(real(E3(:,:,ii))));
end
fprintf(' performance = %6.3f\n',gopt3);

Gcl4 = lft(G,Kpi4);
E3 = eig(Gcl4);
fprintf('Gain from paper:\n')
for ii = 1:size(Pset,2)
    fprintf(' min(real(eig_cl)) = %d,\tmax(real(eig_cl)) = %d\n',...
        min(real(E3(:,:,ii))),max(real(E3(:,:,ii))));
end

figure
step(Gcl1,Gcl2,Gcl3)
legend('Mattei','Crucius','Systune')



return
     






