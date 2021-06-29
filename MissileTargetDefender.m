%% Pursuer-Evader-Defender Game
%  By- Rajeev Shobhit Voleti
%
%       M = missile
%       T = target
%       D = defender  
%       
%-------------------------------------------------------------------------------
clear;clc;clf

%-------------------------------------------------------------------------------
%                           PLAYER DYNAMICS
%
% All players have first order dynamics
%
%            a(s)         1          
%           ------  =  ----------    
%           a_c(s)    tau*s + 1    
% 
%-------------------------------------------------------------------------------

% Number of X states for each player
nx = 3;                 

% Player Time Constants, tau
tauM = 0.1;
tauT = 0.2;
tauD = 0.05;

% Missile States
AM = [0, 1,       0;
      0, 0,       1;
      0, 0, -1/tauM];
BM = [0; 0;  1/tauM];

% Target States
AT = [0, 1,       0;
      0, 0,       1;
      0, 0, -1/tauT];
BT = [0; 0;  1/tauT];

% Defender States
AD = [0, 1,       0;
      0, 0,       1;
      0, 0, -1/tauD];
BD = [0; 0;  1/tauD];

% Initial Conditions
xT0 = [0;  50; 0];
xM0 = [0; 100; 0];
xD0 = [0;  48; 0];

%-------------------------------------------------------------------------------
%                   COST FUNCTION
%-------------------------------------------------------------------------------

% Time Information
tf1  = 3;               % final time for Defender-Missile Interception
tf2  = 5;               % final time for Missile-Target Interception
dt   = 0.01;
tvec = 0:dt:tf2;
tvecSize   = numel(tvec);

tvec_tf1 = 0:dt:tf1;
timeSwitch = numel(tvec_tf1);       %Element number of tf1 in tvec

% Cost Function Data
gMT = 10e6;                                
gDM = 10e6;
GMT = diag([gMT, 0]);
GDM = diag([0, gDM]);
RM  = 1;
RT  = 5/3;
RD  = 3/5;

%-------------------------------------------------------------------------------
%                   ZERO-EFFORT MISS
%-------------------------------------------------------------------------------

% number of ZEM states
nz = 2;                                         

% STATE TRANSITION MATRICES (see Local Functions)
% z_MT(tf2)
phiM_tf2 = stateTransitionMatrix(AM, tf2, dt);
phiT_tf2 = stateTransitionMatrix(AT, tf2, dt);

% z_DM(tf1)
phiM_tf1 = stateTransitionMatrix(AM, tf1, dt);
phiD_tf1 = stateTransitionMatrix(AD, tf1, dt);

% CALCULATE dZdt INPUT MATRICES: BZ, CZ, DZ (see Local Functions)
bZhist_tf2 = zInputScalar(phiM_tf2, BM);
bZhist_tf1 = zInputScalar(phiM_tf1, BM);
cZhist = zInputScalar(phiT_tf2, BT);
dZhist = zInputScalar(phiD_tf1, BD);

BZhist = zeros(tvecSize,2);
CZhist = zeros(tvecSize,2);
DZhist = zeros(tvecSize,2);

% 0 < t <= tf1
for ct = 1:timeSwitch                                   
   BZhist(ct,:) = [ bZhist_tf2(ct), -bZhist_tf1(ct)];
   CZhist(ct,:) = [    -cZhist(ct),               0];
   DZhist(ct,:) = [              0,      dZhist(ct)];
end

% tf1 < t <= tf2
for ct = (timeSwitch+1):tvecSize                        
   BZhist(ct,:) = [ bZhist_tf2(ct), 0];
   CZhist(ct,:) = [    -cZhist(ct), 0];
   DZhist(ct,:) = [              0, 0];
end

%-------------------------------------------------------------------------------
%           DIFFERENTIAL RICATTI EQUATION
%
%       where:
%               P(tf2)  = GMT
%               P(tf1-) = P(tf1+) - GDM
%
%-------------------------------------------------------------------------------
GMTvec = reshape(GMT, [1 numel(GMT)]);
GDMvec = reshape(GDM, [1 numel(GDM)]);

BZvec = flipud(BZhist);
CZvec = flipud(CZhist);
DZvec = flipud(DZhist);

P0vec = GMTvec;
Phist = zeros(tvecSize, numel(P0vec));
Phist(1,:) = P0vec;

% Element of tf1 when going backward in time, (i.e., tf2 is element 1)
reverseSwitch = tvecSize - timeSwitch + 1;  

for ct = 1:tvecSize-1
    
    [Ttemp, Ptemp] = ode45(@RicDE, [tvec(ct) tvec(ct+1)], Phist(ct,:), ...
                     [], BZvec(ct,:), CZvec(ct,:), DZvec(ct,:), RM, RT, RD);
    
    Phist(ct+1,:) = Ptemp(end,:);
    
    % Impulse at tf1 
    if (ct + 1) == reverseSwitch
        Phist(ct+1,:) = Phist(ct+1,:) - GDMvec; 
    end
    
end

% Make Phist be ordered for t = 0 -> tf2
Phist = flipud(Phist);

%-------------------------------------------------------------------------------
%               Z-STATE HISTORIES & CONTROL INPUT HISTORIES
%-------------------------------------------------------------------------------
Zhist  = zeros(tvecSize, nz);
aMhist = zeros(tvecSize,  1);
aThist = zeros(tvecSize,  1);
aDhist = zeros(tvecSize,  1);

% Gains from Solution to DRE
GainMhist = zeros(tvecSize,  nz);
GainThist = zeros(tvecSize,  nz);
GainDhist = zeros(tvecSize,  nz);

% Effective Navigation Gains
NMhist = zeros(tvecSize, 1);
NThist = zeros(tvecSize, 1);
NDhist = zeros(tvecSize, 1);

%Effective Navigation Feed-Forward Gains
FMhist = zeros(tvecSize, 1);
FThist = zeros(tvecSize, 1);
FDhist = zeros(tvecSize, 1);

% Initial Phi(tf) matrices
phiM_tf2_0 = reshape(phiM_tf2(1,:), [nx nx]);
phiT_tf2_0 = reshape(phiT_tf2(1,:), [nx nx]);
phiM_tf1_0 = reshape(phiM_tf1(1,:), [nx nx]);
phiD_tf1_0 = reshape(phiD_tf1(1,:), [nx nx]);

% Initial ZEM Conditions
D = [1, zeros(1, nx-1)];
ZMT_0 = D*(phiM_tf2_0*xM0 - phiT_tf2_0*xT0);
ZDM_0 = D*(phiD_tf1_0*xD0 - phiM_tf1_0*xM0);
Zhist(1,:) = [ZMT_0 ZDM_0];

for ct = 1:tvecSize
    
    t = ct*dt;
    
    Ptemp = reshape(Phist(ct,:), [nz nz]);
    
    GainMhist(ct,:) = -RM\BZhist(ct,:)*Ptemp;
    GainThist(ct,:) =  RT\CZhist(ct,:)*Ptemp;
    GainDhist(ct,:) =  RD\DZhist(ct,:)*Ptemp;
    
    aMhist(ct) = GainMhist(ct,:)*Zhist(ct,:)';
    aThist(ct) = GainThist(ct,:)*Zhist(ct,:)';
    aDhist(ct) = GainDhist(ct,:)*Zhist(ct,:)';
    
    NMhist(ct) = (tf2 - t)^2*GainMhist(ct,1);
    NThist(ct) = (tf2 - t)^2*GainThist(ct,1);
    NDhist(ct) = (tf1 - t)^2*GainDhist(ct,2);
    FMhist(ct) = (tf1 - t)^2*GainMhist(ct,2);
    FThist(ct) = (tf1 - t)^2*GainThist(ct,2);
    FDhist(ct) = (tf1 - t)^2*GainDhist(ct,1);
    
    if ct < tvecSize
    
        [Ttemp, Ztemp] = ode45(@zSysDE, [tvec(ct) tvec(ct+1)], Zhist(ct,:), ...
                     [], BZhist(ct,:), CZhist(ct,:), DZhist(ct,:), ...
                     aMhist(ct), aThist(ct), aDhist(ct));
    
        Zhist(ct+1,:) = Ztemp(end, :);
    
    end
    
end

NMhist = -flipud(NMhist);
NThist = -flipud(NThist);
NDhist = -flipud(NDhist);
FMhist = -flipud(FMhist);
FThist = -flipud(FThist);
FDhist = -flipud(FDhist);


%-------------------------------------------------------------------------------
%               STATE HISTORIES
%-------------------------------------------------------------------------------
XMhist = zeros(tvecSize, nx);
XThist = zeros(tvecSize, nx);
XDhist = zeros(tvecSize, nx);

XMhist(1,:) = xM0';
XThist(1,:) = xT0';
XDhist(1,:) = xD0';

for ct = 1:tvecSize-1
   
    [TtempM, XMtemp] = ode45(@xSysDE, [tvec(ct) tvec(ct+1)], XMhist(ct,:), ...
                        [], AM, BM, aMhist(ct));
    [TtempT, XTtemp] = ode45(@xSysDE, [tvec(ct) tvec(ct+1)], XThist(ct,:), ...
                        [], AT, BT, aThist(ct));
    [TtempD, XDtemp] = ode45(@xSysDE, [tvec(ct) tvec(ct+1)], XDhist(ct,:), ...
                        [], AD, BD, aDhist(ct));
                
    XMhist(ct+1,:) = XMtemp(end,:);
    XThist(ct+1,:) = XTtemp(end,:);
    XDhist(ct+1,:) = XDtemp(end,:);
    
end

%-------------------------------------------------------------------------------
%               Simulate with Initial Ranges
%------------------------------------------------------------------------------- 
range_M0 = 5000;
range_T0 = 0;
range_D0 = range_T0;

% Closing Velocities
V_MT = (range_M0 - range_T0)/tf2;
V_M = -V_MT/2;
V_T =  V_MT/2;

V_D = (range_M0 - range_D0)/tf1 + V_M;

range_M = range_M0 + V_M*tvec;
range_T = range_T0 + V_T*tvec;
range_D = range_D0 + V_D*tvec_tf1;

%-------------------------------------------------------------------------------
%               DISPLACEMENTS BETWEEN PLAYERS
%-------------------------------------------------------------------------------

YMThist = XMhist(:,1) - XThist(:,1);
YDMhist = XMhist(:,1) - XDhist(:,1);

%-------------------------------------------------------------------------------
%               PLOTS
%-------------------------------------------------------------------------------

figure(1)
plot(tvec, Zhist);
title('Zero Effort Miss Terminal Projection');
legend('Z_{MT}', 'Z_{DM}')
xlabel('Time (s)');
ylabel('Zero-Effort-Miss (m)') 

figure(2)
subplot(2,1,1),plot(tvec, aMhist, tvec, aThist, tvec, aDhist), ...
    title('Commanded Acceleration'), ...
    legend({'Missile', 'Target', 'Defender'},'Location','southeast'), ...
    ylim([-70 5]),ylabel('a_c (m/s^2)')

subplot(2,1,2),plot(tvec,XMhist(:,3),tvec,XThist(:,3),tvec,XDhist(:,3)), ...
    title('Actual Accelerations'), ...
    legend({'Missile', 'Target', 'Defender'},'Location','southeast'), ...
    ylim([-70 5]),xlabel('Time (s)'),ylabel('a (m/s^2)')
    

figure(3)
hold on
plot(tvec,XMhist(:,1),tvec,XThist(:,1))
plot(tvec_tf1,XDhist(1:timeSwitch,1)), ...
    title('Displacement Perpendicular to Reference Direction'), ...
    legend({'Missile','Target','Defender'},'Location', 'northwest')
    ylabel('y (m)')
    
figure(4)
hold on
plot(range_M,XMhist(:,1));
plot(range_T,XThist(:,1));
plot(range_D,XDhist(1:timeSwitch,1));
ylabel('Player Perpendicular Separation (m)');
xlabel('Range (m)');
title('Trajectories of the Target, Missile. and Defender');
legend({'Missile','Target','Defender'},'Location', 'northwest');
hold off

figure(5)
plot(tvec_tf1,YMThist(1:timeSwitch),tvec_tf1,YDMhist(1:timeSwitch)), ...
    title('Displacement Difference Between Players'), ...
    legend({'\Delta Y_{MT}', '\Delta Y_{DM}'},'Location','northwest'), ...
    xlabel('Time (s)'), ylabel('\Delta y (m)')

figure(6)
sgtitle({'Effective Navigation Gains and', ...
          'Effective Navigation Feed-forward Gains'})
subplot(3,2,1),plot(tvec,NThist), ylim([0 30]), grid on, ...
    ylabel('N_T'), title('Target gain to \sigma_{TM}')

subplot(3,2,2),plot(tvec,FThist), ylim([-2 2]), grid on, ...
    ylabel('\Gamma_T'), title('Target gain to \sigma_{DM}')

subplot(3,2,3),plot(tvec,NMhist), ylim([0 30]), grid on, ...
    ylabel('N_M'), title('Missile gain to \sigma_{TM}')

subplot(3,2,4),plot(tvec,FMhist), ylim([-5 15]), xlim([0 5]), grid on, ...
    ylabel('\Gamma_M'), title('Missile gain to \sigma_{DM}')

subplot(3,2,5),plot(tvec,FDhist), ylim([0 10]), grid on, ...
    ylabel('\Gamma_D'), title('Defender gain to \sigma_{TM}'), ...
    xlabel('time-to-go (sec)')

subplot(3,2,6),plot(tvec,NDhist), ylim([0 30]), grid on, ...
    ylabel('N_D'), title('Defender gain to \sigma_{DM}'), ...
    xlabel('time-to-go (sec)')




%-------------------------------------------------------------------------------
%               LOCAL FUNCTIONS
%-------------------------------------------------------------------------------

function phi = stateTransitionMatrix(A, tf, dt)
% Computes State Transition matrix for state matrix A over tvec = t0:dt:tf
% Output in form where each row is corresponding time containing all 
% elements. Must reshape each row to (n-by-n) matrix for use later.
    
    tvec = 0:dt:tf;
    phi  = zeros(numel(tvec), numel(A));
    
    for i = 1:numel(tvec)
        t = tvec(i);
        phiMat = expm(A*(tf-t));
        phi(i,:) = reshape(phiMat, [1 numel(phiMat)]);
    end
    
end

function zInputScalarHist = zInputScalar(phi, G)
% Computes the Input Matrices(BZ, CZ, or DZ) for Zdot.
% G = B*u from Xdot = A*x + B*u

    nx = size(G,1);
    D  = [1, zeros(1, nx-1)];
    zInputScalarHist = zeros(size(phi,1), 1);
    
    for i = 1:size(phi,1)
        phiMat = reshape(phi(i,:), [nx nx]); 
        inputScalar = D*phiMat*G;
        zInputScalarHist(i) = inputScalar;
    end
    
end

function PdVec = RicDE(~, Pvec, BzVec, CzVec, DzVec, RM, RT, RD)
    
    nx = sqrt(numel(Pvec));
    P  = reshape(Pvec,  [nx nx]);
    B  = reshape(BzVec, [nx  1]);
    C  = reshape(CzVec, [nx  1]);
    D  = reshape(DzVec, [nx  1]);
    
    
    Pdot  = -P*(B*(RM\B') - C*(RT\C') - D*(RD\D'))*P;
    PdVec = reshape(Pdot, [numel(Pdot) 1]);

end

function ZdVec = zSysDE(~, ~, BZvec, CZvec, DZvec, aMc, aTc, aDc)
    
    BZ = BZvec';
    CZ = CZvec';
    DZ = DZvec';
    ZdVec = BZ*aMc + CZ*aTc + DZ*aDc;
end

function XdVec = xSysDE(~, X, A, B, u)
    XdVec = A*X + B*u;
end










  