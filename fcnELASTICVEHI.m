function [SURF, INPU, VEHI, OUTP, COND, FLAG] = fcnELASTICVEHI(SURF, INPU, VEHI, OUTP, COND, FLAG)

% This function computes the dynamic response of an elastic vehicle using
% Houbolt's recurrence matrix solution method

b = INPU.vecSPAN/2; % Wing half span (Why do they use b for half-span and not b/2? Who knows)

SURF.idxFLEX = find(SURF.vecDVEWING == find(FLAG.vecFLEXIBLE == 1));

% Determine non-dimensional node distribution (lambda) across wing
SURF.vecLAMBDA = [0; SURF.matCENTER(SURF.vecDVELE(SURF.idxFLEX) == 1,2)];
idx1 = 1:length(SURF.matCENTER(SURF.vecDVELE(SURF.idxFLEX) == 1,2));
idx2 = 2:length(SURF.matCENTER(SURF.vecDVELE(SURF.idxFLEX) == 1,2))+1;

SURF.vecLAMBDA = (SURF.vecLAMBDA(idx2,1) - SURF.vecLAMBDA(idx1,1))./b;

n = length(SURF.vecLAMBDA); % Number of stations across half span

% SURF.vecLAMBDA = [0.09; 0.18; 0.17; 0.16; 0.16; 0.16];

%% Form H2 matrix (eqn C10)
matH2 = cumsum(triu(repmat(SURF.vecLAMBDA(2:end)',length(SURF.vecLAMBDA(2:end)),1)),2);

%% Form H1 matrix (eqn C6/C7/C8)
% This was a nightmare. If something is broken in here, gg gents, good luck
idx = 2:n;

% Eqn C8
A = [0;(1/4)*(INPU.matEIx(1,1)./INPU.matEIx(idx-1,1)) + (1/12)*(INPU.matEIx(1,1)./INPU.matEIx(idx,1))];
B = [0;(1/3)*(INPU.matEIx(1,1)./INPU.matEIx(idx-1,1)) + (1/6)*(INPU.matEIx(1,1)./INPU.matEIx(idx,1))];
C = [0;(1/12)*(INPU.matEIx(1,1)./INPU.matEIx(idx-1,1)) + (1/12)*(INPU.matEIx(1,1)./INPU.matEIx(idx,1))];
D = [0;(1/6)*(INPU.matEIx(1,1)./INPU.matEIx(idx-1,1)) + (1/3)*(INPU.matEIx(1,1)./INPU.matEIx(idx,1))];

diagH2 = diag(matH2);

% Temp variables for A, B, C and D coefficients from C7
tempA = tril(repmat((diagH2.*diagH2.*A(2:end))',n-1,1));
tempB = repmat(B(2:end),1,length(B)-1)'.*[(diagH2(1:end-1).*matH2(2:end,:))',zeros(n-1,1)];
tempC = [zeros(n-1,1),tril(repmat((diagH2(1:end-1).*diagH2(1:end-1).*C(2:end-1))',n-1,1))];
tempD = repmat(D(2:end),1,length(D)-1)'.*[zeros(n-1,1),(diagH2(1:end-1).*matH2(2:end,:))'];

temp1 = zeros(n-1);
temp1(:,1) = temp1(:,1) + repmat(SURF.vecLAMBDA(1),n-1,1).*matH2(1,:)';

matH1 = temp1 + tempA + tempB + tempC + tempD;

%% Form A matrix --> Bending stiffness matrix

matB = (INPU.matEIx(1)/(b)^3).*inv(matH1*matH2); % Eqn C16

b0i = sum(matB,2); % Eqn C18

b00 = -sum(b0i,1); % Eqn C21

matA = zeros(n);

matA(1,1) = b00;
matA(2:end,1) = b0i;
matA(1,2:end) = b0i;
matA(2:end,2:end) = matB;

p = 100*ones(n,1);

w = matA\p;

% figure;
% plot(SURF.matCENTER(SURF.vecDVELE(SURF.idxFLEX) == 1,2),w)
% axis equal

%% Form B matrix --> Torsional stiffness matrix

j = (2./(b*SURF.vecLAMBDA(idx))).*(1./(1./INPU.matGJt(idx-1,1) + 1./INPU.matGJt(idx,1))); % Eqn C26

% Creating diagonal elements for B matrix in eqn C29
matB_1 = diag([diag(cumsum(triu(repmat(j',n,1),-1),2));j(end)]); % Yup, many functions in one line. Nothing to see here officer

% Off diagonal elements in eqn C29
matB_2 = diag(-j,1);
matB_3 = diag(-j,-1);

matB = matB_1 + matB_2 + matB_3;

q = 0.000001*ones(n,1);

phi = matB\q;

% figure;
% plot(SURF.matCENTER(SURF.vecDVELE(SURF.idxFLEX) == 1,2),phi)

%% Form coupling matrix C --> Eqn 42 and 43

matC = [matA, zeros(n); zeros(n), matB];

P = [p;q]; % Resultant vector aka external loads

% ================ TIMESTEPPING STARTS HERE ==============================
%% Form S matrices from eqn 61/62

% Eqn A3
etaCOEFF = [-2; 5; -4; 1];
eta = etaCOEFF.*repmat(INPU.vecLM(1:end-1)'./(COND.valSDELTIME*COND.valSDELTIME),4,1);

etapCOEFF = -1.*etaCOEFF;
etap = etaCOEFF.*repmat(INPU.vecLM(1:end-1)'.*SURF.vecLSM(1:end-1)'./(COND.valSDELTIME*COND.valSDELTIME),4,1);

% Eqn A4
nu = etap;

nupCOEFF = etaCOEFF;
nup = nupCOEFF.*repmat(INPU.vecJT(1:end-1)'./(COND.valSDELTIME*COND.valSDELTIME),4,1);

% Eqn 61/62
S0 = [diag(eta(1,:)), diag(etap(1,:)); diag(nu(1,:)), diag(nup(1,:))];

S1 = [diag(eta(2,:)), diag(etap(2,:)); diag(nu(2,:)), diag(nup(2,:))];

S2 = [diag(eta(3,:)), diag(etap(3,:)); diag(nu(3,:)), diag(nup(3,:))];

S3 = [diag(eta(4,:)), diag(etap(4,:)); diag(nu(4,:)), diag(nup(4,:))];

%% Initial response from initial conditions

matD = matC - S0; % Eqn 65 --> Influence matrix to update structure solution at time "n"

matLOAD = [100*ones(n,1); 10*ones(n,1)];

for valSTRUCTTIME = 1:20000

if valSTRUCTTIME <= 3
    
    OUTP.matELASTIC(:,1) = -(matD + S2 + 8*S3)\matLOAD; % Response at t = -1*eps --> Eqn 70
    OUTP.matELASTIC(:,2) = zeros(2*n,1); % Response at t = 0
    OUTP.matELASTIC(:,3) = (matD + S2 + 8*S3)\matLOAD; % Response at t = 1*eps --> Eqn 71
    
else

    vecQn = S1*OUTP.matELASTIC(:,valSTRUCTTIME-1) + S2*OUTP.matELASTIC(:,valSTRUCTTIME-2) + S3*OUTP.matELASTIC(:,valSTRUCTTIME-3); % Eqn 66 --> Resultant vector at time "n"
    OUTP.matELASTIC(:,valSTRUCTTIME) = matD\vecQn; % Elastic response of fuselage and bending/torsion of wing at t = n*eps --> Eqn 78
    
end

end
end