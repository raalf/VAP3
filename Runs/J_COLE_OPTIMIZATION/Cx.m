% Calculate CX
clear,clc
load('CxSweep.mat')

vinf = 77.2;
rho = 1.007;
weightN = 13344.6648; % 3000 lb in N
S    = OUTP(1).valAREA; % ref. platform area
CL   = weightN./(0.5*rho*vinf.^2*S);

ALPHA = interp1([OUTP.vecCL_AVG],[OUTP.vecVEHALPHA],CL);
CDvap = interp1([OUTP.vecVEHALPHA],[OUTP.vecCD_AVG],ALPHA);

CDnac = zeros(size(OUTP,2),1);
for i =1:size(OUTP,2) 
    % V_inf = 80 here because that is what was run in VAP
    CDnac(i,:) = OUTP(i).NACELLE.vecNACDRAG/(0.5*rho*(80^2)*S); 
end

CDnac = interp1([OUTP.vecVEHALPHA],CDnac,ALPHA);

LD = 14;
CD = CL./(LD); % Calculate CD with Borer L/D Data
CX = CD - CDvap - CDnac;
fprintf('CX without props and without nacelles: 0.033\n')
fprintf('CX with props and without nacelles: %.4f\n',CD-CDvap)
fprintf('CX with props and with nacelles: %.4f\n',CX)