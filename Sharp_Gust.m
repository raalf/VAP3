clc
clear
warning off

% filename = 'inputs/X57_Cruise.vap';
% filename = 'inputs/TMotor.vap';
% filename = 'inputs/Goland_Wing.vap';
filename = 'inputs/Kussner.vap';
% alpha_seq = -12:2:8;

% for i = 1:length(alpha_seq)
    
%     VAP_IN.vecVEHALPHA = alpha_seq(i);
%     m = 10;
%     VAP_IN.valDELTIME = (1/m);
%     VAP_IN.valMAXTIME = ceil((9.5/VAP_IN.valDELTIME))+3;
    VAP_IN = [];
    [OUTP, COND, INPU, FLAG, MISC, SURF, VEHI, VISC, WAKE] = fcnVAP_MAIN(filename, VAP_IN);
    CL = OUTP.vecCL(end);
%     CM = OUTP.vecVEHCM
    CD = OUTP.vecCD(end);
%     CDi(i) = OUTP.vecCDI(end);
    
% end

AR = INPU.vecSPAN*INPU.vecSPAN/INPU.vecAREA;

% [OUTP.vecCL] = fcnUNSTEADYWRAPPER(SURF.nfree,[COND.vecVEHVINF,0,0],SURF.vecDVEHVCRD,INPU.vecAREA,SURF.gammaold,...
%     COND.valMAXTIME,COND.valGUSTSTART,COND.valDELTIME,INPU.vecSYM,COND.vecVEHBETA,SURF.en_t);

OUTP.vecCL = OUTP.vecCL.*((AR+2)/AR);

% figure(420)
% hold on
% plot(1:COND.valMAXTIME,OUTP.norm_percent,'-.k','linewidth',1.5)
% grid on
% xlabel('Time step')
% ylabel('Normal Flow Percentage')

figure(6969)
hold on
plot((COND.valGUSTSTART:COND.valMAXTIME).*COND.valDELTIME - COND.valGUSTSTART*COND.valDELTIME,OUTP.vecCL(COND.valGUSTSTART:end),'-k','linewidth',1.5)
% plot((1:COND.valMAXTIME).*COND.valDELTIME,OUTP.vecCL,'-b','linewidth',1.5)
box on
grid on
ylabel('C_l')
xlabel('Time (s)')

Kussner_Function

plot(t,Cl,'--r','linewidth',1.5)

save('Sharp_Edge_m10.mat')