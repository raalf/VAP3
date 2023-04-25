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

[ledves, ~, ~] = find(SURF.vecDVELE > 0 & SURF.vecDVEWING > 0);
lepanels = SURF.vecDVEPANEL(ledves);

isCurWing = SURF.vecDVEWING(ledves) == 1;

idxdve = uint16(ledves(isCurWing));
idxpanel = lepanels(isCurWing);


m = INPU.vecM(idxpanel);
m = m(1);

% Matrix of how much we need to add to an index to get the next chordwise element
% It is done this way because n can be different for each panel. Unlike in the wake,
% we can't just add a constant value to get to the same spanwise location in the next
% row of elements
tempm = repmat(INPU.vecN(idxpanel), 1, m).*repmat([0:m-1],length(idxpanel~=0),1);

rows = repmat(idxdve,1,m) + uint16(tempm);

nforce = (SURF.nfree(rows(5,:),:) + SURF.nind(rows(5,:),:));
vecAREA = sum(sum(SURF.vecDVEAREA(rows(5,:)),2),1);
qinf = 0.5*COND.vecVEHVINF*COND.vecVEHVINF;
CL = sum(nforce,1)./(qinf*vecAREA);

[SURF.nfree] = fcnUNSTEADYWRAPPER(SURF.nfree,[COND.vecVEHVINF,0,0],SURF.vecDVEHVCRD,INPU.vecAREA,SURF.GammaInt,...
    COND.valMAXTIME,COND.valGUSTSTART,COND.valDELTIME,INPU.vecSYM,COND.vecVEHBETA,SURF.en_t);

nforce_fwd = (SURF.nfree(rows(5,:),:) + SURF.nind(rows(5,:),:));
CL_fwd = sum(nforce_fwd,1)./(qinf*vecAREA);

% OUTP.vecCL = OUTP.vecCL.*((AR+2)/AR);

% figure(420)
% hold on
% plot(1:COND.valMAXTIME,OUTP.norm_percent,'-.k','linewidth',1.5)
% grid on
% xlabel('Time step')
% ylabel('Normal Flow Percentage')

figure(6969)
hold on
% Kussner_Function
% plot(s,Clk,'-k','linewidth',1.5)
s2 = 2*((COND.valGUSTSTART:COND.valMAXTIME).*COND.valDELTIME - COND.valGUSTSTART*COND.valDELTIME);
plot(s2./2,CL(COND.valGUSTSTART:end),'-.b','linewidth',1.5)
% plot(s2,CL_fwd(COND.valGUSTSTART:end),'--r','linewidth',1.5)
% plot(s2,CL(COND.valGUSTSTART:end),':m','linewidth',1.5)
% plot((1:COND.valMAXTIME).*COND.valDELTIME,OUTP.vecCL,'-b','linewidth',1.5)
box on
grid on
grid minor
ylabel('C_l')
% xlabel('Distance Traversed Through Gust in Semi-Chords (b)')
% ylim([0 0.4])


% save('Sharp_Edge_m30.mat')