function [valCL,valCDI] = fcnFORCES(matCOEFF,vecK,matDVE,valNELE,matCENTER,matVLST,vecUINF,vecDVELESWP,vecDVEMCSWP,vecDVEHVSPN,vecDVEROLL,vecDVEPITCH,vecDVEYAW,vecDVELE,vecDVETE,matADJE,...
    valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP,vecSYM,vecDVETESWP,valAREA,valBETA,vecDVEWING)
%% Forces package
%place any force functions in here and add a description.

% INPUT:

% OUTPUT:
% valCL - Lift force coefficient


%% Element normal forces, lift forces and side forces (freestream and induced)
 [nfree,nind,liftfree,liftind,sidefree,sideind] = fcnDVENFORCE(matCOEFF,vecK,matDVE,valNELE,matCENTER,matVLST,vecUINF,vecDVELESWP,vecDVEMCSWP,vecDVEHVSPN,vecDVEROLL,vecDVEPITCH,vecDVEYAW,vecDVELE,matADJE,...
    valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP, vecSYM, vecDVETESWP);

%% TE Element induced drag forces
%% Induced Drag force
[inddrag] =fcnDVEINDDRAG(matCOEFF,matDVE,matVLST,vecUINF,vecDVEHVSPN,vecDVETE,valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, ...
    valWSIZE, valTIMESTEP,vecSYM,vecDVEWING );
%% Sum up element forces to generate total wing forces
[valCL, valCLF, valCLI, valCY, valCYF, valCYI,valCDI]= fcnWINGNFORCE(liftfree,liftind,sidefree,sideind,inddrag,vecUINF,valAREA,vecSYM,valBETA);


