function [valCL] = fcnFORCES(matCOEFF,vecK,matDVE,valNELE,matCENTER,matVLST,vecUINF,vecDVELESWP,vecDVEMCSWP,vecDVEHVSPN,vecDVEROLL,vecDVEPITCH,vecDVEYAW,vecDVELE,matADJE,...
    valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP,vecSYM,vecDVETESWP,valAREA,valBETA)
%% Forces package
%place any force functions in here and add a description.

% INPUT:

% OUTPUT:
% valCL - Lift force coefficient


%% Element normal forces, lift forces and side forces (freestream and induced)
 [nfree,nind,liftfree,liftind,sidefree,sideind] = fcnDVENFORCE(matCOEFF,vecK,matDVE,valNELE,matCENTER,matVLST,vecUINF,vecDVELESWP,vecDVEMCSWP,vecDVEHVSPN,vecDVEROLL,vecDVEPITCH,vecDVEYAW,vecDVELE,matADJE,...
    valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP, vecSYM, vecDVETESWP);

%% Sum up element forces to generate CL and CY
[valCL, valCLF, valCLI, valCY, valCYF, valCYI]= fcnWINGNFORCE(nfree,nind,liftfree,liftind,sidefree,sideind,vecUINF,valAREA,vecSYM,valBETA);