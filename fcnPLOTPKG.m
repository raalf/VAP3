function [] = fcnPLOTPKG(flagVERBOSE, flagPLOTWAKEVEL, flagCIRCPLOT, flagPLOTUINF, valNELE, matDVE, matVLST, ....
    matCENTER, matFUSEGEOM, ...
    valWNELE, matWDVE, matWVLST, matWCENTER, matWDVEMP, matWDVEMPIND, matUINF, ...
    vecDVEROLL, vecDVEPITCH, vecDVEYAW, matCOEFF)

[hFig2] = fcnPLOTBODY(flagVERBOSE, valNELE, matDVE, matVLST, matCENTER, []);
[hFig2] = fcnPLOTWAKE(flagVERBOSE, hFig2, valWNELE, matWDVE, matWVLST, matWCENTER);
[hLogo] = fcnPLOTLOGO(0.97,0.03,14,'k','none');

if flagPLOTWAKEVEL == 1
    try
        quiver3(matWDVEMP(:,1),matWDVEMP(:,2),matWDVEMP(:,3),matWDVEMPIND(:,1),matWDVEMPIND(:,2),matWDVEMPIND(:,3));
    end
end
if flagPLOTUINF == 1
    try
        %             quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3),matVEHUINF(:,1),matVEHUINF(:,2),matVEHUINF(:,3),'g');
        %             quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3),matROTORUINF(:,1),matROTORUINF(:,2),matROTORUINF(:,3),'c');
        quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3),matUINF(:,1),matUINF(:,2),matUINF(:,3),'g');
    end
end

if flagCIRCPLOT == 1
    fcnPLOTCIRC(valNELE, matDVE, matVLST, matCENTER, vecDVEROLL, vecDVEPITCH, vecDVEYAW, matCOEFF, 1e2)
end

end

