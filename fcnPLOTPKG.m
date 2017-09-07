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
    for i = 1:valNELE
        endpoint_left = (sum(matVLST([matDVE(i,4); matDVE(i,1)],:),1)./2) - matCENTER(i,:);
        endpoint_right = (sum(matVLST([matDVE(i,2); matDVE(i,3)],:),1)./2) - matCENTER(i,:);
        
        tt = fcnGLOBSTAR([endpoint_left; endpoint_right],repmat(vecDVEROLL(i),2,1), repmat(vecDVEPITCH(i),2,1), repmat(vecDVEYAW(i),2,1));
        
        etas = linspace(tt(1,2),tt(2,2))';
        circ = matCOEFF(i,3).*etas.^2 + matCOEFF(i,2).*etas + matCOEFF(i,1);
        
        pt_loc = [linspace(tt(1,1),tt(2,1))' etas circ];
        
        len = size(circ,1);
        circ_glob = fcnSTARGLOB(pt_loc, repmat(vecDVEROLL(i),len,1), repmat(vecDVEPITCH(i),len,1), repmat(vecDVEYAW(i),len,1));
        circ_glob = circ_glob + matCENTER(i,:);
        hold on
        plot3(circ_glob(:,1), circ_glob(:,2), circ_glob(:,3),'-m','LineWidth',3)
        hold off
        
    end
end


end

