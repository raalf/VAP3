function [] = fcnPLOTPKG(FLAG, SURF, VISC, WAKE, COND)

if FLAG.PREVIEW == 1; FLAG.RELAX = 0; end

[hFig2] = fcnPLOTBODY(FLAG.VERBOSE, SURF.valNELE, SURF.matDVE, SURF.matVLST, SURF.matCENTER, []);
[hFig2] = fcnPLOTWAKE(FLAG.VERBOSE, hFig2, WAKE.valWNELE, WAKE.matWDVE, WAKE.matWVLST, WAKE.matWCENTER, WAKE.vecWDVESURFACE);
% [hLogo] = fcnPLOTLOGO(0.97,0.03,14,'k','none');

if FLAG.PLOTWAKEVEL == 1 && COND.valMAXTIME > 0 && FLAG.RELAX == 1
    try
        quiver3(WAKE.matWDVEMP(:,1),WAKE.matWDVEMP(:,2),WAKE.matWDVEMP(:,3),WAKE.matWDVEMPIND(:,1),WAKE.matWDVEMPIND(:,2),WAKE.matWDVEMPIND(:,3));
    end
end

if FLAG.PLOTUINF == 1
    try
        %             quiver3(SURF.matCENTER(:,1),SURF.matCENTER(:,2),SURF.matCENTER(:,3),matVEHUINF(:,1),matVEHUINF(:,2),matVEHUINF(:,3),'g');
        %             quiver3(SURF.matCENTER(:,1),SURF.matCENTER(:,2),SURF.matCENTER(:,3),matROTORUINF(:,1),matROTORUINF(:,2),matROTORUINF(:,3),'c');
        quiver3(SURF.matCENTER(:,1),SURF.matCENTER(:,2),SURF.matCENTER(:,3),SURF.matUINF(:,1),SURF.matUINF(:,2),SURF.matUINF(:,3),'g');
    end
end

if FLAG.CIRCPLOT == 1
    fcnPLOTCIRC(SURF.valNELE, SURF.matDVE, SURF.matVLST, SURF.matCENTER, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW, SURF.matCOEFF, 1e2)
end

end

