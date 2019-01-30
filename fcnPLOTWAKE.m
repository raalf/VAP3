function [hFig2] = fcnPLOTWAKE(verbose, hFig2, valWNELE, matWDVE, matWVLST, matWCENTER, vecWDVESURFACE, wdve_sym)

col = single(vecWDVESURFACE)./max(single(vecWDVESURFACE));
% col = vecWDVESURFACE;
hold on
colormap parula

patch('Faces',matWDVE,'Vertices',matWVLST,'FaceVertexCData',col,'FaceColor','flat','EdgeAlpha',0.6,'FaceAlpha',0.4);
patch('Faces',matWDVE(wdve_sym,:),'Vertices',[matWVLST(:,1) matWVLST(:,2).*-1 matWVLST(:,3)],'FaceVertexCData',col(wdve_sym),'FaceColor','flat','EdgeAlpha',0.6,'FaceAlpha',0.4);

if verbose == 1
    for ii = 1:valWNELE
        str = sprintf('%d',ii);
        text(matWCENTER(ii,1),matWCENTER(ii,2),matWCENTER(ii,3),str,'Color','k','FontSize',15)
    end
end

hold off   
end