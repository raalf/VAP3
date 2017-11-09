function [hFig2] = fcnPLOTWAKE(verbose, hFig2, valWNELE, matWDVE, matWVLST, matWCENTER, vecWDVESURFACE)

fig = hFig2;

col = vecWDVESURFACE./max(vecWDVESURFACE);
% col = vecWDVESURFACE;
hold on
colormap parula
patch('Faces',matWDVE,'Vertices',matWVLST,'FaceVertexCData',col,'FaceColor','flat','EdgeAlpha',0.6,'FaceAlpha',0.4);


% for i = 1:max(vecWDVESURFACE)
% patch('Faces',matWDVE,'Vertices',matWVLST,'FaceColor','b','EdgeColor','b','FaceAlpha',0.5);
% end

if verbose == 1
    for ii = 1:valWNELE
        str = sprintf('%d',ii);
        text(matWCENTER(ii,1),matWCENTER(ii,2),matWCENTER(ii,3),str,'Color','k','FontSize',15)
    end
    
%     for ii = 1:length(matWVLST(:,1))
%         str = sprintf('%d',ii);
%         text(matWVLST(ii,1),matWVLST(ii,2),matWVLST(ii,3),str,'Color','g','FontSize',15);
%     end
    
end