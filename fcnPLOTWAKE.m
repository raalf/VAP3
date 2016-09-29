function [hFig2] = fcnPLOTWAKE(verbose, hFig2, valWNELE, matWDVE, matWVLST, matWCENTER)

fig = hFig2;

patch('Faces',matWDVE,'Vertices',matWVLST,'FaceColor','b')
hold on

alpha(0.5)

if verbose == 1
    for ii = 1:valWNELE
        str = sprintf('%d',ii);
        text(matWCENTER(ii,1),matWCENTER(ii,2),matWCENTER(ii,3),str,'Color','k','FontSize',15);
    end
    
%     for ii = 1:length(matWVLST(:,1))
%         str = sprintf('%d',ii);
%         text(matWVLST(ii,1),matWVLST(ii,2),matWVLST(ii,3),str,'Color','g','FontSize',15);
%     end
    
end