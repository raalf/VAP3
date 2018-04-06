function [hFig2] = fcnPLOTBODY(verbose, valNELE, matDVE, matVLST, matCENTER, matFDVE, matFVLST, sym)

hFig2 = figure(3);
clf(3);

patch('Faces',matDVE,'Vertices',matVLST,'FaceColor',[255 90 90]./255)
if sym == true
    patch('Faces',matDVE,'Vertices',[matVLST(:,1) matVLST(:,2).*-1 matVLST(:,3)],'FaceColor',[255 90 90]./255)
end
hold on


% alpha(0.5)

if verbose == 1
    for ii = 1:valNELE
        str = sprintf('%d',ii);
        text(matCENTER(ii,1),matCENTER(ii,2),matCENTER(ii,3),str,'Color','k','FontSize',15);
    end
    
    %     for ii = 1:length(matVLST(:,1))
    %         str = sprintf('%d',ii);
    %         text(matVLST(ii,1),matVLST(ii,2),matVLST(ii,3),str,'Color','g','FontSize',15);
    %     end
    
end

if ~isempty(matFDVE)
    sz = size(matFDVE);
    
    if length(sz) < 4
        patch('Faces',matFDVE,'Vertices',matFVLST,'FaceColor',[255 90 90]./255)
    else
        for i = 1:sz(4)
            patch('Faces',matFDVE(:,:,i),'Vertices',matFVLST(:,:,i),'FaceColor',[255 90 90]./255)
        end
    end
    
end

hold off

box on
grid on
axis equal
axis tight

xlabel('X-Dir','FontSize',15);
ylabel('Y-Dir','FontSize',15);
zlabel('Z-Dir','FontSize',15);