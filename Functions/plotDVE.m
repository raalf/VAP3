function [  ] = plotDVE(Output, Temp, FW)
%PLOTDVE Summary of this function goes here
%   Detailed explanation goes here
wake = 1;
names = 0;

hFig1 = figure(1);
clf(1);

for panel = 1:length(Output.(Temp.AI)(end).TimestepData)
    
    %     DVE = FW.Panels(panel).DVE;
    DVE = Output.(Temp.AI)(end).TimestepData(panel).DVE;
    xo = DVE.xo;
    xsi = DVE.xsi;
    eta = DVE.eta;
    phi_LE = DVE.phiLE;
    phi_TE = DVE.phiTE;
    nu = DVE.roll;
    epsilon = DVE.pitch;
    psi = DVE.yaw;
    
    
    hold on
    % DVE Normal Vector
    %         quiver3(xo(:,1),xo(:,2),xo(:,3),DVE.norm(:,1),DVE.norm(:,2),DVE.norm(:,3),0,'b');
    
    for n = 1:length(DVE.Index);
        %computing left-leading edge point in local ref. frame
        tempA(1) = -xsi(n) - eta(n)*tand(phi_LE(n));
        tempA(2) = -eta(n);
        tempA(3) = 0;
        
        tempAA = star_glob(tempA,nu(n),epsilon(n),psi(n));
        x1 = tempAA+reshape(xo(n,:),1,3,1);
        
        % 		computing left-trailing edge point in local ref. frame
        tempA(1) = xsi(n) - eta(n)*tand(phi_TE(n));
        tempA(2) = -eta(n);
        tempA(3) = 0;
        
        tempAA = star_glob(tempA,nu(n),epsilon(n),psi(n));
        x2 = tempAA+reshape(xo(n,:),1,3,1);
        
        %computing right-trailing edge point in local ref. frame
        tempA(1) = xsi(n) + eta(n)*tand(phi_TE(n));
        tempA(2) = eta(n);
        tempA(3) = 0;
        
        tempAA = star_glob(tempA,nu(n),epsilon(n),psi(n));
        x3 = tempAA+reshape(xo(n,:),1,3,1);
        
        %computing right-leading edge point in local ref. frame
        tempA(1) = -xsi(n) + eta(n)*tand(phi_LE(n));
        tempA(2) = eta(n);
        tempA(3) = 0;
        
        tempAA = star_glob(tempA,nu(n),epsilon(n),psi(n));
        x4 = tempAA+reshape(xo(n,:),1,3,1);
        
        fillX = [x1(1) x2(1) x3(1) x4(1) x1(1)];
        fillY = [x1(2) x2(2) x3(2) x4(2) x1(2)];
        fillZ = [x1(3) x2(3) x3(3) x4(3) x1(3)];
        
%         fill3(fillX,fillY,fillZ,'r','EdgeColor','k','LineWidth',2)
        fill3(fillX,fillY,fillZ,'r','FaceAlpha',0.5,'EdgeColor','k','LineWidth',1)
        
        if names == 1
            for ii = 1:length(DVE.Index)
                str = sprintf('%d',DVE.Index(ii));
                text(DVE.xo(ii,1),DVE.xo(ii,2),DVE.xo(ii,3),str);
            end
        end
    end
    
    clear x1 x2 x3 x4 tempA tempAA
    axis equal
    grid on
    
    if wake == 1
        for jj = 2:length(Output.(Temp.AI)(end).TimestepData(panel).WakeDVE)
            
            %             DVE = Output.Panels(panel).WakeDVE(jj);
            DVE = Output.(Temp.AI)(end).TimestepData(panel).WakeDVE(jj);
            xo = DVE.xo;
            xsi = DVE.xsi;
            eta = DVE.eta;
            phi_LE = DVE.phiLE;
            phi_TE = DVE.phiTE;
            nu = DVE.roll;
            epsilon = DVE.pitch;
            psi = DVE.yaw;
            xleft = DVE.xleft;
            xright = DVE.xright;
            
%             scatter3(xleft(:,1), xleft(:,2), xleft(:,3),10,'r','filled');
%             scatter3(xright(:,1), xright(:,2), xright(:,3), 20,'g');
            
            for n = 1:length(DVE.Index);
                %computing left-leading edge point in local ref. frame
                tempA(1) = -xsi(n) - eta(n)*tand(phi_LE(n));
                tempA(2) = -eta(n);
                tempA(3) = 0;
                
                tempAA = star_glob(tempA,nu(n),epsilon(n),psi(n));
                x1 = tempAA+reshape(xo(n,:),1,3,1);
                
                % 		computing left-trailing edge point in local ref. frame
                tempA(1) = xsi(n) - eta(n)*tand(phi_TE(n));
                tempA(2) = -eta(n);
                tempA(3) = 0;
                
                tempAA = star_glob(tempA,nu(n),epsilon(n),psi(n));
                x2 = tempAA+reshape(xo(n,:),1,3,1);
                
                %computing right-trailing edge point in local ref. frame
                tempA(1) = xsi(n) + eta(n)*tand(phi_TE(n));
                tempA(2) = eta(n);
                tempA(3) = 0;
                
                tempAA = star_glob(tempA,nu(n),epsilon(n),psi(n));
                x3 = tempAA+reshape(xo(n,:),1,3,1);
                
                %computing right-leading edge point in local ref. frame
                tempA(1) = -xsi(n) + eta(n)*tand(phi_LE(n));
                tempA(2) = eta(n);
                tempA(3) = 0;
                
                tempAA = star_glob(tempA,nu(n),epsilon(n),psi(n));
                x4 = tempAA+reshape(xo(n,:),1,3,1);
                
                fillX = [x1(1) x2(1) x3(1) x4(1) x1(1)];
                fillY = [x1(2) x2(2) x3(2) x4(2) x1(2)];
                fillZ = [x1(3) x2(3) x3(3) x4(3) x1(3)];
                
                fill3(fillX,fillY,fillZ,'b','FaceAlpha',0.25,'EdgeColor','b')
%                 fill3(fillX,fillY,fillZ,'b','EdgeColor','b')
                if names == 1
                    for ii = 1:length(DVE.Index)
                        str = sprintf('%d',DVE.Index(ii));
                        text(DVE.xo(ii,1),DVE.xo(ii,2),DVE.xo(ii,3),str);
                    end
                end
            end
        end
    end
    
end

% Temp.lasttimestep = 7;
% view([-37 19])
% zlim([0.35 .7]);
% xlim([-0.55 -0.05]);
% % xlim([-Temp.lasttimestep*FW.Deltime+0.2 -Temp.lasttimestep*FW.Deltime+2]);
% ylim([7.2 7.6]);

set(hFig1, 'Renderer', 'painters');
% set(hFig1, 'Position',[500 500 600 500])

xlabel('x (m)','FontSize',15)
ylabel('y (m)','FontSize',15)
zlabel('z (m)','FontSize',15)
% axis tight
hold off


end

