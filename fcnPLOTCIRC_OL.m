function [] = fcnPLOTCIRC_OL(valNELE, matDVE, matVLST, matCENTER, vecDVEROLL, vecDVEPITCH, vecDVEYAW, matCOEFF)

for i = 1:valNELE
    corners = fcnGLOBSTAR(matVLST(matDVE(i,:),:) - matCENTER(i,:), repmat(vecDVEROLL(i),4,1), repmat(vecDVEPITCH(i),4,1), repmat(vecDVEYAW(i),4,1));
    points = polygrid(corners(:,1), corners(:,2), 100);
    
    % points(:,2) is eta in local, points(:,1) is xsi
    circ = matCOEFF(i,3).*points(:,2).^2 + matCOEFF(i,2).*points(:,2) + matCOEFF(i,5).*points(:,1).^2 + matCOEFF(i,4).*points(:,1) + matCOEFF(i,1);
        
    len = size(circ,1);
    circ_glob = fcnSTARGLOB([points circ], repmat(vecDVEROLL(i),len,1), repmat(vecDVEPITCH(i),len,1), repmat(vecDVEYAW(i),len,1));
    circ_glob = circ_glob + matCENTER(i,:);
    hold on
    scatter3(circ_glob(:,1), circ_glob(:,2), circ_glob(:,3),'xk')
    hold off
    
end


function [inPoints] = polygrid( xv, yv, ppa)

	N = sqrt(ppa);
%Find the bounding rectangle
	lower_x = min(xv);
	higher_x = max(xv);

	lower_y = min(yv);
	higher_y = max(yv);
%Create a grid of points within the bounding rectangle
	inc_x = 1/N;
	inc_y = 1/N;
	
	interval_x = lower_x:inc_x:higher_x;
	interval_y = lower_y:inc_y:higher_y;
	[bigGridX, bigGridY] = meshgrid(interval_x, interval_y);
	
%Filter grid to get only points in polygon
	in = inpolygon(bigGridX(:), bigGridY(:), xv, yv);
%Return the co-ordinates of the points that are in the polygon
	inPoints = [bigGridX(in), bigGridY(in)];

end

end

