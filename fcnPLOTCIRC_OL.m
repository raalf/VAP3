function [] = fcnPLOTCIRC_OL(valNELE, matDVE, matVLST, matCENTER, vecDVEROLL, vecDVEPITCH, vecDVEYAW, matCOEFF)

for i = 1:valNELE
    for j = 1:3
        % Spanwise
        if j == 1
            endpoint_left = (sum(matVLST([matDVE(i,4); matDVE(i,1)],:),1)./2) - matCENTER(i,:);
            endpoint_right = (sum(matVLST([matDVE(i,2); matDVE(i,3)],:),1)./2) - matCENTER(i,:);
        elseif j == 2
            endpoint_left = matVLST(matDVE(i,1),:) - matCENTER(i,:);
            endpoint_right = matVLST(matDVE(i,2),:) - matCENTER(i,:);
        elseif j == 3
            endpoint_left = matVLST(matDVE(i,4),:) - matCENTER(i,:);
            endpoint_right = matVLST(matDVE(i,3),:) - matCENTER(i,:);
        end
        
        tt = fcnGLOBSTAR([endpoint_left; endpoint_right],repmat(vecDVEROLL(i),2,1), repmat(vecDVEPITCH(i),2,1), repmat(vecDVEYAW(i),2,1));
        
        etas = linspace(tt(1,2),tt(2,2))';
        circ = matCOEFF(i,3).*etas.^2 + matCOEFF(i,2).*etas + matCOEFF(i,1);
        
        pt_loc = [linspace(tt(1,1),tt(2,1))' etas circ];
        
        len = size(circ,1);
        circ_glob = fcnSTARGLOB(pt_loc, repmat(vecDVEROLL(i),len,1), repmat(vecDVEPITCH(i),len,1), repmat(vecDVEYAW(i),len,1));
        circ_glob = circ_glob + matCENTER(i,:);
        hold on
        plot3(circ_glob(:,1), circ_glob(:,2), circ_glob(:,3),'-k','LineWidth',2)
        hold off
        
        %Chordwise
        
        if j == 1
            endpoint_front = (sum(matVLST([matDVE(i,1); matDVE(i,2)],:),1)./2) - matCENTER(i,:);
            endpoint_back = (sum(matVLST([matDVE(i,3); matDVE(i,4)],:),1)./2) - matCENTER(i,:);
        elseif j == 2
            endpoint_front = matVLST(matDVE(i,1),:) - matCENTER(i,:);
            endpoint_back = matVLST(matDVE(i,4),:) - matCENTER(i,:);
        elseif j == 3
            endpoint_front = matVLST(matDVE(i,2),:) - matCENTER(i,:);
            endpoint_back = matVLST(matDVE(i,3),:) - matCENTER(i,:);
        end
        
        tt = fcnGLOBSTAR([endpoint_front; endpoint_back],repmat(vecDVEROLL(i),2,1), repmat(vecDVEPITCH(i),2,1), repmat(vecDVEYAW(i),2,1));
        
        xsis = linspace(tt(1,1),tt(2,1))';
        circ = matCOEFF(i,5).*xsis.^2 + matCOEFF(i,4).*xsis + matCOEFF(i,1);
        
        pt_loc = [xsis linspace(tt(1,2),tt(2,2))' circ];
        
        len = size(circ,1);
        circ_glob = fcnSTARGLOB(pt_loc, repmat(vecDVEROLL(i),len,1), repmat(vecDVEPITCH(i),len,1), repmat(vecDVEYAW(i),len,1));
        circ_glob = circ_glob + matCENTER(i,:);
        hold on
        plot3(circ_glob(:,1), circ_glob(:,2), circ_glob(:,3),'-k','LineWidth',2)
        hold off
    end
end

end

