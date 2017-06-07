


for i = 1:m
    for j = 1:n
        %computing left-leading edge point in local ref. frame
        tempA(1) = -xsi(i,j) - eta(i,j)*tand(phiLE(i,j));
        tempA(2) = -eta;
        tempA(3) = 0;
        
        tempAA = star_glob3D(tempA,nu(i,j),epsilon(i,j),psi(i,j));
        x1 = tempAA+CP(i,j,:);
        
        % 		computing left-trailing edge point in local ref. frame
        tempA(1) = xsi(i,j) - eta(i,j)*tand(phiTE(i,j));
        tempA(2) = -eta(i,j);
        tempA(3) = 0;
        
        tempAA = star_glob3D(tempA,nu(i,j),epsilon(i,j),psi(i,j));
        x2 = tempAA+CP(i,j,:);
        
        %computing right-trailing edge point in local ref. frame
        tempA(1) = xsi(i,j) + eta(i,j)*tand(phiTE(i,j));
        tempA(2) = eta(i,j);
        tempA(3) = 0;
        
        tempAA = star_glob3d(tempA,nu(i,j),epsilon(i,j),psi(i,j));
        x3 = tempAA+CP(i,j,:);
        
        %computing right-leading edge point in local ref. frame
        tempA(1) = -xsi(i,j) + eta(i,j)*tand(phiLE(i,j));
        tempA(2) = eta(i,j);
        tempA(3) = 0;
        
        tempAA = star_glob3d(tempA,nu(i,j),epsilon(i,j),psi(i,j));
        x4 = tempAA+CP(i,j,:);
               
        
        
    end
end
