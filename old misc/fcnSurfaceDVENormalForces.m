function [N_force] = fcnSurfaceDVENormalForces(FW,Temp,Aircraft)

jj = 1;
for i = 1:Aircraft.General.Panels
    
    for j = 1:length(FW.Panels(i).DVE.Index)
        
        if j <= FW.Panels(i).n % if we are in the first row
            % vector of bound vortex along leading edge in this case
            tempA(1) = tand(FW.Panels(i).DVE.phiLE(j));
            tempA(2) = 1;
            tempA(3) = 0;
            % Transforming to local ref frame
            S = star_glob(tempA,FW.Panels(i).DVE.roll(j),FW.Panels(i).DVE.pitch(j),FW.Panels(i).DVE.yaw(j));
            chordwise_column = j;
        else
            chordwise_column = mod(j,FW.Panels(i).n); % We use this to determine which column of elements we are in (for indexing later with U1, Uo, U2)
            if chordwise_column == 0 
                chordwise_column = chordwise_column + FW.Panels(i).n;
            end
            
            tempA(1) = tand(FW.Panels(i).DVE.phiMID(j));
            tempA(2) = 1;
            tempA(3) = 0;
            
            S = star_glob(tempA,FW.Panels(i).DVE.roll(j),FW.Panels(i).DVE.pitch(j),FW.Panels(i).DVE.yaw(j));
        end
        
        %fprintf('S: %f %f %f\n', S(1), S(2), S(3));
        
        eta = FW.Panels(i).DVE.eta(j);
        eta8 = eta*0.8;
        
        tempA = cross(Temp.u,S);
        UxS = norm(tempA);
        %fprintf('UxS: %f\n', UxS);
        
        eN = tempA.*(1/UxS);
        
        %fprintf('eN: %f %f %f \n', eN(1), eN(2), eN(3));
        
        %the lift direction  eL=Ux[0,1,0]/|Ux[0,1,0]|
        tempS = sqrt(Temp.u(1)^2 + Temp.u(2)^2 + Temp.u(3)^2);
        eL(1) = -Temp.u(3)/tempS;
        eL(2) = 0;
        eL(3) = Temp.u(1)/tempS;
        
        %fprintf('eL: %f %f %f \n', eL(1), eL(2), eL(3));
        
        %the side force direction eS=UxeL/|UxeL|
        tempA = cross(eL,Temp.u);
        tempS = 1/norm(tempA);
        eS = tempA.*tempS;
        
        %         fprintf('eS: %f %f %f \n', eS(1), eS(2), eS(3));
        
        if j <= FW.Panels(i).n % if we are in the first row
            A = FW.Panels(i).DVE.A(j);
            B = FW.Panels(i).DVE.B(j);
            C = FW.Panels(i).DVE.C(j);
        else
            A = FW.Panels(i).DVE.A(j) - FW.Panels(i).DVE.A(j-FW.Panels(i).n);
            B = FW.Panels(i).DVE.B(j) - FW.Panels(i).DVE.B(j-FW.Panels(i).n);
            C = FW.Panels(i).DVE.C(j) - FW.Panels(i).DVE.C(j-FW.Panels(i).n);
        end
        
        % fprintf('A: %f B: %f C: %f \n', A, B, C);
        
        
        % normal force due to free stream
        N_free = (A*2*eta + C/3*2*eta*eta*eta)*UxS;
        
        if j <= FW.Panels(i).n % if we are in the first row
            
            %computing the leading edge midpoint
            tempA(1) = tempA(1) - FW.Panels(i).DVE.xsi(j);
            tempA(2) = 0;
            tempA(3) = 0;
            %transforming into local reference frame
            tempAA = star_glob(tempA,FW.Panels(i).DVE.roll(j),FW.Panels(i).DVE.pitch(j),FW.Panels(i).DVE.yaw(j));
            
            xoLE = FW.Panels(i).DVE.xo(j,:) + tempAA;
            
            %computing the ind. velocity at left (1) edge of bound vortex
            X = S.*-eta8;
            tempA = xoLE + X;
            w1 = fcnDVEInducedVelocity(FW, Aircraft, Temp, tempA);
            
            %computing the ind. velocity at center (0) of bound vortex
            wo = fcnDVEInducedVelocity(FW, Aircraft, Temp, xoLE);
            
            %//computing the ind. velocity at right (2) edge of bound vortex
            %//vector from mid point of elementary wing bound vortex
            %//to edge (2) of bound vortex; x1=xo-x and x2=xo+x
            %//due to singular behavior of velocity at element edge
            %//velocity is computed 0.1eta away from edge (hence factor 0.8)
            
            X = S.*eta8;
            tempA = xoLE+X;
            
            w2 = fcnDVEInducedVelocity(FW, Aircraft, Temp, tempA);
            
            if FW.m > 1 % if multiple lifting lines, compute and store velocities at half chord location for averaging later
                xoLE = FW.Panels(i).DVE.xo(j,:);
                
                tempA(1) = tand(FW.Panels(i).DVE.phiMID(j));
                tempA(2) = 1;
                tempA(3) = 0;
                
                tempAA = star_glob(tempA,FW.Panels(i).DVE.roll(j),FW.Panels(i).DVE.pitch(j),FW.Panels(i).DVE.yaw(j));
                
                %computing the ind. velocity at left (1) edge of bound vortex
                tempA = tempAA*-eta8; % WHAT IS GOING ON HEREEEEEEEE
                tempA = xoLE + X;
                U1(chordwise_column,:) = fcnDVEInducedVelocity(FW, Aircraft, Temp, tempA); % IS THIS COUNTER J RIGHT? HE HAS "K" WHICH IS SPANWISE ELEMENT #
                
                %computing the ind. velocity at center (0) of bound vortex
                Uo(chordwise_column,:) = fcnDVEInducedVelocity(FW, Aircraft, Temp, xoLE);
                
                %//computing the ind. velocity at right (2) edge of bound vortex
                %//vector from mid point of elementary wing bound vortex
                %//to edge (2) of bound vortex; x1=xo-x and x2=xo+x
                %//due to singular behavior of velocity at element edge
                %//velocity is computed 0.1eta away from edge (hence factor 0.8)
                X = tempAA*eta8;
                tempA = xoLE+X;
                U2(chordwise_column,:) = fcnDVEInducedVelocity(FW, Aircraft, Temp, tempA);
                
            end
            
        else
            %//case of multiple lifting lines along the span
            %//the induced velocity at the lifting line is averaged with the
            %//velocities at mid chord locations of the DVES upstream and
            %//downstream of the bound vortex. Otherwise, the singularity of
            %//the bound vortex and the discontinuity of the bound vortex sheet
            %//of a wing with twist causes trouble.  G.B. 1/24/06
            
            xoLE = FW.Panels(i).DVE.xo(j,:);
            
            %computing the ind. velocity at left (1) edge of bound vortex
            X = S*-eta8;
            tempA = xoLE+X;
            tempAA = fcnDVEInducedVelocity(FW, Aircraft, Temp, tempA);
            
            %averaging velocity and reassigning velocity to temporary variable
            w1 = (U1(chordwise_column,:)+tempAA)./2;
            
            %computing the ind. velocity at center (0) of bound vortex
            U1(chordwise_column,:) = tempAA;
            
            tempAA = fcnDVEInducedVelocity(FW, Aircraft, Temp, xoLE);
            
            %averaging velocity and reassigning velocity to temporary variable
            
            wo = (Uo(chordwise_column,:)+tempAA)./2;
            Uo(chordwise_column,:) = tempAA;
            
            %//computing the ind. velocity at right (2) edge of bound vortex
            %//vector from mid point of elementary wing bound vortex
            %//to edge (2) of bound vortex; x1=xo-x and x2=xo+x
            %//due to singular behavior of velocity at element edge
            %//velocity is computed 0.1eta away from edge (hence factor 0.8)
            
            X = S*eta8;
            tempA = xoLE+X;
            
            tempAA = fcnDVEInducedVelocity(FW, Aircraft, Temp, tempA);
            
            %averaging velocity and reassigning velocity to temporary variable
            w2 = (U2(chordwise_column,:)+tempAA)./2;
            U2(chordwise_column,:) = tempAA;
            
        end
        
        %//Integration of induced forces with Simpson's Rule
        %//Integration requires overhanging edges!!
        %//See also KHH linees 2953 - 2967, A23SIM
        
        %//Kutta-Joukowski at left (1), center, right (2) edge
        tempA = cross(w1,S);
        gamma1 = A-B*eta8+C*eta8*eta8;
        R1 = tempA.*gamma1;
        
        tempA = cross(wo,S);
        gammao = A;
        Ro = tempA.*gammao;
        
        tempA = cross(w2,S);
        gamma2 = A+B*eta8+C*eta8*eta8;
        R2 = tempA.*gamma2;
        
        %//The resultierende induced force of element l is
        %//determined by numerically integrating forces acros element
        %//using Simpson's Rule with overhaning parts
        R = (R1 + 4*Ro+R2).*eta8/3;
        
        R = R + ((7.*R1 - 8.*Ro + 7.*R2).*(eta-eta8)./3);
        
        % normal force
        N_force(jj,5) = N_free;
        N_force(jj,6) = dot(R,eN);
        
        %lift force is the normal force in the x-z plane
        N_force(jj,1) = N_free*sqrt(eN(1)*eN(1) + eN(3)*eN(3));
        
        if eN(3) < 0
            N_force(jj,1) = N_force(jj,1)*-1;
        end
        
        N_force(jj,2) = dot(R,eL);
        
        % side force is force in y-dir
        
        N_force(jj,3) = N_free*eN(2);
        N_force(jj,4) = dot(R,eS);
        
        %fprintf('N_force: %f %f %f %f %f %f \n', N_force(jj,1), N_force(jj,2), N_force(jj,3), N_force(jj,4), N_force(jj,5), N_force(jj,6));
        jj = jj + 1;
        
    end
end

end