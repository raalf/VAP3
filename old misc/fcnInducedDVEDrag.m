function [ Output, D_force ] = fcnInducedDVEDrag(FW, Aircraft, Temp,Output)

%copy from freewake 2015beta.
% W.J.M.B, COURTYARD STUDY LOUNGE, STUDENT CENTER, UNIVERSITY OF CALIFORNIA
% - IRVINE, IRVINE, CALIFORNIA, USA, 92697. 2016-01-23
%{
//This function is the DVE expansion of the function Induced_Eppler_Drag
//that computes the induced drag at the trailing edge, where
//the spanwise bound vorticity has been collapsed into a single vortex.
//The method is discussed more thoroughly in:
//Eppler and Schmid-Goeller, "A Method to Calculate the Influence of
//Vortex Roll-Up on the Induced Drag of Wings," Finite Approximations
//in Fluid Mechanics II, Hirsche, E. H. (ed.), Notes on Numerical Fluid
//Mechanics, Volume 25, Friedr. Vieweg und Sohn, Braunschweig, 1990
//Kutta-Joukowsky is being applied to the trailing edge at three points
//along each trailing edge element.  Similarly to the lift computation,
//Simpson's Rule is used to compute the total drag of each element.

//Function computes forces/density for each spanwise section of the
//wing by applying Kutta-Joukowski's theorem to either of its edges and
//its center in order to get the local forces. The total spanwise section
//forces are determined by integrating with Simpson's rule (with overhang)
//Note: the velocities are not computed directly at the edges, but
//(1-delta)/2 of the elementary span further towards the center.
//Otherwise, the computed induced velocity would become singular at that
//point due to the next neighboring elementary wing influence.
//
//input:
//	info		- general information
//	panelPtr	- information on panels
//	surfacePtr	- information on surface DVEs
//	wakePtr		- information on wake DVEs
//	rightnow	- current time step
//
//output:
//	CDi			- total drag coefficient
//  D_force 	- local drag force/density along span
%}
% //#############################################################################
% //							FORCE LOOP - START
% //#############################################################################
% //
% //in this loop the induced forces at the trailing edge are computed for three
% //points along the trailing egd of the wing in the region of the surface DVE
% //with the index 'index'.  The forces are integrated over the surface DVE's
% //span in order to get the induced drag contribution of that span location.
CDi = 0;
TEindex = 1; % //initializing wake DVE index, this is a counter for all TE DVEs


%loop over panels
for i = 1:Aircraft.General.Panels
    
    %loop over trailing edge elements of this panel
    for j = 1:FW.Panels(i).n
        
        %forward index to TE row of elements, then add j, j is the spanwise index
        index = (FW.Panels(i).n * FW.m) - FW.Panels(i).n  + j;
        
        %       //drag force direction
        eD = Temp.u;
        
        
        %         //#########################################################################
        A = FW.Panels(i).DVE.A(index);
        B = FW.Panels(i).DVE.B(index);
        C = FW.Panels(i).DVE.C(index);
        
        %         //#########################################################################
        % 		//Computing the three points along the unswept trailing edge,
        % 		//at which Kutta-Joukowsky is applied.
        %
        % 		//The wing-trailing edge is located at the TRAILING EDGE of most aft
        % 		//surface DVE
        xsiTE = FW.Panels(i).DVE.xsi(index);
        
        %         	//The trailing-edge sweep is phiTE of the DVE
        phiTE = FW.Panels(i).DVE.phiTE(index);
        %phi is in deg
        % 		//the left and right points are 20% of half span away from edge
        % 		//in order to stay away from the singularity along the edge of
        % 		//the DVE
        eta  =	FW.Panels(i).DVE.eta(index);
        eta8=eta*.8;	%//0.8 as done for lift computation,
        
        X(2,:) = fcnEdge_Point(FW,i,index,phiTE,-eta8,xsiTE,FW.Panels(i).DVE.xo(index,:));
        X(3,:) = fcnEdge_Point(FW,i,index,phiTE,eta8,xsiTE,FW.Panels(i).DVE.xo(index,:));
        
        X(1,:) = (X(2,:) + X(3,:)) / 2;
        
        %         //computing the normalized vector along the trailing edge
        S(1) = X(3,1) - X(2,1);
        S(2) = X(3,2) - X(2,2);
        S(3) = X(3,3) - X(2,3);
        tempS= 0.5/eta8;%//1/norm2(S) //I don't why, but it works G.B. 5/30/05
        S = S*tempS;
        
        
        %         //initializing induced velocities of current span location
        w_ind(2,1) = 0; w_ind(2,2)=0; w_ind(2,3) = 0;
        w_ind(1,1) = 0; w_ind(1,2)=0; w_ind(1,3) = 0;
        w_ind(3,1) = 0; w_ind(3,2)=0; w_ind(3,3) = 0;
        
        %  //###########################################################################
        %   //						INDUCED VELOCITY LOOP - START
        %   //###########################################################################
        %   //
        %   //In this loop, the velocity are computed that are induced by the entire flow
        %   //field at the three points along the wing-trailing edge in the region of the
        %   //surface DVE 'index'.  This is done for each point at a time.
        
        %We can't use Wake_DVE_Vel_Induction becuase we have to remove the sweep
        %from each element in the wake. Possibly in the future we could just pass a
        %sweep into Wake_DVE_Vel_Induction???????
        
        % //loop over the three points of trailing edge of the induced element
        for k = 1:3
            
            %loop over all panels (These are inducers)
            for p = 1: length(FW.Panels)
                span = 0; %spanwise counter for wake dves  THIS MIGHT HAVE TO GO BEFORE THE PANEL LOOP STARTS!!!!!!
                
                %loop over TE elements of panel p (inducers)
                for  s = (FW.Panels(p).n * FW.m) - FW.Panels(p).n  + 1 : (FW.Panels(p).n * FW.m) - FW.Panels(p).n  + FW.Panels(p).n
                    span = span+1;
                    if FW.Panels(p).Wing == FW.Panels(i).Wing
                        %if the inducer (s) is on the same wing as the
                        %induced (index), we do the following...
                        
                        
                        %                         //################################################################
                        % 		  	//New method of moving the TE points of the index (induced DVE) points.
                        % 				//17 Oct 2014. Bill B
                        % 			//This method moves the points in the freestream direction into the plane
                        % 				//passing through the TE of the S (inducer) having freestream direction
                        % 				//as the normal
                        
                        % //determine vector on inducer from control point to TE (call it tempB)
                        tempB(1) = FW.Panels(p).DVE.xsi(s);
                        tempB(2) = 0;
                        tempB(3) = 0;
                        
                        %                         //make temp B global, call it tempA
                        tempA = star_glob(tempB,FW.Panels(p).DVE.roll(s),FW.Panels(p).DVE.pitch(s),FW.Panels(p).DVE.yaw(s));
                        
                        %                           //add tempA to location of control point, so we are left with location
                        % 				//of the inducer's (S) TE in global coords. call it tempTE
                        
                        %is this the same as averaging the TEcoordL and TEcoordR which we already
                        %have saved?????
                        tempTE=FW.Panels(p).DVE.xo(s,:)+tempA;
                        
                        %                 //vector from TE of S(inducer) to point k on TE of index (induced). call it delX
                        delX(1) = X(k,1) - tempTE(1);
                        delX(2) = X(k,2) - tempTE(2);
                        delX(3) = X(k,3) - tempTE(3);
                        
                        %                 //delX projected into the freestream direction (magnitude)
                        tempS=dot(delX,Temp.u);
                        
                        %                 //vector from TE of s(inducer) to TE of index (induced) projected into the freestream direction (with direction) to make tempB
                        tempB = Temp.u * tempS;
                        
                        %                 // (X[k] - tempB) should be global origin to new point of interest
                        Xstar(1) = X(k,1) - tempB(1);
                        Xstar(2) = X(k,2) - tempB(2);
                        Xstar(3) = X(k,3) - tempB(3);
                        
                    else
                        
                        %                         //DVE 's' (the inducer) and DVE 'index' (the induced one) are
                        % 			//of different wigns
                        
                        Xstar(1) = X(k,1);
                        Xstar(2) = X(k,2);
                        Xstar(3) = X(k,3);
                    end
                    %                       //###################################################################
                    % 		  //computing the velocity induced by the freshly shed wake in the TE
                    % 		  //The very beginning of the shed strip of a vortex sheet belongs to
                    % 		  //the first (most recently pooped out) row of wake DVEs.
                    %
                    % 			//assigning temporary DVE that induces on trailing edge
                    % 			//as Schmid-Goeller discusses in his dissertation,
                    % 			//it has no sweep and belongs to a spanwise strip of wake
                    % 			//elements that starts at the point of interest
                    
                    tempDVE.xo 	 = FW.Panels(p).WakeDVE(Temp.timestep+1).xo(span,:);
                    
                    tempDVE.phiLE	 = 0;
                    tempDVE.phiTE 	 = 0;% //wakePtr[rightnow][span].phiTE;
                    
                    tempDVE.nu		 = FW.Panels(p).WakeDVE(Temp.timestep+1).roll(span);
                    tempDVE.epsilon  = FW.Panels(p).WakeDVE(Temp.timestep+1).pitch(span);
                    tempDVE.psi		 = FW.Panels(p).WakeDVE(Temp.timestep+1).yaw(span);
                    
                    tempDVE.eta		 = FW.Panels(p).WakeDVE(Temp.timestep+1).eta(span);
                    tempDVE.xsi		 = FW.Panels(p).WakeDVE(Temp.timestep+1).xsi(span);
                    
                    tempDVE.A		 = FW.Panels(p).WakeDVE(Temp.timestep+1).A(span);
                    tempDVE.B		 = FW.Panels(p).WakeDVE(Temp.timestep+1).B(span);
                    tempDVE.C		 = FW.Panels(p).WakeDVE(Temp.timestep+1).C(span);
                    
                    tempDVE.singfct	 = 0; %//surfacPtr[s].singfct;
                    
                    %                     //		type = 4;  //vortex sheet reaches from 0.5xsi to xsi
                    type = 1;  %//DVE is only a vortex sheet SHOULD THIS BE 2???? DOES IT NEED A FILAMENT AT LE??
                    
                    %                     //computes induced velocity in X[k] due to DVE tempDVE
                    %AT THIS POINT, Freewake uses Single_DVE_Induced_Velocity
                    
                    %WE CANT BECUASE OUR fcnSingle_DVE_Induced_Velocity ISNT SET UP TO TAKE IN
                    %A TEMP VARIABLE, SO WE WILL CALL fcnDVEInfluenceCoeff DIRECTLY HERE...
                    %(SORT OF LIKE REWRITING fcnSingle_DVE_Induced_Velocity HERE...)
                    
                    %nu, eps, phile and phite go in as deg
                    [a3, b3, c3] = fcnDVEInduction(Temp, Xstar, tempDVE.xo, tempDVE.nu, tempDVE.epsilon,tempDVE.phiLE, tempDVE.phiTE, tempDVE.psi, tempDVE.eta, tempDVE.xsi, type, tempDVE.singfct);
                    
                    a = a3; b = b3; c = c3;
                    
                    
                    if FW.Sym == 1
                        tempA = [ tempDVE.xo(1) -tempDVE.xo(2)  tempDVE.xo(3)];
                        
                        [a3, b3, c3] = fcnDVEInduction(Temp, Xstar, tempA, -tempDVE.nu, tempDVE.epsilon, -tempDVE.phiLE, -tempDVE.phiTE, -tempDVE.psi, tempDVE.eta, tempDVE.xsi, type, tempDVE.singfct);
                        %     disp(b3)
                        a = a + a3;
                        b = b - b3;
                        c = c + c3;
                        
                    end
                    
                    w(1) = tempDVE.A*a(1) + tempDVE.B*b(1) + tempDVE.C*c(1);
                    w(2) = tempDVE.A*a(2) + tempDVE.B*b(2) + tempDVE.C*c(2);
                    w(3) = tempDVE.A*a(3) + tempDVE.B*b(3) + tempDVE.C*c(3);
                    
                    w= w*-1/(4*pi);
                    
                    %add to total ind. velocity on point k
                    w_ind(k,1)= w_ind(k,1)+w(1);
                    w_ind(k,2)= w_ind(k,2)+w(2);
                    w_ind(k,3)= w_ind(k,3)+w(3);
                    
                    %                      //###################################################################
                    % 		  //computing the induced velocities of the remaining wake DVEs of the
                    % 		  //current spanwise location
                    
                    % //loop across wake elements %%%DO WE HAVE TO DO THIS
                    % AGAIN OR CAN WE JUST DO IT ONCE. THE ONLY DIFFERENCE
                    % IS THAT NOW WE DO THE REST OF THE WAKE INSTEAD OF THE
                    % FIRST ROW. WHY GO THROUGH ALL THIS TWICE? SEEMS LIKE
                    % IT WOULD TAKE TIME
                    for time = 1:Temp.timestep ;
                        tempDVE.xo 	 = FW.Panels(p).WakeDVE(time).xo(span,:); %note I dont add 1 to the time here, I've accounted for it in le for loop.
                        
                        tempDVE.phiLE	 = 0;
                        tempDVE.phiTE 	 = 0;% //wakePtr[rightnow][span].phiTE;
                        
                        tempDVE.nu		 = FW.Panels(p).WakeDVE(time).roll(span);
                        tempDVE.epsilon  = FW.Panels(p).WakeDVE(time).pitch(span);
                        tempDVE.psi		 = FW.Panels(p).WakeDVE(time).yaw(span);
                        
                        tempDVE.eta		 = FW.Panels(p).WakeDVE(time).eta(span);
                        tempDVE.xsi		 = FW.Panels(p).WakeDVE(time).xsi(span);
                        
                        tempDVE.A		 = FW.Panels(p).WakeDVE(time).A(span);
                        tempDVE.B		 = FW.Panels(p).WakeDVE(time).B(span);
                        tempDVE.C		 = FW.Panels(p).WakeDVE(time).C(span);
                        
                        tempDVE.singfct	 = FW.Panels(p).WakeDVE(time).singfct(span);
                        
                        %
                        if time~= 1
                            type = 1;
                        else
                            type = 3;
                        end
                        
                        %                     //computes induced velocity in X[k] due to DVE tempDVE
                        %AT THIS POINT, Freewake uses Single_DVE_Induced_Velocity
                        
                        %WE CANT BECUASE OUR fcnSingle_DVE_Induced_Velocity ISNT SET UP TO TAKE IN
                        %A TEMP VARIABLE, SO WE WILL CALL fcnDVEInfluenceCoeff DIRECTLY HERE...
                        %(SORT OF LIKE REWRITING fcnSingle_DVE_Induced_Velocity HERE...)
                        
                        %nu, eps, phile and phite go in as deg
                        [a3, b3, c3] = fcnDVEInduction(Temp, Xstar, tempDVE.xo, tempDVE.nu, tempDVE.epsilon,tempDVE.phiLE, tempDVE.phiTE, tempDVE.psi, tempDVE.eta, tempDVE.xsi, type, tempDVE.singfct);
                        
                        a = a3; b = b3; c = c3;
                        
                        
                        if FW.Sym == 1
                            tempA = [ tempDVE.xo(1) -tempDVE.xo(2)  tempDVE.xo(3)];
                            
                            [a3, b3, c3] = fcnDVEInduction(Temp, Xstar, tempA, -tempDVE.nu, tempDVE.epsilon, -tempDVE.phiLE, -tempDVE.phiTE, -tempDVE.psi, tempDVE.eta, tempDVE.xsi, type, tempDVE.singfct);
                            
                            a = a + a3;
                            b = b - b3;
                            c = c + c3;
                            
                        end
                        
                        
                        w(1) = tempDVE.A*a(1) + tempDVE.B*b(1) + tempDVE.C*c(1);
                        w(2) = tempDVE.A*a(2) + tempDVE.B*b(2) + tempDVE.C*c(2);
                        w(3) = tempDVE.A*a(3) + tempDVE.B*b(3) + tempDVE.C*c(3);
                        
                        w= w*-1/(4*pi);
                        
                        %add to total ind. velocity on point k
                        w_ind(k,1)= w_ind(k,1)+w(1);
                        w_ind(k,2)= w_ind(k,2)+w(2);
                        w_ind(k,3)= w_ind(k,3)+w(3);
                        
                    end %//end loop over time, along a strip in wake
                    
                end%end loop over spanwise position on this panel(inducer)
            end%end loop over panels (inducers)
            
        end %end loop over the three points k (induced)
        
        %   //###########################################################################
        %   //						INDUCED VELOCITY LOOP - END
        %   //###########################################################################
        
        %  //###########################################################################
        %   //				AND NOW: THE FORCE INTEGRATION
        %   //###########################################################################
        
        % 		//Integration of induced forces with Simpson's Rule
        % 		//Integration requires overhanging edges!!
        % 		//See also KHH linees 2953 - 2967, A23SIM
        
        % 		//Kutta-Joukowski at left (2) edge
        tempA = cross(w_ind(2,:),S);			%// w1xS
        gamma1  = A-B*eta8+C*eta8*eta8;		%//gamma1
        R1 = tempA*gamma1;
        
        % 		//Kutta-Joukowski at center (1)
        tempA = cross(w_ind(1,:),S);				%// woxS
        gammao  = A;
        Ro = tempA*gammao;
        
        %  		//Kutta-Joukowski at right (3) edge
        tempA = cross(w_ind(3,:),S);				%// w2xS
        gamma2  = A+B*eta8+C*eta8*eta8;
        R2 = tempA*gamma2;
        
        % //The resultierende induced force of element l is
        % 		//determined by numerically integrating forces acros element
        % 		//using Simpson's Rule with overhaning parts
        
        R(1)  = (R1(1)+4*Ro(1)+R2(1))*eta8/3;			%//Rx
        R(2)  = (R1(2)+4*Ro(2)+R2(2))*eta8/3;			%//Ry
        R(3)  = (R1(3)+4*Ro(3)+R2(3))*eta8/3;			%//Rz
        
        % 		//plus overhanging parts
        R(1) = R(1)+((7*R1(1)-8*Ro(1)+7*R2(1))*(eta-eta8)/3); %//Rx
        R(2) = R(2)+((7*R1(2)-8*Ro(2)+7*R2(2))*(eta-eta8)/3); %//Ry
        R(3) = R(3)+((7*R1(3)-8*Ro(3)+7*R2(3))*(eta-eta8)/3); %//Rz
        
        %         //#########################################################################
        % 		//the DRAG FORCE/density is the induce force in eD direction
        % 	//#########################################################################
        D_force(TEindex) = dot(R,eD);
        
        %         //add all partial drag/lift/side values [force/density]
        CDi = CDi+D_force(TEindex);
        
        TEindex = TEindex+1; %advancing span index of all TE DVEs
    end %end loop over all the induced elements
    
    
end %end loop over all the induced panels

% //#############################################################################
% //							FORCE LOOP - END
% //#############################################################################


% //non-dimensionalize
tempS = 0.5*FW.Uinf*FW.Uinf*Aircraft.Reference.S;
CDi= CDi/tempS;

if FW.Sym == 1
    
    CDi= CDi *2; %//sym. geometry and flow, twice the drag
end

%% HEREITH LIES SOME STUFF ABOUT TOTAL FORCES. MAYBE ITS FOR THE VISCOUS WRAP BUT I DIDNT COPY IT. I WILL IF NEEDED. BILLLBILLBILLBILL


Output.(Temp.AI)(Temp.timestep+2).CDi = CDi;

Output.(Temp.AI)(Temp.timestep+2).e = ((Output.(Temp.AI)(Temp.timestep+2).CL)^2)/(pi*Aircraft.Reference.AR*CDi);
Output.(Temp.AI)(Temp.timestep+2).deltae = abs(Output.(Temp.AI)(Temp.timestep+2).e - Output.(Temp.AI)(Temp.timestep+1).e);

end %end function





