function [ w_wake ] = fcnWake_DVE_Vel_Induction( FW,Aircraft,Temp,P )

%copy of Wake_DVE_Vel_Induction from freewake 2015beta
%finds the velocity induced by the entire wake at point p

%{
//this function computes the velocity in point P that is indunced by the
//DVEs in the wake.
//
//input:
// info		general element information
// P		point of where induced velocities are being computed
// wakePtr	wake DVE information
// timestep	current timestep
//
//output:
// w_wake	induced velocity in point P
%}

% //initializing velocity w_wake
w_wake(1)=0;
w_wake(2)=0;
w_wake(3)=0;


%loop over each panel
for i = 1:Aircraft.General.Panels
    %i = panel number
    %Step 1. First we take care of the special cases, which are the first and last
    %row of elements in the wake.
    
    %loop over each n in the panel
    for j = 1:FW.Panels(i).n
        %j = element index on panel i
        %         //computing induced velocity at P due to most downstream wake DVEs
        % 		//The very first row of wake DVEs consist of semi-infinite vortex
        % 		//sheets.  At the first time step they have a vortex filametn at
        % 		//their leading edge
        
        %last input to fcnSingle_DVE_Induced_Velocity is type 2, for wake
               
        if(Temp.timestep<1)
            w_ind=fcnSingle_DVE_Induced_Velocity(FW,Temp,0,i,j,P,-3,2);
        else
            w_ind=fcnSingle_DVE_Induced_Velocity(FW,Temp,0,i,j,P,3,2);
        end
        
        %       //adding "starting" DVE influence, type -3 or 3
        w_wake(1)  =  w_wake(1) + w_ind(1);
        w_wake(2)  =  w_wake(2) + w_ind(2);
        w_wake(3)  =  w_wake(3) + w_ind(3);
        
        
        %         //computing ind. vel. at P due to wake DVE just aft of trailing edge
        % 		//These DVEs have a vortex filament only at their leading edge
        if(Temp.timestep>=1) 
            
            w_ind = fcnSingle_DVE_Induced_Velocity(FW,Temp,Temp.timestep,i,j,P,2,2);
            
            % 			//adding last timestep DVE influence, type 2
            w_wake(1)  =  w_wake(1) + w_ind(1);
            w_wake(2)  =  w_wake(2) + w_ind(2);
            w_wake(3)  =  w_wake(3) + w_ind(3);
            
        end
        
    end
    
    
    %Step 2.  Now we loop over the rest of the elements inbetween the first
    %and last row of elements, on this panel
    
    %     //loop over wake DVE's, time and span dependent
    %     for(time=1;time<timestep;time++)
    for h = 1:Temp.timestep-1 % going through all timesteps
        
        for j = 1:FW.Panels(i).n % going through all spanwise wakeDVEs per panel
            
            % 		//computing induced velocity of wake DVE [time][span] on point P.
            w_ind = fcnSingle_DVE_Induced_Velocity(FW,Temp,h,i,j,P,1,2);
            
            % 		//adding delta velocities, element type 1.
            w_wake(1)  =  w_wake(1) + w_ind(1);
            w_wake(2)  =  w_wake(2) + w_ind(2);
            w_wake(3)  =  w_wake(3) + w_ind(3);
        end
        
        
    end
end


