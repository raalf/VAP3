function [ w_surface ] = fcnSurface_DVE_Vel_Induction ( FW,Aircraft,Temp,P)

%copy of Surface_DVE_Vel_Induction from freewake 2015beta
%finds the velocity induced by the entire surface at point p

%initialize for superspeed
w_surface(1) = 0;
w_surface(2) = 0;
w_surface(3) = 0;

%loop over panels
for ii = 1:Aircraft.General.Panels % looping through influencing elements (all other DVEs)
    
    
    %loop over all DVEs on panel i
    for jj = 1:length(FW.Panels(ii).DVE.Index)
        
        %DVE type 0, surface type 1
        w_ind = fcnSingle_DVE_Induced_Velocity(FW,Temp,0,ii,jj,P,0,1);
       
        %add delta velocity
        w_surface(1) = w_surface(1) + w_ind(1);
        w_surface(2) = w_surface(2) + w_ind(2);
        w_surface(3) = w_surface(3) + w_ind(3);
        
    end
end

%{
//initializing w_surface
	w_surface[0]=0;
	w_surface[1]=0;
	w_surface[2]=0;


	//loop over surface DVEs that induce velocities
	for(j=0;j<info.noelement;j++)
	{
			//computing induced velocity of j-th surface DVE on point P.
			Single_DVE_Induced_Velocity(info,surfacePtr[j],P,w_ind,0);

			//adding delta velocities
			w_surface[0] += w_ind[0];
			w_surface[1] += w_ind[1];
			w_surface[2] += w_ind[2];
            }
%}