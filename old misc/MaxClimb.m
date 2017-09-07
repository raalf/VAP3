function [wclimb]...
       = MaxClimb(CL,CD,Rrecip,wmaxth,WSroh)

g = 9.81;

wclimbold = -1000; %initializing climb rate
wclimb = -200;%initializing climb rate

r=1/Rrecip; %assigning first turning radius at outer edge of thermal

while(wclimb > wclimbold)
    wclimbold = wclimb;
%     Prof Tom Eq. 87
    Vsc = CD*(CL^-1.5)*sqrt(WSroh)*(1-(WSroh/(r*g*CL))^2)^(-0.75);
    Rnormal = r*Rrecip; % normalized radius (wrt to thermal radius)
    wthermal = wmaxth*(1-Rnormal*Rnormal);%thermal strenght at radius r
    wclimb = wthermal-Vsc; %climb rate    
    r=r-0.1; %moving on to next radius, 2 meters furhter inside
end


%radius of max climb and subsequent sinkrate
%[r+2,57.3*acos(sqrt(1-(WSroh/(r*g*CL))^2))]
end