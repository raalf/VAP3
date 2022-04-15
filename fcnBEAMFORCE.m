function [SURF, OUTP, INPU] = fcnBEAMFORCE(SURF, OUTP, COND, INPU, FLAG, valTIMESTEP)
% Compute the forces and moments applied at each beam node

[ledves, ~, ~] = find(SURF.vecDVELE > 0);
lepanels = SURF.vecDVEPANEL(ledves);

isCurWing = SURF.vecWINGTYPE(ledves) == 1;

idxdve = uint16(ledves(isCurWing));
idxpanel = lepanels(isCurWing);

m = INPU.vecM(idxpanel);

m = m(1);

% Matrix of how much we need to add to an index to get the next chordwise element
% It is done this way because n can be different for each panel. Unlike in the wake,
% we can't just add a constant value to get to the same spanwise location in the next
% row of elements
tempm = repmat(INPU.vecN(idxpanel), 1, m).*repmat([0:m-1],length(idxpanel~=0),1);

rows = repmat(idxdve,1,m) + uint16(tempm);

idx1 = SURF.vecDVELE == 1;

s(idx1,:) =( SURF.matVLST(SURF.matDVE(idx1,2),:) - SURF.matVLST(SURF.matDVE(idx1,1),:) )  ./ repmat((SURF.vecDVEHVSPN(idx1).*2),1,3); % Vector along bound vortex

% vector along midchord for remaining elements
if any(idx1 == 0)
    tempa = zeros(SURF.valNELE,1);
    tempa(idx1==0,1) = tan(SURF.vecDVEMCSWP(idx1==0));
    tempa= [tempa ones(SURF.valNELE,1)  zeros(SURF.valNELE,1)];
    
    % global vectors to find ind. velocitiess
    s(idx1 ==0,:)= fcnSTARGLOB(tempa(idx1==0,:), SURF.vecDVEROLL(idx1==0), SURF.vecDVEPITCH(idx1==0), SURF.vecDVEYAW(idx1==0));
    
    %using the verticies to do this doesn't match the results from above for idx1 ==0
    % s(idx1==0,:) =(( (SURF.matVLST(SURF.matDVE(idx1,2),:) + SURF.matVLST(SURF.matDVE(idx1,3),:))./2 ) -...
    %     ( (SURF.matVLST(SURF.matDVE(idx1,1),:) + SURF.matVLST(SURF.matDVE(idx1,4),:))./2 ) )./ (SURF.vecDVEHVSPN.*2);
    
end

tempb = cross(SURF.matUINF, s, 2);

uxs = sqrt(sum(abs(tempb).^2,2));

A(idx1) = SURF.matCOEFF(idx1,1);
B(idx1) = SURF.matCOEFF(idx1,2);
C(idx1) = SURF.matCOEFF(idx1,3);

idx2 = SURF.vecDVELE == 1; %idx2 since we need to do this even for triangles
dvenum = find(idx2==0); %dvenum in question
idxf = SURF.matADJE((ismember(SURF.matADJE(:,1), dvenum) & SURF.matADJE(:,2) == 1),3); %upstream dve num
A(idx2 ==0) = (SURF.matCOEFF(idx2==0,1)-SURF.matCOEFF(idxf,1));
B(idx2 ==0) = (SURF.matCOEFF(idx2==0,2)-SURF.matCOEFF(idxf,2));
C(idx2 ==0) = (SURF.matCOEFF(idx2==0,3)-SURF.matCOEFF(idxf,3));

% Store effective A, B, and C values for each DVE
SURF.vecDVEA = A';
SURF.vecDVEB = B';
SURF.vecDVEC = C';

liftperspan = COND.valDENSITY.*(A(SURF.vecWINGTYPE == 1)' - B(SURF.vecWINGTYPE == 1)'.*SURF.vecDVEHVSPN(SURF.vecWINGTYPE == 1) +...
    C(SURF.vecWINGTYPE == 1)'.*SURF.vecDVEHVSPN(SURF.vecWINGTYPE == 1).*SURF.vecDVEHVSPN(SURF.vecWINGTYPE == 1)).*uxs(SURF.vecWINGTYPE == 1).*SURF.matNORMDIR(SURF.vecWINGTYPE == 1,:);

% Beam aerodynamic force per unit length
for i = 1:size(rows,1)
    f_aero(i,:) = sum(liftperspan(rows(i,:),:),1);
end
f_aero = [f_aero; [0,0,0]]; % Zero lift at beam tip

%--------------------------------------------------------------------------
delta = SURF.matAEROCNTR - SURF.matBEAMLOC(:,:,valTIMESTEP-1); % Moment arm between aerodynamic center and beam axis
OUTP.DEBUG.delta(:,:,valTIMESTEP) = delta;
mom_aero = cross(delta(1:end-1,:),dot(f_aero(1:end-1,:),SURF.matDVENORM(ledves(isCurWing),:),2).*SURF.matDVENORM(ledves(isCurWing),:))...
    + [zeros(length(delta)-1,1) 0.5.*COND.valDENSITY.*COND.vecVEHVINF.*COND.vecVEHVINF.*SURF.vecMAC.*SURF.vecMAC.*OUTP.vecCMDIST(1:length(delta)-1) zeros(length(delta)-1,1)];

if FLAG.FLIGHTDYN == 1
    f_inertial = INPU.vecLM(1:end-1).*(repmat([0,0,-9.81],INPU.valNSELE-1,1) - SURF.matBEAMACC(valTIMESTEP-1,:)); % Force due to beam inertia and gravity loads
    
    deltamom = SURF.matBEAMCGLOC(:,:,valTIMESTEP-1) - SURF.matBEAMLOC(:,:,valTIMESTEP-1);
    mom_inertial = cross(deltamom(1:end-1,:),f_inertial);
else
    f_inertial = INPU.vecLM(1:end-1).*(repmat([0,0,-9.81],INPU.valNSELE-1,1)); % Force due to beam inertia and gravity loads
    
    deltamom = SURF.matBEAMCGLOC(:,:,valTIMESTEP-1) - SURF.matBEAMLOC(:,:,valTIMESTEP-1);
    mom_inertial = cross(deltamom(1:end-1,:),f_inertial);
end

OUTP.vecBEAMFORCE = f_aero(1:end-1,:) + f_inertial; % Beam force is sum of aerodynamic and inertial loads

OUTP.vecBEAMFORCE = [dot(OUTP.vecBEAMFORCE,SURF.matDVENORM(SURF.vecDVELE(SURF.vecDVEWING == 1) == 1,:),2); 0]; % Beam force normal to beam
OUTP.vecBEAMMOM = [dot(mom_aero + mom_inertial,s(SURF.vecDVELE(SURF.vecDVEWING == 1) == 1,:),2); 0]; % Beam moment along spanwise axis



end