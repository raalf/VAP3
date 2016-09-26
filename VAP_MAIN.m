clc
clear

tic

disp('===========================================================================');
disp('VAP (Based on FreeWake 2015)');
disp('Running Version 2016.09  .                             .');
disp('Includes stall model    //                             \\');
disp('No trim solution       //                               \\');
disp('                      //                                 \\');
disp('                     //                _._                \\');
disp('                  .---.              .//|\\.              .---.');
disp('         ________/ .-. \_________..-~ _.-._ ~-..________ / .-. \_________');
disp('                 \ ~-~ /   /H-     `-=.___.=-''     -H\   \ ~-~ /');
disp('                   ~~~    / H          [H]          H \    ~~~');
disp('                         / _H_         _H_         _H_ \');
disp('                           UUU         UUU         UUU');
disp('===========================================================================');
disp(' ');

%% Best Practices
% 1. Define wing from one wingtip to another in one direction
% 2. When using symmetry, define from symmetry plane outward

%% Reading in geometry

% strFILE = 'VAP input.txt';
% strFILE = 'VAP 2panel.txt';
strFILE = 'VAP 3split.txt';

[flagRELAX, flagSTEADY, valAREA, valSPAN, valCMAC, valWEIGHT, ...
    seqALPHA, seqBETA, valKINV, valDENSITY, valPANELS, matGEOM, vecSYM, ...
    vecAIRFOIL, vecN, vecM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, ...
    valFTURB, valFPWIDTH, valDELTAE, valDELTIME, valMAXTIME, valMINTIME, ...
    valINTERF] = fcnVAPREAD(strFILE);

% strFILE = 'input.txt';
%
% [flagRELAX, flagSTEADY, valAREA, valSPAN, valCMAC, valWEIGHT, ...
%     seqALPHA, seqBETA, valKINV, valDENSITY, valPANELS, matGEOM, vecSYM, ...
%     vecAIRFOIL, vecN, vecM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, ...
%     valFTURB, valFPWIDTH, valDELTAE, valDELTIME, valMAXTIME, valMINTIME, ...
%     valINTERF] = fcnFWREAD(strFILE);

flagPLOT = 1;

%% Discretize geometry into DVEs

[vecDVECTLPT, vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, vecDVEMCSWP, vecDVETESWP, ...
    vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVEAREA, vecDVENORM, ...
    matVLST, matDVE, valNELE, matADJE, vecDVESYM, vecDVETIP] = fcnGENERATEDVES(valPANELS, matGEOM, vecSYM, vecN, vecM);


%% Plotting Wing

if flagPLOT == 1
    [hFig2] = fcnPLOTBODY(0, valNELE, matDVE, matVLST, vecDVECTLPT, vecDVENORM);
    [hLogo] = fcnPLOTLOGO(0.97,0.03,14,'k','none');
end

%% Add boundary conditions to D-Matrix

all_DVEs = [1:valNELE]';

% Vorticity and circulation equations between DVEs ----------------------------------------------------------------------

idx1 = matADJE(:,2) == 2; % Borders on 2nd local edge (right side)
idx2 = matADJE(:,2) == 4; % Borders on 4th local edge (left side)

len = length(idx1(idx1>0,1));

% Mapping equations evaluated on edge 2 to equations evaluated on edge 4
% This returns a vector where we can see for edge 2 equations, which is the
% matching edge 4 equation
% This will work well for split vorticity but it will need modification for split circulation later on
% [~, LB] = ismembertol([matADJE(idx2,3) matADJE(idx2,1)],[matADJE(idx1,1) matADJE(idx1,3)],'ByRows',true,'OutputAllIndices',true);
[~, LB] = ismembertol([matADJE(idx1,1) matADJE(idx1,3)],[matADJE(idx2,3) matADJE(idx2,1)],'ByRows',true,'OutputAllIndices',true);
eqn_num = cell2mat(LB);

idx3 = matADJE(idx1,1);
idx4 = matADJE(idx2,1);
idx5 = idx4(eqn_num); % Corresponding edge 4 DVE for every edge 2 DVE

% dgamma1 = dgamma2
% B + 2*C*eta = B - 2*C*eta
% !!!!ACROSS A SPLIT, DGAMMA IS CONSTANT FOR ALL PANELS!!!!

% dgamma = B + 2*C*eta
dgamma1 = [zeros(len,1) ones(len,1) 2.*vecDVEHVSPN(idx3)];

% dgamma2 = B - 2*C*eta
% Multiplied by -1 as we bring it to the other side of the equal sign
dgamma2 = [zeros(len,1) ones(len,1) -2.*vecDVEHVSPN(idx5)].*-1;

% Getting appropriate row and column numbers to assign the above
% dgamma1 and dgamma2 into the D-matrix
col1 = reshape([repmat((idx3.*3)-2,1,3) + repmat([0:2], len,1)]',[],1);
col2 = reshape([repmat((idx5.*3)-2,1,3) + repmat([0:2], len,1)]',[],1);
rows = reshape(repmat([1:len]',1,3)',[],1);

vort = zeros(len,valNELE*3);

vort(sub2ind(size(vort),rows,col1)) = reshape(dgamma1',[],1);
vort(sub2ind(size(vort),rows,col2)) = reshape(dgamma2',[],1);

% gamma1 = gamma2
% A + B*eta + C*eta^2 = A - B*eta + C*eta^2
% For split:
% gamma12 + gamma22 + ... + gamma2N = gamma41 + gamma42 + ... + gamma4N
% where the left side of the equal side is DVEs with edge 2 at the split
% and the right side is DVEs with edge 4 at the split

% Firstly, the splits are dealt with

% DVE numbers with local edge 2 at the split, along with the degree of the split
[idx13, idx14] = histcounts(idx3, all_DVEs);
dsplit2 = idx14(idx13 > 1);
sdegree = idx13(idx13 > 1);

% Finding the corresponding DVEs with local edge 4 at the split
[idx26,~] = find(repmat(matADJE(idx2,3),1,length(dsplit2)) == repmat(dsplit2',length(matADJE(idx2,3)),1));

% Reassigning dsplit2 so it lines up with dsplit4
dsplit2 = matADJE(idx2,3);
dsplit2 = dsplit2(idx26);
dsplit4 = idx4(idx26);

circ = zeros(len,valNELE*3);
count = 1;

udsplit2 = unique(dsplit2);
for i = 1:length(udsplit2)
    
    idx = find(dsplit2 == udsplit2(i));
    len2 = length(idx);
    % DVE on left of split
    gamma1s = [1 vecDVEHVSPN(udsplit2(i)) vecDVEHVSPN(udsplit2(i)).^2];
    % DVEs on right of split
    gamma2s = [ones(len2,1) -vecDVEHVSPN(dsplit4(idx)) vecDVEHVSPN(dsplit4(idx)).^2].*-1;
    
    col1 = reshape([repmat((udsplit2(i).*3)-2,1,3) + repmat([0:2], 1, 1)]',[],1);
    col2 = reshape([repmat((dsplit4(idx).*3)-2,1,3) + repmat([0:2], len2, 1)]',[],1);
    
    circ(sub2ind(size(circ),repmat(count,length(col1),1),col1)) = reshape(gamma1s',[],1);
    circ(sub2ind(size(circ),repmat(count,length(col2),1),col2)) = reshape(gamma2s',[],1);
    
    count = count + 1;
end

% Now we move on to regular circulation between two neighbouring panels
% gamma1 = gamma2
% A + B*eta + C*eta^2 = A - B*eta + C*eta^2

% DVE numbers with local edge 2 neighbouring another DVE
d2202 = idx14(idx13 == 1);

% Finding the DVE numbers with local edge 4 neighbouring the above corresponding DVEs
[idx27,~] = find(repmat(matADJE(idx1,1),1,length(d2202)) == repmat(d2202',length(matADJE(idx1,1)),1));
d2204 = matADJE(matADJE(:,2) == 2,3);
d2204 = d2204(idx27);

len3 = length(d2204);

gamma1 = [ones(len3,1) vecDVEHVSPN(d2202) vecDVEHVSPN(d2202).^2];
gamma2 = [ones(len3,1) -vecDVEHVSPN(d2204) vecDVEHVSPN(d2204).^2].*-1;

% Getting appropriate row and column numbers to assign the above
% gamma1 and gamma2 into the D-matrix
col1 = reshape([repmat((d2202.*3)-2,1,3) + repmat([0:2], len3,1)]',[],1);
col2 = reshape([repmat((d2204.*3)-2,1,3) + repmat([0:2], len3,1)]',[],1);
rows = reshape(repmat([1:len3]',1,3)',[],1) + count-1;

circ(sub2ind(size(circ),rows,col1)) = reshape(gamma1',[],1);
circ(sub2ind(size(circ),rows,col2)) = reshape(gamma2',[],1);

circ(~any(circ,2),:) = [];

% END vorticity and circulation equations between DVEs ----------------------------------------------------------------

% Vorticity equations at symmetry -------------------------------------------------------------------------------------
if ~isempty(vecDVESYM) == 1
    
    % DVE number where we have a symmetry condition
    [idx10, ~] = find(vecDVESYM);
    len = length(idx10);
    locedge = vecDVESYM(idx10);
    
    % dgamma_sym = 0
    % B + 2*C*eta = 0 for right edge (local edge 2), though I doubt this one will be used
    % B - 2*C*eta = 0 for left edge (local edge 4)
    
    dgamma_sym(locedge == 2,:) = [zeros(length(idx10(locedge == 2)),1) ones(length(idx10(locedge == 2)),1) 2.*vecDVEHVSPN(idx10(locedge == 2))];
    dgamma_sym(locedge == 4,:) = [zeros(length(idx10(locedge == 4)),1) ones(length(idx10(locedge == 4)),1) -2.*vecDVEHVSPN(idx10(locedge == 4))];
    
    % Getting appropriate row and column numbers to assign the above
    % dgamma1 and dgamma2 into the D-matrix
    col3 = reshape([repmat((idx10.*3)-2,1,3) + repmat([0:2], len,1)]',[],1);
    rows = reshape(repmat([1:len]',1,3)',[],1);
    
    vort_sym = zeros(len,valNELE*3);
    
    vort_sym(sub2ind(size(vort_sym),rows,col3)) = reshape(dgamma_sym',[],1);
    
end
% END vorticity equations at symmetry -----------------------------------------------------------------------------------

% Circulation equations at wingtip --------------------------------------------------------------------------------------

idx21 = matADJE(:,2) == 2; % Borders on 2nd local edge (right side)
idx22 = matADJE(:,2) == 4; % Borders on 4th local edge (left side)

len = length(idx21(idx21>0,1));

% Finding DVEs with nothing bordering on local edge 2
% This is done by finding which DVEs are absent from matADJE for edge 2
[~, lb2] = ismembertol(all_DVEs,unique(matADJE(idx21,1),'rows'),'ByRows',true,'OutputAllIndices',true);
tip_locedge2 = all_DVEs(cell2mat(lb2) == 0);

% Removing entries in tip_locedge2 that are symmetry conditions
if ~isempty(vecDVESYM) == 1
    temp20 = idx10(locedge == 4); % Symmetry edges with local edge 4
    [~, LB3] = ismembertol(temp20,tip_locedge2,'ByRows',true,'OutputAllIndices',true);
    tip_locedge2(nonzeros(cell2mat(LB3))) = [];
end

% test2 = unique(matADJE(idx22,1:2),'rows')

locedge = vecSYM(idx10);
% The input defines local edge number as 1 or 2, so here it is
% transformed to 2 or 4
locedge(locedge == 1) = 4;

% dgamma_sym = 0
% B + 2*C*eta = 0 for right edge (local edge 2), though I doubt this one will be used
% B - 2*C*eta = 0 for left edge (local edge 4)

dgamma_sym(locedge == 2,:) = [zeros(length(idx10(locedge == 2)),1) ones(length(idx10(locedge == 2)),1) 2.*vecDVEHVSPN(idx10(locedge == 2))];
dgamma_sym(locedge == 4,:) = [zeros(length(idx10(locedge == 4)),1) ones(length(idx10(locedge == 4)),1) -2.*vecDVEHVSPN(idx10(locedge == 4))];

% Getting appropriate row and column numbers to assign the above
% dgamma1 and dgamma2 into the D-matrix
col3 = reshape([repmat((idx10.*3)-2,1,3) + repmat([0:2], len,1)]',[],1);
rows = reshape(repmat([1:len]',1,3)',[],1);

vort_sym = zeros(len,valNELE*3);

vort_sym(sub2ind(size(vort_sym),rows,col3)) = reshape(dgamma_sym',[],1);

% END Circulation equations at wingtip ----------------------------------------------------------------------------------
%}
%% Add kinematic conditions to D-Matrix

%% Create D-Resultant, solve D-Matrix

%% Timestep to solution
%   Move wing
%   Generate new wake elements
%   Create W-Matrix and W-Resultant
%   Solve W-Matrix
%   Relaxation procedure (Relax, create W-Matrix and W-Resultant, solve W-Matrix)
%   Calculate surface normal forces
%   Calculate DVE normal forces
%   Calculate induced drag
%   Calculate cn, cl, cy, cdi
%   Calculate viscous effects

%% Viscous wrapper

toc

% whos