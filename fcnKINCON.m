function [matD, collpts] = fcnKINCON(matD, SURF, INPU, FLAG, valTIMESTEP)

% Flow tangency is to be enforced at all control points on the surface HDVEs
% In the D-Matrix, dot (a,b,c) of our influencing HDVE with the normal of the point we are influencing on

%% Adding king kong conditions to bottom 1/3 of D-matrix

tedvenum = [];

% Points we are influencing:
fpg = SURF.matCENTER;

% Move trailing edge element collocation point from matCENTER location to
% the TE
% fpg_te = SURF.matCENTER(SURF.vecDVETE == 3,:);
% fpg_te = fcnGLOBSTAR(fpg_te, SURF.vecDVEROLL(SURF.vecDVETE == 3), SURF.vecDVEPITCH(SURF.vecDVETE == 3), SURF.vecDVEYAW(SURF.vecDVETE == 3));
% 
% lambda = 0.9;
% fpg_te = fpg_te + [lambda*SURF.vecDVEHVCRD(SURF.vecDVETE == 3),zeros(length(find(SURF.vecDVETE == 3)),2)];
% fpg_te = fcnSTARGLOB(fpg_te, SURF.vecDVEROLL(SURF.vecDVETE == 3), SURF.vecDVEPITCH(SURF.vecDVETE == 3), SURF.vecDVEYAW(SURF.vecDVETE == 3));
% fpg(SURF.vecDVETE == 3,:) = fpg_te;
collpts = fpg; % Save new collocation points to output

% List of DVEs we are influencing from (each one for each of the above fieldpoints)
len = length(fpg(:,1));
dvenum = reshape(repmat(1:SURF.valNELE,len,1),[],1);

% Identify TE DVEs and determine what rows in dvenum they correspond to
tedves = find(SURF.vecDVETE > 0);
for i = 1:length(tedves)
    tedvenum = [tedvenum; find(dvenum == tedves(i))];
end    

fpg = repmat(fpg,SURF.valNELE,1);

% DVE type 0 is a surface element
dvetype = zeros(length(dvenum),1);

% if valTIMESTEP < 1
if FLAG.STEADY == 1
    dvetype(tedvenum) = 2; % Trailing edge DVEs are type 2 (no filament at TE) --> GB change from 8/10/21; translated from C++ nonsense to a real language by MM
end
% end

%set singfct to zero temporarily. Why? Not gonna do this, we have NaN CL because of this. T.D.K 2017-04-26
% [a, b, c] = fcnDVEINF(dvenum, dvetype, fpg, zeros(size(SURF.vecK,1),1), SURF.matDVE, SURF.matVLST, SURF.vecDVEHVSPN, SURF.vecDVEHVCRD,SURF.vecDVEROLL, SURF.vecDVEPITCH, vecDVEYAW, SURF.vecDVELESWP, SURF.vecDVETESWP, vecSYM);
[a, b, c] = fcnDVEINF(dvenum, dvetype, fpg, SURF.vecK, SURF.matDVE, SURF.matVLST, SURF.vecDVEHVSPN, SURF.vecDVEHVCRD, SURF.vecDVEROLL, SURF.vecDVEPITCH, SURF.vecDVEYAW, SURF.vecDVELESWP, SURF.vecDVETESWP, INPU.vecSYM, FLAG.GPU);

% List of normals we are to dot the above with
normals = repmat(SURF.matDVENORM,SURF.valNELE,1); % Repeated so we can dot all at once

% Dotting a, b, c with the normals of the field points
temp60 = [dot(a,normals,2) dot(b,normals,2) dot(c,normals,2)];

% Reshaping and inserting into the bottom of the D-Matrix
rows = [1:len]';

king_kong = zeros(len, SURF.valNELE*3);
king_kong(rows,:) = reshape(permute(reshape(temp60',3,[],SURF.valNELE),[2 1 3]),[],3*SURF.valNELE,1);

matD = [matD; king_kong];

end

