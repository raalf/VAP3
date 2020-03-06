% Airfoil Data Generator for VAP3.1
% 
% Input:	airfoil.dat or coordinates
%           Alpha Sweep     (Default: 150e3 200e3 250e3 500e3 1e6)
%           Re Sweep        (Default: -10:0.5:20)
% 
% Output:	polar 3D matrix
%           polar meshgrid
% 



airfoilFilename = 'RAALF-B8';

% determines if the airfoilFilename is a .mat or .dat
% if it is a .mat, load the coordinates matrix
% load ('RAALF-B8');
% if it is a .dat, import the coordinates to a MATLAB matrix
airfoil_name = 'NACA0003';
coord = dlmread([airfoil_name, '.dat'],'',1,0);




AlphaRange = -10:0.5:20;
% if not specified, default values are [-10:0.5:20];

ReRange    = 150000:150000:5e6;
% ReRange = 50000:25000:1e6;
% ReRange = 1e6:1e6:6e6;
% ReRange = [150e3 200e3 250e3 500e3 1e6];
% if not specified, default values are [150e3 200e3 250e3 500e3 1e6];

% preallocate polar output from xfoil
pol = nan(length(AlphaRange),9,length(ReRange));
%%
% Verify Alpha and Re lists
AlphaRange = sort(unique(AlphaRange));
countAlpha = length(AlphaRange(:));
ReRange    = sort(unique(ReRange));
countRe    = length(ReRange(:));

%%
% check if aux_files folder exist


% par for / for
parfor n = 1:length(ReRange)
    
    tempPol = nan(countAlpha, 9);
    Re = ReRange(n);
    
    % store xfoil output in tempPol
    xfoilPol = xfoil(coord, AlphaRange, Re, 10, [1 1], false);
    
    % identify missing alphas by comparing the tempPol with AlphaRange
    [idxPol,idxtempPol] = ismember(AlphaRange,xfoilPol(:,1));
    
    % store only valid xfoil results in the correct array order for
    % meshgrid generation down the line
    tempPol(idxPol,:) = xfoilPol(idxtempPol(idxtempPol~=0),:); 
    
    pol(:,:,n) = tempPol;
end



%%
% Deal with results
meshgridAlpha = reshape(pol(:,1,:),countAlpha,countRe);
meshgridCl    = reshape(pol(:,2,:),countAlpha,countRe);
meshgridCd    = reshape(pol(:,3,:),countAlpha,countRe);
meshgridCm    = reshape(pol(:,5,:),countAlpha,countRe);
meshgridRe    = reshape(pol(:,8,:),countAlpha,countRe);


% detech CLmax at each Re
[Clmax,b] = max(meshgridCl,[],1);
idxClmax = sub2ind([countAlpha, countRe],b,1:countRe);

ClmaxRe = meshgridRe(idxClmax);
ClmaxAlpha = meshgridAlpha(idxClmax);
ClmaxCd = meshgridCd(idxClmax);


%%
clf
scatter3(meshgridCd(:),meshgridCl(:),meshgridRe(:),10,meshgridAlpha(:),'filled');
colorbar
hold on
plot3(ClmaxCd,Clmax,ClmaxRe,'k-o','LineWidth',1)
hold off



%%
% clf
% surf(meshgridCd,meshgridCl,meshgridRe,meshgridAlpha)
% hold on
% plot3(ClmaxCd,Clmax,ClmaxRe,'r-o','LineWidth',1)
% hold off
% 
% 
% %%
% clc
% clf
% idxgood=~(isnan(meshgridCl) | isnan(meshgridRe) | isnan(meshgridCd));
% [Y1,Y2] = meshgrid(0:0.1:2, 0:0.5e6:5e6);
% Vq = interpn(meshgridCl(idxgood),meshgridRe(idxgood),meshgridCd(idxgood),Y1,Y2,'cubic');
% % mesh(Y1,Y2,Vq)
% % title('Refined Grid');
% 
% 
% 
% 
% %%
% idxNotNan = ~(isnan(ClmaxRe)|isnan(Clmax)|isnan(CLmaxCd));
% a = interp1(ClmaxRe(idxNotNan),Clmax(idxNotNan),1.5e6,'linear');
% a = interp1(ClmaxRe(idxNotNan),Clmax(idxNotNan),1.5e6,'linear');
% hold on
% 
% 
% 
% hold off












































