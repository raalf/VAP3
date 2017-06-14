function [matNEWGEOM, valPANELS, vecN, vecM, vecSYM] = fcnMULTIBLADE(valPANELS, valNUMB, vecN, vecM, vecSYM, matGEOM)
% This function creates the geometry for multiple blades.
% **THIS WILL NOT WORK FOR THIS SYSTEM AND IS ONLY HERE FOR REFERENCE**

% Angle of each blade from original (+ve CCW & in rad)
tempTHETA = (0:2*pi/valNUMB:2*pi)';
tempTHETA(valNUMB+1) = [];

% Create rotation matrix
tempTOP = [cos(tempTHETA), -sin(tempTHETA), zeros([valNUMB,1])];
tempMID = [sin(tempTHETA), cos(tempTHETA), zeros([valNUMB,1])];
tempLOW = [zeros([valNUMB, 1]), zeros([valNUMB, 1]), ones([valNUMB, 1])];
tempROTATE = reshape(tempTOP',[1, 3, valNUMB]);
tempROTATE(2,:,:) = reshape(tempMID',[1, 3, valNUMB]);
tempROTATE(3,:,:) = reshape(tempLOW',[1, 3, valNUMB]);

% Rearrange matGEOM to only x,y,z components in a n by 3 matrix
tempGEOM = (reshape((permute(matGEOM(:,1:3,:),[2,1,3])),[3,0.2*numel(matGEOM)]))';

% Rotate the x,y,z coordinates.
for i = 1:valNUMB
    matNEWGEOM(:,:,i) = (tempROTATE(:,:,i)*tempGEOM')';
end

% Converty new geometry into a 2D matrix
matNEWGEOM = (reshape((permute(matNEWGEOM,[2,1,3])),[3,numel(matNEWGEOM)/3]))';

% Add the chord values and angles to new geometry
matNEWGEOM(:,4) = repmat((reshape((matGEOM(:,4,:)),[1,0.2*numel(matGEOM)]))',valNUMB,1);
matNEWGEOM(:,5) = repmat((reshape((matGEOM(:,5,:)),[1,0.2*numel(matGEOM)]))',valNUMB,1);

% Return matrix into original form
reshape(matNEWGEOM',[2,5,0.1*numel(matNEWGEOM)])
temp = reshape(reshape(matNEWGEOM',[1,numel(matNEWGEOM)]),[1,10,0.1*numel(matNEWGEOM)]);
matNEWGEOM = permute(reshape(permute(temp,[2,1,3]),[5,2,0.1*numel(matNEWGEOM)]),[2,1,3]);

% Increase number of panels
valPANELS = valPANELS*valNUMB;

% Duplicated symetry and chord/spanwise element numbers for each blade
vecSYM = repmat(vecSYM,valNUMB,1);
vecN = repmat(vecN,valNUMB,1);
vecM = repmat(vecM,valNUMB,1);

% Plot to ensure rotation is as expected
figure(1)
clf(1)
hold on
for i = 1:0.1*numel(matNEWGEOM)
    scatter3(matNEWGEOM(:,1,i),matNEWGEOM(:,2,i),matNEWGEOM(:,3,i))
end
axis equal
xlabel('X-Dir','FontSize',15);
ylabel('Y-Dir','FontSize',15);
zlabel('Z-Dir','FontSize',15);
box on
grid on
axis equal
axis tight
hold off


end

