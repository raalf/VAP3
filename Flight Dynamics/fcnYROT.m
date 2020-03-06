function [rotfpg, fpg2] = fcnYROT(deltaEPS,fpg,rotpt)
% Function to perform rotation of points about y-axis and specified point

% Rotate points about rotpt based on angle deltaEPS
tempfpg = fpg - rotpt;

ROT = [cos(deltaEPS) 0 sin(deltaEPS); 0 1 0; -sin(deltaEPS) 0 cos(deltaEPS)];

fpg2 = (ROT*tempfpg')' + rotpt; % New points

rotfpg = fpg2 - fpg; % Change in position of points from original location

end