% Symmetry condition was handled incorrectly, casuing non-sym edge being
% overwritted. This problem is now fixed. (Jan18,2016) -Alton
% VAP2.0 - Rewrite by Alton (Sep 2016)

% INPUT 
% valPANELS - loop counter
% vecSYM
clc

%% create Edges matrix, this holds the information of all the edges.
% Next part of the code will determine which edges are touching each other
% This is equvalent to VAP1.0 panelBoundary line 6-24
edges = [reshape(repmat(1:valPANELS,2,1),[],1),...
    repmat([1;2],3,1),...
    reshape(permute(matGEOM,[1,3,2]),6,[],1)];

%% Take out edge with symmetry condition
% Edge index if sym is true
% This is equvalent to VAP1.0 panelBoundary line 27-36
symEdgeIdx = find(vecSYM)*2+vecSYM(vecSYM~=0)-2;% This point out the row number of symm edge

% Delete symmetry edge with NaNs
edges = edges(~ismember(1:6,symEdgeIdx),:);

%% This is equvalent to VAP1.0 panelBoundary line 44-68

joint_num = 1;

for n = 1:length(edges(:,1)) % loop number of all edges

    if isnan(edges(n,1)) == 0 % Skip Line with Nan, they have already been processed
        
        FindEdge = edges(n,3:7); % edge we want to find duplicates of
        idx = ismember(edges(:,3:7),FindEdge,'rows'); % find duplicates and store their indexes
        result = edges(idx,:);  % result matrix with all the edges touching the same joint
        % Save result to Joint
        cellJointPanel{joint_num,1} = result(:,1);
        cellJointEdge{joint_num,1} = result(:,2);
        
        %Overwrite current and result edges on the same joint with NaN
        % so they won't get repeated
        edges(n,:) = deal(NaN);
        edges(idx,:) = deal(NaN);
        
        joint_num = joint_num+1;
    end   
end


%%  Detect which panels are connected to form a wing
%   Results are output in FW.Panels().Wing
%   Update on Jan19, 2016 -Alton

PanelCount = length(cellJointPanel);
LinkedPanel = zeros(PanelCount,PanelCount);
for n = 1:PanelCount;
    P = cellJointPanel{n};
    [C,ia,ib] = intersect(LinkedPanel,P);
    % First Panel
    if isempty(C) == 1
        wing = sum(LinkedPanel(:,1)~=0) + 1; % Count Existing wing in matrix
        LinkedPanel(wing,1:length(P)) = P;
    else
        new = P(P~=C);  %New connected Panel
        begin = sum(LinkedPanel(ib,:)~=0); %Begin writing panel number
        LinkedPanel(ib,begin+1:begin+length(new)) = new;
    end   
end

for n = 1:PanelCount;
     % By finding the position of the panel index in LinkedPanel
     % This will return index of the panel. 
     % mod funtion will return it's row number aka wing index
    wing = mod(find(LinkedPanel==n),PanelCount);
    
    if wing == 0    % Avoid Wing = 0
        wing = wing+PanelCount;
    end
    
wing
end


































