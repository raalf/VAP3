

numPanel = Aircraft.General.Panels;

Edges = zeros(numPanel*2,7);
Symmetry = zeros(numPanel,1);

%
for n = 1:numPanel %loop over number of panels
    Panel = Aircraft.Surface(n).Panel;
    %                   panel root/tip         X          Y   Z   Chord   Twist
    %                         (1 or 2)
    Edges((n*2)-1,:) = [n,       1, Panel.X(1),Panel.Y(1),Panel.Z(1),Panel.Chord(1),Panel.Twist(1)];
    Edges((n*2),:)   = [n,       2, Panel.X(2),Panel.Y(2),Panel.Z(2),Panel.Chord(2),Panel.Twist(2)];
    
    if isempty(Aircraft.Surface(n).Sym) == 1
        Symmetry(n) = 0;
    else
        Symmetry(n) = Aircraft.Surface(n).Sym;
    end
end
%

% Take out edge with symmetry condition
% Edge index if sym is true
%   eg            [1,2,3,4,5]         -1   +[1 2 0 0 0].* [1 1 0 0 0]
Sym_edge_idx = ((1:length(Symmetry))'-1+Symmetry).*(Symmetry>0);
% Sym_edge_idx holds info about which panel has symmetry
% and Symmetry holds the info about which edge of the panel has symmetry

% Reduce the variable to non-zeros only
Sym_edge_idx = Sym_edge_idx(Symmetry>0);
Symmetry = Symmetry(Symmetry>0);

% Overwrite symmetry edge with NaNs
for n = 1:length(Sym_edge_idx)
    idx = find(Edges(:,1)==Sym_edge_idx(n) & Edges(:,2)==Symmetry(n) == 1);
    Edges(idx,:) = deal(NaN);
end

%
joint_num = 1;


for n = 1:length(Edges(:,1)) % loop number of all edges

    if isnan(Edges(n,1)) == 0 % Skip Line with Nan, they have already been processed
        
        FindEdge = Edges(n,3:7); % edge we want to find duplicates of
        idx = ismember(Edges(:,3:7),FindEdge,'rows'); % find duplicates and store their indexes
        result = Edges(idx,:);  % result matrix with all the edges touching the same joint
        % Save result to Joint
        Joint(joint_num).Panel = result(:,1);
        Joint(joint_num).Edge = result(:,2);
        
        %Overwrite current and result edges on the same joint with NaN
        % so they won't get repeated
        Edges(n,:) = deal(NaN);
        Edges(idx,:) = deal(NaN);
        
        joint_num = joint_num+1;
    end
    
    
end



% Output Joint to FW
FW.Joint = Joint;