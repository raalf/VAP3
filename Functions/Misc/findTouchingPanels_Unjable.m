[Aircraft, FW] = generateDVEs(Aircraft, FW)

numPanel = Aircraft.General.Panels;

Edges = zeros(numPanel*2,7);

for n = 1:numPanel %go through each panel
    Panel = Aircraft.Surface(n).Panel;
    %                   panel root/tip         X          Y   Z   Chord   Twist
    %                         (1 or 2)
    Edges((n*2)-1,:) = [n,       1, Panel.X(1),Panel.Y(1),Panel.Z(1),Panel.Chord(1),Panel.Twist(1)];
    Edges((n*2),:)   = [n,       2, Panel.X(2),Panel.Y(2),Panel.Z(2),Panel.Chord(2),Panel.Twist(2)];
end


L = zeros(numPanel,1);
R = zeros(numPanel,1);

for n = 1:length(Edges(:,1))
    
    % Get Target Panel Number
    panel = Edges(n,1);
    
    % Isolate Edge we want to search for matches
    FindEdge = Edges(n,3:7);
    
    % Search the Edge from all edges
    idx = ismember(Edges(:,3:7),FindEdge,'rows');
    MatchedPanel = Edges(idx,1);
    % Remove the panel result if it's the same panel
    MatchedPanel = MatchedPanel(MatchedPanel~=panel);
    
    % Storing Reults in matrix L and matrix R
    if isempty(MatchedPanel) ~= 1
        A = length(MatchedPanel);
        if Edges(n,2) == 1
            L(panel,1:A) = MatchedPanel';
        elseif Edges(n,2) == 2
            R(panel,1:A) = MatchedPanel';
        end
    end
end




end