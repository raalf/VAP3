function [Aircraft, FW] = generateSingfct(Aircraft, FW, ii)

n = FW.Panels(ii).n;
% m = FW.Panels(ii).m;
m = FW.m;

%% singularity factor for the wing elements
%the singfct is set equal for all wake elements of a specific wing
%singfct is 1% of the tip element half-span(eta)

%ii is the panel number we are on

%now we need the tip panel for this wing
for k = 1:size(FW.Joint,2)
    if size(FW.Joint(k).Panel) == 1
        tippanel = FW.Joint(k).Panel;
        %does this match the current panel's wing?
        if FW.Panels(ii).Wing == FW.Panels(tippanel).Wing
            break
        end
        
    end
end

tempS = 0.01*FW.Panels(tippanel).DVE.eta(end);
FW.Panels(ii).DVE.singfct(1:n*m) = tempS;
clear i indx tempS


end

