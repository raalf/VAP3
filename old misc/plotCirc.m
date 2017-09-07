function [ ] = plotCirc(FW, Aircraft, Output, Temp, style, leg)

hFig2 = figure(2);
% clf(2);
count = 1;

hold on
if FW.m == 1
    for i = 1:Aircraft.General.Panels
        for j = 1:length(Output.(Temp.AI)(end).TimestepData(i).DVE.Index)
            
            eta_range = [-FW.Panels(i).DVE.eta(j):FW.Panels(i).DVE.eta(j)/500:FW.Panels(i).DVE.eta(j)];
            eta_range_plot = FW.Panels(i).DVE.LECoordL(j,2):norm(FW.Panels(i).DVE.LECoordL(j,2)-FW.Panels(i).DVE.LECoordR(j,2))/1000:FW.Panels(i).DVE.LECoordR(j,2);
            
            circ = Output.(Temp.AI)(end).TimestepData(i).DVE.A(j) + Output.(Temp.AI)(end).TimestepData(i).DVE.B(j).*eta_range + Output.(Temp.AI)(end).TimestepData(i).DVE.C(j).*(eta_range.^2);
            %             vort = FW.Panels(i).DVE.B(j) + 2.*FW.Panels(i).DVE.C(j).*(eta_range);
            plot(eta_range_plot, circ, style);
            if isempty(leg) == 0
                legappend(leg)
            end
            grid on
            
            count = count + 1;
            
        end
    end
end
hold off
axis tight

end

