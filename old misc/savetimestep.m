function [Output] = savetimestep(FW,Temp,Output)

%function to save timestep data to struct Output
%saves FW.Panels to output struct
% if Temp.timestep < 0
%     Output.(Temp.AI)(Temp.timestep+2).(strrep(strcat('timestep',num2str(Temp.timestep)),'-','N')) = FW.Panels; %change '-' to 'N'  
    Output.(Temp.AI)(Temp.timestep+2).TimestepData = FW.Panels; %change '-' to 'N'
    Output.(Temp.AI)(Temp.timestep+2).TimestepNumber = (strrep(strcat('timestep',num2str(Temp.timestep)),'-','N')); %change '-' to 'N' 
% else
%     Output.(Temp.AI)(Temp.timestep+2).TimestepNumber.(strcat('timestep',num2str(Temp.timestep))) = FW.Panels; 
% end
end