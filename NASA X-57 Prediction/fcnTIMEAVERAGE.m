function [avg_value] = fcnTIMEAVERAGE(value, valROTORRPM, valDELTIME)

timestep_rpm = ceil(1/((valROTORRPM/60)*valDELTIME));
avg_value = mean(value(end - timestep_rpm:end,:,:),1);

end

