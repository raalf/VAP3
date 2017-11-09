function [avg_value] = fcnTIMEAVERAGE(value, valROTORRPM, valDELTIME)

timestep_rpm = ceil(1/((valROTORRPM/60)*valDELTIME));

try
   avg_value = mean(value(end - timestep_rpm:end,:,:),1);
catch
   avg_value =  mean(value,1);
end

end

