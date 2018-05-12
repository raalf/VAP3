function [avg_value] = fcnTIMEAVERAGE(value, valROTORRPM, valDELTIME)
valROTORRPM = abs(valROTORRPM);
timestep_rpm = ceil(1/((valROTORRPM/60)*valDELTIME));

try
   avg_value = nanmean(value(end - timestep_rpm:end,:,:),1);
catch
   avg_value =  nanmean(value,1);
end

end

