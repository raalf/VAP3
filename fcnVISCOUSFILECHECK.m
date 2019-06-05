function [ FLAG ] = fcnVISCOUSFILECHECK(FLAG, VISC)
%FCNVISCOUSFILECHECK Summary of this function goes here
%   Detailed explanation goes here


if FLAG.VISCOUS == 1
    
    % case 1, no airfoil file specified
    if isempty(VISC.cellAIRFOIL) == 1 %|| all(cellfun(@isnan, VISC.cellAIRFOIL,'UniformOutput',false))
        % no airfoil tag found in the input file
        disp('FLAG.VISCOUS = 0: No airfoil tag found in the input file.')
        FLAG.VISCOUS = 0;

    else
        [uniqueAirfoil,~,~] = unique(VISC.cellAIRFOIL);
        lenAirfoil = length(uniqueAirfoil);
        
        % get file list within airfoils folder
        filelist = dir('airfoils');
        % use regular expression to filter out only the .mat files
        idxDotMat = ~cellfun('isempty',regexp({filelist.name}','.mat$'));
        dotMatFiles = {filelist(idxDotMat).name}';
        
        % add '.mat' after each airfoil name from input,
        % compare with the airfoil .mat files in the folder
        idxMatch = ismember(strcat(uniqueAirfoil,'.mat'),dotMatFiles);
        
        if sum(idxMatch) ~= lenAirfoil
            % not all airfoil found
            FLAG.VISCOUS = 0;
            disp('FLAG.VISCOUS = 0: Airfoil data file missing.');
            disp(VISC.cellAIRFOIL(~idxMatch));
        end
            
        
    end
    
    
else % if FLAG.VISCOUS is anything other than 1,
    % set flag to zero
    FLAG.VISCOUS = 0;
    disp('FLAG.VISCOUS = 0: Unknown viscous flag input.')
end




end

