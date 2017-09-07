function [HB] = HBrunner(Aircraft, Conditions, FW, alpha)
%% This function will write an input file from ASI to Horstmann's 
% clc
% clear

if isdeployed == 0
    addpath('Functions');
end

% alpha = 5;

% load('Aircraft/simplewing.mat');
% load('Aircraft/StandardCirrus.mat');
% load('Aircraft/designOne.mat');
% load('Aircraft/workingsplit.mat');
% alpha = 6;

FW = panelBoundary(Aircraft,FW);
% Conditions.Alpha = alpha;

inp_string = 'test1';
fp = fopen(strcat(inp_string,'.inp'),'w');

%% Header stuff

fprintf(fp, ...
['|======================================================================================================================================|\n', ...
' LIFTING_LINE INPUTFILE\n', ...
'|--------------------------------------------------------------------------------------------------------------------------------------|\n', ...
' VERSION V2.5\n', ...
'|--------------------------------------------------------------------------------------------------------------------------------------|\n', ...
'|**************************************************************************************************************************************|\n', ...
'|--------------------------------------------------------------------------------------------------------------------------------------|\n', ...
'|                                                           NAME_OF_DATASET                                                            |\n', ...
' LIFTING_LINE V2.5, VAP: "SOMETHING"\n']);

%% First block
fprintf(fp, ...
['|--------------------------------------------------------------------------------------------------------------------------------------|\n', ...    
'|   SYMMETRY   | GEO_SCALING  |  TWIST_DISTR | QSI_STDY_ROT | N_PROPELLERS | N_ADD_LINES  |   XML_DATA   |                             |\n', ... 
'              %d  %f                  %d              %d              %d              %d              %d                              \n' 
], FW.Sym, 1, 0, 0, 0, 0, 0  ...
);

fprintf(fp, ...
['|--------------------------------------------------------------------------------------------------------------------------------------|\n', ...    
'|   REF_AREA   |   REF_SPAN   | REF_LEN_CMX  | REF_LEN_CMY  | REF_LEN_CMZ  |  MOM_REF_X   |  MOM_REF_Y   |  MOM_REF_Z   |              |\n', ... 
'   %f      %f       %f       %f       %f       %f       %f       %f                   \n' 
], Aircraft.Reference.S, Aircraft.Reference.b, 0, 0, 0, 0, 0, 0 ...
);

fprintf(fp, ...
['|--------------------------------------------------------------------------------------------------------------------------------------|\n', ...    
'|   N_WINGS    |\n', ... 
'              %d\n' 
], Aircraft.General.Wings ...
);

fprintf(fp, ...
['|--------------------------------------------------------------------------------------------------------------------------------------|\n', ...    
'|   N_PW_SPAN  |  N_PW_CHORD  |\n']);

for i = Aircraft.General.Wings:-1:1 % Going backwards so we match the geometry at the bottom of the input file
    [~, indx] = find([FW.Panels.Wing] == i);
%     fprintf(fp, '              %d              %d\n', i, length(indx)); 
    fprintf(fp, '              %d              %d\n', length(indx), 1);
end

fprintf(fp, ...
['|--------------------------------------------------------------------------------------------------------------------------------------|\n', ...
'|**************************************************************************************************************************************|\n', ...
'|--------------------------------------------------------------------------------------------------------------------------------------|\n', ...
'|   N_ALPHA    |    ALPHA     |    ALPHA     |    ALPHA     |    ALPHA     |    ALPHA     |    ALPHA     |    ALPHA     |      ...      \n', ...
'              %d '], length([Conditions.Alpha]));

if length([Conditions.Alpha]) > 7
   disp('A$AP ROCKY ft. Drake, 2 Chainz, Kendrick Lamar')
end

for i = 1:length([Conditions.Alpha])
   fprintf(fp,' %f      ', Conditions.Alpha(i));
end

fprintf(fp,['\n|--------------------------------------------------------------------------------------------------------------------------------------|\n', ...
'| N_TARGET_CFZ |  TARGET_CFZ  |  TARGET_CFZ  |  TARGET_CFZ  |  TARGET_CFZ  |  TARGET_CFZ  |  TARGET_CFZ  |  TARGET_CFZ  |      ...      \n', ...
'              0 \n', ...
'|--------------------------------------------------------------------------------------------------------------------------------------|\n', ...
'|  COMP_MA_NO  |     BETA     | \n', ...
'   %f       %f\n' ...
], 0, Conditions.Beta);

fprintf(fp, ...
['|--------------------------------------------------------------------------------------------------------------------------------------|\n', ...
'|**************************************************************************************************************************************|\n', ...
'|--------------------------------------------------------------------------------------------------------------------------------------|\n']);

%% Second Block
fprintf(fp, ...
['|   INDEX_PW   |   N_PANELS   |  PANEL_DISTR |FACT_CL_ALPHA | FUSELAGE_INF | ELLIPT_CHORD |                                            |\n', ...
'|    X_PW_1    |    Y_PW_1    |    Z_PW_1    |  CHORD_PW_1  |  TWIST_PW_1  | COUPL_COND_1 | FORKING_NO_1 |                             |\n', ...
'|    X_PW_2    |    Y_PW_2    |    Z_PW_2    |  CHORD_PW_2  |  TWIST_PW_2  | COUPL_COND_2 | FORKING_NO_2 |                             |\n', ...
'|--------------------------------------------------------------------------------------------------------------------------------------|\n' ...
]);

count = 1;
for i = Aircraft.General.Panels:-1:1
                % Panel number,     N_spanwise,     Panel_Distr,    CL alpha corr,  Fuselage,   Ellipt_chord
    fprintf(fp, '             %d             %d              %d  %f                  %d              %d\n', count, 10, 0, 1, 0, 0);
    
    
    % finding BC of outer edge (Edge 2 in VAP usually but edge 1 in Horstmann)
    if Aircraft.Surface(i).Sym == 2
        BC2 = 'F T F F F F'; % symmetry plane
        fork = 0;
    else
        for ii = 1:length(FW.Joint)
            [indx,~] = find(FW.Joint(ii).Panel == i);
            
            if isempty(indx) == 0
                if FW.Joint(ii).Edge(indx) == 2 && length(FW.Joint(ii).Panel) == 1
                    BC2 = 'T F F F F F'; % free end
                    fork = 0;
                    break;
                elseif FW.Joint(ii).Edge(indx) == 2 && length(FW.Joint(ii).Panel) == 2
                    BC2 = 'F F T T F F'; % 220 condition
                    fork = 0;
                    break;
                elseif FW.Joint(ii).Edge(indx) == 2 && length(FW.Joint(ii).Panel) == 3
                    BC2 = 'F F T F T F'; % split 
                    fork = 1;
                    break;
                end
            end
        end
    end 

    fprintf(fp, '  %f       %f       %f      %f       %f        %s              %d\n', ... 
        (Aircraft.Surface(i).Panel.X(2)+(Aircraft.Surface(i).Panel.Chord(2)*0.25)), Aircraft.Surface(i).Panel.Y(2), Aircraft.Surface(i).Panel.Z(2), Aircraft.Surface(i).Panel.Chord(2), Aircraft.Surface(i).Panel.Twist(2), BC2, fork);
    
    clear fork 
    % finding BC of inner edge (Edge 1 in VAP usually but edge 2 in Horstmann)
    if Aircraft.Surface(i).Sym == 1
        BC1 = 'F T F F F F'; % symmetry plane
        fork = 0;
    else
        for ii = 1:length(FW.Joint)
            [indx,~] = find(FW.Joint(ii).Panel == i); % SOMETHING WRONG HEREEEEEE
            
            if isempty(indx) == 0
                if FW.Joint(ii).Edge(indx) == 1 && length(FW.Joint(ii).Panel) == 1
                    BC1 = 'T F F F F F'; % free end
                    fork = 0;
                    break;
                elseif FW.Joint(ii).Edge(indx) == 1 && length(FW.Joint(ii).Panel) == 2
                    BC1 = 'F F T T F F'; % 220 condition
                    fork = 0;
                    break;
                elseif FW.Joint(ii).Edge(indx) == 1 && length(FW.Joint(ii).Panel) == 3
                    BC1 = 'F F T F T F'; % split  
                    fork = 1;
                    break;
                end
            end
        end
    end    

    fprintf(fp, '  %f       %f       %f      %f       %f        %s              %d\n', ... 
   (Aircraft.Surface(i).Panel.X(1)+(Aircraft.Surface(i).Panel.Chord(1)*0.25)), Aircraft.Surface(i).Panel.Y(1), Aircraft.Surface(i).Panel.Z(1), Aircraft.Surface(i).Panel.Chord(1), Aircraft.Surface(i).Panel.Twist(1), BC1, fork);
    
    fprintf(fp, '|--------------------------------------------------------------------------------------------------------------------------------------|\n');   
    count = count + 1;
end 

%% End
fprintf(fp, '|=====================================================END OF LIFTING_LINE INPUTFILE====================================================|\n');
fclose(fp);

%% Run

[~, ~] = system(sprintf('LIFTING_LINE_WINDOWS_64BIT.exe -om:overwrite -pj:%s.inp',inp_string));
% system(sprintf('LIFTING_LINE_WINDOWS_64BIT.exe -om:overwrite -pj:%s.inp',inp_string))

%% Results

fp = fopen(strcat(inp_string,'.lili.V2.5\results\',inp_string,'.out'),'r');

% Lift coefficient
ch = fscanf(fp,'%s', 1);
while strcmp(ch, 'GESAMTAUFTRIEBSBEIWERT') == 0 && feof(fp) == 0
    ch = fscanf(fp,'%s', 1);
end

while strcmp(ch, '=') == 0 && feof(fp) == 0
    ch = fscanf(fp,'%c',1);
end

HB.CL = fscanf(fp,'%f');

% Induced drag coefficient
ch = fscanf(fp,'%s', 1);
while strcmp(ch, 'CWI') == 0 && feof(fp) == 0
    ch = fscanf(fp,'%s', 1);
end

while strcmp(ch, '=') == 0 && feof(fp) == 0
    ch = fscanf(fp,'%c',1);
end

HB.CDi = fscanf(fp,'%f');

HB.e = (HB.CL^2)/(pi*Aircraft.Reference.AR*HB.CDi);

fclose(fp);
clear inp_string indx ii i fp fork count ch BC1 BC2 ans 














