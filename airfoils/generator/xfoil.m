function [pol, foil] = xfoil(coord, alphaRange, ReRange, pause_length, transition, output_foil)
% Created by A.Y. initially, modified by T.D.K for optimization
% Modified filenames to be unique each run so it could be parallelized
% T.D.K. 2018-01-05

% Transition is a 1x2 of US and LS transition points (x/c)
% Output foil is true or false

upper_transition = transition(1);
lower_transition = transition(2);

foil = [];

%% Write airfoil coordinates to file

file_coord= [tempname('aux_files') '.foil'];
foil_name = file_coord(strfind(file_coord,'\')+1:strfind(file_coord,'.')-1);

if exist(file_coord,'file'),  delete(file_coord); end
fid = fopen(file_coord,'w');
if (fid <= 0)
    error([mfilename ':io'],'Unable to create file %s',file_coord);
else
    fprintf(fid,'%s\n',foil_name);
    fprintf(fid,'%9.5f   %9.5f\n',coord');
    fclose(fid);
end

%%
filename = file_coord; % add .dat

exeName = ['xfoil_', foil_name, '.exe'];
exePath = ['aux_files\', exeName];
copyfile('xfoil.exe', exePath);

len_alpha = length(alphaRange);
len_Re = length(ReRange);

pol = nan(len_alpha, 9, len_Re);

for j = 1:len_Re
    
    Re = ReRange(j);
    
    temp_command = [tempname('aux_files'), '.txt'];
    fileID = fopen(temp_command, 'w');
    temp_polar_out = [tempname('aux_files'), '.dat'];
    temp_output = [tempname('aux_files'), '.txt'];
    
    fprintf(fileID,'plop\n');
    fprintf(fileID,'g\n');
    fprintf(fileID,'\n');
    fprintf(fileID,'load\n');
    fprintf(fileID,'%s\n', filename);
    fprintf(fileID,'oper\n');
    fprintf(fileID,'iter\n');
    fprintf(fileID,'15\n');
    fprintf(fileID,'visc\n');
    fprintf(fileID,'%i\n', Re);
    
    % Forcing transition early on (% of chord on upper and lower)
    fprintf(fileID,'vpar\n');
    fprintf(fileID,'xtr\n');
    fprintf(fileID,'%0.2f\n', upper_transition);
    fprintf(fileID,'%0.2f\n', lower_transition);
    fprintf(fileID,'\n');
    
    fprintf(fileID,'pacc\n');
    fprintf(fileID,'%s\n', temp_polar_out);
    fprintf(fileID,'\n');
    
    [file_dump, file_cpwr] = deal(cell(1,length(alphaRange))); % Preallocate cell arrays
    for ii = 1:length(alphaRange)
        alpha = alphaRange(ii);
        %         if alpha < -5 || alpha > 10
        %             fprintf(fileID,'init\n');
        %         end
        fprintf(fileID,'alfa %0.3f\n', alpha);
        
        if output_foil == 1
            file_dump{ii} = sprintf('%s_a%06.3f_dump.dat',foil_name,alpha);
            file_cpwr{ii} = sprintf('%s_a%06.3f_cpwr.dat',foil_name,alpha);
            fprintf(fileID,'dump %s\n',file_dump{ii});
            fprintf(fileID,'cpwr %s\n',file_cpwr{ii});
        end
    end
    
    % Polar output filename
    if output_foil
        file_pwrt = sprintf('aux_files\%s_pwrt.dat',foil_name);
        fprintf(fileID,'pwrt\n%s\n',file_pwrt);
        fprintf(fileID,'plis\n');
    end
    fprintf(fileID,'pacc\n');
    fprintf(fileID,'\nquit\n');
    fclose(fileID);
    
    %     [~,~] = system(sprintf('%s <%s> %s', exePath, temp_command, temp_output));
    %     pause(0.5);
    [~,~] = system(sprintf(['start /b %s <%s> %s'], exePath, temp_command, temp_output));
    pause(pause_length);
    [~,~] = system(sprintf('taskkill /IM %s /F', exeName));
    
    if output_foil == 1
        jj = 0;
        
        for ii = 1:length(alphaRange)
            jj = jj + 1;
            fid = fopen([file_dump{ii}],'r');
            if (fid<=0)
                error([mfilename ':io'],'Unable to read xfoil output file %s',file_dump{ii});
            else
                D = textscan(fid,'%f%f%f%f%f%f%f%f','Delimiter',' ','MultipleDelimsAsOne',true,'CollectOutput',1,'HeaderLines',1);
                fclose(fid);
                delete([file_dump{ii}]);
                
                if ii == 1 % Use first run to determine number of panels (so that NACA airfoils work without vector input)
                    Npanel = length(D{1}); % Number of airfoil panels pulled from the first angle tested
                    % Preallocate Outputs
                    [foil.s, foil.x, foil.y, foil.UeVinf, foil.Dstar, foil.Theta, foil.Cf, foil.H] = deal(zeros(Npanel,length(alphaRange)));
                end
                
                %             % store data
                %             if ((jj>1) && (size(D{1},1)~=length(foil(ind).x)) && sum(abs(foil(ind).x(:,1)-size(D{1},1)))>1e-6 ),
                %                 ind = ind + 1;
                %                 jj = 1;
                %             end;
                foil.s(:,jj) = D{1}(:,1);
                foil.x(:,jj) = D{1}(:,2);
                foil.y(:,jj) = D{1}(:,3);
                foil.UeVinf(:,jj) = D{1}(:,4);
                foil.Dstar(:,jj) = D{1}(:,5);
                foil.Theta(:,jj) = D{1}(:,6);
                foil.Cf(:,jj) = D{1}(:,7);
                foil.H (:,jj)= D{1}(:,8);
            end
            
            foil.alpha(1,jj) = alphaRange(ii);
            
            % Read cp file
            fid = fopen([file_cpwr{ii}],'r');
            if (fid<=0)
                error([mfilename ':io'],'Unable to read xfoil output file %s',file_cpwr{ii});
            else
                C = textscan(fid, '%10f%9f%f', 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,3, 'ReturnOnError', false);
                fclose(fid);
                delete([file_cpwr{ii}]);
                % store data
                if ii == 1 % Use first run to determine number of panels (so that NACA airfoils work without vector input)
                    NCp = length(C{1}); % Number of points Cp is listed for pulled from the first angle tested
                    % Preallocate Outputs
                    [foil.xcp, foil.cp] = deal(zeros(NCp,length(alphaRange)));
                    foil.xcp = C{1}(:,1);
                end
                foil.cp(:,jj) = C{3}(:,1);
            end
        end
        
    end
    %% Polars
    fileID = fopen(temp_polar_out,'r');
    formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';
    dataArray = textscan(fileID, formatSpec, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'HeaderLines' ,12, 'ReturnOnError', false);
    fclose(fileID);
    
    len = size([dataArray{:,1:7}],1);
    pol(1:len, 1:8, j) = [[dataArray{:,1:7}] repmat(Re, len, 1)];
    pol(1:len, 9, j) = (pol(1:len, 2, j).^3)./(pol(1:len, 3, j).^2);
    
    %     matArray = [dataArray{:,1:7}];
    %     heading = {'alpha','CL','CD','CDp','CM','Top_Xtr','Bot_Xtr'};
    %
    %     for n = 1:length(matArray(:,1))
    %         for i = 1:length(matArray(1,:))
    %             Polars.(sprintf('Re_%i',Re'))(n).(heading{i}) = matArray(n,i);
    %             Polars.(sprintf('Re_%i',Re'))(n).('Re') = Re;
    %         end
    %     end
    %
    %
    %     if length(ReRange) == 1
    %         Polars = Polars.(sprintf('Re_%i',Re'));
    %     end
    
    delete(temp_polar_out,temp_command,temp_output);
    
end

delete(file_coord);
delete(exePath);

end