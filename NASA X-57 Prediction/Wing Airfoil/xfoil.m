function [Polars, foil] = xfoil(coord, alphaRange, ReRange, pause_length)

foil = [];
% Put to 1 to output foil (cp curves etc)
output_foil = 0;

%% Write airfoil coordinates to file
wd = 'aux_files/';
delete([wd 'XFOIL_Polar_Output.dat']);
delete([wd 'output.txt']);

foil_name = strcat(wd,mfilename);
file_coord= [foil_name '.foil'];

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
exePath = 'xfoil.exe';

len_alpha = length(alphaRange);
len_Re = length(ReRange);

pol = nan(len_alpha, 9, len_Re);

for j = 1:len_Re
    delete([wd 'XFOIL_Polar_Output.dat']);
    delete([wd 'output.txt']);

    Re = ReRange(j);
    
    fileID = fopen([wd 'XFOIL_Commands.txt'],'w');
    fprintf(fileID,'plop\n');
    fprintf(fileID,'g\n');
    fprintf(fileID,'\n');
    fprintf(fileID,'load\n');
    fprintf(fileID,'%s\n', filename);
    fprintf(fileID,'oper\n');
    %     fprintf(fileID,'iter\n');
    %     fprintf(fileID,'10\n');
    fprintf(fileID,'visc\n');
    fprintf(fileID,'%i\n', Re);
    fprintf(fileID,'pacc\n');
    fprintf(fileID,'%s\n', [wd 'XFOIL_Polar_Output.dat']);
    fprintf(fileID,'\n');
    
    [file_dump, file_cpwr] = deal(cell(1,length(alphaRange))); % Preallocate cell arrays
    for ii = 1:length(alphaRange)
        alpha = alphaRange(ii);
        if alpha < -5 || alpha > 1
            fprintf(fileID,'init\n');
        end
        fprintf(fileID,'alfa %0.3f\n', alpha);
        
        if output_foil == 1
            file_dump{ii} = sprintf('%s_a%06.3f_dump.dat',foil_name,alpha);
            file_cpwr{ii} = sprintf('%s_a%06.3f_cpwr.dat',foil_name,alpha);
            fprintf(fileID,'dump %s\n',file_dump{ii});
            fprintf(fileID,'cpwr %s\n',file_cpwr{ii});
        end
    end
    
    % Polar output filename
    if output_foil == 1
        file_pwrt = sprintf('%s_pwrt.dat',foil_name);
        fprintf(fileID,'pwrt\n%s\n',file_pwrt);
        fprintf(fileID,'plis\n');
    end
    fprintf(fileID,'pacc\n');
    fprintf(fileID,'\nquit\n');
    fclose(fileID);
    
%     [~,~] = system(sprintf(['%s <' wd 'XFOIL_Commands.txt> aux_files/output.txt'], exePath));
%     pause(3); % The & at the end of the command means MATLAB doesn't wait for XFOIL to run, so we need to

    [~,~] = system(sprintf(['start /b %s <' wd 'XFOIL_Commands.txt> ' wd 'output.txt'], exePath));
    pause(pause_length);
    [~,~] = system(sprintf('taskkill /IM %s /F', exePath));

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
    fileID = fopen([wd 'XFOIL_Polar_Output.dat'],'r');
    formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';
    dataArray = textscan(fileID, formatSpec, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'HeaderLines' ,12, 'ReturnOnError', false);
    fclose(fileID);
    
%     len = size([dataArray{:,1:7}],1);
%     pol(1:len, 1:8, j) = [[dataArray{:,1:7}] repmat(Re, len, 1)];
%     pol(1:len, 9, j) = (pol(1:len, 2, j).^3)./(pol(1:len, 3, j).^2);
    
        matArray = [dataArray{:,1:7}];
        heading = {'alpha','CL','CD','CDp','CM','Top_Xtr','Bot_Xtr'};
    
        for n = 1:length(matArray(:,1))
            for i = 1:length(matArray(1,:))
                Polars.(sprintf('Re_%i',Re'))(n).(heading{i}) = matArray(n,i);
                Polars.(sprintf('Re_%i',Re'))(n).('Re') = Re;
            end
        end
    
    
        if length(ReRange) == 1
            Polars = Polars.(sprintf('Re_%i',Re'));
        end
       
end

end