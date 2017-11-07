clear
clc

AirfoilName = 'MAXWELL';
load('x57_airfoil.mat')

Polars = xfoil(coord, -2:0.5:15, [150000:100000:3e6], 10);

fn = fieldnames(Polars);
output = [];
for n = 1:length(fn)
    Polar = Polars.(fn{n});
    output_temp = [[Polar.alpha]',[Polar.CL]',[Polar.CD]',[Polar.Re]',[Polar.CM]'];
    [~,idx] = sort([Polar.CL]');
    output_temp = output_temp(idx,:);
    output = [output;output_temp];
end

%
fileID = fopen('freewake_airfoil.dat','w');

fprintf(fileID, '%s,        free transition row= %i\n', AirfoilName, length(output(:,1)));
for n = 1:length(output(:,1));
    fprintf(fileID,'%.4f %.4f %.4f %i %.4f \n',...
        output(n,1),output(n,2),output(n,3),output(n,4),output(n,5));
end
fclose(fileID);
