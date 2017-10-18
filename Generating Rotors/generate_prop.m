clc
clear

load('NASA_MIL_Propeller.mat');

sections = 10;
root = 2;
tip = 11.5;
section_loc = root:(tip-root)/sections:tip;

[twist, ~] = fcnCREATEFIT(Beta_vs_r(:,1), Beta_vs_r(:,2));
[le, ~] = fcnCREATEFIT(LE_vs_r(:,1), LE_vs_r(:,2));
[te, ~] = fcnCREATEFIT(TE_vs_r(:,1), TE_vs_r(:,2));

twist = twist(section_loc);
chord = (le(section_loc) - te(section_loc)).*0.0254;
le = [-le(section_loc) section_loc' zeros(length(section_loc),1)].*0.0254;

%%
% chord = chord.*2.5;
% le = le.*2.5;

symmetry = 0;
N = 1;

filename = 'inputs/NASA_MIL_Propeller_Large.vapgeom';
fp = fopen(filename,'w');

fprintf(fp,'<panel>\n');
fprintf(fp,'\t<symmetry>%d</symmetry>\n',symmetry);
fprintf(fp,'\t<N>%d</N>\n',N);

for i = 1:sections
    fprintf(fp,'\t<section>\n');
    fprintf(fp,'\t\t<x>%f</x>\n',le(i,1));
    fprintf(fp,'\t\t<y>%f</y>\n',le(i,2));
    fprintf(fp,'\t\t<z>%f</z>\n',le(i,3));
    fprintf(fp,'\t\t<chord>%f</chord>\n',chord(i));
    fprintf(fp,'\t\t<twist>%f</twist>\n',twist(i));
    fprintf(fp,'\t</section>\n');
end

fprintf(fp,'</panel>\n');

fclose(fp);