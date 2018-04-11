function str = fcnQMIL2VAP(filename, n, airfoil)

% Skipping everything except the geometry at the bottom
fp = fopen(filename);
start = 1;
ch = fgetl(fp);
while(isempty(ch) || ch(1)~='#')
ch = fgetl(fp);
start = start + 1;
end
fclose(fp);

% Reading in the geometry
qmil_geom = dlmread(filename, '', start, 0);

section_y = linspace(qmil_geom(1,1), qmil_geom(end,1), n + 1)';
section_z = section_y.*0;
section_chord = interp1(qmil_geom(:,1), qmil_geom(:,2), section_y);
section_twist = interp1(qmil_geom(:,1), qmil_geom(:,3), section_y);
section_x = -section_chord./4;

% Building the string for the VAP3 XML geometry
str = sprintf('<panel>\n\t<airfoil>%s</airfoil>\n\t<symmetry>0</symmetry>\n\t<N>1</N>\n', airfoil);
for i = 1:(n + 1)
    str_tmp = sprintf('\t<section>\n\t\t<x>%.6f</x>\n\t\t<y>%.6f</y>\n\t\t<z>%.6f</z>\n\t\t<chord>%0.6f</chord>\n\t\t<twist>%.6f</twist>\n\t</section>\n',section_x(i), section_y(i), section_z(i), section_chord(i), section_twist(i));
    str = [str, str_tmp];
end
str = [str, sprintf('</panel>\n')];

end