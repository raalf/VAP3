function qmil_filename = fcnQMILCREATE(blades, thrust, vinf, rpm, diam)
cl0 = 0;
cla = 6.2832;
clmin = -0.8;
clmax = 1.2;
cd0 = 0.01;
cd2u = 0.008;
cd2l = 0.006;
clcd0 = 0.4;
reref = 150000;
reexp = -0.5;
xides = [0 0.5 1];
cldes = [0.6 0.5 0.4];
hub_radius = 0.15;

str_1 = sprintf('temp prop\n\n%d !Nblades\n\n%.4f %.4f ! CL0 CL_A\n%.4f %.4f ! CLmin CLmax\n\n%.4f %.4f %.4f %.4f ! CD0 CD2u CD2l CLCD0\n%.4f %.4f ! Reref REexp\n\n%.4f %.4f %.4f ! XIdes (r over R locations where design cl is specified)\n%.4f %.4f %.4f ! CLdes (specified cl)\n\n',...
    blades, cl0, cla, clmin, clmax, cd0, cd2u, cd2l, clcd0, reref, reexp, xides(1), xides(2), xides(3), cldes(1), cldes(2), cldes(3));

str_2 = sprintf('%.4f ! hub radius(m)\n%.4f ! tip radius (m)\n%.4f ! speed (m/s)\n%.4f ! rpm\n\n', hub_radius, diam/2, vinf, rpm);
str_3 = sprintf('%.4f ! Thrust(N)\n%.4f ! Power(W)\n\n%.4f %.4f ! Ldes KQdes\n\n%d ! Nout number of output stations (optional)', thrust, 0, 0, 0, 30);

str = [str_1, str_2, str_3];

qmil_filename = [tempname('aux_files'), '.prop'];

fp = fopen(qmil_filename, 'w');
fprintf(fp, str);
fclose(fp);

end