clc
clear

cd C:\Users\travi\OneDrive\Desktop\GIT\VAP3\Runs\J_COLE_OPTIMIZATION\Analysis

z(1,:) = [76 74 72 68 66 65 63 59 57 56 53 155 2250 -2 482 0 1 8 700 0 1 18 900 0 1 27 1100 0 1 35 1300 0 1 51 1600 0 1]; % 76907.3
z(2,:) = [76 74 72 70 66 64 62 62 57 55 52 142 2209 -13 249 0 1 -3 476 1 0 18 888 1 1 35 1309 -2 1 52 1607 1 1 61 1801 1 1]; % 72133.5
z(3,:) = [75 74 72 68 66 64 61 59 57 53 53 110 2100 -15 200 -1 1 -8 351 0 0 -2 483 0 1 16 900 -0 1 33 1300 1 1 52 1601 -0 1]; % 73381.7

rotors = [1, 2, 3];

legend_entry = {'Baseline Design', 'Design 1', 'Design 2', 'Design 3'};

linestyles = {'--';'-.';'-';':'};
markers = {'o';'x';'s';'^';'*';'d';'v';'>';'<';'p';'h'};
colors = {'k';'b';'r';'m';'c';'g'};

cd ./../
for i = 1:size(z,1) + 1
    if i == size(z,1) + 1
        [out(1,i), Design(i).ITER, Design(i).ITEROUTP] = fcnBASELINE_OBJ();
    else
        [out(1,i), Design(i).ITER, Design(i).ITEROUTP] = fcnOBJECTIVE(z(i,:), 11, 6, 4);
    end
end
cd Analysis/
% save('matlab.mat');

load('matlab.mat');

%% Drag Bar Graph
x = categorical({'Total Drag', 'Induced Drag', 'Profile Drag'});

baseline = [Design(end).ITEROUTP(end).OUTP.vecCD_AVG; ...
    Design(end).ITEROUTP(end).OUTP.vecCDI_AVG; ...
    Design(end).ITEROUTP(end).OUTP.vecCDP_AVG];

designs = [Design(1).ITEROUTP(end).OUTP.vecCD_AVG, Design(2).ITEROUTP(end).OUTP.vecCD_AVG, Design(3).ITEROUTP(end).OUTP.vecCD_AVG; ...
    Design(1).ITEROUTP(end).OUTP.vecCDI_AVG, Design(2).ITEROUTP(end).OUTP.vecCDI_AVG, Design(3).ITEROUTP(end).OUTP.vecCDI_AVG;...
    Design(1).ITEROUTP(end).OUTP.vecCDP_AVG, Design(2).ITEROUTP(end).OUTP.vecCDP_AVG, Design(3).ITEROUTP(end).OUTP.vecCDP_AVG];

hFig200 = figure(200);
clf(200);

b = bar( x, 100.*(designs - repmat(baseline,[1,3]))./repmat(baseline,[1,3]) );
b(1).LineStyle = '--';
b(1).LineWidth = 1;
b(2).LineStyle = '-.';
b(2).LineWidth = 1;
b(3).LineStyle = ':';
b(3).LineWidth = 1;
grid minor
box on
axis tight

ylabel('Percent Change from Baseline', 'FontSize', 15)
legend('Design 1', 'Design 2', 'Design 3')

%% Power Bar Graph
x = categorical({'Total Power', 'Induced Power', 'Profile Power'});

baseline = [Design(end).ITEROUTP(end).OUTP.vecCP_AVG; ...
    Design(end).ITEROUTP(end).OUTP.vecCPI_AVG; ...
    Design(end).ITEROUTP(end).OUTP.vecCPP_AVG];

designs = [mean(Design(1).ITEROUTP(end).OUTP.vecCP_AVG,2), mean(Design(2).ITEROUTP(end).OUTP.vecCP_AVG,2), mean(Design(3).ITEROUTP(end).OUTP.vecCP_AVG,2); ...
    mean(Design(1).ITEROUTP(end).OUTP.vecCPI_AVG,2), mean(Design(2).ITEROUTP(end).OUTP.vecCPI_AVG,2), mean(Design(3).ITEROUTP(end).OUTP.vecCPI_AVG,2);...
    mean(Design(1).ITEROUTP(end).OUTP.vecCPP_AVG,2), mean(Design(2).ITEROUTP(end).OUTP.vecCPP_AVG,2), mean(Design(3).ITEROUTP(end).OUTP.vecCPP_AVG,2)];

hFig201 = figure(201);
clf(201);

b = bar( x, 100.*(designs - repmat(baseline,[1,3]))./repmat(baseline,[1,3]) );
b(1).LineStyle = '--';
b(1).LineWidth = 1;
b(2).LineStyle = '-.';
b(2).LineWidth = 1;
b(3).LineStyle = ':';
b(3).LineWidth = 1;
grid minor
box on
axis tight

ylabel('Percent Change from Baseline', 'FontSize', 15)
legend('Design 1', 'Design 2', 'Design 3')

%% Prop Thrust Distribution
hFig202 = figure(202);
clf(202);
loc = 2.*(Design(end).ITEROUTP(end).OUTP.ROTOR.vecSPANLOC)./Design(end).ITEROUTP(end).OUTP.vecROTDIAM(1);
p = plot(loc, mean(Design(end).ITEROUTP(end).OUTP.ROTOR.vecTHRUSTDIST_AVG,3));
p.LineStyle = linestyles{1,:};
p.Marker = markers{1,:};
p.Color = colors{1,:};

grid minor
box on
axis tight
xlabel('Nondimensionalized Radial Location', 'FontSize', 15);
ylabel('Time-Averaged Thrust (N)', 'FontSize', 15);

hold on
for i = 1:size(z,1)
loc = 2.*(Design(i).ITEROUTP(end).OUTP.ROTOR(1).vecSPANLOC)./Design(i).ITEROUTP(end).OUTP.vecROTDIAM(1);
thrust = mean(mean(cat(1, Design(i).ITEROUTP(end).OUTP.ROTOR(:).vecTHRUSTDIST_AVG),3),1);
p = plot(loc, thrust);
p.LineStyle = linestyles{i+1,:};
p.Marker = markers{i+1,:};
p.Color = colors{i+1,:};
end
lgnd = legend(legend_entry,'Location','NorthWest');
hold off

%% Prop Power Distribution
hFig203 = figure(203);
clf(203);
loc = 2.*(Design(end).ITEROUTP(end).OUTP.ROTOR.vecSPANLOC)./Design(end).ITEROUTP(end).OUTP.vecROTDIAM(1);
p = plot(loc, mean(Design(end).ITEROUTP(end).OUTP.ROTOR.vecTORQUEDIST_AVG,3));
p.LineStyle = linestyles{1,:};
p.Marker = markers{1,:};
p.Color = colors{1,:};
grid minor
box on
axis tight
xlabel('Nondimensionalized Radial Location', 'FontSize', 15);
ylabel('Time-Averaged Torque (Nm)', 'FontSize', 15);

hold on
for i = 1:size(z,1)
loc = 2.*(Design(i).ITEROUTP(end).OUTP.ROTOR(1).vecSPANLOC)./Design(i).ITEROUTP(end).OUTP.vecROTDIAM(1);
torque = mean(mean(cat(1, Design(i).ITEROUTP(end).OUTP.ROTOR(:).vecTORQUEDIST_AVG),3),1);
p = plot(loc, torque);
p.LineStyle = linestyles{i+1,:};
p.Marker = markers{i+1,:};
p.Color = colors{i+1,:};
end
lgnd = legend(legend_entry,'Location','NorthWest');
hold off

%% Lift Distribution
hFig204 = figure(204);
clf(204);
loc = 2.*(Design(1).ITEROUTP(end).OUTP.WING(1).vecSPANLOC)./10;
p = plot(loc, Design(end).ITEROUTP(end).OUTP.WING(1).vecLDIST_AVG);
p.LineStyle = linestyles{1,:};
p.Marker = markers{1,:};
p.Color = colors{1,:};
grid minor
box on
axis tight
xlabel('Nondimensionalized Spanwise Location', 'FontSize', 15);
ylabel('Lift (N)', 'FontSize', 15);

hold on
for i = 1:size(z,1)
loc = 2.*(Design(i).ITEROUTP(end).OUTP.WING(1).vecSPANLOC)./10;
lift = Design(i).ITEROUTP(end).OUTP.WING(:).vecLDIST_AVG;
p = plot(loc, lift);
p.LineStyle = linestyles{i+1,:};
p.Marker = markers{i+1,:};
p.Color = colors{i+1,:};
end
lgnd = legend(legend_entry,'Location','NorthWest');
hold off

%% Drag Distribution
hFig205 = figure(205);
clf(205);
loc = 2.*(Design(1).ITEROUTP(end).OUTP.WING(1).vecSPANLOC)./10;
drag = Design(end).ITEROUTP(end).OUTP.WING(1).vecDIDIST_AVG(11:end) + Design(end).ITEROUTP(end).OUTP.WING(1).vecDPDIST_AVG;
p = plot(loc, drag);
p.LineStyle = linestyles{1,:};
p.Marker = markers{1,:};
p.Color = colors{1,:};
grid minor
box on
axis tight
xlabel('Nondimensionalized Spanwise Location', 'FontSize', 15);
ylabel('Drag (N)', 'FontSize', 15);

hold on
for i = 1:size(z,1)
    loc = 2.*(Design(i).ITEROUTP(end).OUTP.WING(1).vecSPANLOC)./10;
    drag = Design(i).ITEROUTP(end).OUTP.WING(1).vecDIDIST_AVG(2:2:end) + Design(i).ITEROUTP(end).OUTP.WING(1).vecDPDIST_AVG;
    p = plot(loc, drag);
    p.LineStyle = linestyles{i+1,:};
    p.Marker = markers{i+1,:};
    p.Color = colors{i+1,:};
end
lgnd = legend(legend_entry,'Location','NorthWest');
hold off

