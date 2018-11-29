clc
clear

cd C:\Users\travi\OneDrive\Desktop\GIT\VAP3\Runs\J_COLE_OPTIMIZATION\Analysis

addpath('../../../')
addpath('../../../airfoils')
addpath('./../../../Runs/J_COLE_OPTIMIZATION/aux_files')
addpath('./../../../Runs/J_COLE_OPTIMIZATION/')

z(1,:) = [75.176 71.3419 69.6033 68.1224 64.8988 64.9348 61.9753 59.1098 59.9463 57.9888 52.4117 159.43 2276.93 464.578 -9.70714 0.886147]; %75529
z(2,:) = [75.445 72.566 70.296 67.249 65.547 63.779 60.086 60.699 58.761 60.124 54.47 100 4982.9 463.55 -9.288 0.70461]; %80359

z(3,:) = [linspace(75.6, 53, 11) 100 4153 463.55 0 1]; %81000
z(4,:) = [linspace(75.6, 53, 11) 100 4153 170.00 0 0.1]; %81500

rotors = [1, 1, 1, 1, 1, 3];
linestyles = {'--';'-.';'-';':'};
markers = {'o';'x';'s';'^';'*';'d';'v';'>';'<';'p';'h'};
colors = {'k';'b';'r';'m';'c';'g'};


air_temp = -1; % Celsius at 8000 ft altitude
c = 331.3*sqrt(1 + air_temp/273.15);
max_tip_speed = 0.84*c; % Mach 0.84 max tip speed
min_tip_speed = 0.5*c;

% cd ./../
% for i = 1:size(z,1) + 1
%     if i == size(z,1) + 1
%         [out_all(1,i), Design_all(i).ITER, Design_all(i).ITEROUTP] = fcnBASELINE_OBJ();
%     else
%         [out_all(1,i), Design_all(i).ITER, Design_all(i).ITEROUTP] = fcnOBJECTIVE(z(i,:), 11, 0, rotors(i), 3, max_tip_speed, min_tip_speed, pwd);
%     end
% end
% cd Analysis/
% % save('matlab.mat');
% 

% load('matlab4.mat');
% Design_all(6) = [];
% rotors = [1, 1, 1, 1, 1, 3];
% i = 6;
% z = [76.3609 73.34 71.08 69.358 66.56 64.9754 62.04 59.78 57.5467 55.4369 53 107.184 2980.88 129.922 -9.82695 0.687749 262.291 -1.80909 0.764347 396.786 -5.68468 0.630234]; % 76519.1;
% cd ./../
% [out_all(1,i), Design_all(i).ITER, Design_all(i).ITEROUTP] = fcnOBJECTIVE(z, 11, 0, rotors(i), 3, max_tip_speed, min_tip_speed, 'G:\GIT\VAP3\Runs\J_COLE_OPTIMIZATION');
% cd Analysis/
% save('matlab5.mat');

load('matlab5.mat');
idx = [5 1 2 6];
len = length(idx);

m_wing = [2 2 2 2 2 1];
n_rotor = [19 19 19 19 19 10];

Design = Design_all;
out = out_all;

legend_entry = {'Baseline', 'Large Propeller', 'Small Propeller', 'Three Propeller'};

Design(6).ITEROUTP(end).OUTP.valDENSITY = .955;

% %% Drag Bar Graph
vinf = 77.2;
% baseline = [Design(idx(1)).ITEROUTP(end).OUTP.vecCD_AVG; ...
%     Design(idx(1)).ITEROUTP(end).OUTP.vecCDI_AVG; ...
%     Design(idx(1)).ITEROUTP(end).OUTP.vecCDP_AVG];
% baseline = baseline.*(0.5.*Design(idx(1)).ITEROUTP(end).OUTP.valDENSITY.*vinf.^2.*(Design(idx(1)).ITEROUTP(end).OUTP.valAREA));
% 
% for i = 2:len
%     designs(:,i-1) = [Design(idx(i)).ITEROUTP(end).OUTP.vecCD_AVG; Design(idx(i)).ITEROUTP(end).OUTP.vecCDI_AVG; Design(idx(i)).ITEROUTP(end).OUTP.vecCDP_AVG];
%     designs(:,i-1) = designs(:,i-1).*(0.5.*Design(idx(i)).ITEROUTP(end).OUTP.valDENSITY.*vinf.^2.*(Design(idx(i)).ITEROUTP(end).OUTP.valAREA));
% end
% 
% drag_table = array2table([baseline designs]','VariableNames',{'D', 'D_i', 'D_p'},'RowNames',legend_entry)
% 
% %% Power Bar Graph
% x = categorical({'Total Power', 'Induced Power', 'Profile Power'});
% 
% baseline = [mean(Design(idx(1)).ITEROUTP(end).OUTP.vecCP_AVG); ...
%     mean(Design(idx(1)).ITEROUTP(end).OUTP.vecCPI_AVG); ...
%     mean(Design(idx(1)).ITEROUTP(end).OUTP.vecCPP_AVG)];
% baseline = baseline.*(Design(idx(1)).ITEROUTP(end).OUTP.valDENSITY.*(abs(Design(idx(1)).ITEROUTP(end).OUTP.vecROTORRPM(1)./60).^3).*(Design(idx(1)).ITEROUTP(end).OUTP.vecROTDIAM(1).^5));
% 
% designs = [];
% for i = 2:len
%     designs(:,i-1) = [mean(Design(idx(i)).ITEROUTP(end).OUTP.vecCP_AVG,2); mean(Design(idx(i)).ITEROUTP(end).OUTP.vecCPI_AVG,2); mean(Design(idx(i)).ITEROUTP(end).OUTP.vecCPP_AVG,2)];
%     designs(:,i-1) = designs(:,i-1).*rotors(idx(i)).*(Design(idx(i)).ITEROUTP(end).OUTP.valDENSITY.*(abs(Design(idx(i)).ITEROUTP(end).OUTP.vecROTORRPM(1)./60).^3).*(Design(idx(i)).ITEROUTP(end).OUTP.vecROTDIAM(1).^5));
% end
% 
% power_table = array2table([baseline.*2 designs.*2]','VariableNames',{'P', 'P_i', 'P_p'},'RowNames',legend_entry)
% 
% 
% % hFig201 = figure(201);
% % clf(201);
% % 
% % b = bar( x, 100.*(designs - repmat(baseline,[1,size(z,1)]))./repmat(baseline,[1,size(z,1)]) );
% % b(1).LineStyle = '--';
% % b(1).LineWidth = 1;
% % % b(2).LineStyle = '-.';
% % % b(2).LineWidth = 1;
% % % b(3).LineStyle = ':';
% % % b(3).LineWidth = 1;
% % grid minor
% % box on
% % axis tight
% % 
% % ylabel('Percent Change from Baseline', 'FontSize', 15)
% % % legend('Design 1', 'Design 2', 'Design 3')
% 
% %% Prop Thrust Distribution
% hFig202 = figure(202);
% clf(202);
% loc = 2.*(Design(idx(1)).ITEROUTP(end).OUTP.ROTOR.vecSPANLOC)./Design(idx(1)).ITEROUTP(end).OUTP.vecROTDIAM(1);
% p = plot(loc, mean(reshape(Design(idx(1)).ITEROUTP(end).OUTP.ROTOR.vecTHRUSTDIST_AVG, 1, 19, 3),3)./mean(diff(loc)));
% p.LineStyle = linestyles{1,:};
% p.Marker = markers{1,:};
% p.Color = colors{1,:};
% 
% grid minor
% box on
% % axis tight
% xlabel('Nondimensionalized Radial Location', 'FontSize', 12);
% ylabel('Time-Averaged Thrust Distribution (N/m)', 'FontSize', 12);
% 
% hold on
% for i = 2:len
%     loc = 2.*(Design(idx(i)).ITEROUTP(end).OUTP.ROTOR(1).vecSPANLOC)./Design(idx(i)).ITEROUTP(end).OUTP.vecROTDIAM(1);
%     thrust = mean(mean(reshape(cat(1, [Design(idx(i)).ITEROUTP(end).OUTP.ROTOR(:).vecTHRUSTDIST_AVG]), 1, n_rotor(idx(i)), rotors(idx(i))*3),3),1)./mean(diff(loc));
%     p = plot(loc, thrust);
%     p.LineStyle = linestyles{i,:};
%     p.Marker = markers{i,:};
%     p.Color = colors{i,:};
% end
% lgnd = legend(legend_entry,'Location','NorthWest');
% hold off
% 
% % Figure handle, figure name, figure path, [WIDTH HEIGHT] in inches
% fcnFIG2LATEX(hFig202, 'prop_thrust_dist.pdf', 'figures', [7 5]);
% 
% %% Prop Power Distribution
% hFig203 = figure(203);
% clf(203);
% loc = 2.*(Design(idx(1)).ITEROUTP(end).OUTP.ROTOR.vecSPANLOC)./Design(idx(1)).ITEROUTP(end).OUTP.vecROTDIAM(1);
% p = plot(loc, mean(reshape(Design(idx(1)).ITEROUTP(end).OUTP.ROTOR.vecTORQUEDIST_AVG, 1, 19, 3),3)./mean(diff(loc)));
% p.LineStyle = linestyles{1,:};
% p.Marker = markers{1,:};
% p.Color = colors{1,:};
% grid minor
% box on
% % axis tight
% xlabel('Nondimensionalized Radial Location', 'FontSize', 12);
% ylabel('Time-Averaged Torque Distribution (Nm/m)', 'FontSize', 12);
% 
% hold on
% for i = 2:len
%     loc = 2.*(Design(idx(i)).ITEROUTP(end).OUTP.ROTOR(1).vecSPANLOC)./Design(idx(i)).ITEROUTP(end).OUTP.vecROTDIAM(1);
%     torque = mean(mean(reshape(cat(1, [Design(idx(i)).ITEROUTP(end).OUTP.ROTOR(:).vecTORQUEDIST_AVG]), 1, n_rotor(idx(i)), rotors(idx(i))*3),3),1)./mean(diff(loc));
%     p = plot(loc, torque);
%     p.LineStyle = linestyles{i,:};
%     p.Marker = markers{i,:};
%     p.Color = colors{i,:};
% end
% lgnd = legend(legend_entry,'Location','NorthWest');
% hold off
% 
% % Figure handle, figure name, figure path, [WIDTH HEIGHT] in inches
% fcnFIG2LATEX(hFig203, 'prop_power_dist.pdf', 'figures', [7 5]);
% 
% %% Lift Distribution
% hFig204 = figure(204);
% clf(204);
% loc = [Design(6).ITEROUTP(end).OUTP.WING(1).vecSPANLOC_PROJ; 4.82]./4.82;
% lift = [Design(idx(1)).ITEROUTP(end).OUTP.WING(1).vecLDIST_AVG 0]./(4.82/20);
% 
% p = plot(loc, lift);
% p.LineStyle = linestyles{1,:};
% p.Marker = markers{1,:};
% p.Color = colors{1,:};
% grid minor
% box on
% % axis tight
% xlabel('Nondimensionalized Spanwise Location', 'FontSize', 12);
% ylabel('Lift Distribution (N/m)', 'FontSize', 12);
% 
% hold on
% for i = 2:len
% %     loc = 2.*(Design(i).ITEROUTP(end).OUTP.WING(1).vecSPANLOC)./10;
%     lift = [Design(idx(i)).ITEROUTP(end).OUTP.WING(1).vecLDIST_AVG 0]./(4.82/20);
%     p = plot(loc, lift);
%     p.LineStyle = linestyles{i,:};
%     p.Marker = markers{i,:};
%     p.Color = colors{i,:};
%  end
% lgnd = legend(legend_entry,'Location','SouthWest');
% hold off
% 
% % Figure handle, figure name, figure path, [WIDTH HEIGHT] in inches
% fcnFIG2LATEX(hFig204, 'lift_dist.pdf', 'figures', [7 5]);
% 
% %% Drag Distribution
% hFig205 = figure(205);
% clf(205);
% loc = [Design(6).ITEROUTP(end).OUTP.WING(1).vecSPANLOC_PROJ; 4.82]./4.82;
% drag = [Design(idx(1)).ITEROUTP(end).OUTP.WING(1).vecDIDIST_AVG(21:end) + Design(idx(1)).ITEROUTP(end).OUTP.WING(1).vecDPDIST_AVG 0]./(4.82/20);
% p = plot(loc, drag);
% p.LineStyle = linestyles{1,:};
% p.Marker = markers{1,:};
% p.Color = colors{1,:};
% grid minor
% box on
% % axis tight
% xlabel('Nondimensionalized Spanwise Location', 'FontSize', 12);
% ylabel('Drag Distribution (N/m)', 'FontSize', 12);
% 
% hold on
% for i = 2:len
% %     loc = 2.*(Design(i).ITEROUTP(end).OUTP.WING(1).vecSPANLOC)./10;
%     drag = [[nonzeros(Design(idx(i)).ITEROUTP(end).OUTP.WING(1).vecDIDIST_AVG)]' + Design(idx(i)).ITEROUTP(end).OUTP.WING(1).vecDPDIST_AVG 0]./(4.82/20);
%     p = plot(loc, drag);
%     p.LineStyle = linestyles{i,:};
%     p.Marker = markers{i,:};
%     p.Color = colors{i,:};
% end
% lgnd = legend(legend_entry,'Location','NorthWest');
% hold off
% 
% % Figure handle, figure name, figure path, [WIDTH HEIGHT] in inches
% fcnFIG2LATEX(hFig205, 'drag_dist.pdf', 'figures', [7 7]);
% 
% % %% Drag Distribution
% % hFig808 = figure(808);
% % clf(808);
% % loc = [Design(6).ITEROUTP(end).OUTP.WING(1).vecSPANLOC_PROJ; 4.82]./4.82;
% % drag = [Design(idx(1)).ITEROUTP(end).OUTP.WING(1).vecDPDIST_AVG 0];
% % p = plot(loc, drag);
% % p.LineStyle = linestyles{1,:};
% % p.Marker = markers{1,:};
% % p.Color = colors{1,:};
% % grid minor
% % box on
% % axis tight
% % xlabel('Spanwise Location (m)', 'FontSize', 15);
% % ylabel('Drag Distribution (N/m)', 'FontSize', 15);
% % 
% % hold on
% % for i = 2:len
% % %     loc = 2.*(Design(i).ITEROUTP(end).OUTP.WING(1).vecSPANLOC)./10;
% %     drag = [Design(idx(i)).ITEROUTP(end).OUTP.WING(1).vecDPDIST_AVG 0];
% %     p = plot(loc, drag);
% %     p.LineStyle = linestyles{i,:};
% %     p.Marker = markers{i,:};
% %     p.Color = colors{i,:};
% % end
% % lgnd = legend(legend_entry,'Location','NorthWest');
% % hold off
% % 
% % % Figure handle, figure name, figure path, [WIDTH HEIGHT] in inches
% % fcnFIG2LATEX(hFig808, 'drag_dist.pdf', 'figures', [7 7]);
% 
% %% Fitness vs Generation
% hFig206 = figure(206);
% clf(206);
% 
% i = [1, 3, 5];
% for j = 1:length(i)
%     maxes = [];
%     A = dlmread(['case_', num2str(i(j)), '.txt'],'');
%     
%     subplot(length(i),1,j);
%     
%     tm1 = 50;
%     if i(j) == 3; tm1 = 7; end
%     for jj = 1:tm1:size(A,1)
%         if jj+tm1 > size(A,1)
%             maxes = [maxes; min(A(jj,1):A(end,1))];
%         else
%             maxes = [maxes; min(A(jj,1):A(jj+tm1,1))];  
%         end            
%     end
%     scatter(1:length(maxes), sort(maxes,'descend'), 'xk');
%     
%     hold on
%     plot(1:length(maxes), repmat(out(idx(1)), length(maxes),1),'r','LineWidth',2)
%     hold off
%     
%     xlabel('Generation','FontSize',10);
%     ylabel('Best Fitness (W)','FontSize',10);
%     box on
%     grid minor
%     title(['Case ', num2str(i(j))],'FontSize',12)
% %     axis tight
%     
% end
% % A = dlmread('fitness.txt','',2,0);
% 
% % Figure handle, figure name, figure path, [WIDTH HEIGHT] in inches
% fcnFIG2LATEX(hFig206, 'fitness_vs_generation.pdf', 'figures', [7 5]);
% 
% %%

idx = [3 4];
len = length(idx);
designs = [];
%% Drag Bar Graph
for i = 1:len
    designs(:,i) = [Design(idx(i)).ITEROUTP(end).OUTP.vecCD_AVG; Design(idx(i)).ITEROUTP(end).OUTP.vecCDI_AVG; Design(idx(i)).ITEROUTP(end).OUTP.vecCDP_AVG];
    designs(:,i) = designs(:,i).*(0.5.*Design(idx(i)).ITEROUTP(end).OUTP.valDENSITY.*vinf.^2.*(Design(idx(i)).ITEROUTP(end).OUTP.valAREA));
    d_tot(i,1) = designs(1,i);
end

legend_entry = {'2y_p/b = 0.96, Inboard up', '2y_p/b = 0.35, Inboard down'};
drag_table_small_props = array2table([designs]','VariableNames',{'D', 'D_i', 'D_p'},'RowNames',legend_entry)

%% Power
designs = [];
for i = 1:len
    designs(:,i) = [mean(Design(idx(i)).ITEROUTP(end).OUTP.vecCP_AVG,2); mean(Design(idx(i)).ITEROUTP(end).OUTP.vecCPI_AVG,2); mean(Design(idx(i)).ITEROUTP(end).OUTP.vecCPP_AVG,2)];
    designs(:,i) = designs(:,i).*rotors(idx(i)).*(Design(idx(i)).ITEROUTP(end).OUTP.valDENSITY.*(abs(Design(idx(i)).ITEROUTP(end).OUTP.vecROTORRPM(1)./60).^3).*(Design(idx(i)).ITEROUTP(end).OUTP.vecROTDIAM(1).^5));
    eff(:,i) = (vinf.*(d_tot(i)+500))./...
        (mean(Design(idx(i)).ITEROUTP(end).OUTP.vecCP_AVG,2).*rotors(idx(i)).*(Design(idx(i)).ITEROUTP(end).OUTP.valDENSITY.*(abs(Design(idx(i)).ITEROUTP(end).OUTP.vecROTORRPM(1)./60).^3).*(Design(idx(i)).ITEROUTP(end).OUTP.vecROTDIAM(1).^5)).*2);
end

power_table_small_props = array2table([designs.*2]','VariableNames',{'P', 'P_i', 'P_p'},'RowNames',legend_entry)

%% Prop Thrust Distribution
hFig210 = figure(210);
clf(210);
% loc = 2.*(Design(idx(1)).ITEROUTP(end).OUTP.ROTOR.vecSPANLOC)./Design(idx(1)).ITEROUTP(end).OUTP.vecROTDIAM(1);
% p = plot(loc, mean(reshape(Design(idx(1)).ITEROUTP(end).OUTP.ROTOR.vecTHRUSTDIST_AVG, 1, 19, 3),3));
% p.LineStyle = linestyles{1,:};
% p.Marker = markers{1,:};
% p.Color = colors{1,:};

grid minor
box on
% axis tight
xlabel('Nondimensionalized Radial Location', 'FontSize', 12);
ylabel('Time-Averaged Thrust Distribution (N/m)', 'FontSize', 12);

hold on
for i = 1:len
    loc = 2.*(Design(idx(i)).ITEROUTP(end).OUTP.ROTOR(1).vecSPANLOC)./Design(idx(i)).ITEROUTP(end).OUTP.vecROTDIAM(1);
    thrust = mean(mean(reshape(cat(1, [Design(idx(i)).ITEROUTP(end).OUTP.ROTOR(:).vecTHRUSTDIST_AVG]), 1, n_rotor(idx(i)), rotors(idx(i))*3),3),1)./mean(diff(loc));
    p = plot(loc, thrust);
    p.LineStyle = linestyles{i,:};
    p.Marker = markers{i,:};
    p.Color = colors{i,:};
end
lgnd = legend(legend_entry,'Location','SouthEast');
hold off

% Figure handle, figure name, figure path, [WIDTH HEIGHT] in inches
fcnFIG2LATEX(hFig210, 'small_inboard_v_outboard_thrust.pdf', 'figures', [7 5]);

%% Prop Power Distribution
hFig211 = figure(211);
clf(211);
% loc = 2.*(Design(idx(1)).ITEROUTP(end).OUTP.ROTOR.vecSPANLOC)./Design(idx(1)).ITEROUTP(end).OUTP.vecROTDIAM(1);
% p = plot(loc, mean(reshape(Design(idx(1)).ITEROUTP(end).OUTP.ROTOR.vecTORQUEDIST_AVG, 1, 19, 3),3));
% p.LineStyle = linestyles{1,:};
% p.Marker = markers{1,:};
% p.Color = colors{1,:};
grid minor
box on
% axis tight
xlabel('Nondimensionalized Radial Location', 'FontSize', 12);
ylabel('Time-Averaged Torque Distribution (Nm/m)', 'FontSize', 12);

hold on
for i = 1:len
    loc = 2.*(Design(idx(i)).ITEROUTP(end).OUTP.ROTOR(1).vecSPANLOC)./Design(idx(i)).ITEROUTP(end).OUTP.vecROTDIAM(1);
    torque = mean(mean(reshape(cat(1, [Design(idx(i)).ITEROUTP(end).OUTP.ROTOR(:).vecTORQUEDIST_AVG]), 1, n_rotor(idx(i)), rotors(idx(i))*3),3),1)./mean(diff(loc));
    p = plot(loc, torque);
    p.LineStyle = linestyles{i,:};
    p.Marker = markers{i,:};
    p.Color = colors{i,:};
end
lgnd = legend(legend_entry,'Location','SouthEast');
hold off

% Figure handle, figure name, figure path, [WIDTH HEIGHT] in inches
fcnFIG2LATEX(hFig211, 'small_inboard_v_outboard_power.pdf', 'figures', [7 5]);

%% Lift Distribution
hFig214 = figure(214);
clf(214);
loc = [Design(6).ITEROUTP(end).OUTP.WING(1).vecSPANLOC_PROJ; 4.82]./(4.82);
% p = plot(loc, [Design(idx(1)).ITEROUTP(end).OUTP.WING(1).vecLDIST_AVG 0]);
% p.LineStyle = linestyles{1,:};
% p.Marker = markers{1,:};
% p.Color = colors{1,:};
grid minor
box on
% axis tight
xlabel('Nondimensionalized Spanwise Location', 'FontSize', 12);
ylabel('Lift Distribution (N/m)', 'FontSize', 12);

hold on
for i = 1:len
%     loc = 2.*(Design(i).ITEROUTP(end).OUTP.WING(1).vecSPANLOC)./10;
    lift = [Design(idx(i)).ITEROUTP(end).OUTP.WING(1).vecLDIST_AVG 0]./(4.82/20);
    p = plot(loc, lift);
    p.LineStyle = linestyles{i,:};
    p.Marker = markers{i,:};
    p.Color = colors{i,:};
end
lgnd = legend(legend_entry,'Location','SouthWest');
hold off

% Figure handle, figure name, figure path, [WIDTH HEIGHT] in inches
fcnFIG2LATEX(hFig214, 'lift_dist_small_props.pdf', 'figures', [7 5]);

%% Drag Distribution
hFig215 = figure(215);
clf(215);
loc = [Design(6).ITEROUTP(end).OUTP.WING(1).vecSPANLOC_PROJ; 4.82]./(4.82);
% drag = [Design(idx(1)).ITEROUTP(end).OUTP.WING(1).vecDIDIST_AVG(21:end) + Design(idx(1)).ITEROUTP(end).OUTP.WING(1).vecDPDIST_AVG 0]./(4.82/20);
% p = plot(loc, drag);
% p.LineStyle = linestyles{1,:};
% p.Marker = markers{1,:};
% p.Color = colors{1,:};
grid minor
box on
% axis tight
xlabel('Nondimensionalized Spanwise Location', 'FontSize', 12);
ylabel('Drag Distribution (N/m)', 'FontSize', 12);

hold on
for i = 1:len
%     loc = 2.*(Design(i).ITEROUTP(end).OUTP.WING(1).vecSPANLOC)./10;
    drag = [[nonzeros(Design(idx(i)).ITEROUTP(end).OUTP.WING(1).vecDIDIST_AVG)]' + Design(idx(i)).ITEROUTP(end).OUTP.WING(1).vecDPDIST_AVG 0]./(4.82/20);
    p = plot(loc, drag);
    p.LineStyle = linestyles{i,:};
    p.Marker = markers{i,:};
    p.Color = colors{i,:};
end
lgnd = legend(legend_entry,'Location','NorthEast');
hold off

% Figure handle, figure name, figure path, [WIDTH HEIGHT] in inches
fcnFIG2LATEX(hFig215, 'drag_dist_small_props.pdf', 'figures', [7 7]);

