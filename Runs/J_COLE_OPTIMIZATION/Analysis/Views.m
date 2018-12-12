clc
clear

% cd G:\GIT\VAP3\Runs\J_COLE_OPTIMIZATION\Analysis

figure(300);
clf(300);

addpath('../../../')
addpath('../../../airfoils')
addpath('./../../../Runs/J_COLE_OPTIMIZATION/aux_files')
addpath('./../../../Runs/J_COLE_OPTIMIZATION/')
addpath('./../../../Runs/J_COLE_OPTIMIZATION/Analysis')


z{1} = [76.3609 73.34 71.08 69.358 66.56 64.9754 62.04 59.78 57.5467 55.4369 53 107.184 2980.88 129.922 -9.82695 0.687749 262.291 -1.80909 0.764347 396.786 -5.68468 0.630234]; % 76519;

% z{2} = [75.176 71.3419 69.6033 68.1224 64.8988 64.9348 61.9753 59.1098 59.9463 57.9888 52.4117 159.43 2276.93 464.578 -9.70714 0.886147];
% z{3} = [75.445 72.566 70.296 67.249 65.547 63.779 60.086 60.699 58.761 60.124 54.47 100 4982.9 463.55 -9.288 0.70461];
rotors = [3 1 1 1];

air_temp = -1; % Celsius at 8000 ft altitude
c = 331.3*sqrt(1 + air_temp/273.15);
max_tip_speed = 0.84*c; % Mach 0.84 max tip speed
min_tip_speed = 0.5*c;
cases = [1 5 3];

for i = 1:4
    
    if i <= 3
        make_vap_file(z{i}, i, 11, 0, rotors(i), 3, max_tip_speed, min_tip_speed)
    end
    
    cd '../../../'
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = 0;
    VAP_IN.vecCOLLECTIVE = repmat(0, rotors(i), 1);
    VAP_IN.vecVEHVINF = 77.2;
    VAP_IN.valMAXTIME = 1
    VAP_IN.valSTARTFORCES = 1
    VAP_IN.valDELTIME = (1/60)/(2250/60);
    if i <= 3
        OUTP = fcnVAP_MAIN(sprintf('Design_%d.vap',i), VAP_IN);
    else
        OUTP = fcnVAP_MAIN('X57_BASELINE.vap', VAP_IN);
    end
    cd 'Runs/J_COLE_OPTIMIZATION/Analysis/'
    
    hFig3 = gcf;
    hFig_view = findobj('Parent',hFig3,'Type','axes');
    hFig_temp = figure(300+i);
    clf(300+i)
    hSub1 = subplot(2,1,1);
    copyobj(get(hFig_view(1),'Children'),hSub1);
    gcf;
    view([-90 0])
    
    title('Front View','FontSize',15)
    xlabel('Y-Direction (m)','FontSize',10);
    zlabel('Z-Direction (m)','FontSize',10);
%     title(titles{i},'FontSize',15);
    hold on
    %     circular_arrow(hFig_temp, 1.1*(z(i,12)/2)./100, z(i,15:-1:14)./100, 90, 70, 2, 'k', 1000);
    hold off
    set(gca,'XTick',(-1:0.5:1))
    box off
    grid on
    axis image
    
    hSub2 = subplot(2,1,2);
    copyobj(get(hFig_view(1),'Children'),hSub2);
    gcf;
    view([90 90])
    box off
    grid on
    axis image
    
    title('Top View','FontSize',15)
    xlabel('X-Direction (m)','FontSize',10);
    ylabel('Y-Direction (m)','FontSize',10);
    
    if i <= 3
        fcnFIG2LATEX(gcf, ['View_Case_',num2str(cases(i)),'.pdf'],[9 4]);
    else
        fcnFIG2LATEX(gcf, ['View_Baseline.pdf'],[9 4]);
    end
end

% len = size(z,1) + 1;
% for i = 1:len
%     cd '../../../'
%     VAP_IN = [];
%     VAP_IN.vecVEHALPHA = 0;
%     VAP_IN.vecCOLLECTIVE = repmat(0, num_props(i), 1);
%     VAP_IN.vecVEHVINF = 77.2;
%     VAP_IN.valMAXTIME = 0
%     VAP_IN.valSTARTFORCES = 0
%     VAP_IN.valDELTIME = (1/60)/(2250/60);
%     if i <= len - 1
%         OUTP = fcnVAP_MAIN(sprintf('Design_%d.vap',i), VAP_IN);
%     else
%         OUTP = fcnVAP_MAIN('X57_BASELINE.vap', VAP_IN);
%     end
%     cd 'Runs/J_COLE_OPTIMIZATION/Analysis/'
%
%     hFig4 = gcf;
%     hFig_view = findobj('Parent',hFig4,'Type','axes');
%     hFig_temp = figure(301);
% %     clf(301)
%     hSub1 = subplot(len,1,(len + 1) - i);
%     copyobj(get(hFig_view(1),'Children'),hSub1);
%     view([90 90])
%
% %     hold on
% % %     circular_arrow(hFig_temp, 1.1*(z(i,12)/2)./100, z(i,15:-1:14)./100, 90, 70, 2, 'k', 1000);
% %     hold off
% %     set(gca,'XTick',(-1:0.5:1))
% %     box off
% %     grid on
% %     axis image
% %
% %     hSub2 = subplot(2,1,2);
% %     copyobj(get(hFig_view(1),'Children'),hSub2);
% %     view([90 90])
%
%     box off
%     grid on
%     axis image
% end
%
%
% for i = 1:len
%     try clf(302); end;
%     cd '../../../'
%     VAP_IN = [];
%     VAP_IN.vecVEHALPHA = 0;
%     VAP_IN.vecCOLLECTIVE = repmat(0, num_props(i), 1);
%     VAP_IN.vecVEHVINF = 77.2;
%     VAP_IN.valMAXTIME = 0
%     VAP_IN.valSTARTFORCES = 0
%     VAP_IN.valDELTIME = (1/60)/(2250/60);
%     if i <= len - 1
%         OUTP = fcnVAP_MAIN(sprintf('Design_%d.vap',i), VAP_IN);
%     else
%         OUTP = fcnVAP_MAIN('X57_BASELINE.vap', VAP_IN);
%     end
%     cd 'Runs/J_COLE_OPTIMIZATION/Analysis/'
%
%     hFig5 = gcf;
%     hFig_view = findobj('Parent',hFig5,'Type','axes');
%     hFig_temp = figure(302);
% %     clf(301)
%     hSub1 = axes;
%     copyobj(get(hFig_view(1),'Children'), hSub1);
%     view([90 90])
%
% %     hold on
% % %     circular_arrow(hFig_temp, 1.1*(z(i,12)/2)./100, z(i,15:-1:14)./100, 90, 70, 2, 'k', 1000);
% %     hold off
% %     set(gca,'XTick',(-1:0.5:1))
% %     box off
% %     grid on
% %     axis image
% %
% %     hSub2 = subplot(2,1,2);
% %     copyobj(get(hFig_view(1),'Children'),hSub2);
% %     view([90 90])
%
%     box off
%     grid on
%     axis image
% end