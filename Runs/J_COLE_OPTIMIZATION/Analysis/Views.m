clc
clear

cd C:\Users\travi\OneDrive\Desktop\GIT\VAP3\Runs\J_COLE_OPTIMIZATION\Analysis

% z(:,12) = diameter
% z(:,14:16) = loc
z(1,:) = [76 74 72 68 66 65 63 59 57 56 53 155 2250 -2 482 0 1 8 700 0 1 18 900 0 1 27 1100 0 1 35 1300 0 1 51 1600 0 1]; % 76907.3
z(2,:) = [76 74 72 70 66 64 62 62 57 55 52 142 2209 -13 249 0 1 -3 476 1 0 18 888 1 1 35 1309 -2 1 52 1607 1 1 61 1801 1 1]; % 72133.5
z(3,:) = [75 74 72 68 66 64 61 59 57 53 53 110 2100 -15 200 -1 1 -8 351 0 0 -2 483 0 1 16 900 -0 1 33 1300 1 1 52 1601 -0 1]; % 73381.7

num_props = [1 2 3 1];

for i = 1:4
    cd '../../../'
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = 0;
    VAP_IN.vecCOLLECTIVE = repmat(0, num_props(i), 1);
    VAP_IN.vecVEHVINF = 77.2;
    VAP_IN.valMAXTIME = 0
    VAP_IN.valSTARTFORCES = 0
    VAP_IN.valDELTIME = (1/60)/(2250/60);
    if i <= 3
        OUTP = fcnVAP_MAIN(sprintf('Design_%d.vap',i), VAP_IN);
    else
        OUTP = fcnVAP_MAIN('X57_BASELINE.vap', VAP_IN);
    end
    cd 'Runs/J_COLE_OPTIMIZATION/Analysis/'
    
    hFig3 = gcf;
    hFig_view = findobj('Parent',hFig3,'Type','axes');
    hFig_temp = figure(300);
    clf(300)
    hSub1 = subplot(2,1,1);
    copyobj(get(hFig_view(1),'Children'),hSub1);
    view([-90 0])
    
    hold on
    circular_arrow(hFig_temp, 1.1*(z(i,12)/2)./100, z(i,15:-1:14)./100, 90, 70, 2, 'k', 1000);
    hold off
    set(gca,'XTick',(-1:0.5:1))
    box off
    grid on
    axis image
    
    hSub2 = subplot(2,1,2);
    copyobj(get(hFig_view(1),'Children'),hSub2);
    view([90 90])
    
    box off
    grid on
    axis image
end