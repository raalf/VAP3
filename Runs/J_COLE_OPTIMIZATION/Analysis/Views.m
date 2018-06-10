clc
clear

addpath('../../../')
addpath('../../../airfoils')
addpath('./../../../Runs/J_COLE_OPTIMIZATION/aux_files')
addpath('./../../../Runs/J_COLE_OPTIMIZATION/')

% z(:,12) = diameter
% z(:,14:16) = loc
z(1,:) = [76 74 69 69 66 64 63 60 55 55 54 zeros(1, 11) 160 2229 -5 451 -4 1 11 701 2 1 17 894 1 1];

num_props = [1 1];

for i = 1:size(z,1) + 1
    cd '../../../'
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = 0;
    VAP_IN.vecCOLLECTIVE = repmat(0, num_props(i), 1);
    VAP_IN.vecVEHVINF = 77.2;
    VAP_IN.valMAXTIME = 0
    VAP_IN.valSTARTFORCES = 0
    VAP_IN.valDELTIME = (1/60)/(2250/60);
    if i <= 1
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
%     circular_arrow(hFig_temp, 1.1*(z(i,12)/2)./100, z(i,15:-1:14)./100, 90, 70, 2, 'k', 1000);
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