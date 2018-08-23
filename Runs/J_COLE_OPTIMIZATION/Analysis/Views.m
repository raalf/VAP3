clc
clear

cd G:\GIT\VAP3\Runs\J_COLE_OPTIMIZATION\Analysis

try clf(301); end

addpath('../../../')
addpath('../../../airfoils')
addpath('./../../../Runs/J_COLE_OPTIMIZATION/aux_files')
addpath('./../../../Runs/J_COLE_OPTIMIZATION/')
addpath('./../../../Runs/J_COLE_OPTIMIZATION/Analysis')

% z(:,12) = diameter
% z(:,14:16) = loc
z(1,:) = [71.3941 68.7533 68.6304 65.9588 64.5861 63.2590 61.8437 61.4030 61.1291 60.0904 59.0279 0 0 0 0 0 0 0 0 0 0 0 158.186 2226.86 450.979 -12.5485 0.770671]; % 74748.4
% z(1,:) = [85.9495 83.1285 80.1618 77.3493 75.7004 63.8025 60.9942 59.5311 58.3303 57.1046 52.2952 10.054 17.625	17.9299	20.7694	31.5342	33.8535	38.4692	44.6215	49.6259	64.9211	67.2667	159.808	2097.24	144.366	-13.7057 0.605461]; %75836.2

rotors = [1 1];

for i = 1:size(z,1) + 1

    if i <= size(z,1)
    make_vap_file(z(i,:), i, 11, rotors(i), 3)
    end
    
    cd '../../../'
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = 0;
    VAP_IN.vecCOLLECTIVE = repmat(0, rotors(i), 1);
    VAP_IN.vecVEHVINF = 77.2;
    VAP_IN.valMAXTIME = 0
    VAP_IN.valSTARTFORCES = 0
    VAP_IN.valDELTIME = (1/60)/(2250/60);
    if i <= size(z,1)
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

len = size(z,1) + 1;
for i = 1:len
    cd '../../../'
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = 0;
    VAP_IN.vecCOLLECTIVE = repmat(0, num_props(i), 1);
    VAP_IN.vecVEHVINF = 77.2;
    VAP_IN.valMAXTIME = 0
    VAP_IN.valSTARTFORCES = 0
    VAP_IN.valDELTIME = (1/60)/(2250/60);
    if i <= len - 1
        OUTP = fcnVAP_MAIN(sprintf('Design_%d.vap',i), VAP_IN);
    else
        OUTP = fcnVAP_MAIN('X57_BASELINE.vap', VAP_IN);
    end
    cd 'Runs/J_COLE_OPTIMIZATION/Analysis/'
    
    hFig4 = gcf;
    hFig_view = findobj('Parent',hFig4,'Type','axes');
    hFig_temp = figure(301);
%     clf(301)
    hSub1 = subplot(len,1,(len + 1) - i);
    copyobj(get(hFig_view(1),'Children'),hSub1);
    view([90 90])
    
%     hold on
% %     circular_arrow(hFig_temp, 1.1*(z(i,12)/2)./100, z(i,15:-1:14)./100, 90, 70, 2, 'k', 1000);
%     hold off
%     set(gca,'XTick',(-1:0.5:1))
%     box off
%     grid on
%     axis image
%     
%     hSub2 = subplot(2,1,2);
%     copyobj(get(hFig_view(1),'Children'),hSub2);
%     view([90 90])
    
    box off
    grid on
    axis image
end


for i = 1:len
    try clf(302); end;
    cd '../../../'
    VAP_IN = [];
    VAP_IN.vecVEHALPHA = 0;
    VAP_IN.vecCOLLECTIVE = repmat(0, num_props(i), 1);
    VAP_IN.vecVEHVINF = 77.2;
    VAP_IN.valMAXTIME = 0
    VAP_IN.valSTARTFORCES = 0
    VAP_IN.valDELTIME = (1/60)/(2250/60);
    if i <= len - 1
        OUTP = fcnVAP_MAIN(sprintf('Design_%d.vap',i), VAP_IN);
    else
        OUTP = fcnVAP_MAIN('X57_BASELINE.vap', VAP_IN);
    end
    cd 'Runs/J_COLE_OPTIMIZATION/Analysis/'
    
    hFig5 = gcf;
    hFig_view = findobj('Parent',hFig5,'Type','axes');
    hFig_temp = figure(302);
%     clf(301)
    hSub1 = axes;
    copyobj(get(hFig_view(1),'Children'), hSub1);
    view([90 90])
    
%     hold on
% %     circular_arrow(hFig_temp, 1.1*(z(i,12)/2)./100, z(i,15:-1:14)./100, 90, 70, 2, 'k', 1000);
%     hold off
%     set(gca,'XTick',(-1:0.5:1))
%     box off
%     grid on
%     axis image
%     
%     hSub2 = subplot(2,1,2);
%     copyobj(get(hFig_view(1),'Children'),hSub2);
%     view([90 90])
    
    box off
    grid on
    axis image
end