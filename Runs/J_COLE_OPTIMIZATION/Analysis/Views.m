clc
clear

cd G:\GIT\VAP3\Runs\J_COLE_OPTIMIZATION\Analysis

try clf(301); end

addpath('../../../')
addpath('../../../airfoils')
addpath('./../../../Runs/J_COLE_OPTIMIZATION/aux_files')
addpath('./../../../Runs/J_COLE_OPTIMIZATION/')
addpath('./../../../Runs/J_COLE_OPTIMIZATION/Analysis')

% z(1,:) = [85.5402	83.9746	76.0336	72.5489	71.8403	71.1037	70.6373	62.761	59.9215	56.108	52.742	125.063	2730.78	116.264	-4.25819	0.647225	262.327	5.79578	0.86871	408.39	-11.7387	0.633426];
% rotors = [3];
% z(1,:) = [95.4245	92.8958	92.1974	94.3823	87.104	80.0929	78.1465	66.699	58.2872	57.833	55.6887	160	2370.87	434.46	-8.43005	1];
% rotors = 1;

z(1,:) = [95.445	92.7811	92.1621	94.3397	87.2277	80.1124	78.1725	66.7732	58.3955	57.8899	55.8115	160	2346.6	435.491	-8.41503	1];
rotors = 1;

for i = 1:size(z,1) + 1

    if i <= size(z,1)
    make_vap_file(z(i,:), i, 11, 0, rotors(i), 3)
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