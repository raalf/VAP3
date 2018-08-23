clc
clear

cd G:\GIT\VAP3\Runs\J_COLE_OPTIMIZATION\Analysis

try clf(301); end

addpath('../../../')
addpath('../../../airfoils')
addpath('./../../../Runs/J_COLE_OPTIMIZATION/aux_files')
addpath('./../../../Runs/J_COLE_OPTIMIZATION/')
addpath('./../../../Runs/J_COLE_OPTIMIZATION/Analysis')

z(1,:) = [96.5209203911987,93.9926180988485,93.9042645278558,93.5531270131122,88.7110422982123,86.6438164842058,80.8364050348930,79.1773155264245,78.0742363801537,65.2577667383824,52.6907600591495,2.50283349085077,3.56067448180554,6.63610548525079,11.9593033289778,12.1554365561683,15.4154947001996,33.4634950385750,36.1568526174222,58.4279985750835,65.1574769811730,70.0532521654860,139.619967982518,4956.09473672779,101.739854217312,11.2698248840502,0.877676813272480,262.359822198856,-8.88350665193784,0.132397495015192,422.979790180400,-10.6811043298787,0.640425530678599]
rotors = [3,5];

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