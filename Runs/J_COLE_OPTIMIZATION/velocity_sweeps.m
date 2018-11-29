clc
clear
cd G:\GIT\VAP3\Runs\J_COLE_OPTIMIZATION

% N_prop = 1
% file = 'Design_1' % Case 3
% rotor_diam = 1.594300;
% rotor_rpm = 2276.93;

% % N_prop = 3
% % file = 'Design_3' % Case 3
% % rotor_diam = 1.133130;
% % rotor_rpm = 2963.080000;
% 
N_prop = 1
file = 'Baseline' % Case 3
rotor_diam = 1.524;
rotor_rpm = 2250;


% %% Wing Sweep
% 
% wing_sweep_filename = ['Runs/J_COLE_OPTIMIZATION/',file, '/', file, '_wing.vap'];
% 
% cd '../../'
% seqALPHA = [0:11];
% parfor i = 1:length(seqALPHA)
%     VAP_IN = [];
%     VAP_IN.vecVEHALPHA = seqALPHA(i);
%     VAP_IN.valDELTIME = .25;
%     VAP_IN.valSTARTFORCES = 30;
%     VAP_IN.valMAXTIME = 30;
% %                     VAP_IN.valSTARTFORCES = 5
% %                     VAP_IN.valMAXTIME = 10
%     WING_SWEEP(i) = fcnVAP_MAIN(wing_sweep_filename, VAP_IN);
%     %     view([90 90]);
% end
% cd 'Runs/J_COLE_OPTIMIZATION/'
% 
% sweep_vinf = [WING_SWEEP.vecVINF];
% sweep_vinf = sweep_vinf(end,:);
% 
% %% Propeller Sweep at appropriate J (fix RPM)
% J = 77.2/((rotor_rpm/60)*rotor_diam);
% sweep_rpm = (sweep_vinf./(J.*rotor_diam)).*60;
% prop_sweep_filename = ['Runs/J_COLE_OPTIMIZATION/',file, '/', file, '_rotor.vap'];
% 
% cd '../../'
% for jj = 1:length(sweep_rpm)
%     vecCOLLECTIVE = [-7:2:7];
%     parfor i = 1:length(vecCOLLECTIVE)
%         VAP_IN = [];
%         VAP_IN.vecCOLLECTIVE = vecCOLLECTIVE(i);
%         VAP_IN.vecVEHALPHA = 0;
%         VAP_IN.vecVEHVINF = sweep_vinf(jj);
%         VAP_IN.vecROTORRPM = sweep_rpm(jj);
%         VAP_IN.valSTARTFORCES = 100;
%         VAP_IN.valMAXTIME = 100;
% %                         VAP_IN.valSTARTFORCES = 5
% %                         VAP_IN.valMAXTIME = 10
%         VAP_IN.valDELTIME = (1/60)/(sweep_rpm(jj)/60);
%         PROP_SWEEP(i,jj) = fcnVAP_MAIN(prop_sweep_filename, VAP_IN);
%         %     view([90 90]);
%     end
% end
% cd 'Runs/J_COLE_OPTIMIZATION/'
% 
% save([file, '/', file, '_one.mat']);
% 
% load([file, '/', file, '_one.mat'])
% 
% %% Trim
% rho = 1.007;
% weightN = 13344.6648; % 3000 lb in N
% 
% parfor i = 1:length(sweep_vinf)
% % for i = 1:length(sweep_vinf)
% 
%     vinf = sweep_vinf(i);
%     
%     ROTDIAM  = PROP_SWEEP(1,i).vecROTDIAM(1); % should be 1.524
%     ROTORRPM = PROP_SWEEP(1,i).vecROTORRPM(1); %should be 2250
%     rps      = ROTORRPM/60;
%     
%     % Calculate CT required to maintain S&L flight
%     % Calculate CL at VINF and S&L flight
%     S    = WING_SWEEP(i).valAREA; % ref. platform area
%     CL   = weightN./(0.5*rho*vinf.^2*S);
% 
%     % interpolate alpha to maintain steady level flight at VINF
%     % using wing only data
%     ALPHA = interp1([WING_SWEEP.vecCLv_AVG],[WING_SWEEP.vecVEHALPHA],CL,'linear','extrap');
% 
%     % LD = interp1(borer(:,1),borer(:,2),vinf*1.94384); % get L/D from Borer Data
%     LD = 14;
%     CD = CL./(LD); % Calculate CD with Borer L/D Data
%     D  = 0.5*rho*vinf.^2.*CD*S; % Calulate drag force in Newton
%     thrust(i)  = (D./cosd(ALPHA))/(2*N_prop); % Calculate Thrust force required from EACH PROP
%     CT = thrust(i)/(rho*rps^2*ROTDIAM^4);
%     
%     % use scatteredInterpolant to avoid meshgrid
%     propCT   = [PROP_SWEEP(:,i).vecCT_AVG]';
%     propColl = [PROP_SWEEP(:,i).vecCOLLECTIVE]';
%     propVINF = repmat(vinf,size(propCT));
%     
%     % meshgrids of 3d vinf,ct,collective data
%     propVINF = reshape(propVINF,[],length(unique(propVINF)));
%     propCT   = reshape(propCT,[],length(unique(propVINF)));
%     propColl = reshape(propColl,[],length(unique(propVINF)));
%     
%     [~,maxCTidx] = max(propCT);
%     propCT(maxCTidx+1:end) = nan;
%     
%     % idx = propColl<=0; % quick way to hack off the partially stalled propeller
%     % idx = ~isnan(propCT);
%     % F = scatteredInterpolant(propVINF(idx), propCT(idx), propColl(idx),'linear','nearest');
% 
%     vecCOLLECTIVE = interp1(propCT(~isnan(propCT)), propColl(~isnan(propCT)), CT, 'linear', 'extrap');
%     
%     vecCOLLECTIVE = vecCOLLECTIVE(~isnan(vecCOLLECTIVE));
%     CT = CT(~isnan(vecCOLLECTIVE));
%     vinf = vinf(~isnan(vecCOLLECTIVE));
%     
%     ITER_tmp.maxIter = 5;
%     ITER_tmp.numCase = length(vinf);
%     
%     ITER_tmp.Iteration = repmat((1:ITER_tmp.maxIter)',1,ITER_tmp.numCase);
%     ITER_tmp.CL        = nan(1, ITER_tmp.numCase);
%     ITER_tmp.CT        = nan(1, ITER_tmp.numCase);
%     ITER_tmp.CD        = nan(1, ITER_tmp.numCase);
%     ITER_tmp.AOA       = nan(1, ITER_tmp.numCase);
%     ITER_tmp.CLTV      = nan(1, ITER_tmp.numCase);
%     
%     seqALPHA = ALPHA;
%     
%     TRIMMED = false;
%     
%     
%     for n = 1:ITER_tmp.maxIter
%         
%         if n == 2
%             dCL1 = CL - ITER_tmp.CL(1,:);
%             dCT1 = CT - ITER_tmp.CT(1,:);
%             % New sets of AOA input for 2nd iteration in order to hit the targeted CL
%             seqALPHA = interp1([WING_SWEEP.vecCLv_AVG],[WING_SWEEP.vecVEHALPHA],CL + dCL1, 'linear', 'extrap');
%             % New sets of collective pitch input for 2nd iteration in order to hit the targeted CT
%             vecCOLLECTIVE = interp1(propCT(~isnan(propCT)), propColl(~isnan(propCT)), CT + dCT1, 'linear', 'extrap');
%         elseif n > 2
%             % New sets of AOA input for the next iteration in order to hit the targeted CL
%             seqALPHA = interp1(ITER_tmp.CL(1:n-1),ITER_tmp.AOA(1:n-1),CL,'linear','extrap');
%             
%             % New sets of collective pitch input input for the next iteration in order to hit the targeted CT
%             vecCOLLECTIVE = interp1(ITER_tmp.CT(1:n-1),ITER_tmp.CLTV(1:n-1),CT,'linear','extrap');
%             
%         end
%         
%         ITER_tmp.AOA(n)  = seqALPHA;
%         ITER_tmp.CLTV(n) = vecCOLLECTIVE;
%         
%         cd '../../'
%         vap_filename = ['Runs/J_COLE_OPTIMIZATION/',file, '/', file, '.vap'];
%         VAP_IN = [];
%         VAP_IN.vecVEHALPHA = seqALPHA;
%         VAP_IN.vecCOLLECTIVE = repmat(vecCOLLECTIVE, N_prop, 1);
%         VAP_IN.vecVEHVINF = vinf;
%         VAP_IN.valMAXTIME = 160;
%         VAP_IN.valSTARTFORCES = VAP_IN.valMAXTIME-20;
% %                         VAP_IN.valMAXTIME = 10
% %                         VAP_IN.valSTARTFORCES = 6
%         VAP_IN.valDELTIME = (1/60)/(sweep_rpm(i)/60);
%         OUTP = fcnVAP_MAIN(vap_filename, VAP_IN);
%         cd 'Runs/J_COLE_OPTIMIZATION/'
%         
%         % Write results
%         %         CL_star = OUTP.vecCL_AVG + (OUTP.vecCT_AVG
%         ITER_tmp.CL(n,:) = OUTP.vecCL_AVG + (2.*(sind(ALPHA).*(sum(OUTP.vecCT_AVG).*((OUTP.vecROTORRPM(1)/60).^2).*(OUTP.vecROTDIAM(1).^4).*rho))./(rho.*S.*vinf.^2));
%         ITER_tmp.CT(n,:) = nanmean(OUTP.vecCT_AVG);
%         ITER_tmp.CD(n,:) = [OUTP.vecCD_AVG];
%         %         CDtemp = [OUTP.vecCD_AVG];
%         %         CDtemp(isnan([OUTP.vecCLv])) = nan;
%         %         ITER.CD(n,:) = nanmean(CDtemp,1);
%         
% %         ITEROUTP(i,n).OUTP = OUTP;
%         
%         dCT = abs((ITER_tmp.CT(end) - CT)./CT);
%         dCL = abs((ITER_tmp.CL(end) - CL)./CL);
%         if dCT <= 0.02 && dCL <= 0.02
%             break;
%         end
%         
%     end
%     
%     ITER(i) = ITER_tmp;
%     ITEROUTP(i).OUTP = OUTP;
%     
%     ITER_tmp = struct;
%     
% end   
% save([file, '/', file, '_sweep.mat']);
load([file, '/', file, '_sweep.mat']);

hFig700 = figure(700);
clf(700)
rho = .957;
len = size(sweep_vinf,2);
for i = 1:len
   aoa(i,1) = ITER(i).AOA(end);
   power(i,1) = sum(2.*ITEROUTP(i).OUTP.vecCP_AVG).*((sweep_rpm(i)/60).^3).*(rotor_diam.^5).*rho;
end
plot(sweep_vinf, power, '-ok');
xlabel('Flight Speed (m/s)','FontSize',15);
ylabel('Trimmed Power Required (W)','FontSize',15);
grid minor
box on
axis tight

N_prop = 1
file = 'Design_1' % Case 3
rotor_diam = 1.594300;
rotor_rpm = 2276.93;
clear power
load([file, '/', file, '_sweep.mat']);
rho = .957;
len = size(sweep_vinf,2);
for i = 1:len
   aoa(i,1) = ITER(i).AOA(end);
   power(i,1) = sum(2.*ITEROUTP(i).OUTP.vecCP_AVG).*((sweep_rpm(i)/60).^3).*(rotor_diam.^5).*rho;
end
hold on
plot(sweep_vinf, power, '--*b');
hold off

    
