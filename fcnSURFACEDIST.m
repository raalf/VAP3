function [OUTP, SURF] = fcnSURFACEDIST(COND, INPU, SURF, OUTP, valTIMESTEP)

% Function to calculate the force distributions on each lifting surface.
% This takes several components/logic from fcnDVENFORCE
ledves = find(SURF.vecDVELE == 1);

[matROWS] = fcnDVEROW(ledves, SURF, INPU);

matROWS = cellfun(@transpose, matROWS, 'UniformOutput', false); % Matrix of DVE numbers and their location on the wing (vecN x vecM large)

for i = 1:max(SURF.vecDVEWING)
    
    liftfree_edge = [];
    
    idx1 = SURF.vecDVELE == 1; %index of LE vectors (will be the same)
    % vector of bound vortex along leading edge for LE row
    % tempa(idx1,1) = tan(SURF.vecDVELESWP(idx1));
    s = zeros(SURF.valNELE,3);
    s(idx1,:) =( SURF.matVLST(SURF.matDVE(idx1,2),:) -SURF.matVLST(SURF.matDVE(idx1,1),:) )  ./ repmat((SURF.vecDVEHVSPN(idx1).*2),1,3); %?? why is the S vector non-dim. to the span?
    
    % vector along midchord for remaining elements
    if any(idx1 == 0)
        tempa = zeros(SURF.valNELE,1);
        tempa(idx1==0,1) = tan(SURF.vecDVEMCSWP(idx1==0));
        tempa= [tempa ones(SURF.valNELE,1)  zeros(SURF.valNELE,1)];
        
        % global vectors to find ind. velocitiess
        s(idx1 ==0,:)= fcnSTARGLOB(tempa(idx1==0,:), SURF.vecDVEROLL(idx1==0), SURF.vecDVEPITCH(idx1==0), SURF.vecDVEYAW(idx1==0));
        
    end
    
    tempb = cross(SURF.matUINF, s, 2);
    
    uxs = sqrt(sum(abs(tempb).^2,2));
    
    en = tempb.*repmat((1./uxs),1,3);
    
    % if first row, A=A, B=B, C=C
    A = zeros(1,SURF.valNELE);
    B = zeros(1,SURF.valNELE);
    C = zeros(1,SURF.valNELE);
    
    A(idx1) = SURF.matCOEFF(idx1,1);
    B(idx1) = SURF.matCOEFF(idx1,2);
    C(idx1) = SURF.matCOEFF(idx1,3);
    % if any other row, A= A-Aupstream, B= B-Bupstream, C= C-Cupstream
    
    idx2 = SURF.vecDVELE == 1; %idx2 since we need to do this even for triangles
    dvenum = find(idx2==0); %dvenum in question
    idxf = SURF.matADJE((ismember(SURF.matADJE(:,1), dvenum) & SURF.matADJE(:,2) == 1),3); %upstream dve num
    A(idx2 ==0) = (SURF.matCOEFF(idx2==0,1)-SURF.matCOEFF(idxf,1));
    B(idx2 ==0) = (SURF.matCOEFF(idx2==0,2)-SURF.matCOEFF(idxf,2));
    C(idx2 ==0) = (SURF.matCOEFF(idx2==0,3)-SURF.matCOEFF(idxf,3));
    
    % Loop through each element and get lift at 100 spanwise points based
    % on local A, B and C
    nspnele = 100; % Number of local spanwise elements to use for force distribution
    for j = 1:max(max(matROWS{i}))
        span_loc = linspace(-SURF.vecDVEHVSPN(j),SURF.vecDVEHVSPN(j),nspnele)';
        SURF.Gamma(:,j,valTIMESTEP) = A(j) + B(j).*span_loc + C(j).*span_loc.*span_loc;
        nfree(:,j) = (A(j) + B(j).*span_loc + C(j).*span_loc.*span_loc).*uxs(j);
        
        if valTIMESTEP > 1
            
            nfree(:,j) = nfree(:,j) + 2*SURF.vecDVEHVCRD(j).*(SURF.Gamma(:,j,valTIMESTEP) - SURF.Gamma(:,j,valTIMESTEP-1))./COND.valDELTIME;
            
        end
        
        en_loc = repmat(en(j,:),nspnele,1); % Make local normal direction that is nspnele long to match the nspnele points for the force dist.
        liftfree(:,j) = nfree(:,j).*sqrt(en_loc(:,1).*en_loc(:,1) + en_loc(:,3).*en_loc(:,3)); %does this work with beta?
        liftfree(en_loc(j,3)<0) = -liftfree(en_loc(j,3)<0);
    end
    
    % Rearrange lift distribution to 3D matrix
    for k = 1:size(matROWS{i},1)
        temp_lift(:,k,:) = liftfree(:,matROWS{i}(k,:));
    end
    
    liftfree = reshape(sum(temp_lift,2),[nspnele*size(matROWS{i},2) 1]);
    
%     liftfree = nfree.*sqrt(en(matROWS{i},1).*en(matROWS{i},1) + en(matROWS{i},3).*en(matROWS{i},3)); %does this work with beta?
%         
%     nfree = ((A(:,matROWS{i}') - B(:,matROWS{i}').*SURF.vecDVEHVSPN(matROWS{i},:)' + C(:,matROWS{i}').*SURF.vecDVEHVSPN(matROWS{i},:)'.*SURF.vecDVEHVSPN(matROWS{i},:)').*uxs(matROWS{i},:)')'; % Normal force dist./density at each DVE left edge
%     nfree_tip = ((A(matROWS{i}(:,end)) + B(matROWS{i}(:,end)).*SURF.vecDVEHVSPN(matROWS{i}(:,end))' + C(matROWS{i}(:,end)).*SURF.vecDVEHVSPN(matROWS{i}(:,end))'.*SURF.vecDVEHVSPN(matROWS{i}(:,end))').*uxs(matROWS{i}(:,end))')'; % Normal force/density at wing tip
%     
% %     nfree_mid = (A(:,matROWS{i}').*uxs(matROWS{i},:)')';
%     
%     % Convert normal force to lift direction for all DVEs
%     liftfree = nfree.*sqrt(en(matROWS{i},1).*en(matROWS{i},1) + en(matROWS{i},3).*en(matROWS{i},3)); %does this work with beta?
%     liftfree(en(matROWS{i},3)<0) = -liftfree(en(matROWS{i},3)<0);
%     
% %     liftfree_mid = nfree_mid.*sqrt(en(matROWS{i},1).*en(matROWS{i},1) + en(matROWS{i},3).*en(matROWS{i},3)); %does this work with beta?
% %     liftfree_mid(en(matROWS{i},3)<0) = -liftfree_mid(en(matROWS{i},3)<0);
%     
%     % Convert normal force to lift direction for tip elements
%     liftfree_tip = nfree_tip.*sqrt(en(matROWS{i}(:,end),1).*en(matROWS{i}(:,end),1) + en(matROWS{i}(:,end),3).*en(matROWS{i}(:,end),3));
%     liftfree_tip(en(matROWS{i}(:,end),3)<0) = -liftfree_tip(en(matROWS{i}(:,end),3)<0);
%     
%     % Rearrange lift distribution to be same size as matROWS. This shows
%     % the lift distribution at each chordwise lifting-line
%     if i > 1
%         idx_multi = matROWS{i} - min(min(matROWS{i})) + 1; % Special case to adjust index for multiple wings. Seems sketchy but also seems to work so nothing to see here
%         liftfree_edge(1:size(matROWS{i},1),1:size(matROWS{i},2)) = liftfree(idx_multi);
%         
% %         liftfree_mid(1:size(matROWS{i},1),1:size(matROWS{i},2)) = liftfree_mid(idx_multi);
%     else
%         liftfree_edge(1:size(matROWS{i},1),1:size(matROWS{i},2)) = liftfree(matROWS{i});
%         
% %         liftfree_mid(1:size(matROWS{i},1),1:size(matROWS{i},2)) = liftfree_mid(matROWS{i});
%     end
%     
%     % Add on tip lift to distribution (obviously should be very very close
%     % to 0)
%     liftfree_edge = [liftfree_edge, liftfree_tip];
%     
%     liftfree_edge = liftfree_edge.*COND.valDENSITY; % Convert from force/density to actual force
%     
% %     liftfree_mid = liftfree_mid.*COND.valDENSITY; % Convert from force/density to actual force
%     
%     idxtol = abs(liftfree_edge) < 1e-10; % Check if lift is below tolerance. This will essentially set tip forces to 0 instead of a small positive or negative value. This helps for plotting
%     liftfree_edge(idxtol) = 0;
    
%     OUTP.vecLIFTDIST{i}(:,valTIMESTEP) = sum(liftfree_edge,1); % Sum lift contribution at each lifting line
    OUTP.vecLIFTDIST{i}(:,valTIMESTEP) = liftfree.*COND.valDENSITY; % Sum lift contribution at each lifting line
    SURF.vecLIFTDISTCOORD{i} = linspace(0,INPU.vecSPAN(i)/2,nspnele*size(matROWS{i},2))';
%     SURF.vecLIFTDISTCOORD{i} = [SURF.matVLST(SURF.matDVE(matROWS{i}(1,:),1),2); SURF.matVLST(SURF.matDVE(matROWS{i}(1,end),2),2)]; % Spanwise coordinate for plotting lift distribution. This probably won't work for rotors since the y coordinate will rotate. Need to think of something else to try
    
%     OUTP.vecLIFTDIST_MID{i}(:,valTIMESTEP) = sum(liftfree_mid,1); % Sum lift contribution at each lifting line
    
end

end

