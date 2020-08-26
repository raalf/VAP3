function [OUTP, SURF] = fcnSURFACEDIST(COND, INPU, SURF, OUTP, FLAG, valTIMESTEP)

% Function to calculate the force distributions on each lifting surface.
% This takes several components/logic from fcnDVENFORCE
ledves = find(SURF.vecDVELE == 1);

[matROWS] = fcnDVEROW(ledves, SURF, INPU, unique(SURF.vecWINGTYPE));

matROWS = cellfun(@transpose, matROWS, 'UniformOutput', false); % Matrix of DVE numbers and their location on the wing (vecN x vecM large)

for i = 1:1
    
    liftfree_edge = [];
    
    nfree = [];
    liftfree = [];
    temp_lift = [];
    
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
    
    % Loop through each element and get lift at several spanwise points based
    % on local A, B and C
    SURF.nspnele = 10; % Number of local spanwise elements to use for force distribution
    SURF.span_glob = [];
    int = [];
    matROWSt = matROWS{i}';
    for j = min(min(matROWS{i})):max(max(matROWS{i}))
        span_loc = linspace(-SURF.vecDVEHVSPN(j)*0.95,SURF.vecDVEHVSPN(j)*0.95,SURF.nspnele)'; % Span location on DVE in local coordinates. Excluding DVE edges to avoid overlapping points
        span_loc_temp = span_loc + SURF.vecDVEHVSPN(j);
        span_glob = fcnSTARGLOB([zeros(length(span_loc),1),span_loc_temp,zeros(length(span_loc),1)], SURF.vecDVEROLL(j), SURF.vecDVEPITCH(j), SURF.vecDVEYAW(j));
%         SURF.span_glob = [SURF.span_glob; linspace(SURF.matVLST(SURF.matDVE(matROWSt(j),1),2),SURF.matVLST(SURF.matDVE(matROWSt(j),2),2),SURF.nspnele)'];
        SURF.span_glob = [SURF.span_glob; SURF.matVLST(SURF.matDVE(matROWSt(j),1),:)+span_glob];
        SURF.Gamma(:,j,valTIMESTEP) = A(j) + B(j).*span_loc + C(j).*span_loc.*span_loc;
        nfree(:,j-min(min(matROWS{i}))+1) = (A(j) + B(j).*span_loc + C(j).*span_loc.*span_loc).*repmat(uxs(j),SURF.nspnele,1).*(span_loc(end)-span_loc(1))/(SURF.nspnele-1);
%         nfree(:,j-min(min(matROWS{i}))+1) = (2*A(j).*SURF.vecDVEHVSPN(j) + (2/3)*C(j).*SURF.vecDVEHVSPN(j).*SURF.vecDVEHVSPN(j).*SURF.vecDVEHVSPN(j)).*uxs(j);
        
        if valTIMESTEP > 1
            
%             nfree(:,j-min(min(matROWS{i}))+1) = nfree(:,j-min(min(matROWS{i}))+1) + 0*2*SURF.vecDVEHVCRD(j).*(SURF.Gamma(:,j,valTIMESTEP) - SURF.Gamma(:,j,valTIMESTEP-1))./COND.valDELTIME;
            
        end
        
%         en_loc = repmat(en(j-min(min(matROWS{i}))+1,:),nspnele,1); % Make local normal direction that is nspnele long to match the nspnele points for the force dist.
%         liftfree(:,j) = nfree(:,j-min(min(matROWS{i}))+1).*sqrt(en_loc(:,1).*en_loc(:,1) + en_loc(:,3).*en_loc(:,3)); %does this work with beta?
%         liftfree(en_loc(1:nspnele,3)<0) = -liftfree(en_loc(1:nspnele,3)<0);
    end
    
    SURF.span_glob = uniquetol(SURF.span_glob(:,2),1e-3);
    % Rearrange lift distribution to 3D matrix
    for k = 1:size(matROWS{i},1)
%         temp_lift(:,k,:) = liftfree(:,matROWS{i}(k,:));
    end
    
%     liftfree = reshape(sum(temp_lift,2),[nspnele*size(matROWS{i},2) 1]);

%     OUTP.vecLIFTDIST{i}(:,valTIMESTEP) = liftfree.*COND.valDENSITY; % Sum lift contribution at each lifting line
%     SURF.vecLIFTDISTCOORD{i} = linspace(0,SURF.matVLST(SURF.matDVE(matROWS{i}(1,end),2),2),nspnele*size(matROWS{i},2))';

    % Calculate normal force distribution at each chordwise lifting line
    SURF.ndist_chord{i} = [];
    for k = 1:size(matROWS{i},1)
        SURF.ndist_chord{i} = [SURF.ndist_chord{i}; reshape(nfree(:,matROWS{i}(k,:)),[1 size(nfree(:,matROWS{i}(k,:)),1)*size(nfree(:,matROWS{i}(k,:)),2)])];
    end
    
    SURF.ndist{i} = sum(SURF.ndist_chord{i},1);
    
%     SURF.ndist{i} = flip(uniquetol(SURF.ndist{i},1e-8));
%     SURF.span_glob = unique(SURF.span_glob);
    
    nfree = reshape(nfree,[SURF.nspnele*max(max(matROWS{i})),1]);
    SURF.vecLIFTDISTCOORD{i} = SURF.matCENTER(matROWS{i}(1,1:end),2);
%     SURF.vecLIFTDISTCOORD{i} = linspace(SURF.vecDVEHVSPN(1),INPU.vecSPAN/2,nspnele*size(matROWS{i},2))';
    
    if FLAG.STRUCTURE == 1
    % Calculate pitching moment about elastic axis due to lift
        q_inf = 0.5*COND.valDENSITY*COND.vecVEHVINF*COND.vecVEHVINF;
        delX = SURF.matAEROCNTR - SURF.matEA;

%         OUTP.vecMOMDIST = cross(delX,[zeros(size(SURF.vecSTRUCTSPNDIST,1),2), interp1(SURF.vecLIFTDISTCOORD{1},OUTP.vecLIFTDIST{1}(:,end),SURF.vecSTRUCTSPNDIST,'linear','extrap')]);
%         OUTP.vecMOMDIST = OUTP.vecMOMDIST + [zeros(size(SURF.vecSTRUCTSPNDIST,1),1), q_inf*interp1(SURF.vecSPANDIST,SURF.vecMAC,SURF.vecSTRUCTSPNDIST,'linear','extrap').*...
%             interp1(SURF.vecSPANDIST,SURF.vecMAC,SURF.vecSTRUCTSPNDIST,'linear','extrap').*interp1(SURF.vecSPANDIST,OUTP.vecCMDIST(SURF.vecWINGTYPE(SURF.vecDVELE == 1) == 1),SURF.vecSTRUCTSPNDIST,'linear','extrap'), zeros(size(SURF.vecSTRUCTSPNDIST,1),1)];
    end
       
end

end

