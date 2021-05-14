%%
halewingspan = halewingpoints(:,end);
i = 1:4:size(halewingspan,1);
halewingspan = halewingspan(i);
halewingspan(end) = [];

halewinggamma(:,3) = -halewinggamma(:,3);
gamma = halewinggamma(:,3);
gamma = reshape(gamma,[3,64]);
for j = 1:64
for i = 2:-1:1
gamma(i+1,j) = gamma(i+1,j) - gamma(i,j);
end
end
gamma_dist = sum(gamma,1)';

span = 16;
% for i = 1:3
    figure(3)
    hold on
    plot(halewingspan,gamma_dist,'s','linewidth',1.5)
    grid on
    box on
    xlabel('Y Location (m)')
    ylabel('Circulation (m^2/s)')
% end

[ledves, ~, ~] = find(SURF.vecDVELE > 0);
[tedves, ~, ~] = find(SURF.vecDVETE > 0);
lepanels = SURF.vecDVEPANEL(ledves);
isCurWing = SURF.vecWINGTYPE(ledves) == 1;
idxdve = uint16(ledves(isCurWing));
idxpanel = lepanels(isCurWing);
m = INPU.vecM(idxpanel);
m = m(1);
% Matrix of how much we need to add to an index to get the next chordwise element
% It is done this way because n can be different for each panel. Unlike in the wake,
% we can't just add a constant value to get to the same spanwise location in the next
% row of elements
tempm = repmat(INPU.vecN(idxpanel), 1, m).*repmat([0:m-1],length(idxpanel~=0),1);
rows = repmat(idxdve,1,m) + uint16(tempm);

idx1 = SURF.vecDVELE == 1; %index of LE vectors (will be the same)

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

circ = A(rows)';
circ_dist = sum(A(rows)',1);

% for i = 1:3
    figure(3)
    hold on
    plot(SURF.matCENTER(1:15,2),circ_dist,'-','linewidth',1.5)
    grid on
    box on
    xlabel('Y Location (m)')
    ylabel('Circulation (m^2/s)')
% end