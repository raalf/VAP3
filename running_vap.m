clc
clear

filename = {['inputs/J_COLE_BASELINE_SYM.vap'];['inputs/J_COLE_BASELINE_SYM_CLOCKWISE.vap']};

for i = 1:2
    [vecCL(:,:,i), vecCDI(:,:,i), vecE(:,:,i), vecCT(:,:,i)] = fcnVAP_MAIN(char(filename(i)));
end
%%
hFig10 = figure(10);
plot(vecE(:,:,1),'-ok');
hold on
plot(vecE(:,:,2),'--xr');
grid minor
box on
axis tight