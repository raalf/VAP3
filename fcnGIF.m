function [hFig3] = fcnGIF(valTIMESTEP, FLAG, SURF, VISC, WAKE, COND, INPU, case_num)

if exist('GIF/','dir') ~= 7; mkdir('GIF'); end

fcnPLOTPKG(valTIMESTEP, FLAG, SURF, VISC, WAKE, COND, INPU)
% view([33 22])
% view([-45 20])
view([-20 20])
set(gcf, 'Position', get(0, 'Screensize'));
hFig3 = figure(3);
zlim([-30 30]);
ylim([-30 30]);
xlim([-120 15]);
frame = getframe(hFig3);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);

% Write to the GIF File

gif_str = ['GIF/output_',num2str(case_num),'.gif'];

if valTIMESTEP == 1
    imwrite(imind,cm, gif_str,'gif', 'Loopcount',inf);
else
    imwrite(imind,cm, gif_str,'gif','WriteMode','append', 'DelayTime', 0.05);
end