function [hFig3] = fcnGIF(flagVERBOSE, valTIMESTEP, valNELE, matDVE, matVLST, matCENTER, matFUSEGEOM,...
                    valWNELE, matWDVE, matWVLST, matWCENTER, vecWDVESURFACE, case_num)

[hFig3] = fcnPLOTBODY(flagVERBOSE, valNELE, matDVE, matVLST, matCENTER, []);
[hFig3] = fcnPLOTWAKE(flagVERBOSE, hFig3, valWNELE, matWDVE, matWVLST, matWCENTER, vecWDVESURFACE);
% view([33 22])
view([-45 20])

frame = getframe(hFig3);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);

% Write to the GIF File

gif_str = ['GIF/output_',num2str(case_num),'.gif'];

if valTIMESTEP == 1
    imwrite(imind,cm, gif_str,'gif', 'Loopcount',inf);
else
    imwrite(imind,cm, gif_str,'gif','WriteMode','append');
end