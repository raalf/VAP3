function [hFig3] = fcnGIF(flagVERBOSE, valTIMESTEP, valNELE, matDVE, matVLST, matCENTER, matFUSEGEOM,...
                    valWNELE, matWDVE, matWVLST, matWCENTER)

[hFig3] = fcnPLOTBODY(flagVERBOSE, valNELE, matDVE, matVLST, matCENTER, matFUSEGEOM);
[hFig3] = fcnPLOTWAKE(flagVERBOSE, hFig3, valWNELE, matWDVE, matWVLST, matWCENTER);
% view([33 22])
view([-45 20])

frame = getframe(hFig3);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);

% Write to the GIF File

if valTIMESTEP == 1
    imwrite(imind,cm,'GIF/output.gif','gif', 'Loopcount',inf);
else
    imwrite(imind,cm,'GIF/output.gif','gif','WriteMode','append');
end