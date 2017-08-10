if flagPRINT == 1 && valTIMESTEP == 1
    
    header1 = ['             '];
    header2 = [' ',sprintf('TIMESTEP'),'    '];
    header3 = ['------------'];
    for j = 1:valVEHICLES
        header1 = [header1,sprintf('VEHICLE %d',j),'                '];
        header2 = [header2,[sprintf('CL'),'          ',sprintf('CDI'),'          ']];
        header3 = [header3,['-------------------------']];
    end
    
    disp(header1);
    disp(header2);
    disp(header3);
end

txtout = ['  ', sprintf('%4d',valTIMESTEP),'       '];
for j = 1:valVEHICLES
    txtout = [txtout, sprintf('%0.5f',vecCL(valTIMESTEP,j)), '     ', sprintf('%0.5f',vecCDI(valTIMESTEP,j)), '      '];
end
disp(txtout)