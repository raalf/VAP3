function [] = fcnFLIGHTDYNPRINTOUT(flagPRINT, valTIMESTEP, valVEHICLES, vecCL, vecCDI, matDEFGLOB, vecCG, perturb, i)
if flagPRINT == 1 && valTIMESTEP == 1
    rotor_count = 1;
    
    header1 = [];
    header2 = sprintf('TIMESTEP');
    header3 = ['------------'];
    for j = 1:valVEHICLES
        header1 = [header1, sprintf('\t\t\tVEHICLE %d\t',j)];
        header2 = [header2, [sprintf('\tCL'), sprintf('\t\t\tCDI'), sprintf('\t\t\tTip Def. (m)'), sprintf('\t\tVeh. Z Loc. (m)'), sprintf('\t\tPitch Angle (deg)')]];
        header3 = [header3, ['-------------------------------------------------------------------------']];      

    end
    
    fprintf([header1 '\n']);
    fprintf([header2 '\n']);
    fprintf([header3 '\n']);
end

txtout = ['\t', sprintf('%4d',valTIMESTEP)];

for j = 1:valVEHICLES
%     txtout = [txtout, sprintf('\t%0.4f',vecCL(valTIMESTEP,j,i)), sprintf('\t\t%0.4f',vecCDI(valTIMESTEP,j,i)), sprintf('\t\t%0.4f',matDEFGLOB(valTIMESTEP,end)), sprintf('\t\t\t\t%0.4f',vecCG(3)), sprintf('\t\t\t\t%0.4f', (180/pi)*perturb(valTIMESTEP,4))];

    txtout = [txtout, sprintf('\t%0.4f',vecCL(valTIMESTEP,j,i)), sprintf('\t\t%0.4f',vecCDI(valTIMESTEP,j,i)), sprintf('\t\t%0.4f',0), sprintf('\t\t\t\t%0.4f',vecCG(3)), sprintf('\t\t\t\t%0.4f', (180/pi)*perturb(valTIMESTEP,4))];

end

fprintf([txtout '\n']);