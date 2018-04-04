function [] = fcnPRINTOUT(flagPRINT, valTIMESTEP, valVEHICLES, vecCL, vecCDI, vecCT, vecROTORJ, vecROTORVEH,i)
if flagPRINT == 1 && valTIMESTEP == 1
    rotor_count = 1;
    
    header1 = [];
    header2 = sprintf('TIMESTEP');
    header3 = ['------------'];
    for j = 1:valVEHICLES
        header1 = [header1, sprintf('\t\t\tVEHICLE %d\t',j)];
        header2 = [header2, [sprintf('\tCL'), sprintf('\t\t\tCDI')]];
        header3 = [header3, ['-------------------------']];
        
        for jj = 1:length(nonzeros(vecROTORVEH == j))
            header1 = [header1, sprintf('\t\tJ = %0.3f', vecROTORJ(i,jj))];
            header2 = [header2, sprintf('\t\t\tCT_%d', rotor_count)];
            rotor_count = rotor_count + 1;
            header3 = [header3, ['---------------']];
        end
        

    end
    
    fprintf([header1 '\n']);
    fprintf([header2 '\n']);
    fprintf([header3 '\n']);
end

txtout = ['\t', sprintf('%4d',valTIMESTEP)];

for j = 1:valVEHICLES
    txtout = [txtout, sprintf('\t%0.4f',vecCL(valTIMESTEP,j,i)), sprintf('\t\t%0.4f',vecCDI(valTIMESTEP,j,i))];
%     , sprintf('\t\t%0.4f',vecCT(valTIMESTEP,j))

        for jj = 1:length(nonzeros(vecROTORVEH == j))

            txtout = [txtout sprintf('\t\t%0.4f',vecCT(valTIMESTEP,jj,i))];
            
        end

end

fprintf([txtout '\n']);