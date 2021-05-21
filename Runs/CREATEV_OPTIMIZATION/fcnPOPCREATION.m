function [pop] = fcnPOPCREATION(nvars,seed_size,N_bendstiff,N_torstiff,N_elasticaxis,N_massaxis,ub,lb,A,b)

totalPopSize = seed_size;
pop = zeros(totalPopSize,nvars);

i = 1;

while i <= totalPopSize
    
    % Creates random stiffness values between minimum and 1,000,000 Nm2
    minEI = lb(1);
    EI(1) = minEI + (1e6-minEI)*rand(1);
    
    % Making sure EI follows relative change constraints
    for j = 2:N_bendstiff
        EI(j) = ((1.05-0.95)*rand(1)+0.95)*EI(j-1);
    end
    
    minGJ = lb(N_bendstiff+1);
    GJ(1) = minGJ + (1e6-minGJ)*rand(1);
    
    % Making sure GJ follows relative change constraints
    for j = 2:N_torstiff
        GJ(j) = ((1.05-0.95)*rand(1)+0.95)*GJ(j-1);
    end
    
    % Places elastic and mass axis randomly between upper and lower bounds
    % of constraints
    minEA = lb(N_bendstiff+N_torstiff+1);
    maxEA = ub(N_bendstiff+N_torstiff+1);
    EA(1) = minEA + (maxEA - minEA)*rand(1);
    
    % Making sure EA follows relative change constraints
    for j = 2:N_elasticaxis
        EA(j) = ((1.025-0.975)*rand(1)+0.975)*EA(j-1);
    end
    
    minCG = lb(N_bendstiff+N_torstiff+N_elasticaxis+1);
    maxCG = ub(N_bendstiff+N_torstiff+N_elasticaxis+1);
    CG(1) = minCG + (maxCG - minCG)*rand(1);
    
    % Making sure CG follows relative change constraints
    for j = 2:N_elasticaxis
        CG(j) = ((1.025-0.975)*rand(1)+0.975)*CG(j-1);
    end
    
    pop(i,:) = [EI,GJ,EA,CG];
    
    if all(A*pop(i,:)' - b <= 0)
        i = i + 1;
    else
        pop(i,:) = pop(i,:).*0;
    end
    
end