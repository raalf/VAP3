function [pop] = fcnPOPCREATION(nvars,seed_size,N_bendstiff,N_torstiff,N_elasticaxis,N_massaxis,ub,lb,A,b)

totalPopSize = seed_size;
pop = zeros(totalPopSize,nvars);

i = 1;

while i <= totalPopSize
    
    % Creates random stiffness values between minimum and 1,000,000 Nm2
    minEI = lb(1);
    EI = repmat(minEI + (1e6-minEI)*rand(1),1,N_bendstiff);
    
    minGJ = lb(N_bendstiff+1);
    GJ = repmat(minGJ + (1e6-minGJ)*rand(1),1,N_torstiff);
    
    % Places elastic and mass axis randomly between upper and lower bounds
    % of constraints
    minEA = lb(N_bendstiff+N_torstiff+1);
    maxEA = ub(N_bendstiff+N_torstiff+1);
    EA = repmat(minEA + (maxEA - minEA)*rand(1),1,N_elasticaxis);    
    
    minCG = lb(N_bendstiff+N_torstiff+N_elasticaxis+1);
    maxCG = ub(N_bendstiff+N_torstiff+N_elasticaxis+1);
    CG = repmat(minCG + (maxCG - minCG)*rand(1),1,N_massaxis);
    
    pop(i,:) = [EI,GJ,EA,CG];
    
    if all(A*pop(i,:)' - b <= 0)
        i = i + 1;
    else
        pop(i,:) = pop(i,:).*0;
    end
    
end