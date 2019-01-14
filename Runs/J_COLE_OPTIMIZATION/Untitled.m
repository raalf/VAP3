clc
clear

A = dlmread('opthistory.txt');

idx = A > 1.5e+4;

out = [];
row = 0;
column = 0;
for i = 1:size(idx,1)
    for j = 1:size(idx,2)
        if idx(i,j) == true
            row = row + 1;
            column = 1;
        end
        out(row,column) = A(i,j);
        column = column + 1;
    end
end

if (any(any(out(:,24:end)))) == false
    disp('All good')
end
out = out(:,1:23);

case1 = [];
case3 = [];
case5 = [];

for i = 1:size(out,1)
    len = sum(any(out(i,:),1));
    
    if len > 17
        case3(end+1,:) = out(i,:);
    elseif out(i,17) == 1 && out(i,13) == 160
        case1(end+1,:) = out(i,1:17);
    else
        case5(end+1,:) = out(i,1:17);
    end
    
end

dlmwrite('case1.txt',case1,' ');
dlmwrite('case3.txt',case3,' ');
dlmwrite('case5.txt',case5,' ');

% row = 1;
% while true
%     if any(A(row,17:end))
%        A(row + 2:end + 1,:) = A(row + 1:end,:);
%        A(row + 1,:) = [A(row,17:end) zeros(1,16)];
%        A(row,:) = [A(row,1:16) zeros(1,(t1-1)*16)];
%        row = row + 1;
%     else
%        row = row + 1;
%     end
%
%     if row > size(A,1)
%        break;
%     end
% end
%
% out = A(:,1:16);
% out(~any(out,2),:) = [];