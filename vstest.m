clc
clear

le = [0 1 0; 0 2 0];

gran = 0.5;
x = -5:gran:5;
y =-5:gran:5;
z = -0.5:gran:0.5;
[X,Y,Z] = meshgrid(x,y,z);
fp_0 = [X(:) Y(:) Z(:)];
len = size(fp_0, 1);

hspan = le(2,2) - le(1,2);
hchord = 1;
phi = 0;
k = 0.1;

hspan = repmat(hspan,len,1);
hchord = repmat(hchord,len,1);
phi = repmat(phi,len,1);
k = repmat(k,len,1);

matCOEFF = [0 1 0];

fp1 = fp_0 - mean(le,1);
[a, b, c] = fcnVSIND(hspan, hchord, phi, fp1, k);
D = [zeros(len,3) b c];
D = reshape(reshape(D', 1, 9, []), 3, 3, len);
v = permute(sum(D.*repmat(reshape(matCOEFF',1,3,[]),3,1,1),2),[2 1 3]);
v = reshape(permute(v,[3 1 2]),[],3,1)./(-4*pi);

hFig2 = figure(2);
clf(2);

plot3(le(:,1), le(:,2), le(:,3), '-or','linewidth',2);
hold on
quiver3(fp_0(:,1), fp_0(:,2), fp_0(:,3), v(:,1), v(:,2), v(:,3),'r')


xlabel('X-Dir','FontSize',15);
ylabel('Y-Dir','FontSize',15);
zlabel('Z-Dir','FontSize',15);

axis equal
grid on
box on

%%
le2 = [-le(:,2) le(:,1) le(:,3)] + [2 2 0];

fp2 = fp_0 - mean(le2,1);
fp2 = [-fp2(:,2) fp2(:,1) fp2(:,3)];
[a, b, c] = fcnVSIND(hspan, hchord, phi, fp2, k);
b = [b(:,2) -b(:,1) b(:,3)];
c = [c(:,2) -c(:,1) c(:,3)];

D = [zeros(len,3) b c];
D = reshape(reshape(D', 1, 9, []), 3, 3, len);
v = permute(sum(D.*repmat(reshape(matCOEFF',1,3,[]),3,1,1),2),[2 1 3]);
v = reshape(permute(v,[3 1 2]),[],3,1)./(-4*pi);

quiver3(fp_0(:,1), fp_0(:,2), fp_0(:,3), v(:,1), v(:,2), v(:,3),'g')
plot3(le2(:,1), le2(:,2), le2(:,3), '-og','linewidth',2);




hold off
