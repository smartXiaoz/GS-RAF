clear all;
close all;

i = 1;
j = 1;
for ii = 2 : 1 : 10
    for jj = 0.1 : 0.1 : 5
        Err_GAF = F_GS_RAF(ii, jj);
        ERR_S(i, j) = mean(Err_GAF(:,1000:end));
        j = j + 1;
    end
    i = i + 1;
    j = 1;
    disp(ii)
end

jj = 2 : 1 : 10;
ii = 0.1 : 0.2 : 5;
[x,y] = meshgrid(ii,jj);
surf(x',y',ERR_S');
colormap jet
set(gca,'fontsize',24);
xlabel('\xi','FontSize',24);
ylabel('L','FontSize',24);
zlabel('MSD','FontSize',24);