clear all;
close all;

LL = 2:1:50;
LEN = size(LL, 2);
parfor ii = 1 : LEN
    L = LL(ii);
%     [Err_GAF0, ERR_S0(ii)] = F_GS_RAF_L(L, 0.1);
%     [Err_GAF1, ERR_S1(ii)] = F_GS_RAF_L(L, 0.5);
    [Err_GAF2, ERR_S2(ii)] = F_GS_RAF_L(L, 3);
%     [Err_GAF3, ERR_S3(ii)] = F_GS_RAF_L(L, 5);
%     [Err_GAF4, ERR_S4(ii)] = F_GS_RAF_L(L, 10);
    disp(ii)
end

figure
hold on
box on
wid = 3;
MarkerSize = 6;
% plot(ERR_S0, '-r', 'LineWidth', wid, 'Marker', 's', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'r')
% plot(ERR_S1, '-g', 'LineWidth', wid, 'Marker', 'd', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'g')
plot(ERR_S2, '-b', 'LineWidth', wid, 'Marker', 'h', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'b')
% plot(ERR_S3, '-c', 'LineWidth', wid, 'Marker', 'p', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'c')
% plot(ERR_S4, '-m', 'LineWidth', wid, 'Marker', '*', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'm')
% h=legend('GS-RAF(\lambda=0.995, \xi=0.1)','GS-RAF(\lambda=0.995, \xi=0.5)','GS-RAF(\lambda=0.995, \xi=2)', 'GS-RAF(\lambda=0.995, \xi=5)', 'GS-RAF(\lambda=0.995, \xi=10)');
% set(h,'FontName','Times New Roman','FontSize',24,'FontWeight','normal');
h=legend('GS-RAF(\lambda=0.995, \xi=3)');
set(h,'FontName','Times New Roman','FontSize',24,'FontWeight','normal');

set(gca,'fontsize',24);
set(gca,'fontsize',24);
xlabel('The node number L','FontSize',24);
ylabel('MSDs(dB)','FontSize',24);

