%% Í¼ÂË²¨RLS
clear all; close all;

%% MEE-AF-V10
L = 2500;
p = 10;
q = 10;
lamda1 = 0.995;%RLS
lamda2 = 0.99;%RMC
lamda3 = 0.995;%RMEE
lamda4 = 0.995;%GSAF
sigma1= 2;%%RMC
sigma2 = 0.5;%%RMEE

tic
% addpath(toolboxdir('stats'))
% addpath(genpath('Statistics and Machine Learning Toolbox'))
A_g = ones(q);
D_g = diag(sum(A_g));
L_g = D_g - A_g - (q-1)*eye(q);
LL_g = D_g - A_g;
for mm = 1:50
    vv = randn(1,L) * 0.1;
%     vv = rand(1,L) -0.5;
    v1=randn(1,L)*0.1; v2=randn(1,L)*10; v3=randn(1,L)*10;
    v4=randn(1,L)*0.1+0.5;
    v5=randn(1,L)*0.1-0.5;
    rp=rand(1,L);
    %    vv = v1 + (rp>0.95).*v2;
%     vv = (rp<=0.9).*v1 + (rp>0.9 & rp<=0.95).*v2 + (rp>0.95).*v3;
    vv = (rp<=0.98).*v1 + (rp>0.98).*v3;
%     vv = (rp<=0.48).*v4 + (rp>0.52).*v5 + (rp>0.48 & rp<0.52).*v3;

    %      vv = exp(-vv.^2/2/sigma1^2).*vv;
%     vv=randn(1,L)*0.1;
%     vv = stable(1.5,0,1,0,p,L);
%      vv = alpha_stable_noise(1.1,0.1,0,0,L);
%    wo1 = randn(p,1);wo2 = randn(p,1);wo3 = randn(p,1);
      wo1 = [0.5,0.4,0.3, 0.2,0.1, 0.1,0.2,0.3, 0.4,0.5]';
%     wo = [ kron(wo1, ones(1,L/3)) kron(wo2, ones(1,L/3)) kron(wo3, ones(1,L/3))];
     wo = [ kron(wo1, ones(1,L)) ];
  %   wo = [ kron(randn(p,1), ones(1,L/2)),  kron(randn(p,1), ones(1,L/2)) ];
    uu = randn(p,L);
%     path(path,'e:\work\speech enhancement\dbs');
%     [s1,FS,NBITS]=wavread('sp02.wav');
%     s1 = s1*1;
%     
%     for ii = 1 : p
%         uu(ii, :) = s1(3000+ii : 3000+ii+L-1)*10;
%     end
    for ii = 1 : L
        dd(ii) = wo(:,ii)' * uu(:,ii) + vv(ii);
    end
    
    % dd = wo' * uu + vv;
    
    w_LMS = zeros(p,1);   
    w_M_RLS1 =w_LMS;%%MEE
    w_M_RLS2 =w_LMS;%%MEE
    w_M_RLS3 =w_LMS;%%MEE
    w_M_RLS4 =w_LMS;%%MEE
    w_M_RLS5 =w_LMS;%%MEE
    w_M_RLS6 =w_LMS;%%MEE
    w_M_RLS7 =w_LMS;%%MEE
    w_RLS = w_LMS;
    w_C_RLS = w_LMS;
   
    
    %% RLS MEE, algorithm 7
    m = p;
%     L = q;
    PL = eye(m,m)*0.1; 
    eL = zeros(q, 1);
    uL = zeros(m,q);
    PHI = zeros(q, q);
    for ii = q : L
        Err_RLS_MEE7(mm,ii) = (wo(:,ii)  - w_M_RLS7)' * (wo(:,ii)  - w_M_RLS7);
        for jj = 1 : q
            eL(jj) = dd(ii - q + jj) - w_M_RLS7' * uu(: , ii - q + jj);
            uL(:,jj) = uu(: , ii - q + jj);
        end              
        L_1 = L_g;
        L_2 = L_g;
        L_3 = L_g;
         L_4 = L_g;
        %% ÍØÆË±ä»¯
        u = eL / median(abs(eL - median(eL)));
        for jj = 1 : q
            if abs(u(jj)) > 1
                L_1(jj,:)=0;
                L_1(:,jj)=0;
            end
            if abs(u(jj)) > 3
                L_2(jj,:)=0;
                L_2(:,jj)=0;
            end
             
            if abs(u(jj)) > 5
                L_3(jj,:)=0;
                L_3(:,jj)=0;
            end
            if abs(u(jj)) > 10
                L_4(jj,:)=0;
                L_4(:,jj)=0;
            end
        end      
        for jj = 1 : q          
             L_1(jj,jj)=-sum(L_1(jj,:));
              L_2(jj,jj)=-sum(L_2(jj,:));
               L_3(jj,jj)=-sum(L_3(jj,:));
               L_4(jj,jj)=-sum(L_4(jj,:));
        end   
       
%         Smooth_pre(mm,ii) = u'*LL_g*u;
%         Smooth_pos(mm,ii) = u'*L_1*u;
        PHI = L_1  + 1/1000*eye(q);
        Smooth_pre(mm,ii) = eL'*LL_g*eL;
        Smooth_pos1(mm,ii) = eL'*L_1*eL;
        Smooth_pos2(mm,ii) = eL'*L_2*eL;
        Smooth_pos3(mm,ii) = eL'*L_3*eL;
        Smooth_pos4(mm,ii) = eL'*L_4*eL;
        KL =PL * uL * inv ( lamda4*inv(PHI)+uL'*PL*uL);
        PL = 1/lamda4 * (PL - KL * uL' * PL);
        w_M_RLS7 = w_M_RLS7 + KL * eL;        
    end
    
    disp(mm)
end
toc
wid = 1.5;
figure,hold on;
box on
plot(10*log10(mean(Smooth_pre)),'LineWidth', wid,'Color', 'c');
% plot(10*log10(mean(Smooth_pos1)),'LineWidth', wid,'Color', 'g');
plot(10*log10(mean(Smooth_pos2)),'LineWidth', wid,'Color', 'b');
plot(10*log10(mean(Smooth_pos3)),'LineWidth', wid,'Color', 'r');
plot(10*log10(mean(Smooth_pos4)),'LineWidth', wid,'Color', 'g');
legend('without topological transformation ','topological transformation (\xi=3)','topological transformation (\xi=5)','topological transformation (\xi=10)');
xlabel('Iterations');ylabel('Signal smoothness(dB)');
