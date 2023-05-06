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
sigma2 = 0.2;%%RMEE

tic
% addpath(toolboxdir('stats'))
% addpath(genpath('Statistics and Machine Learning Toolbox'))
A_g = ones(q);
D_g = diag(sum(A_g));
L_g = D_g - A_g - (q-1)*eye(q);

for mm = 1:10
%     vv = chi2rnd(1,[1, LEN]);
%     vv = trnd(5,[1, L]) * 0.1;



%% case1
    vv = randn(1,L) * 0.1;
    v1=randn(1,L)*0.1; v2=randn(1,L)*10; 
    rp=rand(1,L);

    vv = (rp<=0.98).*v1 + (rp>0.98).*v2;
    uu = randn(p,L);
    wo1 = [0.5,0.4,0.3, 0.2,0.1, 0.1,0.2,0.3, 0.4,0.5]';
    wo = [ kron(wo1, ones(1,L)) ];
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
    % RLS_MSE
    tic 
    for zz = 1 : 100
    m=p;
    Pn = eye(m)*1;
    for ii = 1 : L
        Err_RLS(mm,ii) = (wo(:,ii)  - w_RLS)' * (wo(:,ii)  - w_RLS);
        dn = dd(ii);
        un = uu(:,ii);
        en = dn - w_RLS' * un;
        kn = Pn * un / ( lamda1+ un' * Pn * un );
        Pn = 1/lamda1 * ( Pn - kn * un' * Pn);
        w_RLS = w_RLS +kn * en;
    end
    end
    toc
     %% RLS-MCC
    tic 
    for zz = 1 : 100  
    Pn = eye(m)*1;
    for ii = 1 : L
        Err_MCC_RLS(mm,ii) = (wo(:,ii)  - w_C_RLS)' * (wo(:,ii)  - w_C_RLS);
        dn = dd(ii);
        un = uu(:,ii);
        en = dn - w_C_RLS' * un;
         kn = Pn * un / ( 1+ un' * Pn * un );
        kn = Pn * un / ( exp(en^2/2/sigma1^2)*lamda2 + un' * Pn * un );
        
        Pn = 1/lamda2 * ( Pn - kn * un' * Pn);
        w_C_RLS = w_C_RLS +kn * en;
    end
    end
    toc
    %% RLS MEE, algorithm 1
    tic
    for zz  = 1 : 100
    PL = eye(p)*1; 
    eL = zeros(q, 1); 
    PHI = zeros(q, q);
    
    for ii = q : L
        Err_RLS_MEE1(mm,ii) = (wo(:,ii)  - w_M_RLS1)' * (wo(:,ii)  - w_M_RLS1);
        for jj = 1 : q
            eL(jj) = dd(ii - q + jj) - w_M_RLS1' * uu(: , ii - q + jj);
        end
        u0 = uu(:,ii);
        e0 = eL(q);
        phi0 = 0;
        for kk = 2:q
            ek = eL(q-kk+1);
            phi0 = phi0 + lamda3^(kk-1) * exp(- (e0-ek).^2/sigma2^2/2);  
        end
        KL = PL * u0 * inv ( lamda3^2 / phi0 + u0' * PL * u0 );
        PL = 1/lamda3^2 * (PL - KL * u0' * PL);
        w_M_RLS1 = w_M_RLS1 + KL * e0;        
    end
    end
    toc
   
    
    %% RLS MEE, algorithm 7
    m = p;
%     L = q;
    PL = eye(m,m)*0.1; 
    eL = zeros(q, 1);
    uL = zeros(m,q);
    PHI = zeros(q, q);
    Err_GS_RAF = zeros(100,2500);
    tic
    for zz = 1 : 100
    for ii = q : L
        Err_GS_RAF(mm,ii) = (wo(:,ii)  - w_M_RLS7)' * (wo(:,ii)  - w_M_RLS7);
        eL = dd(ii - 9:ii)' - uu(: , ii - 9:ii)' * w_M_RLS7;
        uL = uu(: , ii - 9:ii);         
        L_1 = L_g;
        %% ÍØÆË±ä»¯
        u = eL / median(abs(eL - median(eL)));
        for jj = 1 : q
            if abs(u(jj)) > 4
                L_1(jj,:)=0;
                L_1(:,jj)=0;
            end
        end      
        for jj = 1 : q          
                L_1(jj,jj)=-sum(L_1(jj,:));
        end          
        PHI = L_1  + 1/1000*eye(q);
        KL =PL * uL * inv ( lamda4*inv(PHI)+uL'*PL*uL);
        PL = 1/lamda4 * (PL - KL * uL' * PL);
        w_M_RLS7 = w_M_RLS7 + KL * eL;        
    end
    end
    toc
     
    disp(mm)
end
toc
wid = 1.5;
figure,hold on;
box on
plot(10* log10(mean(Err_RLS)),'-o','LineWidth', wid,'Color', 'c','MarkerFaceColor', 'c', 'MarkerIndices', 1:200:L);
plot(10* log10(mean(Err_MCC_RLS)),'-s','LineWidth', wid,'Color', 'g','MarkerFaceColor', 'g', 'MarkerIndices', 1:200:L);
plot(10* log10(mean(Err_RLS_MEE1)),'-d','LineWidth', wid,'Color', 'b','MarkerFaceColor', 'b', 'MarkerIndices', 1:200:L);
plot(10* log10(mean(Err_GS_RAF)),'-^','LineWidth', wid,'Color', 'r','MarkerFaceColor', 'r', 'MarkerIndices', 1:200:L);


legend('RLS (\lambda=0.995)','RMC (\lambda=0.99, \sigma=2)','RMEE (\lambda=0.995, \sigma=0.2)','GS-RAF (\lambda=0.995)');
xlabel('Iterations');ylabel('MSD');

RLS = 10* log10(mean(mean(Err_RLS(:,1000:end))))
RMC = 10* log10(mean(mean(Err_MCC_RLS(:,1000:end))))
RMEE = 10* log10(mean(mean(Err_RLS_MEE1(:,1000:end))))
GS_RAF = 10* log10(mean(mean(Err_GS_RAF(:,1000:end))))
