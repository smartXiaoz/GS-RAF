clear all; close all
%% 邻接矩阵A；度矩阵D；拉普拉斯矩阵L；
%% 都是对称矩阵
n = 5; %% 图的节点的个数；
A = [ 0 1 0 1 1;
         1 0 1 0 0;
         0 1 0 0 1;
         1 0 0 0 1;
         1 0 1 1 0;];
     %% 节点23标号互换
 A2 = [ 0 0 1 1 1;
         0 0 1 0 1;
          0 0 0 1 1;
        1 0 0 0 1;
         1 1 0 1 0;];
          %% 节点24标号互换
 A3 = [ 0 1 0 1 1;
          1 0 0 0 1;
           0 0 0 1 1;
         1 0 1 0 0;
         1 1 1 0 0;];
%  A = ones(n);
 diag_D = sum(A);
 D = diag(diag_D);
 L = D - A;
 
 [Ua,Va]=eig(A);
 %% 对A做特征值分解后，其最大特征值，对应的特征向量，全是正的。
 %% U(:,end)' = 0.5299    0.3578    0.3578    0.4271    0.5299
 %% 对应节点的度为：3        2              2               2           3
 %% 物理意义：度大的对应的值大，节点4虽然度为2，但它与节点一、四相连；
%% A 图中2和3标号换一下，最大特征值对应的特征向量还是一样，其他特征向量数值相同，但符号有变化
%% A 图中2和3标号换一下，特征向量矩阵数值上24互换，正负号有变化。

[Ul, Vl]=eig(L);
%% 对L做特征值分解后，其最小特征值为0。
%% We denote the entire spectrum：特征值由小到大
%% 谱的定义：特征值由小到大的排列就是一个图的谱。
%% 频率的相似性：最小特征值为0的特征向量为相同值1/sqrt(N);
%% 频率的相似性：频率越大，过零点越多。
%% 度相同的节点，特征向量对应的数值相同。如：23、15节点在每个特征向量上对应的值相同

%% 计算过零点数值
for ii = 1 : n
    ul = Ul(:, ii);
    Num_cross_zero(ii) = 0;
    for kk = 1 : n
        for jj = 1 : n
            if A(kk,jj) > 0
                Num_cross_zero(ii) = Num_cross_zero(ii) + (ul(kk) * ul(jj) <0) ;
            end
        end
    end
end
figure,plot(Num_cross_zero)


%% 图中每个节点有一个值，构成一个向量f(i), i = 1,..., n, n为图中节点数。
%% f' * L * f >= 0.
f = randn(n,1);
f'*L*f;
% f = rand(n,1);f'*L*f

%% 图傅里叶变换：f_hat = U' * f
%% 图傅里叶逆变换：f = U * f_hat
f = randn(n,1);
f(end) = 5;
f_gf = Ul * inv(eye(n) + 0.1* Vl) * Ul' * f

%% 基于多项式的图滤波
k = 2;
Als = zeros(n,k);
for ii = 1 : k
    Als(:, ii) = L^(ii-1) * f;
end

bls = f - 0.1* L * f;
Theta = pinv(Als) * bls;
f_gfk = 0;
for ii = 1 : k
    f_gfk = f_gfk + Theta(ii) * L^(ii-1) * f;
end
f'*L*f
f_gf'*L*f_gf
f_gfk'*L*f_gfk
f, f_gf, f_gfk
