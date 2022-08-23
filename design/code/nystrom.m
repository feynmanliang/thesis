% Nystrom plots

% bodyfat: sgm = 2000,1000,500
% eunite2001: sgm = 10,5,2.5,1
n=200;
T = 100;
R = 7;
sparsity = 1;
experiments = false;
from_file = false;
dataset = 2;
data_dir = '../../tail/data'; %'../../Microsoft/big_data'; %
expand_dims=0;
kernelize = 0;
data_files = {'bodyfat','eunite2001','mg','triazines','housing','space_ga','abalone','mpg',...
              'MSD','a9a/a9a.t','covtype/bet_stages/covtype.t',...
              'w8a/w8a.t','mushrooms.txt','phishing.txt',...
              'madelon.t','splice_scale.txt','sonar_scale.txt'};%'rcv1/rcv1_train.binary',
% data_files = {'cadata','MSD',...
%     'housing','space_ga','cpusmall','cpusmall_scale','abalone',...
%     'abalone_scale','bodyfat','bodyfat_scale','mg','mg_scale',...
%     'housing_scale','triazines','pyrim','mpg','mpg_scale'};
lower = [];

if from_file
    [~, X] = libsvmread(sprintf('%s/%s',data_dir,data_files{dataset}));
    [n,d] = size(X);
    if expand_dims        
        Xe = zeros(d^2,n);
        for i=1:n
            A = X(i,:)'*X(i,:);
            Xe(:,i) = A(:);
        end
        X = Xe'./svds(Xe,1);
    end
    [n,d] = size(X);
    if kernelize
        sgm = 250;
        C = 3*sgm;
        eta = ceil((d+1)/2);
        K=zeros(n,n);
        for i=1:n
            for j=1:n
                ll = norm(X(i,:)-X(j,:));
                K(i,j) =  exp(-(ll/sgm)^2);%max((1-ll/C)^eta,0) *                 
            end
        end
        X = chol(K)';
    end

else
    if dataset == 1
        X = randn(n,n);
    elseif dataset == 2
        [U,~,V] = svd(randn(n));
        l = ones(n,1);
        l(30:n) = 0.01;
        X = U*diag(l)*V';
    elseif dataset == 3
        [U,~,V] = svd(randn(n));
        l = ones(n,1);
        l(21:n) = 0.01;
        l(41:n) = 0.0001;
        X = U*diag(l)*V';   
        lower = [20,17;40,17]; 
    elseif dataset == 4
        [U,~,V] = svd(randn(n));
        l = 1./(1:n);
    %    l(ceil(n/2):n) = 1/ceil(n/2);
        X = U*diag(l)*V';    
    end
    k_max = n/4;    
end

if from_file && ~experiments
    lbds = svds(X,d).^2;
    lbds = lbds/lbds(1);
    lbds(lbds<1e-30)=0;
    if length(lbds)>1000
        lbds = lbds(1:1000);
    end
else
    L = decompose_kernel(full(X*X'));
    lbds = sort(L.D,'descend');
end
n = sum(lbds>0);    
k_max = min(50,n-1);
k_list = 1:k_max;
k_num = length(k_list);
uratio_mean = ones(k_num,1);
uratio_std = zeros(k_num,1);
ratio_mean = ones(k_num,1);
ratio_std = zeros(k_num,1);
ratio_sympoly = ones(k_num,1);
% theory_ratio = {};
% theory_k = {};
sympoly = elem_sympoly(lbds,length(lbds));

if experiments || from_file
    for i=k_list
        k = k_list(i);
        if k>=n, break; end
        Ek = sum(lbds(k+1:n));
        if experiments
            for r=1:R
                [k,r]
                T_temp = T*2^r;
                ratio_list = zeros(T_temp,1);
                uratio_list = zeros(T_temp,1);
                for t=1:T_temp
                    fail = true;
                    while fail
                        try
                            S = sample_dpp(L,k);
                            S2 = randsample(n,k);
                            fail=false;
                        catch err
                        end
                    end
                    ErS = norm(X - X*pinv(X(S,:))*X(S,:),'fro')^2;
                    ErS2 = norm(X - X*pinv(X(S2,:))*X(S2,:),'fro')^2;
                    ratio_list(t) = ErS/Ek;
                    uratio_list(t) = ErS2/Ek;
                end
                ratio_mean(i) = mean(ratio_list);
                ratio_std(i) = 3*std(ratio_list)/sqrt(T_temp);
                uratio_mean(i) = mean(uratio_list);
                uratio_std(i) = 3*std(uratio_list)/sqrt(T_temp);
                if ratio_std(i) < ratio_mean(i)/6, break; end
            end
        end
        ratio_sympoly(i) = (k+1)*(sympoly(k+2,n+1)/sympoly(k+1,n+1))/Ek;
    end
end
kd_list = (1:0.1:k_max)';
kd_num = length(kd_list);
theory_min = 100*kd_list;%1+ kd_list;
theory_curves = {};
theory_k = {};
s_list = [10,15,25];
k_range = 1:kd_num;
for s = 0:sparsity:n-1
    r = sum(lbds(s+1:n))/lbds(s+1);
    kappa = lbds(s+1)/lbds(n);
    ks_ind = (kd_list > s) & (kd_list < s+r);
    ks = kd_list(ks_ind);
    ks_range = k_range(ks_ind);
    ks_num = length(ks);
    if ks_num == 0, continue; end
    ratios_a = zeros(ks_num,1);
    ratios_b = zeros(ks_num,1);
    for i=1:ks_num
        k = ks(i) - s;
        ratios_a(i) = (1 + s/k) * sqrt(1 + 2*k/(r-k));
        ratios_b(i) = (1 + s/k) * kappa;
        theory_min(ks_range(i)) = min([theory_min(ks_range(i)),ratios_a(i)]);
    end
    if sum(s_list==s)>0
        theory_curves{end+1} = ratios_a;
        theory_k{end+1} = ks;
    end
end

set(0,'DefaultAxesFontSize', 13);
set(0,'DefaultTextFontSize', 13);
figure('Color','white','Position', [1000, 1000, 500, 200])

ph_list = [];
legend_list = {};
colors = {[1.0, 0.5, 0], [0, 0.45, 0], [0, 0, 1],[0,0,0.5]};
%colors = {[1.0, 0.5, 0], [0, 0, 1], [0, 0.5, 0],[0,0,0.5]};
types = {':','--','-.'};

% ph_list(end+1) = plot(k_list,1+k_list,':','LineWidth',1,'Color',[0,0,0]);
% legend_list{end+1} = 'Deshpande et al. (2006)';
% hold on
for i=1:length(theory_curves) 
    ph_list(end+1) = plot(theory_k{i},theory_curves{i},types{i},'LineWidth',2,'Color',colors{i});
    legend_list{end+1} = sprintf('\\Phi_{%d}(k)',s_list(i));
    hold on
    ks_min = theory_k{i}(1);
    ks_max = theory_k{i}(end);
    plot(theory_k{i}, i/2+zeros(length(theory_k{i}),1),types{i},'LineWidth',2,'Color',colors{i});
    plot([ks_min,ks_min],i/2+[0,30],':','LineWidth',1,'Color',colors{i});
    plot([ks_max,ks_max],i/2+[0,30],':','LineWidth',1,'Color',colors{i});
end
% ph_list(end+1) = plot(k_list,theory_a(k_list),'--','LineWidth',1,'Color',[1,0,0,]);
% legend_list{end+1} = 'Stable rank bound';
% ph_list(end+1) = plot(k_list,theory_b(k_list),'-.','LineWidth',1,'Color',[1,0,0,]);
% legend_list{end+1} = 'Condition number bound';

ph_list(end+1) = plot(kd_list,theory_min,'-','LineWidth',2,'Color',[1,0,0]);
%plot(kd_list(1:20:kd_num),theory_min(1:20:kd_num),':','Color',[1,0,0,0.5]);
% legend_list{end+1} = 'This paper (upper)';
legend_list{end+1} = sprintf('min \\Phi_s(k)');
hold on
% ph_list(end+1) = plot(k_list,ratio_sympoly,'-.','LineWidth',2,'Color','b');
% legend_list{end+1} = 'Volume sampling';


for i=1:length(lower)
    k = lower(i,1);
    v = lower(i,2);
    plot([k,k], [0,v],'-.','LineWidth',2,'Color',[0,0,0]);
end
if isempty(lower)==false
    ph_list(end+1) = plot([1,n], [0,0],'-.','LineWidth',2,'Color',[0,0,0]);
    legend_list{end+1} = 'This paper (lower)';
end

if experiments
    ph_list(end+1) = errorbar(k_list, ratio_mean,ratio_std,'o-','LineWidth',1,'Color','b');
    legend_list{end+1} = 'Volume sampling';
    hold on
%     ph_list(end+1) = errorbar(k_list,uratio_mean,uratio_std,'s-','LineWidth',1,'Color','r');
%     legend_list{end+1} = 'Uniform';    
end

xlabel('Subset size k');
ylabel('Approximation factor')
xlim([0,k_max]);
ylim([0,30]);%ylim([0,1.1*k_max]);
legend(ph_list,legend_list,'Location','best');
