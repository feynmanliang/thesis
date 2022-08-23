% Nystrom plots

% toy: sgm = 10,5,1
% bodyfat: sgm = 2000,1000,200
% eunite2001: sgm = 10,5,1 or 8,4,1 ??
T = 1000;
R = 2;
sparsity = 1;
dataset = 2;
method = 1;
data_dir = '../../tail/data'; %'../../Microsoft/big_data'; %
kernelize = 1;
data_files = {'bodyfat','eunite2001','mg','triazines','housing','space_ga','abalone','mpg',...
              'MSD','a9a/a9a.t','covtype/bet_stages/covtype.t',...
              'w8a/w8a.t','mushrooms.txt','phishing.txt',...
              'madelon.t','splice_scale.txt','sonar_scale.txt'};
method_names = {'k-DPP','greedy'};
k_max = 40;
k_list = 1:k_max;
k_num = length(k_list);
sgm_list =  [8,4,1]; %[1400,500,200];%
sgm_num = length(sgm_list);
ratio_mean = ones(sgm_num,k_num);
ratio_std = zeros(sgm_num,k_num);
lbd_lists = zeros(sgm_num,k_num);

for sgm_i = 1:sgm_num
    sgm = sgm_list(sgm_i);
    if dataset==0
        n=50;
        d=50;
        [U,~,V] = svd(randn(n,d));
        l = ones(n,1);
        l(21:n) = 1/sgm;
        X = U*diag(l)*V';
    else
        [~, X] = libsvmread(sprintf('%s/%s',data_dir,data_files{dataset}));
        [n,d] = size(X);
    end
    if kernelize
        eta = ceil((d+1)/2);
        K=zeros(n,n);
        for i=1:n
            for j=1:n
                ll = norm(X(i,:)-X(j,:));
                K(i,j) =  exp(-ll^2/(2*sgm^2));
            end
        end
        X = chol(K)';
    end

    L = decompose_kernel(full(X*X'));
    lbds = sort(L.D,'descend');
    lbd_lists(sgm_i,:) = lbds(k_list);

    for i=k_list
        k = k_list(i);
        if k>=n, break; end
        Ek = sum(lbds(k+1:n));
        if method == 1
            for r=1:R
                [k,r]
                T_temp = T*2^r;
                ratio_list = zeros(T_temp,1);
                for t=1:T_temp
                    fail = true;
                    while fail
                        try
                            S = sample_dpp(L,k);
                            fail=false;
                        catch err
                        end
                    end
                    ErS = norm(X - X*pinv(X(S,:))*X(S,:),'fro')^2;
                    ratio_list(t) = ErS/Ek;
                end
                ratio_mean(sgm_i,i) = mean(ratio_list);
                ratio_std(sgm_i,i) = 3*std(ratio_list)/sqrt(T_temp);
                if ratio_std(sgm_i,i) < ratio_mean(sgm_i,i)/6, break; end
            end
        elseif method == 2
            k
            S = [];
            ErS = norm(X,'fro')^2;
            for j=1:k
                r_best = 1;
                for r = 1:n
                    S_temp = [S,r];
                    ErS_temp = norm(X - X*pinv(X(S_temp,:))*X(S_temp,:),'fro')^2;
                    if ErS_temp < ErS
                        ErS = ErS_temp;
                        r_best = r;
                    end
                end
                S = [S,r_best];
            end
            ratio_mean(sgm_i,i) = ErS/Ek;
        end
    end
end



set(0,'DefaultAxesFontSize', 13);
set(0,'DefaultTextFontSize', 13);
figure('Color','white','Position', [1000, 1000, 500, 500])
colors = {[1.0, 0.5, 0], [0, 0.45, 0], [0, 0, 1],[0,0,0.5]};
types = {'-','--','-.',':'};
symbols = {'o-','x-','s-'};

if dataset==0
    data_name = 'toy example';
    parameter_name = 'kappa';
else
    data_name = data_files{dataset};
    parameter_name = 'sigma';
end

subplot(2,1,1);
ph_list = [];
legend_list = {};

for i=1:sgm_num
    ph_list(end+1) = errorbar(k_list, ratio_mean(i,:),ratio_std(i,:),...
        symbols{i},'LineWidth',.5,'CapSize',2,'Color',colors{i});
    legend_list{end+1} = sprintf('\\%s = %.f',parameter_name,sgm_list(i));
    hold on
end

xlabel('Subset size k');
ylabel('Approximation factor')
xlim([0,k_max]);
legend(ph_list,legend_list,'Location','best');
title(sprintf('\\bf Nystr\\"om method (%s)',data_name),'Interpreter','latex');

subplot(2,1,2);
ph_list = [];
legend_list = {};

for i=1:sgm_num
    ph_list(end+1) = semilogy(k_list, lbd_lists(i,:),...
        symbols{i},'LineWidth',.5,'Color',colors{i});
    legend_list{end+1} = sprintf('\\%s = %.f',parameter_name,sgm_list(i));
    hold on
end

xlabel('Index');
ylabel('Eigenvalue')
xlim([0,k_max]);
legend(ph_list,legend_list,'Location','best');
title(sprintf('\\bf Spectral decay (%s)',data_name),'Interpreter','latex');
