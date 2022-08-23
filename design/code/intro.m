d=100;

n_list = unique(ceil(linspace(5,2*d,25)'));
n_list = n_list(abs(n_list-d)>5);
n_max = max(n_list);
n_min = min(n_list);
n_list_big = (1:n_max)';
n_num = length(n_list);
T =50;
%R=5;

% B = randn(d,d);
% [U,~,V] = svd(B,'econ');
%Sgm = eye(d); %Sgm(1,1) = 0.0001; %diag(logspace(-5,0,d));%diag((1:d).^(-2)); %
%B = sqrt(Sgm);
%lbds = eig(Sgm);
%w = [1;zeros(d-1,1)];
%w_squared = norm(w)^2;

rng(2);

num_opts = 3;
show = [true,true,true];
show_norm = true;
plot_type = 'mse';
show_experiments = true;
show_iid = false;
sg = .4;
w_norm = 1;
logkappas = [0,2,5];
mse_mean = zeros(n_num,num_opts);
bias_mean = zeros(n_num,num_opts);
shrinkage_mean = zeros(n_num,num_opts);
mse_error = zeros(n_num,num_opts);
bias_error = zeros(n_num,num_opts);
shrinkage_error = zeros(n_num,num_opts);
mse_count = 1;
theory = zeros(length(n_list_big),num_opts);
theory2 = zeros(length(n_list_big),num_opts);
bias_theory = zeros(length(n_list_big),num_opts);
shrinkage_theory = zeros(length(n_list_big),num_opts);
lambda_theory = zeros(length(n_list_big),num_opts);
for i=1:length(n_list_big)
    n = n_list_big(i);
    for opt=1:num_opts
        Sgm = diag(logspace(0,-logkappas(opt),d));
        Sgm = Sgm * trace(inv(Sgm)) / d;
        lbds = eig(Sgm);
        B = sqrt(Sgm);
        w = w_norm(1)*ones(d,1)/sqrt(d);
        if d <= n
            theory(i,opt) = sg^2*sum(1./lbds);
            bias_theory(i,opt) = 0;
            shrinkage_theory(i,opt) = w_norm;
            lambda_theory(i,opt) = 0;
            if d < n-1
                theory2(i,opt) = sg^2*sum(1./lbds)/(n-d-1); 
            else
                theory2(i,opt) = theory2(i-1,opt);
            end
            if d < n, theory(i,opt) = theory(i,opt) * (1-exp(d-n))/(n-d); end
        else
            lbd = 1e10;
            while sum(lbds./(lbds+lbd)) < n, lbd = lbd/1.1; end
            while sum(lbds./(lbds+lbd)) > n, lbd = lbd*1.001; end
            while sum(lbds./(lbds+lbd)) < n, lbd = lbd/1.00001; end
            gamma = det((Sgm+lbd*eye(d))\Sgm);
            theory(i,opt) = sg^2*sum(1./(lbds+lbd))*(1-gamma)/(d-n) +...
                    (d-n)*w'*inv(Sgm+lbd*eye(d))*w / sum(1./(lbds+lbd));
            bias_theory(i,opt) = norm(w - inv(Sgm+lbd*eye(d))*Sgm*w);
            shrinkage_theory(i,opt) = norm(inv(Sgm+lbd*eye(d))*Sgm*w);
            lambda_theory(i,opt) = lbd;
            if n < d-1
                theory2(i,opt) = sg^2*n/(d-n-1) + norm(w)^2*(d-n)/d;
            else
                theory2(i,opt) = theory2(i-1,opt);
            end
        end
        if sum(n_list == n) > 0
            n
            mse_list = zeros(T,1);
            bias_list = zeros(T,d);
            shrinkage_list = zeros(T,d);
%             mse_mean_list = zeros(R,1);
%             bias_mean_list = zeros(R,1);
%             shrinkage_mean_list = zeros(R,1);
%             for r=1:R
            for t=1:T
                X = randn(n,d)*B;
                Xinv = pinv(X);
                P = Xinv*X;
                Xinv_sq = Xinv'*Xinv;
                mse_list(t) = sg^2*trace(Xinv_sq) + norm(w)^2 - w'*P*w;
                shrinkage_list(t,:) = P*w;
                bias_list(t,:) = w-P*w;
            end
%             mse_mean_list(r) = mean(mse_list);
%             bias_mean_list(r) = norm(mean(bias_list));
%             shrinkage_mean_list(r) = norm(mean(shrinkage_list));
            %btstrp = randsample(T,T,true); 
            
            mse_mean(mse_count,opt) = mean(mse_list);
            mse_error(mse_count,opt) = 3*std(mse_list)/sqrt(T); %2*abs(mean(mse_list) - mean(mse_list(btstrp)));
            bias_mean(mse_count,opt) = norm(mean(bias_list));
            bias_error(mse_count,opt) = 3*mean(std(bias_list))/sqrt(T);  %norm(mean(bias_list) - mean(bias_list(btstrp)));
            shrinkage_mean(mse_count,opt) = norm(mean(shrinkage_list));
            shrinkage_error(mse_count,opt) = 3*mean(std(shrinkage_list))/sqrt(T); 
            %norm(mean(shrinkage_list) - mean(shrinkage_list(btstrp)));
            if opt == num_opts, mse_count = mse_count + 1; end
        end    
    end
end

set(0,'DefaultAxesFontSize', 13);
set(0,'DefaultTextFontSize', 13);
figure('Color','white','Position', [1000, 1000, 400, 300])
width_list = [1.5,1,1,1.5];
type_list = {'-','--','-.',':'};
marker_list = {'o','d','s','.'};
colors = [0,0,0; 1,0,0; 0,0.5,0; 0,0,1]; %[1:5]'*[.2, .1, 0];
if strcmp(plot_type,'bias')
    offset = 0;
    ph_list = zeros(num_opts+offset,1);
    legend_list = cell(num_opts+offset,1);

    for opt=num_opts:-1:1
        if show(opt)
            ph_list(opt+offset) = plot(n_list_big,bias_theory(:,opt),type_list{opt},...
                'LineWidth',width_list(opt),'Color',colors(opt,:));
            if opt==1
                legend_list{opt+offset} = 'isotropic';
            else
                legend_list{opt+offset} = sprintf('\\kappa = 1e%d',logkappas(opt));
            end
            hold on
            if show_experiments
                n_filter = true(n_num,1);
                if opt>1, n_filter = n_list < d; end
                errorbar(n_list(n_filter),bias_mean(n_filter,opt),bias_error(n_filter,opt),'.',...
                    'MarkerSize',10,'CapSize',2,'Color',colors(opt,:));
            end
        end
    end
    ylabel('Bias ||E[X^+y] - w||');
    ylim([0,1]);
    xlim([1,n_max]);
elseif strcmp(plot_type,'shrinkage')
    ph_list = [];
    legend_list = {};
    for opt=1:num_opts
        if show(opt)
            ph_list(end+1) = plot(n_list_big,shrinkage_theory(:,opt),type_list{opt},...
                'LineWidth',width_list(opt),'Color',colors(opt,:));
            if opt==1
                legend_list{end+1} = 'isotropic theory';
            else
                legend_list{end+1} = sprintf('\\kappa = 1e%d theory',logkappas(opt));
            end
            hold on
            if show_experiments
                n_filter = true(n_num,1);
                %if opt>1, n_filter = n_list < d; end
                ph_list(end+1) = errorbar(n_list(n_filter),shrinkage_mean(n_filter,opt),...
                    shrinkage_error(n_filter,opt),marker_list{opt},...
                    'MarkerSize',5,'CapSize',2,'Color',colors(opt,:));
                if opt==1
                    legend_list{end+1} = 'isotropic empirical';
                else
                    legend_list{end+1} = sprintf('\\kappa = 1e%d empirical',logkappas(opt));
                end                
            end
        end
    end
    ylabel('Norm ||E[X^+y]||');
    ylim([0,1]);
    xlim([1,n_max]);    
elseif strcmp(plot_type,'lambda')
    offset = 0;
    ph_list = zeros(num_opts+offset,1);
    legend_list = cell(num_opts+offset,1);
    for opt=num_opts:-1:1
        if show(opt)
            ph_list(opt+offset) = semilogy(n_list_big,lambda_theory(:,opt),type_list{opt},...
                'LineWidth',width_list(opt),'Color',colors(opt,:));
            hold on
            if opt==1
                legend_list{opt+offset} = 'isotropic';
            else
                legend_list{opt+offset} = sprintf('\\kappa = 1e%d',logkappas(opt));
            end
        end
    end
    ylabel('Implicit regularizer \lambda_n');
else    
    offset = 1;
    ph_list = []; %zeros(num_opts+offset,1);
    legend_list = {}; %cell(num_opts+offset,1);
    if show_norm
        plot(n_list_big,w_norm^2*ones(length(n_list_big),1),'-','LineWidth',.5,'Color','k');
        %legend_list{1} = '||w||^2';
        hold on
    end
    if show_iid
        ph_list(end+1) = plot(n_list_big,theory2(:,1),':','LineWidth',2,'Color','k');
        legend_list{end+1} = 'isotropic i.i.d.';
        hold on
    end
    for opt=1:num_opts
        if show(opt)
            ph_list(end+1) = plot(n_list_big,theory(:,opt),type_list{opt},'LineWidth',width_list(opt),'Color',colors(opt,:));
            if opt==1
                legend_list{end+1} = 'isotropic theory';
            else
                legend_list{end+1} = sprintf('\\kappa = 1e%d theory',logkappas(opt));
            end
            hold on
            if show_experiments            
                n_filter = true(n_num,1);
                %if opt>1, n_filter = n_list < d; end
                ph_list(end+1) =errorbar(n_list(n_filter),mse_mean(n_filter,opt),mse_error(n_filter,opt),marker_list{opt},...
                    'MarkerSize',5,'CapSize',2,'Color',colors(opt,:));
                if opt==1
                    legend_list{end+1} = 'isotropic empirical';
                else
                    legend_list{end+1} = sprintf('\\kappa = 1e%d empirical',logkappas(opt));
                end                
            end
        end
    end
    ylabel('MSE');
    %ph_list = ph_list([true,show]);
    %legend_list = legend_list{[true,show]};
    y_min = min(min(theory));
    ylim([0.9*y_min,12*y_min]);
end
xlabel('n');

legend(ph_list,legend_list,'Location','best');
hold off