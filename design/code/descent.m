n=100;

d_list = unique(ceil(logspace(1,3,20)'));
d_list = d_list(abs(d_list-n)>5);
d_max = max(d_list);
d_min = min(d_list);
d_list_big = (d_min:d_max)';
d_num = length(d_list);
T = 20;
show_experiments = false;
show_norm = true;

% B = randn(d,d);
% [U,~,V] = svd(B,'econ');
%Sgm = eye(d); %Sgm(1,1) = 0.0001; %diag(logspace(-5,0,d));%diag((1:d).^(-2)); %
%B = sqrt(Sgm);
%lbds = eig(Sgm);
%w = [1;zeros(d-1,1)];
%w_squared = norm(w)^2;



num_opts = 3;
sg = 1;
w_norm = 1;
logkappas = [0,2,5];
mse_mean = zeros(d_num,num_opts);
mse_error = zeros(d_num,num_opts);
mse_count = 1;
theory = zeros(length(d_list_big),num_opts);
for i=1:length(d_list_big)
    d = d_list_big(i);
    for opt=1:num_opts
        if opt <= 5
            Sgm = diag(logspace(0,-logkappas(opt),d));
        else
            Sgm = diag(logspace(0,-1,d));
        end
        Sgm = Sgm * trace(inv(Sgm)) / d;
        B = sqrt(Sgm);
        lbds = eig(Sgm);
        if opt>5
            w = w_norm(opt-5)*ones(d,1)/sqrt(d); %w_norm * [zeros(d-1,1);1]; %w = w_norm * ones(d,1)/sqrt(d);
        else
            w = w_norm(1)*ones(d,1)/sqrt(d);
        end
%         ws(1,1) = w_norm(1);
%         ws(ceil(d/2),2) = w_norm(2);
%         ws(d,3) = w_norm(3);
        if d <= n
            theory(i,opt) = sg^2*sum(1./lbds);
            if d < n, theory(i,opt) = theory(i,opt) * (1-exp(d-n))/(n-d); end
        else
            lbd = 1e10;
            while sum(lbds./(lbds+lbd)) < n, lbd = lbd/1.1; end
            while sum(lbds./(lbds+lbd)) > n, lbd = lbd*1.001; end
            while sum(lbds./(lbds+lbd)) < n, lbd = lbd/1.00001; end
            sum(lbds./(lbds+lbd))-n;
            gamma = det((Sgm+lbd*eye(d))\Sgm);
            theory(i,opt) = sg^2*sum(1./(lbds+lbd))*(1-gamma)/(d-n) +...
                    (d-n)*w'*inv(Sgm+lbd*eye(d))*w / sum(1./(lbds+lbd));
        end
        if show_experiments && sum(d_list == d) > 0
            d
            mse_list = zeros(T,1);
            for t=1:T
                X = randn(n,d)*B;
                Xinv = pinv(X);
                P = Xinv*X;
                Xinv_sq = Xinv'*Xinv;
                mse_list(t) = sg^2*trace(Xinv_sq) + norm(w)^2 - w'*P*w;
            end
            mse_mean(mse_count,opt) = mean(mse_list);
            mse_error(mse_count,opt) = std(mse_list);
            if opt == num_opts, mse_count = mse_count + 1; end
        end    
    end
end

set(0,'DefaultAxesFontSize', 15);
set(0,'DefaultTextFontSize', 15);
figure('Color','white','Position', [1000, 1000, 400, 225])
ph_list = zeros(num_opts,1);
width_list = [1.5,1,1,1.5];
type_list = {'-','--','-.',':'};
legend_list = cell(num_opts,1);
colors = [0,0,0; 1,0,0; 0,0.5,0; 0,0,1]; %[1:5]'*[.2, .1, 0];
if show_norm
    semilogx(d_list_big./n,w_norm^2*ones(length(d_list_big),1),'-','LineWidth',.5,'Color','k');
    %legend_list{1} = '||w||^2';
    hold on
end
for opt=num_opts:-1:1
    ph_list(opt) = semilogx(d_list_big./n,theory(:,opt),type_list{opt},'LineWidth',width_list(opt),'Color',colors(opt,:));
    if opt==1
        legend_list{opt} = 'isotropic';
    else
        legend_list{opt} = sprintf('\\kappa = 1e%d',logkappas(opt));
    end
    hold on
    if show_experiments
        d_filter = true(d_num,1);
        if opt>1, d_filter = d_list > n; end
        errorbar(d_list(d_filter)./n,mse_mean(d_filter,opt),mse_error(d_filter,opt),'.',...
            'MarkerSize',10,'CapSize',2,'Color',colors(opt,:));
    end
%     legend_list{end+1} = sprintf('%d: i.i.d. (empirical)',opt);
end
xlabel('d/n');
ylabel('MSE');
ylim([theory(1,1),2*theory(end,3)]);
xlim([d_min,d_max]/n);
legend(ph_list,legend_list,'Location','best');
hold off