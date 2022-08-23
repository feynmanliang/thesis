d=50;

n_list = ceil(0.9*d);
n_num = length(n_list);
T = 10000;

% B = randn(d,d);
% [U,~,V] = svd(B,'econ');
Sgm = eye(d); Sgm(1,1) = 0.0001; %diag(logspace(-5,0,d));%diag((1:d).^(-2)); %
B = sqrt(Sgm);
lbds = eig(Sgm);
w = [1;zeros(d-1,1)];
w_squared = norm(w)^2;
sg = 0;

set(0,'DefaultAxesFontSize', 15);
set(0,'DefaultTextFontSize', 15);
figure('Color','white','Position', [500, 500, 500, 400])
ph_list = [];
legend_list = {};
colors = {[1.0, 0.5, 0], [1, 0, 0], [0, 0, 0]};

mse_mean = zeros(n_num,T);
theory = zeros(n_num,1);
lbd = 1e10;
for i=1:n_num
    n = n_list(i)
    while sum(lbds./(lbds+lbd)) < n, lbd = lbd/1.1; end
    while sum(lbds./(lbds+lbd)) > n, lbd = lbd*1.001; end
    while sum(lbds./(lbds+lbd)) < n, lbd = lbd/1.00001; end
    sum(lbds./(lbds+lbd))-n
    gamma = det((Sgm+lbd*eye(d))\Sgm);
    theory(i) = sg^2*sum(1./(lbds+lbd))*(1-gamma)/(d-n) + (d-n)*w'*inv(Sgm+lbd*eye(d))*w / sum(1./(lbds+lbd));
    mse_list = zeros(T,1);
    error_list = zeros(T,1);
    for t=1:T
        X = randn(n,d)*B;
        Xinv = pinv(X);
        mse_list(t) = sg^2*trace(Xinv'*Xinv) + w_squared - w'*Xinv*X*w;
%         w_est = pinv(X)*(X*w+sg*randn(n,1));
%         mse_list(t) = norm(w_est-w)^2;
        if t==1 
            mse_mean(i,t) = mse_list(t); 
        else
            mse_mean(i,t) = mse_mean(i,t-1) + mse_list(t);
        end
        error_list(t) = abs(mse_mean(i,t)/t - theory(i))/theory(i);
    end
    
    ph_list(end+1) = semilogy(1:T,mse_mean(i,:)./(1:T),'LineWidth',1,'Color',colors{i});
    legend_list{end+1} = sprintf('n=%d',n);
    hold on
%     ph_list(end+1) = loglog(1:T,error_list,'LineWidth',1,'Color',colors{i});
%     legend_list{end+1} = sprintf('n=%d',n);
    hold on

end

legend(ph_list,legend_list,'Location','best');
hold off