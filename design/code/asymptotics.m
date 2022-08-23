d_list = floor(logspace(1,2,20));
alpha_list = [0.5,0.75,0.9];
w_list = [1,1,1];
d_num = length(d_list);
alpha_num = length(alpha_list);
T = 100000;
mse_error = zeros(d_num,alpha_num);
mse_std = zeros(d_num,alpha_num);
sg = 0; %1

for i=1:d_num
    d = d_list(i)
    Sgm = diag(logspace(-2,0,d));%diag((d:-1:1).^(-2));%eye(d); Sgm(1,1)=0.1;% eye(d);%
    B = sqrt(Sgm);
    lbds = eig(Sgm);
    %w = zeros(d,1); w(ceil(0.25*d))=1;%[zeros(d-1,1);1]; %[1;zeros(d-1,1)]; %ones(d,1);%zeros(d,1);%  
    
    lbd = 1e10;
    for j=1:alpha_num
        w_ind = max(1,ceil(w_list(j)*d));
        w = zeros(d,1); w(w_ind)=1;
        w_squared = norm(w)^2;        
        alpha = alpha_list(j);
        n = floor(alpha*d)
        while sum(lbds./(lbds+lbd)) < n, lbd = lbd/(1+1e-3); end
        while sum(lbds./(lbds+lbd)) > n, lbd = lbd*(1+1e-6); end
        while sum(lbds./(lbds+lbd)) < n, lbd = lbd/(1+1e-9); end
        sum(lbds./(lbds+lbd))
        gamma = det((Sgm+lbd*eye(d))\Sgm);
        theory = sg^2*sum(1./(lbds+lbd))*(1-gamma)/(d-n) + (d-n)*w'*inv(Sgm+lbd*eye(d))*w / sum(1./(lbds+lbd));
        mse_est = zeros(T,1);
        for t=1:T
            X = randn(n,d)*B;
            Xinv = pinv(X);
            mse_est(t) = sg^2*trace(Xinv'*Xinv) + w_squared - w'*Xinv*X*w;
        end
        mse_std(i,j) = std(mse_est) / (sqrt(T) * theory);
        mse_error(i,j) = (mean(mse_est) - theory) / theory;
    end
end
    
set(0,'DefaultAxesFontSize', 15);
set(0,'DefaultTextFontSize', 15);
figure('Color','white','Position', [500, 500, 500, 400])
ph_list = [];
legend_list = {};
colors = {[0.9, 0.75, 0], [1.0, 0.5, 0], [1, 0, 0], [0, 0, 0]};

for j=1:alpha_num
    ph_list(end+1) = errorbar(d_list, abs(mse_error(:,j)),mse_std(:,j),'LineWidth',2,'Color',colors{j});
    legend_list{end+1} = sprintf('n/d=%.2f, sp=%.2f',alpha_list(j),w_list(j));
    hold on
end
set(gca, 'XScale','log', 'YScale','log')
xlabel('Dimension d')
ylabel('|MSE[X^+y] - theory| / theory');
legend(ph_list,legend_list,'Location','best');
hold off
