d=100;

n_list = 1:1:d-2;
n_num = length(n_list);
r = 1;
T = 1000;

% B = randn(d,d);
% [U,~,V] = svd(B,'econ');
% D = diag((1:d).^(-2));%diag(logspace(-2,0,d));
% B = U*D*V;
% Sgm = B'*B;
% lbds = eig(Sgm);
% w = 10*randn(d,1);

sg = 1;
set(0,'DefaultAxesFontSize', 15);
set(0,'DefaultTextFontSize', 15);
figure('Color','white','Position', [500, 500, 500, 400])
ph_list = [];
legend_list = {};
colors = [1,0.8,0];

mse_std = zeros(n_num,1);
theory = zeros(n_num,1);
lbd = 1e10;
for i=1:n_num
    n = n_list(i);
    while sum(lbds./(lbds+lbd)) < n, lbd = lbd/2; end
    while sum(lbds./(lbds+lbd)) > n, lbd = lbd*1.01; end
    gamma = det((Sgm+lbd*eye(d))\Sgm);
    theory(i) = sg^2*sum(1./(lbds+lbd))/(d-n-1) + (d-n)*w'*inv(Sgm+lbd*eye(d))*w / sum(1./(lbds+lbd));
    mse_list = zeros(T,1);
    for t=1:T
        X = randn(n,d)*B;
        w_est = pinv(X)*(X*w+sg*randn(n,1));
        mse_list(t) = norm(w_est-w)^2;
    end
    mse_mean(i) = mean(mse_list);
    mse_std(i) = std(mse_list)/sqrt(T);
%     rectangle('Position',[n_list(i)-1.5,theory(i)-mse_std(i),1,2*mse_std(i)],...
%         'Curvature',[0.3 0.3],'FaceColor',[.9,.9,.9,0.3],'EdgeColor',[1,1,1,0]);
%     if n>80
%         mid_x = (ridge_x(i-1)+ridge_x(i))/2;
%         mid_y = (ridge_y(i-1)+ridge_y(i))/2;
%         rectangle('Position',[mid_x-mse_std(i),mid_y-mse_std(i),2*mse_std(i),2*mse_std(i)],...
%             'Curvature',[0.3 0.3],'FaceColor',[.9,.9,.9,0.3],'EdgeColor',[1,1,1,0]);
%         if n>90
%             mid_x = (mid_x+ridge_x(i))/2;
%             mid_y = (mid_y+ridge_y(i))/2;
%             rectangle('Position',[mid_x-mse_std(i),mid_y-mse_std(i),2*mse_std(i),2*mse_std(i)],...
%                 'Curvature',[0.3 0.3],'FaceColor',[.9,.9,.9,0.3],'EdgeColor',[1,1,1,0]);
%         end
%     end
end
hold on
ph_list(end+1) = plot(n_list,theory,'LineWidth',1,'Color',[0,0,0,0.5]);
legend_list{end+1} = 'theory';

for i=1:n_num
    n = n_list(i);
    ph = scatter(n_list(i),mse_mean(i),16,colors*0.5,'x');
    if i==1 || i==n_num || i==round(n_num/2)
        ph_list(end+1) = ph; 
        legend_list{end+1} = sprintf('n=%d',n);
    end    
end

% xlim([0,x_max]);
% ylim([0,0.5]);
xlim([1,d]);
legend(ph_list,legend_list,'Location','best');
hold off