d=100;

n_list = 1:1:d;
n_num = length(n_list);
r = 1;
T = 1000;

% B = randn(d,d);
% [U,~,V] = svd(B,'econ');
% D = diag(logspace(-2,0,d));
% B = U*D*V;
% Sgm = B'*B;
% lbds = eig(Sgm);
% w = randn(d,1);

set(0,'DefaultAxesFontSize', 15);
set(0,'DefaultTextFontSize', 15);
figure('Color','white','Position', [500, 500, 500, 400])
ph_list = [];
legend_list = {};

x_max = 0;
w_mean = zeros(2,n_num);
w_std = zeros(n_num,1);
ridge_x = zeros(n_num,1);
ridge_y = zeros(n_num,1);
colors = [1,0.8,0];
lbd = 1e10;
for i=1:n_num
    n = n_list(i);
    while sum(lbds./(lbds+lbd)) < n, lbd = lbd/2; end
    while sum(lbds./(lbds+lbd)) > n, lbd = lbd*1.01; end
    w_ridge = (Sgm+lbd*eye(d))\Sgm * w;

    for j=1:r   
        x_list = zeros(T,1);
        y_list = zeros(T,1);
        for t=1:T
            X = randn(n,d)*B;
            w_est = pinv(X)*X*w;
            x_list(t) = w_est(1);
            y_list(t) = w_est(2);
            %w_min(:,j) = w_min(:,j) + pinv(X)*X*w;           
        end
        w_mean(1,i) = mean(x_list);
        w_mean(2,i) = mean(y_list);
        w_std(i) = (std(x_list)+std(y_list))/sqrt(T);
        %w_min(:,j) = w_min(:,j)/T;
        x_max = max(x_max, abs(w_mean(1,i)));
        x_max = max(x_max, abs(w_mean(2,i)));
    end
    n_eff = sum(lbds./(lbds+lbd));
    ridge_x(i) = w_ridge(1);
    ridge_y(i) = w_ridge(2);
    %ph_list(end+1) = plot(w_ridge(1),w_ridge(2),40,colors*n_eff/d,'filled');
    %legend_list{end+1} = sprintf('n=%d',round(n_eff));
    %hold on
    rectangle('Position',[w_ridge(1)-w_std(i),w_ridge(2)-w_std(i),2*w_std(i),2*w_std(i)],...
        'Curvature',[0.3 0.3],'FaceColor',[.9,.9,.9,0.3],'EdgeColor',[1,1,1,0]);
    if n>80
        mid_x = (ridge_x(i-1)+ridge_x(i))/2;
        mid_y = (ridge_y(i-1)+ridge_y(i))/2;
        rectangle('Position',[mid_x-w_std(i),mid_y-w_std(i),2*w_std(i),2*w_std(i)],...
            'Curvature',[0.3 0.3],'FaceColor',[.9,.9,.9,0.3],'EdgeColor',[1,1,1,0]);
        if n>90
            mid_x = (mid_x+ridge_x(i))/2;
            mid_y = (mid_y+ridge_y(i))/2;
            rectangle('Position',[mid_x-w_std(i),mid_y-w_std(i),2*w_std(i),2*w_std(i)],...
                'Curvature',[0.3 0.3],'FaceColor',[.9,.9,.9,0.3],'EdgeColor',[1,1,1,0]);
        end
    end
end
hold on
ph_list(end+1) = plot(ridge_x,ridge_y,'LineWidth',2,'Color','black');
legend_list{end+1} = sprintf('ridge',round(n_eff));

for i=1:n_num
    n = n_list(i);
    ph = scatter(w_mean(1,i),w_mean(2,i),16,colors*(d-n)/d,'x');
    if i==1 || i==n_num || i==round(n_num/2)
        ph_list(end+1) = ph; 
        legend_list{end+1} = sprintf('n=%d',n);
    end    
end

xlim([0,x_max]);
ylim([0,0.5]);
legend(ph_list,legend_list,'Location','best');
hold off