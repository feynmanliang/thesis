d_list = [20,25,30,35,40];
alpha_list = 0.2;
d_num = length(d_list);
alpha_num = length(alpha_list);
T = 40000;
R = 4;
proj_error = zeros(d_num,alpha_num);

proj_std = zeros(d_num,alpha_num);


for i=1:d_num
    d = d_list(i);
    Sgm = eye(d); Sgm(1,1)=100000000;%diag((d:-1:1).^(-2));% eye(d);%=diag(logspace(-2,0,d));%
    B = sqrt(Sgm);
    lbds = eig(Sgm);
    w = zeros(d,1);
    w(1) = 1;
    
    lbd = 1e10;
    for j=1:alpha_num
        alpha = alpha_list(j);
        n = floor(alpha*d);
        fprintf('d=%d, n=%d\n',d,n);
        while sum(lbds./(lbds+lbd)) < n, lbd = lbd/(1+1e-3); end
        while sum(lbds./(lbds+lbd)) > n, lbd = lbd*(1+1e-6); end
        while sum(lbds./(lbds+lbd)) < n, lbd = lbd/(1+1e-9); end
        %sum(lbds./(lbds+lbd))
        SqSgmlbd = sqrt(eye(d) + Sgm/lbd);    
        theory = w'*inv(eye(d) + Sgm/lbd)*w;
        proj_global = 0;
        %proj_global = zeros(d,d);
        cnt = 0;
        TT = T;
        for r=1:R
            fprintf('...sampling T=%d...\n',TT);
            btstrp_l = randsample(TT,TT,true); 
            btstrp_c = zeros(TT,1);
            for t=1:TT
                btstrp_c(btstrp_l(t)) = btstrp_c(btstrp_l(t)) + 1;
            end
            proj_mean = 0;
            btstrp_mean = 0;
            %proj_mean = zeros(d,d);
            %btstrp_mean = zeros(d,d);
            for t=1:TT
                cnt = cnt + 1;
                X = randn(n,d)*B;                
                M = eye(d) - pinv(X)*X;
                Bias = w'*M*w;
                proj_mean = proj_mean + Bias; %M
                proj_global = proj_global + Bias; %M
                btstrp_mean = btstrp_mean + btstrp_c(t) * Bias; %M
            end
            btstrp_mean = btstrp_mean / TT;
            proj_mean = proj_mean / TT;
            proj_std(i,j) = abs(proj_mean - btstrp_mean) / theory;
            proj_error(i,j) = abs((proj_mean / theory) - 1); 
            %proj_std(i,j) = norm(SqSgmlbd*(proj_mean-btstrp_mean)*SqSgmlbd);
            %proj_error(i,j) = norm(SqSgmlbd*proj_mean*SqSgmlbd - eye(d));
            if proj_std(i,j) < proj_error(i,j)/5
                break;
            else
                TT = TT * 2;
            end
        end
        proj_global = proj_global / cnt;
        proj_error(i,j) = abs((proj_global / theory) - 1); 
        %proj_error(i,j) = norm(SqSgmlbd*proj_global*SqSgmlbd - eye(d));
        fprintf('...done.\n');
    end
end
    
set(0,'DefaultAxesFontSize', 15);
set(0,'DefaultTextFontSize', 15);
figure('Color','white','Position', [500, 500, 500, 400])
ph_list = [];
legend_list = {};
colors = {[1.0, 0.5, 0], [1, 0, 0], [0, 0, 0]};

for j=1:alpha_num
    ph_list(end+1) = errorbar(d_list, proj_error(:,j),proj_std(:,j),'LineWidth',2,'Color',colors{j});
    legend_list{end+1} = sprintf('n/d=%.2f',alpha_list(j));
    hold on
end
set(gca, 'XScale','log', 'YScale','log')
xlabel('Dimension d')
ylabel('Normalized matrix norm difference');
legend(ph_list,legend_list,'Location','best');
hold off
