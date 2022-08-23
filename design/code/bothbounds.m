
r = 30; 

all_eps =  0.01:0.001:10;
num_eps = length(all_eps);
all_ub = zeros(num_eps,1); 
all_lb = zeros(num_eps,1);
all_k = zeros(num_eps,1);


for i = 1:num_eps
    eps = all_eps(i)
    k = eps/(1+eps) * r;
    all_k(i) = k; 
    all_ub(i) = 1 + 3*eps;
    all_lb(i) = 1 + eps * (k-1)/ (2*k*k);
end   
sel_k = all_k > 1;
set(0,'DefaultAxesFontSize', 15);
set(0,'DefaultTextFontSize', 15);
figure('Color','white','Position', [500, 500, 500, 400])


%    ph_list(end+1) = semilogx(ks,ratios,'LineWidth',1,'Color',colors{delta});
%    legend_list{end+1} = sprintf('\\delta=%.3f',delta_list(delta));
    hold on
    %semilogx(1:d,theory_ratios,'--','LineWidth',1,'Color',colors{delta});    
    %scatter(sr,0,20,colors{delta},'LineWidth',2,'MarkerEdgeAlpha',0.5);
    plot (all_k(sel_k), all_lb(sel_k), '--r', all_k(sel_k), all_ub(sel_k), '--b');

    x2 = [all_k(sel_k), fliplr(all_k(sel_k))];
inBetween = [all_lb(sel_k), fliplr( all_ub(sel_k))];
%patch(x2, inBetween, 'g');

%legend(ph_list,legend_list,'Location','best');
xlim([0,r])
%title(decay_list{decay});
hold off