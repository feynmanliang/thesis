d=100;
m=20;
decay = 3; %
decay_list = {'exponential decay', 'one-step drop','mixed','linear decay','polynomial'};

alpha_list = logspace(-10,1,100)';
alpha_num = length(alpha_list);
ratios = zeros(alpha_num,1);
ks = zeros(alpha_num,1);
delta_list = [0.01];
theory_ratios = zeros(d,1);


set(0,'DefaultAxesFontSize', 15);
set(0,'DefaultTextFontSize', 15);
figure('Color','white','Position', [500, 500, 500, 400])
ph_list = [];
legend_list = {};
colors = {[0 0 1],[1 0 0],[0 0.8 0]};
for delta=1:length(delta_list)
    if decay==1        
        lbd = 0.9.^(1:d)';
    elseif decay==2
        lbd = [ones(1,m),delta_list(delta)*ones(1,d-m)]';
    elseif decay==3
        lbd = [ones(1,m),0.9.^(m+1:2*m),0.1*0.9^(2*m)*ones(1,d-2*m)]';
    elseif decay==4
        m=20;
        lbd = [ones(1,m),1e-1*ones(1,m),1e-2*ones(1,m),5e-4*ones(1,2*m)]';
    elseif decay==5
        lbd = (1:d).^(-1/4);
    end
    sr_list = [];
    next = 1;
    while next < d        
        sr_list(end+1) = next-1 + sum(lbd(next:end))/lbd(next);
        next = ceil(sr_list(end));
    end
    for i=1:alpha_num
        alpha  = alpha_list(i);
        ks(i) = sum(lbd./(lbd+alpha));
        kceil = ceil(ks(i));
        Ek = sum(lbd(kceil:end));
        if kceil > 1
             Ek = Ek + (kceil-ks(i))*lbd(kceil-1);
        end        
        ratios(i) = alpha*ks(i)/Ek;
    end
    for k=1:floor(sr)-1
        theory_ratios(k) = 1 + k/(sr-k);
    end
    for k=floor(sr):d
        theory_ratios(k)=k+1;
    end

    ph_list(end+1) = plot(ks,ratios,'LineWidth',1,'Color',colors{delta});
    yylim = max(ratios);
    legend_list{end+1} = sprintf('\\delta=%.3f',delta_list(delta));
    hold on
    %semilogx(1:d,theory_ratios,'--','LineWidth',1,'Color',colors{delta});    
    %scatter(sr,0,20,colors{delta},'LineWidth',2,'MarkerEdgeAlpha',0.5);
    for sr=sr_list(1:end-1)
        line([sr,sr],[0,yylim],'Color','red','LineWidth',2,'LineStyle',':');
    end
    %line([sr*2/3,sr*2/3],[2,0],'Color',colors{delta},'LineWidth',1,'LineStyle',':');
end


%legend(ph_list,legend_list,'Location','best');
xlim([1,d])
ylim([0,yylim]);
ylabel('Approximation ratio');
xlabel('Value of k');
title('');
%title(decay_list{decay});
hold off