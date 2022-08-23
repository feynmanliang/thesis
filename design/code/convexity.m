% nystrom symmetric polynomials

n=51;
T = 10;
R = 7;
sparsity = 1;
dataset = 5;
lower = [];

if dataset == 1
    X = randn(n,n);
elseif dataset == 2
    [U,~,V] = svd(randn(n));
    l = ones(n,1);
    l(floor(n/2)+1:n) = 0.01;
    X = U*diag(l)*V';
elseif dataset == 3
    [U,~,V] = svd(randn(n));
    l = ones(n,1);
    l(21:n) = 0.01;
    l(40:n) = 0.001;
    X = U*diag(l)*V';   
    %lower = [19,0.8*21;59,0.8*51]; 
elseif dataset == 4
    [U,~,V] = svd(randn(n));
    l = 1./(1:n).^2;
%    l(ceil(n/2):n) = 1/ceil(n/2);
    X = U*diag(l)*V';    
elseif dataset==5
    l=[ones(10,1);.1*ones(2,1)];
    l = l/sum(l);
    X = diag(sqrt(l));
    [n,~]=size(X);
end

k_list = 1:n-1;
k_num = length(k_list);
L = decompose_kernel(X*X');
lbds = sort(L.D,'descend');
theory = ones(k_num,1);
poly = elem_sympoly(lbds,n);


for i=k_list
    k = k_list(i);
    if k>=n, break; end
    %Ek = sum(lbds(k+1:n));
    theory(i) = ((k+1)/n)*poly(k+2,n+1)/poly(k+1,n+1);
end


% set(0,'DefaultAxesFontSize', 13);
% set(0,'DefaultTextFontSize', 13);
% figure('Color','white','Position', [1000, 1000, 500, 400])

ph_list = [];
legend_list = {};
%colors = {[1.0, 0.5, 0], [1, 0, 0], [0, 0.5, 0],[0,0,0.5]};
%types = {'-','--','-.',':'};

ph_list(end+1) = plot(k_list,theory,'-','LineWidth',2,'Color',[1,0,0,0.5]);
legend_list{end+1} = 'Expected error';
