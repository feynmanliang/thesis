d=10;
options = optimset('PlotFcns',@optimplotfval,'MaxIter',5000);
eps = 0.01;

lbd = sort(abs(randn(1,d)))

for i=1:2
    [lbd,val] = fminsearch(@bound,lbd,options);
    i, -val %, sort(abs(lbd),'descend')
    lbd = lbd + eps*randn(1,d);
end