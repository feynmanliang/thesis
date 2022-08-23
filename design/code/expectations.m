d = 5;
n = 5;
X = randn(n,d);
[~,D,V] = svd(X,'econ');
D = diag(1:d);%eye(d);%diag(logspace(-1,0,d));
X = U*D*V;

Einv1 = 0;
Einv2 = 0;
Ek = 0;
Zval = 0;

Z = det(eye(d) + X'*X);

for k=0:d
    S = nchoosek(1:n,k);
    for i=1:nchoosek(n,k)
        XS = X(S(i,:),:);
        PS = det(XS*XS')/Z;
        Einv1 = Einv1 + PS * trace(pinv(XS'*XS)); 
        Einv2 = Einv2 + PS * trace(X*pinv(XS'*XS)*X');
        Ek = Ek + PS * k;
        Zval = Zval + PS;
    end
end

Zval

Ek
Thk = trace(X'*X*inv(eye(d)+X'*X))
Diffk = Ek - Thk

Einv1
Th1 = trace(inv(eye(n)+X*X')) - (n-d)*det(X'*X)/Z
Diff1 = Einv1 - Th1

Einv2