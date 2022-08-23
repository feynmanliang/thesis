function val = bound( lbd , c)
    if nargin == 1, c = sqrt(3); end
    lbd = abs(lbd);
    lbd = sort(lbd,'descend');
    d = length(lbd);
    r = sum(lbd) / lbd(1);
    Er = sum(lbd(1:d>=r/c));
    %if ceil(r/c) >= 2, Er = Er + (ceil(r/c)-(r/c))*lbd(ceil(r/c)-1); end
    alpha = c^2 * Er / r;
    k = sum(lbd./(lbd+alpha));
    val = -k*c/r;
    %val = - alpha*k/(c*Er); 