function I = GQ(func,d,n)
%GQ Compute the Gaussian quadrature rule for region Omega=[0,1]^d
%
%CALL:  I = GQ(func,d,n)
%  I = Gaussian quadrature approximation
%  func = function of interest
%  d = dimension of the region
%  n = number of abscissas along one dimension
assert(d==1|d==2,"Invalid dimension")

[x,w] = qrule(n);
x = 1/2*(x+1);
w = 1/2*w;

if d==1
    % fill in here
    I = sum(w .* func(x));
else
    [y,v] = qrule(n);
    y = 1/2*(y+1);
    v = 1/2*v;
    
    % fill in here    
    I=0;
    for i=1:n
        for j=1:n
            I = I + func(x(i), y(j)) * w(i) * v(j);
        end
    end
    
end

end