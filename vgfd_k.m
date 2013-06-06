function k = vgfd_k(x)
global theta nu sigma
lambdap=1/(sqrt(theta^2*nu^2/4+sigma^2*nu/2)+theta*nu/2);
lambdan=1/(sqrt(theta^2*nu^2/4+sigma^2*nu/2)-theta*nu/2);
[m n] = size(x);
k = zeros(m,n);
for j = 1:n
    for i = 1:m
        if x(i,j) > 0
            k(i,j) = 1/nu*exp(-lambdap*abs(x(i,j)))/abs(x(i,j));
        elseif x(i,j) < 0
            k(i,j) = 1/nu*exp(-lambdan*abs(x(i,j)))/abs(x(i,j));
        else
            k(i,j) = 0;
        end
    end
end