function J = vgfd_J(x,w,M)
dx = x(2) - x(1);
J = zeros(size(w));
N = length(w);
for i=1:N
    for m=-M:M
        if m == -M || m == M || m == -1 || m == 1
            rho = 0.5;
        else
            rho = 1;
        end
        if (i+m) < 1 || (i+m) > N
            J(i) = J(i) + vgfd_psi(x(i)+m*dx)*vgfd_k(m*dx)*rho;
        else
            J(i) = J(i) + w(i+m)*vgfd_k(m*dx)*rho;
        end        
    end
end
J = dx*J;
