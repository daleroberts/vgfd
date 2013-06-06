function z = vgfd_psi(x)
global K
z = max(exp(x) - K,0);
