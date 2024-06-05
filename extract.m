function [rho, kperp, drho] = extract(variable)
rho = variable.out.shot57558.d19.rhodif;
kperp = variable.out.shot57558.d19.k_perp;
drho = variable.out.shot57558.d19.drho;
end