function K_SI=calcul_K(V,rhoprot)
V=(7.91e+4)*10^-24   %cm^3
Na=6.02e+23  
rel=0.28e-12 %cm
M=40370       %g/mol
rhosolvant=344e+21;  %cm-3
rhoprot=384e+21; %cm-3

delta_rho=rhoprot-rhosolvant;
K_SI=((V*Na/M*rel*delta_rho)^2)/Na %g{-2}*mol*cm{2}
end