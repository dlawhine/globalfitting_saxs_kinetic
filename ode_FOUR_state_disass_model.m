function dcdt = ode_FOUR_state_disass_model(t,c,kfwd,kback,alphaS,alphaB,betaS,betaB)

%   C <-> betaB*Ib+gammaB*D   kfwd(1), kback(1)
%   Ib <-> betaS*Is+gammaS*D  kfwd(2), kback(2)
%   Is <-> alphaS*D    kfwd(3), kback(3)

%c(1)=dimer,c(2)=intermediaireSMALL,c(3)=intermediaireBIG,c(4)=capsid

gammaB=90-alphaB*betaB;
gammaS=alphaB-alphaS*betaS;
mu_gammaB=1;
mu_gammaS=1;
mu_alphaS=1;
mu_betaS=2;
mu_betaB=2;



dcdt = zeros(4,1);
dcdt(1) = alphaS*(kfwd(3)*c(2)-kback(3)*c(1)^mu_alphaS)+...
    gammaS*(kfwd(2)*c(3)-kback(2)*c(1)^mu_gammaS*c(2)^mu_betaS)+...
    gammaB*(kfwd(1)*c(4)-kback(1)*c(1)^mu_gammaB*c(3)^mu_betaB);
dcdt(2) = betaS*(kfwd(2)*c(3)-kback(2)*c(1)^mu_gammaS*c(2)^mu_betaS)+...
    (kback(3)*c(1)^mu_alphaS-kfwd(3)*c(2));
dcdt(3) = betaB*(kfwd(1)*c(4)-kback(1)*c(1)^mu_gammaB*c(3)^mu_betaB)+...
    (kback(2)*c(1)^mu_gammaS*c(2)^mu_betaS-kfwd(2)*c(3));
dcdt(4) = (-kfwd(1)*c(4)+kback(1)*c(1)^mu_gammaB*c(3)^mu_betaB);

end