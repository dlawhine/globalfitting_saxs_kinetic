function Res = fit_Bf4species_Pintextracted(x,q,I,IErr,t,par,K,param,matrix,P,method)

    i = find(par == 0);
    par(i) = x;
    c0 = [0;0;0;0.5*2.9/param.Mcapsid*1e+6];  kfwd_log = [par(1);par(3);par(5)]; kback_log = [par(2);par(4);par(6)];
    alphaS= par(7);alphaB=par(8);betaS=par(9);betaB=par(10);

    kfwd=10.^kfwd_log;kback=10.^kback_log;
    
C = FOURstate_Disassembly(t,param,matrix,c0,kfwd,kback,alphaS,alphaB,betaS,betaB);
%     i = find(P(1,:) ~= 0); j = find(P(1,:) == 0);
%     Ir = I-P(:,i)*C(i,:);  Cr = C(j,:);
%     Pr = fit_basis_spectra(Ir,Cr,method);
%     B=P; B(:,j) = Pr;
%     
    B=P;
    Res = normmat((I-B*C)./IErr);
end


