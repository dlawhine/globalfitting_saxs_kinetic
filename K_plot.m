function [K,K_umol,KMdimer2]=K_plot(Nconf,Confparam,param,mode,struct_spectra)
q=Confparam(Nconf).q;
t=Confparam(Nconf).time; 
I=Confparam(Nconf).Intensities;
c=param.c;


TR=strcmp(mode,'automatic');
if TR==1
K=2.7785e-13;
K_umol=K;
KMdimer2=[];
else

%% Nconf=1
if Nconf==1
K=2.7785e-13;
K_umol=K;
KMdimer2=[];
end

%% Nconf=2
if Nconf==2
figure
plot(q.^2,log(I(:,1)/c(1)),'r.-');
title('log(I(t=0))=f(q^{2}) Dimer')
ylabel('log(I(t=0)) Dimer')
xlabel('q^{2} in A^{-1}')


% On calcule KM� � partir de la valeur de I0 lue sur le graphe
log_I0 = input('Value of log(I0) read in the graph [by default : -7.7]: ');        
if isempty(log_I0)
    log_I0 = -7.7;
end

I0=exp(log_I0);
KMdimer2=I0;
disp('K=K_umol in unit(I)*g^{-2}*mol^{2}*µmol^{-1}*L')
K_umol=KMdimer2/((param.Mdimer)^2);
K=K_umol;
end

end
end
