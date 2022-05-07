clear all
clc
close all
time_start_cpu=cputime;
time_start_real=tic;
mode_Confparam = 'automatic'; % 'manual' or 'automatic'
mode_K_plot='automatic';
mode_svd='automatic';

%%
%% Reading Experimental Data : I,q,Ierr  
%% Choosing configuration (see Configurations_avec_esrf.m)
%% Preprocessing : Selecting Dimer spectrum, Removing outliers
%%
load('sc3852_data3_ierrsample00243_00392_Buffer00233_00242_absunits_qAngstroms.mat')
Confparam=struct('q',{},'tConflength',{},'time',{},'Intensities',{},'concentration_gL',{},'tlength',{});
[Confparam,Nconf]=Configurations_avec_esrf(mode_Confparam,q_sample00243_00392_Buffer00233_00242,I_sample00243_00392_Buffer00233_00242,tsample00243_00392_Buffer00233_00242,IErr_sample00243_00392_Buffer00233_00242);
str = sprintf('Configuration is: Configuration %d',Nconf);
disp(str);
q=Confparam(Nconf).q;q=q(1:end);
t=Confparam(Nconf).time; t([1,3:40,48,50,52])=[];
I=Confparam(Nconf).Intensities;I=I(1:end,:); Pdimer=mean(I(:,end-4:end),2);Pdimer=Pdimer/3.8e-2;  I(:,[1,3:40,48,50,52])=[];
IErr=Confparam(Nconf).IntensitiesERR; IErr=IErr(1:end,:);IErr(:,[1,3:40,48,50,52])=[];

%% Fixing stoechiometric coefficients
alphaB=35;
betaB=2;
betaS=1;

%% Fixing parameters
param.Mdimer=40.37*10^3;%40.6kDa
param.Mcapsid=90*param.Mdimer;
param.c_gL=2.9;
c=param.c_gL/param.Mdimer*10^6;param.c=c;
disp('-----------')
struct_spectra=struct('q',q);
Volume=(7.91e+4)*10^-24; %cm^3
rhoprot=384e+21; %cm-3  
K_SI=calcul_K(Volume,rhoprot);
disp(sprintf('K calculated with V=7.91e+4A^{3} and rhoprot=384e+21cm{-3} (g^{-2}*mol*cm^{2})    = %5.3g',K_SI));

%% Calculation of the parameter K if needed, otherwise we fix it with appropriate value
[K,K_umol,KMdimer2]=K_plot(Nconf,Confparam,param,mode_K_plot,struct_spectra); 
K=2.9e-13;

%% SVD analysis                
[Ip,U,S,V,RankTrunc,Rank]=svd_analysis(mode_svd,q,I,IErr);

%% Form factors : Dimer, Big intermediate, Capsid
Bdim=K*(param.Mdimer)^2*Pdimer;
P=zeros(length(q),4);
P(:,1)=Bdim;
load('BigINT_formfactor.mat')
Bbigint=interp1(q_bigint_early_70mermasse_Rouge_betaB2,Pbigint_early_70mermasse_Rouge_betaB2,q);
P(:,3)=Bbigint*1.035;Bbigint=P(:,3)*(alphaB/35)^2;
load('Capsid_SVD_2015_03_05.mat')
I_capsidSVD=interp1(q_capsidSVD,I_capsidSVD,q);
Bcaps=I_capsidSVD*K*param.Mcapsid^2/(0.5*2.9)*4.3;
P(:,4)=Bcaps;

%% Other Parameters : Stoechiometric coefficients
alphaSstart=9:18;
field='betaS_numbers';
value={0;[10 24 30 45];[3 6 8 12 15 18 30];[1 2 5 8 10 12 15 20];...
    [2 5 8 10 12 14 16 18];[2 5 8 10 12 15];[2 5 8 10 12];[2 5 8 10];[betaS];...
    [betaS];[betaS];[betaS];...% 10,11,12
    [betaS];[betaS];[betaS];...%13,14,15
    [betaS];[betaS];[betaS];...%16,17,18
    [1];[1];[1];...%19-->21
    [1];[1];[1];[1];...%22-->25
    [1];[1];[1];[1];[1]};%26-->30
betaS_start=struct(field,value);

for i_index=1:length(alphaSstart)
alphaS=alphaSstart(i_index);
betaS_temp=betaS_start(alphaS).betaS_numbers;

%% Form factor : small intermediate
load('SmallInt_late_50uM_P12.mat')
Bsmallint=interp1(q_smallint_late_50um_P12,Psmallint_late_50um_P12,q,'spline');
P(:,2)=Bsmallint/0.0830222*K*param.Mdimer^2*alphaS^2;Bsmallint=P(:,2);

for i=1:length(betaS_temp)
betaS=betaS_temp(i);

%%
%% Global fitting
%%

%% Initialisation of parameters
param.C0_dimerinit=0.5*2.9/param.Mdimer*1e+6;
matrix=zeros(4,length(t));
matrix(1,:)=param.C0_dimerinit*ones(1,length(t));

    kfwd1max=0.2;
    kfwd1min=0.0918;
    kback1max=1e-6;
    kback1min=1e-9;
    kfwd2max=1e-1;
    kfwd2min=1e-3;
    kback2max=1e-6;
    kback2min=1e-9;
    kfwd3max=3e-3;
    kfwd3min=1.5e-3;
    kback3max=1e-6;
    kback3min=1e-9;
    lb1=[kfwd1min;kback1min;kfwd2min;kback2min;kfwd3min;kback3min]; lb1=log10(lb1);
    ub1=[kfwd1max;kback1max;kfwd2max;kback2max;kfwd3max;kback3max];ub1=log10(ub1);
    lb=[lb1];ub=[ub1];
    
    Number_x0=3;
    x0_four_state=zeros(size(lb,1),Number_x0);r=zeros(size(lb));
    for j=1:Number_x0
    for i=1:size(lb)
    r(i) = lb(i) + (ub(i)-lb(i)).*rand(1,1);
    end
    x0_four_state(:,j)=r;
    end
    
    lb1=[kfwd1min;kback1min;kfwd2min;kback2min;kfwd3min;kback3min]; lb1=log10(lb1);
    ub1=[kfwd1max;kback1max;kfwd2max;kback2max;kfwd3max;kback3max];ub1=log10(ub1);
    lb=[lb1];ub=[ub1];
    
    kfwd1_log=zeros(1,Number_x0);
    kback1_log=zeros(1,Number_x0);
    kfwd2_log=zeros(1,Number_x0);
    kback2_log=zeros(1,Number_x0);
    kfwd3_log=zeros(1,Number_x0);
    kback3_log=zeros(1,Number_x0);
    Chi2temp=zeros(1,Number_x0);
    RFactortemp=zeros(1,Number_x0);
    FEVAL=zeros(1,Number_x0);
    EXITFLAG=zeros(1,Number_x0);

    for k=1:Number_x0
        FixPar_four_state=[0;0;0;0;0;0;alphaS;alphaB;betaS;betaB];
    options = optimoptions(@fmincon,'Algorithm','sqp');
    [x1,FEVAL(k),EXITFLAG(k),OUTPUT] = fmincon(@(x1) fit_Bf4species_Pintextracted(x1,q,I,IErr,t,FixPar_four_state,K,param,matrix,P,'global'),...
        x0_four_state(:,k),[],[],[],[],lb,ub,[],options);
    iPar = find(FixPar_four_state == 0);   FixPar_four_state(iPar) = x1;
    c0 = [0;0;0;0.5*2.9/param.Mcapsid*1e+6];  
    kfwd_temp = [FixPar_four_state(1);FixPar_four_state(3);FixPar_four_state(5)];
    kback_temp = [FixPar_four_state(2);FixPar_four_state(4);FixPar_four_state(6)];
    alphaS=FixPar_four_state(7);alphaB=FixPar_four_state(8);
    betaS=FixPar_four_state(9);betaB=FixPar_four_state(10);
    
    kfwd1_log(k)=kfwd_temp(1);kfwd2_log(k)=kfwd_temp(2);kfwd3_log(k)=kfwd_temp(3);
    kback1_log(k)=kback_temp(1);kback2_log(k)=kback_temp(2);kback3_log(k)=kback_temp(3);
    kfwd(1)=10.^kfwd1_log(k);kfwd(2)=10.^kfwd2_log(k);kfwd(3)=10.^kfwd3_log(k);
    kback(1)=10.^kback1_log(k);kback(2)=10.^kback2_log(k);kback(3)=10.^kback3_log(k);

%%
%% The optimal parameters have been found. Now we compute matrices and figures of merits.
%%     
C = FOURstate_Disassembly(t,param,matrix,c0,kfwd,kback,alphaS,alphaB,betaS,betaB);

B=P;
    Chi2temp(k) = normmat((I-B*C)./IErr)^2/(numel(I)-1);
    RFactortemp(k) = sum(sum(abs(I-B*C)))/sum(sum(abs(I)));
  
    kfwd1(k)=kfwd(1);kfwd2(k)=kfwd(2);kfwd3(k)=kfwd(3);
    kback1(k)=kback(1);kback2(k)=kback(2);kback3(k)=kback(3);
    end
    disp(sprintf('Chi-square    = %5.3g',Chi2temp));
    disp(sprintf('R-Factor      = %5.3g',RFactortemp));
    
index_min_n=find(Chi2temp==min(Chi2temp),1);
Chi2=Chi2temp(index_min_n);RFactor=RFactortemp(index_min_n);
kfwd1_fin=kfwd1(index_min_n);kfwd2_fin=kfwd2(index_min_n);kfwd3_fin=kfwd3(index_min_n);
kback1_fin=kback1(index_min_n);kback2_fin=kback2(index_min_n);kback3_fin=kback3(index_min_n);
clear kfwd kback
kfwd=[kfwd1_fin;kfwd2_fin;kfwd3_fin]
kback=[kback1_fin;kback2_fin;kback3_fin]


PintBig=Bbigint;
PintSmall=Bsmallint;
B=horzcat(Bdim,PintSmall,PintBig,Bcaps);

%%
%% We save what we need.
%%

save(['Output_corrMasseSpectresfixes_kfwdetkbacklibre_2015_03_qmin4,2_qmax0,1_alphaB',num2str(alphaB),'betaB',num2str(betaB),'alphaS',num2str(alphaS),'betaS',num2str(betaS),'.mat'],...
    'K','t','IErr','I','q','kfwd','kback','kfwd1_fin','kfwd2_fin','kfwd3_fin','kback1_fin','kback2_fin','kback3_fin','B','Chi2','RFactor')
dataIntSmall=horzcat(q,B(:,2));
save(['PintSmall_corrMasseSpectresfixes_kfwdetkbacklibre_qmin4,2_2015_03_alphaB',num2str(alphaB),'betaB',num2str(betaB),...
    'alphaS',num2str(alphaS),'betaS',num2str(betaS),'.txt'],'dataIntSmall','-ascii');       
dataIntBig=horzcat(q,B(:,3));
save(['PintBig_corrMasseSpectresfixes_kfwdetkbacklibre_qmin4,2_2015_03_alphaB',num2str(alphaB),'betaB',num2str(betaB),...
    'alphaS',num2str(alphaS),'betaS',num2str(betaS),'.txt'],'dataIntBig','-ascii');       
dataCapsid=horzcat(q,B(:,4));
save(['PcapsSVD_corrMasseSpectresfixes_kfwdetkbacklibre_qmin4,2_2015_03_alphaB',num2str(alphaB),'betaB',num2str(betaB),...
    'alphaS',num2str(alphaS),'betaS',num2str(betaS),'.txt'],'dataCapsid','-ascii');       

     
close all
end
end
running_time_cpu=cputime-time_start_cpu;
running_time_real=toc(time_start_real)
save('running_time.mat','running_time_real')