function [Confparam,Nconf,TR]=Configurations_avec_esrf(mode,q_sc3852,I_sc3852,t_sc3852,IErr_sc3852)
%% List of the different configurations


disp('-----')
disp('Configuration 1: ')
disp('q=q(1)-->q(end):4e-3A-1-->0.4A-1')
q_1=q_sc3852(1:705);I_1=I_sc3852(1:705,:);IErr_1=IErr_sc3852(1:705,:);
Confparam(1).q=q_1;
Confparam(1).time=t_sc3852;
Confparam(1).tConflength=length(t_sc3852);
Confparam(1).concentration_gL=2.9;
Confparam(1).Intensities=I_1;
Confparam(1).IntensitiesERR=IErr_1;

disp('-----')
disp('Configuration 2: ')
disp('q=q(1)-->q(350):4e-3A-1-->0.2A-1')
q_2=q_sc3852(1:350);I_2=I_sc3852(1:350,:);IErr_2=IErr_sc3852(1:350,:);
Confparam(2).q=q_2;
Confparam(2).time=t_sc3852;
Confparam(2).tConflength=length(t_sc3852);
Confparam(2).concentration_gL=2.9;
Confparam(2).Intensities=I_2;
Confparam(2).IntensitiesERR=IErr_2;

disp('-----')
disp('Configuration 3: ')
disp('q=q(1)-->q(350):4e-3A-1-->0.1A-1')
q_2=q_sc3852(1:175);I_2=I_sc3852(1:175,:);IErr_2=IErr_sc3852(1:175,:);
Confparam(3).q=q_2;
Confparam(3).time=t_sc3852;
Confparam(3).tConflength=length(t_sc3852);
Confparam(3).concentration_gL=2.9;
Confparam(3).Intensities=I_2;
Confparam(3).IntensitiesERR=IErr_2;
%% Choice of Nconf according to mode (automatic or manual)
disp('-----')
TR=strcmp(mode,'manual');
if TR==1
Nconf = input('Configuration [by default : 3]: ');
  if isempty(Nconf)
    Nconf = 3;
  end
end

if TR==0 
    Nconf = 3;  
end


end
