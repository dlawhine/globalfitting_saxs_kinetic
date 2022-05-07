function B=construction_B(K,param,n,Bdim,Bcaps,FF)


B(:,1)=Bdim;
B(:,3)=Bcaps;

B(:,2)=K*(n*param.Mdimer)^2*FF;


end