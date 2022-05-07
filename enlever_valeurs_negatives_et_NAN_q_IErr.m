function [qout,IErrout,messageerreur]=enlever_valeurs_negatives_et_NAN_q_IErr(q,IErr)
Lcolumn=size(IErr,2);
Lrow=size(IErr,1);
%Valeurs n�gatives
for j=1:Lcolumn
for i=1:Lrow
    if IErr(i,j)<0
        IErr(i,j)=0;
    end
end
end


%Valeurs NaN
k=1;
irepere=zeros(1);   
for j=1:Lcolumn
for i=1:Lrow
    if isnan(I(i,1))==1 
        irepere(k)=i;
        k=k+1; 
    end
end
    tailleirepere(j)=length(irepere);
end
flag=0;

for j=2:Lcolumn
    if tailleirepere(j)~=tailleirepere(j-1)+tailleirepere(1)
    flag=1;       
    end
end

if flag==1
    afficher='Erreur 1? : Pb au niveau du nombre de Nan';
else
    afficher='Erreur 1 ? : no error';
end
indicelignesNAN=irepere(1:tailleirepere(1));


in=IErr;
notnan = sum(arrayfun(@(x) ~isnan(x),in(:,1:end)),2); %ftp://ftp-developpez.com/grin/intro%20matlab%20v1.3.3.pdf
out = in(notnan~=0,:);


%% On pense � enlever les valeurs de q et IErr dont l'intensite est NaN.
for i=1:Lrow
    for k=1:length(indicelignesNAN)
    if i==indicelignesNAN(k)
        q(i)=NaN;
    end
    end
end

in2=q;
in2=in2';
notnan2 = sum(arrayfun(@(x) ~isnan(x),in2(:,1:end)),2); %ftp://ftp-developpez.com/grin/intro%20matlab%20v1.3.3.pdf
out2 = in2(notnan2~=0,:);

IErrout=out;
qout=out2;
messageerreur=afficher;

end