function [Ip,U,S,V,RankTrunc,Rank]=svd_analysis(mode,q,I,IErr)
TR=strcmp(mode,'automatic');

if TR==1
    RankTrunc =10;
else

    RankTrunc = input('Rank of the truncated matrix [10]: ');
    if isempty(RankTrunc)
        RankTrunc = 10;
    end
    if RankTrunc > 10
        disp('Rank limited to 10!');
        RankTrunc = 10;
    end
end  

disp('SVD analysis');
disp('------------');
Noise = mean(mean(IErr));
disp(sprintf('Mean error = %5.3g',Noise));
Nq = size(I,1); Nt = size(I,2);
[U,S,V] = svd(I);
S = diag(S);
    % Display singular values and vectors information 
    disp('Singular values:');
    disp(S(1:5)');
    CorrU = diag(U(1:(Nq-1),:)'*U(2:Nq,:));
    CorrV = diag(V(1:(Nt-1),:)'*V(2:Nt,:));
    disp('Autocorrelations of singular vectors:');
    disp('          U           V');
    disp([CorrU(1:5) CorrV(1:5)]);   
% Residual computation and display    
Rank = length(S); ResErr = fliplr(sqrt(cumsum(fliplr(S'.^2))));
NoiseLev = Noise*sqrt((Nq-(1:Rank)).*(Nt-(1:Rank)));

figure
subplot(2,1,1);
semilogy(S,'o');
xlim([0 10]);
ylabel('Singular value');
subplot(2,1,2);
semilogy(1:Rank,ResErr,'b+-',1:Rank,NoiseLev,'r');
xlim([0 10]);
xlabel('Channel');
ylabel('Residual');

    
    [U,S,V]=svds(I,RankTrunc);
    Up = U; Vp = (S*V')';
    Ip = Up*Vp';
    if ~exist('Bf_3species') 
    Bf_3species = zeros(length(q),3);
    end
    i = find(Bf_3species(1,:) ~= 0);
    Cr = mldivide(Up,Bf_3species(:,i)); Cf = zeros(size(Cr,1),size(Bf_3species,2));
    Cf(:,i) = Cr;
%     for k = 1:RankTrunc
%     subplot(RankTrunc,2,2*k-1); semilogx(q,U(:,k),'b');
%     subplot(RankTrunc,2,2*k);   plot(t,V(:,k),'b');
%     end
end

