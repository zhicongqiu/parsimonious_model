function [BIC_incomplete]= ...
	 calculate_GEM_incomplete(Priors,Pxi,Psi,V,N,K,M)

%calculate GEM BIC based on Graham, 2006
%number of free parameters = 2 for NB
%%CC-term
F_k = sum(V,2);
CC = (M-1)/2*log(N);
for i=1:K
    if F_k(i)==0
       CC = CC+log(N);
    elseif F_k(i)==M
      CC = CC+M*log(N);
    elseif F_k(i)<M
      CC = CC+(1+F_k(i))*log(N)+M*log(2);
    end
end

BIC_incomplete = zeros(N,1);

for i=1:M
  sum_numerator = zeros(N,1);
  for j=1:K
    if V(j,i)==1
      %component specific
      sum_numerator = sum_numerator.*Pxi(j).P(:,i);
    else
      %shared
      sum_numerator = sum_numerator.*Psi(:,j);
    end
  end
  sum_numerator = sum_numerator.*Priors(i);
  BIC_incomplete = BIC_incomplete+sum_numerator;
end

BIC_incomplete(find(BIC_incomplete<realmin))=realmin;
BIC_incomplete = CC-sum(log(BIC_incomplete));
