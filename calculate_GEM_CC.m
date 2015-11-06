function [CC]= ...
	 calculate_GEM_CC(V,N,K,M)

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
