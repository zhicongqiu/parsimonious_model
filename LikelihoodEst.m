function [Pzi] = LikelihoodEst(Priors,Mu,Mu_s,Sigma,Sigma_s,V,N,K,M)

global Data;
  %compute Pz in each dimension and in each component
  Pzi = zeros(N,M);
  sum_numerator = zeros(N,M);
  for i=1:M
    for j=1:K
	if V(j,i)==1
	  %component specific
	  sum_numerator(:,i) = ...
	  sum_numerator(:,i)...
	  -log(sqrt(2*pi)*Sigma(j,i))...
	  -(Data(:,j)-Mu(j,i)).^2./(2*Sigma(j,i)^2);
	else
	  %shared
	  sum_numerator(:,i) = ...
	  sum_numerator(:,i)...
	  -log(sqrt(2*pi)*Sigma_s(j))...
	  -(Data(:,j)-Mu_s(j)).^2./(2*Sigma_s(j)^2);
	end

    end
  end
  temp_mult = max(sum_numerator')';
  sum_numerator = ...
  bsxfun(@plus,sum_numerator,-temp_mult);
  %disp(min(sum_numerator)')
  sum_numerator = ...
  bsxfun(@plus,sum_numerator',log(Priors'))';
  temp = exp(sum_numerator);
  temp(temp<realmin)=realmin;
  Px = sum(temp,2);
  for i=1:M
    Pzi(:,i) = temp(:,i)./Px;
  end

