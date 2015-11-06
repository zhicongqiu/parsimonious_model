function [LIKK]= ...
	 calculate_GEM_complete_lik_parallelK(Data_k,Mu,Mu_s,...
					      Sigma,Sigma_s,V,N,k,M)

%calculate GEM BIC based on Graham, 2006
%fprintf('eliminate 2nd dimension\n');
%size(Pxi)
global Pzi;
%size(Psi)
if sum(V)==M
  temp=0;
  %PP = 0;
else
  %PP = -log(sqrt(2*pi)*Sigma_s)...
  %     -(Data_k-Mu_s).^2./(2*Sigma_s^2);
  temp = bsxfun(@times,-log(sqrt(2*pi)*Sigma_s)...
		       -(Data_k-Mu_s).^2./(2*Sigma_s^2),ones(N,M)).*Pzi*(1-V');
end

%P = bsxfun(@minus,-bsxfun(@rdivide,bsxfun(@minus,Data_k,Mu).^2,...
%			  2*Sigma.^2),log(sqrt(2*pi).*Sigma));
LIKK = ...
-sum(bsxfun(@minus,-bsxfun(@rdivide,bsxfun(@minus,Data_k,Mu).^2,...
			   2*Sigma.^2),log(sqrt(2*pi).*Sigma)).*Pzi*V'+temp);

%clear P;
clear temp;
