function [Priors Mu Mu_s Sigma Sigma_s] = ...
	 GEM_Max(V,share_exist,N,K,M);

global Data Pzi;
Mu = zeros(K,M);
Sigma = zeros(K,M);
Mu_s = zeros(1,K);
Sigma_s = zeros(1,K);
%% M-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E = sum(Pzi,1);

E(E<realmin)=realmin;
%component specific parameters
%Update the priors
Priors = E./N;
%Update the centers
Mu = bsxfun(@rdivide,Data'*Pzi,E);
for i=1:M
    temp = bsxfun(@minus,Data,Mu(:,i)');
    Sigma(:,i)=sqrt(Pzi(:,i)'*(temp.^2)./E(i));
end
Sigma(Sigma<1e-5) = 1e-5;

%Update shared parameters
Pzi_s = Pzi*(1-V');
Up = sum(Data.*Pzi_s,1);
Down = sum(Pzi_s,1);
Down(Down<realmin) = realmin;

Mu_s = Up./Down;

Up = sum(Pzi_s.*bsxfun(@minus,Data,Mu_s).^2,1);
Sigma_s = sqrt(Up./Down);
Sigma_s(Sigma_s<1e-5)=1e-5;


