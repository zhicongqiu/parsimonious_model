function [error] = ...
calculate_error(Priors,Mu,Mu_s,Sigma,Sigma_s,...
			 share_exist,V,N,K,num_comp_least)

global Data Label;% Pxi;
[Pzi] = ...
LikelihoodEst(Priors,Mu,Mu_s,Sigma,...
     Sigma_s,V,N,K,num_comp_least);
[Pzi_max Pzi_hard] = max(Pzi');

owner = zeros(1,length(Priors));
freq = zeros(1,length(Priors));
num_total = zeros(1,length(Priors));
error = 0;
for i=1:length(Priors)
  temp = Label(find(Pzi_hard==i));
  if ~isempty(temp)
    [owner(i) freq(i)]= mode(temp);
    num_total(i) = length(temp);
    error = error+num_total(i)-freq(i);
  end
end
  
%the number of unique elements
%[freq element]=hist(Pzi_hard,unique(Pzi_hard));
error = error/N;
%{
ERROR.error = error;
ERROR.owner = owner;
ERROR.num = num_total;
ERROR.freq = freq;
%}
