function [Prior0 Mu0 Sigma0 V share_exist] = ...
	 delete_strategy(Priors,Mu,Mu_s,Sigma,Sigma_s,V,share_exist)
%need to delete one component
%strategy 1: delete component with the least mass

%strategy 1
fprintf('call delete strategy\n')
[a index] = min(Priors);
fprintf('%d th component deleted.\n',index)
Prior0 = [Priors(1:index-1) Priors(index+1:end)];
Mu0 = [Mu(:,1:index-1) Mu(:,index+1:end)];
Sigma0 = [Sigma(:,1:index-1) Sigma(:,index+1:end)];
V = [V(:,1:index-1) V(:,index+1:end)];
%re-estimate share_exist
M = size(Mu0,2);
for i=1:size(V,1)
  if sum(V(i,:))<M
    share_exist(i)=1;
  else
    share_exist(i)=0;
  end
end
%reevaluate the priors
Prior0 = Prior0./sum(Prior0);

%Other strategies later


