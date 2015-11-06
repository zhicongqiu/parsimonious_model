function [V_least share_exist_least BIC_least] = ...
	 configureSwitch(Data_k,Mu,Mu_s,Sigma,Sigma_s,...
			 share_exist,V_trial,j,N,M,updated)

global Pzi;
%configure switch on the kth feature
BIC_least = Inf;
case_num=0;
%
%BIC_least = BIC_old;
V_least = V_trial;
share_exist_least = share_exist;

if updated==0
  BIC_least = ...
  calculate_GEM_complete_lik_parallelK(Data_k,Mu,Mu_s,...
				       Sigma,Sigma_s,V_trial,N,j,M);
  
else
  %assume P=2 (naive bayes)
  %V_trial = V;
  %case 1: all switches are off
  V_trial(:) = 0;
  %recalculate input likelihood
  BIC_new = log(N) + ...
	    calculate_GEM_complete_lik_parallelK(Data_k,Mu,Mu_s,...
						 Sigma,Sigma_s,...
						 V_trial,N,j,M);
  %fprintf('case 1 BIC is %f\n',BIC_new);
  if BIC_new<BIC_least
    case_num=1;
  V_least = V_trial;
  BIC_least = BIC_new;
  share_exist_least = 1;
  end
  
  %case 2: all switches are on
  V_trial(:) = 1;
  %recalculate input likelihood
  BIC_new = M*log(N) + ...
	    calculate_GEM_complete_lik_parallelK(Data_k,Mu,Mu_s,...
						 Sigma,Sigma_s,...
						 V_trial,N,j,M);
  %fprintf('case 2 BIC is %f\n',BIC_new);
  
  if BIC_new<BIC_least
    case_num=2;
    V_least = V_trial;
    BIC_least = BIC_new;
    share_exist_least = 0;
  end

  %case 3:
  %configure switch based on the updated parameters
  sum_share = 0;
  for i=1:M
    %try each switch configuration, make sure Psi exists!!
      temp1 = -log(sqrt(2*pi)*Sigma_s)...
	     -(Data_k-Mu_s).^2./(2*Sigma_s^2);
      temp2 = -log(sqrt(2*pi)*Sigma(i))...
	      -(Data_k-Mu(i)).^2./(2*Sigma(i)^2);
      %if Pzi_c(:,i)'*(log(Pxi.P(:,i))-log(Psi))>log(N)
      if Pzi(:,i)'*(temp2-temp1)>log(N)
	V_trial(i)=1;
      else
	V_trial(i)=0;
	sum_share=sum_share+1;
      end
  end

  %recalculate input likelihood
  BIC_new = log(N)*(1+sum(V_trial)) + M*log(2) + ...
	    calculate_GEM_complete_lik_parallelK(Data_k,Mu,Mu_s,...
						 Sigma,Sigma_s,...
						 V_trial,N,j,M);
  %fprintf('case 3 BIC is %f\n',BIC_new);
  
  if BIC_new<BIC_least
    case_num = 3;
    BIC_least = BIC_new;
    V_least = V_trial;
    if sum_share>=1
      share_exist_least=1;
    else
      share_exist_least=0;
    end
  end
end

