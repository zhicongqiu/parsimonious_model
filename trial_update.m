function [Mu_s_least Sigma_s_least updated] = ...
	 trial_update(Data_k,Mu,Mu_s,Sigma,Sigma_s,share_exist,V_trial,N,M)

global Pzi;
%trial update shared parameters
BIC_least = Inf;
updated = 0;

%fprintf('begin...\n');
if share_exist==0
  %find the components with the closest distance

  temp_dis_least = inf;
  for j=1:M
    for k=j+1:M
	temp_dis = abs(Mu(j)-Mu(k))
	if temp_dis<temp_dis_least
	   temp_dis_least = temp_dis;
	   index_1 = j;
	   index_2 = k;
	end
    end
  end
  V_trial(index_1)=0;V_trial(index_2)=0;share_exist_trial = 1;
  updated = 1;
  %Update shared parameters, using the previous Pzi???
  Up = 0; Down = 0;
  %Update the centers using the subset of sample belong to common
  %Update shared parameters
  %Pzi_s = Pzi*(1-V_trial');%cancelled for uncessary memory usgage
  Up = sum(Data_k.*(Pzi*(1-V_trial')));
  Down = sum((Pzi*(1-V_trial')));
  if Down<realmin
    Down = realmin;
  end
  Mu_s_least = Up/Down;
  Up = sum((Pzi*(1-V_trial')).*(Data_k-Mu_s_least).^2);
  Sigma_s_least = sqrt(Up/Down);
  Sigma_s_least(Sigma_s_least<1e-5)=1e-5;
	    
  %{  
  for j=1:M
    %trial-set mu_s and sigma_s to any of the component    
    Mu_s=Mu(j);
    Sigma_s=Sigma(j);
    [V_trial share_exist_trial BIC_new] = ...
    configureSwitch(Data_k,Mu,Mu_s,Sigma,Sigma_s,...
		    share_exist,V_trial,i,N,M,1);    
    %fprintf('end configure\n');
    if BIC_new<BIC_least
      %remember this setting
      BIC_least = BIC_new;
      V_least = V_trial;
      share_exist_least = share_exist_trial;
      if share_exist_least==1
	updated = 1;
	%Update shared parameters, using the previous Pzi???
	Up = 0; Down = 0;
	%Update the centers using the subset of sample belong to common
	%Update shared parameters
	Pzi_s = Pzi*(1-V_least');
	Up = sum(Data_k.*Pzi_s);
	Down = sum(Pzi_s);
	if Down<realmin
	  Down = realmin;
	end
	Mu_s_least = Up/Down;
	Up = sum(Pzi_s.*(Data_k-Mu_s_least).^2);
	Sigma_s_least = sqrt(Up/Down);
        Sigma_s_least(Sigma_s_least<1e-5)=1e-5;
      end
    end
  end
  if updated == 0
    %randomly set to the first component
    Mu_s_least = Mu(1);
    Sigma_s_least = Sigma(1);
  end
%}
%use the updated values
else
  updated = 1;
  Mu_s_least = Mu_s;
  Sigma_s_least = Sigma_s;
end
