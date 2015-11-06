%function [Priors Mu Mu_s Sigma Sigma_s ...
%	 V share_exist BIC_old wrong Log_incomplete Log_complete] ...
%	 = GEM(Data, Priors, Mu, Mu_s, Sigma, Sigma_s, V,share_exist)


%indicate wrongnes
wrong = 0;
%% Criterion to stop the EM iterative update
M = size(Sigma,2);

nbStep = 0;
%% E-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parallelizable...
[Pzi] = ...
LikelihoodEst(Priors,Mu,Mu_s,Sigma,Sigma_s,V,N,K,M);

%Pzi_c = Pzi;
%BIC_old
[CC]= calculate_GEM_CC(V,N,K,M);
Idx = 1:K;
[LIKK]= ...
ndpar_arrayfun(nproc_temp,@calculate_GEM_complete_lik_parallelK,...
	       Data,Mu,Mu_s,Sigma,Sigma_s,V,N,Idx,M,...
	       'Vectorized', true, 'ChunksPerProc',Chunks,...
	       'VerboseLevel', 1,'IdxDimensions', ...
	       [2 1 2 1 2 1 0 2 0],'CatDimensions', [2]);
BIC_old = CC+sum(LIKK)-sum(Pzi*log(Priors)');
LIKK_old = BIC_old-CC;

fprintf('evaluate BIC_old=%f and LIKK=%f.\n',BIC_old,LIKK_old);

while nbStep<1e4
  %checkpoint for nan number
  if isnan(BIC_old)==true %|| isnan(incomplete_old)==true
    fprintf('complete or incomplete is nan, stop immediately.\n');
    wrong=1;
    break;
  end
  nbStep = nbStep+1;

  fprintf('%dth iteration\n',nbStep)

  %Compute cumulated posterior probability
  %Maximization evaluation
  [Priors_new Mu_new Mu_s_new Sigma_new Sigma_s_new] = ...
  GEM_Max(V,share_exist,N,K,M);

  %trial set Mu_s and Sigma_s if it is not set
  %parallelizable
  updated = zeros(1,K);
  fprintf('begin trial updates...\n');

  [Mu_s_new Sigma_s_new updated] = ...
  ndpar_arrayfun(nproc_temp,@trial_update,Data,...
		 Mu,Mu_s_new,Sigma,Sigma_s_new,share_exist,V,N,M,...
		 'Vectorized', true, 'ChunksPerProc', Chunks,...
		 'VerboseLevel', 0,'IdxDimensions', ...
		 [2 1 2 1 2 2 1 0 0],'CatDimensions', [2 2 2]);
  

  fprintf('configuring switches\n')

  %configure switches, and return complete BIC_new, using the updated
  %parameters
  %parallelizable
  %BIC_newK not used. If used, modified accordingly
  [V share_exist BIC_newK] = ...
  ndpar_arrayfun(nproc_temp,@configureSwitch,Data,Mu_new,Mu_s_new,...
		 Sigma_new,Sigma_s_new,...
		 share_exist,V,Idx,N,M,updated, ...
		 'Vectorized', true, 'ChunksPerProc', Chunks,...
		 'VerboseLevel', 0,'IdxDimensions', ...
		 [2 1 2 1 2 2 1 2 0 0 2],'CatDimensions', [1 2 2]);
 
  [CC]= calculate_GEM_CC(V,N,K,M);
  LIKK=zeros(1,K);
  [LIKK]= ...
  ndpar_arrayfun(nproc_temp,@calculate_GEM_complete_lik_parallelK,...
		 Data,Mu_new,Mu_s_new,Sigma_new,Sigma_s_new,V,N,Idx,M,...
		 'Vectorized', true, 'ChunksPerProc',Chunks,...
		 'VerboseLevel', 0,'IdxDimensions', ...
		 [2 1 2 1 2 1 0 2 0],'CatDimensions', [2]);
  BIC_new = CC+sum(LIKK)-sum(Pzi*log(Priors_new)');
  LIKK_new = BIC_new-CC;

  temp = abs(BIC_old-BIC_new);
  if BIC_new>BIC_old && temp>1e-4
    wrong=1; 
    fprintf('wrong after switches:\n')
    fprintf('update complete log likelihood from %f to %f.\n',...
    	    BIC_old,BIC_new)
    fprintf('diff is %f.\n',temp*1e6)
    break;  
  else
    %stop criterion: relative difference of BIC
    if abs(BIC_new-BIC_old)<learning_threshold
      fprintf('current in-BIC: %f, from %f, done!\n',...
	      BIC_new,BIC_old)
      PostEnt = 0;
      for ii=1:N
	  for jj=1:M
	    if Pzi(ii,jj)~=0
	      PostEnt = PostEnt+Pzi(ii,jj)*log(Pzi(ii,jj));
	    end
	  end
      end      
      BIC_incomplete = BIC_old-PostEnt;
      break;
    end
    fprintf('current in-BIC: %f, log=%f,from %f,%f. update parameters\n',...
	    BIC_new,LIKK_new,BIC_old,LIKK_old)
    fprintf('log_new-log_old=%f\n',LIKK_new-LIKK_old)
    fprintf('number of on-switches: %d\n',sum(sum(V)))
    Priors = Priors_new;
    Mu = Mu_new;
    Mu_s = Mu_s_new;
    Sigma = Sigma_new;
    Sigma_s = Sigma_s_new;

    %reevaluate Pzi and Px using new switches, E-step
    [Pzi] = LikelihoodEst(Priors,Mu,Mu_s,Sigma,Sigma_s,V,N,K,M);
    [CC]= calculate_GEM_CC(V,N,K,M);
    [LIKK]= ...
    ndpar_arrayfun(nproc_temp,@calculate_GEM_complete_lik_parallelK,...
		   Data,Mu,Mu_s,Sigma,Sigma_s,V,N,Idx,M,...
		   'Vectorized', true, 'ChunksPerProc',Chunks,...
		   'VerboseLevel', 1,'IdxDimensions', ...
		   [2 1 2 1 2 1 0 2 0],'CatDimensions', [2]);
    BIC_old = CC+sum(LIKK)-sum(Pzi*log(Priors)');
    LIKK_old = BIC_old-CC;
  end
end
