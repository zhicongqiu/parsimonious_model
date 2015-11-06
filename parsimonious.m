function [MODEL,METRICS] = parsimonious(Data_whole,M_max)

%parsimonious model based on Graham, 2006
%input:
%Data_whole: N*(K+1) metrix, each row is a sample and each column is a
%feature, last column is label
%M_max: maximum component considered
%output:
%MODEL: shared models obtained
%METRICS: measure of BIC, error rate and k_means error at M_max

%octave specific command, load statistics and parallel toolbox
%pkg load statistics
%pkg load ndpar
%generate_random_num;

global Pzi Data Label ;
Data = Data_whole(:,1:end-1);
[N K] = size(Data);
Label = Data_whole(:,end);
Log_incomplete = zeros(M_max-1,1);
Log_complete = zeros(M_max-1,1);
BIC_least = Inf;
%parallel setting with the actual number of processors
nproc_temp = nproc;
Chunks = K;
%Criterion to stop the EM iterative update
learning_threshold = 1e-4;

error = zeros(1,M_max);
%time stamp
tic();
for i=M_max:-1:2

  if i==M_max %initialize parameters
    [Priors Mu Sigma] = GEM_init_kmeans(M_max);
    Mu_s = zeros(1,K);
    Sigma_s = 1e-6*ones(1,K);
    %begin by setting share_exist to zeros, and all Vs to ones
    share_exist = zeros(1,K);
    V = ones(K,M_max);
    [error(1)] = ...
    calculate_error(Priors,Mu,Mu_s,...
		    Sigma,Sigma_s,share_exist,...
		    V,N,K,M_max);
    fprintf('initial k-mean error:%f\n',error(1));
  else
    %specify the deletion of a single component
    [Priors Mu Sigma V share_exist]  = ...
    delete_strategy(Priors,Mu,Mu_s,Sigma,Sigma_s,V,share_exist);
   
  end
  %clear functions;
  fprintf('begine GEM with number of component = %d\n',i)
  GEM;
  Log_incomplete(M_max-i+1) = BIC_incomplete;
  Log_complete(M_max-i+1) = BIC_old;
  %calculate error in this iteration
  [error(M_max-i+2)] = ...
  calculate_error(Priors,Mu,Mu_s,...
		  Sigma,Sigma_s,share_exist,...
		  V,N,K,i);
  fprintf('error at comp=%d is %f...\n',i,error(M_max-i+2));
  if wrong==1
    fprintf('sth is wrong, stop immediately!\n')
    break;
  end
  if BIC_incomplete<=BIC_least
    %favor smaller number of components when tied
    fprintf('BIC new min at num_comp=%d, update parameters\n',i)
    BIC_least = BIC_incomplete;
    
    share_exist_least = share_exist;
    V_least = V;
    Priors_least = Priors;
    Mu_least = Mu;
    Mu_s_least = Mu_s;
    Sigma_least = Sigma;
    Sigma_s_least = Sigma_s;
    num_comp_least = i;
  end
end
MODEL.SWITCH = V_least;
MODEL.Priors = Priors_least;
MODEL.Mu = Mu_least;
MODEL.Mu_s = Mu_s_least;
MODEL.Sigma = Sigma_least;
MODEL.Sigma_s = Sigma_s_least;
METRICS.Log_incomplete = Log_incomplete;
METRICS.ERROR = [error(2:end)];
METRICS.km_error = error(1);
toc();
