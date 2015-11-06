function [Priors, Mu, Sigma] = GEM_init_kmeans(nbStates)

% Modified to evaluate GEM-NB
% Sigma is the variance matrix, each column coresponds to a given
% compoenents, dim is K*M
% This function initializes the parameters of a Gaussian Mixture Model 
% (GMM) by using k-means clustering algorithm.
%
% Author:	Sylvain Calinon, 2009
%			http://programming-by-demonstration.org
%
% Inputs -----------------------------------------------------------------
%   o Data:     D x N array representing N datapoints of D dimensions.
%   o nbStates: Number K of GMM components.
% Outputs ----------------------------------------------------------------
%   o Priors:   1 x K array representing the prior probabilities of the
%               K GMM components.
%   o Mu:       D x K array representing the centers of the K GMM components.
%   o Sigma:    D x D x K array representing the covariance matrices of the 
%               K GMM components.
% Comments ---------------------------------------------------------------
%   o This function uses the 'kmeans' function from the MATLAB Statistics 
%     toolbox. If you are using a version of the 'netlab' toolbox that also
%     uses a function named 'kmeans', please rename the netlab function to
%     'kmeans_netlab.m' to avoid conflicts. 
global Data;
[nbData nbVar] = size(Data);

%Use of the 'kmeans' function from the MATLAB Statistics toolbox
[Data_id, Centers] = kmeans(Data, nbStates , 'emptyaction', 'singleton'); % !!!!!!!!! modified by Fatih
Mu = Centers';
for i=1:nbStates
  idtmp = find(Data_id==i);
  Priors(i) = length(idtmp);
  Sigma(:,i) = sqrt(var(Data(idtmp,:))');
  %Sigma(:,:,i) = cov([Data(:,idtmp) Data(:,idtmp)]');
  %Add a tiny variance to avoid numerical instability
  %Sigma(:,:,i) = Sigma(:,:,i) + 1E-5.*diag(ones(nbVar,1));
end
Sigma(Sigma<1e-5)=1e-5;
Priors = Priors ./ sum(Priors);


