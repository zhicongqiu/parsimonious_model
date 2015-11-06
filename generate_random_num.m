%generate random Gaussian distributed vectors
%{
mu1=[];mu2=[];mu3=[];mu4=[];
sigma1=[];sigma3=[];
for i=1:1000
  mu1 = [mu1 0 -4 0 0 0];
  mu2 = [mu2 0 0 0 0 0];
  mu3 = [mu3 -3 4 0 0 0];
  mu4 = [mu4 3 4 0 0 0];
  sigma1 = [sigma1 sqrt(3) 1 1 1 1];
  sigma3 = [sigma3 1 1 1 1 1];
end
%}

mu1 = [0 -4 0 0 0];
mu2 = [0 0 0 0 0];
mu3 = [-3 4 0 0 0];
mu4 = [3 4 0 0 0];
sigma1 = [sqrt(3) 1 1 1 1];
sigma3 = [1 1 1 1 1];
sigma2 = sigma1;
sigma4 = sigma3;

Data_whole = zeros(2000,6);
for i=1:500
  Data_whole(i,1:end-1) = normrnd(mu1',sigma1')';
end
Data_whole(1:500,end) = 1;

for i=501:1000
  Data_whole(i,1:end-1) = normrnd(mu2',sigma2')';
end
Data_whole(501:1000,end) = 2;

for i=1001:1500
  Data_whole(i,1:end-1) = normrnd(mu3',sigma3')';
end
Data_whole(1001:1500,end) = 3;

for i=1501:2000
  Data_whole(i,1:end-1) = normrnd(mu4',sigma4')';
end
Data_whole(1501:2000,end) = 4;

