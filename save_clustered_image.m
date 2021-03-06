 function [] = save_clustered_image(GMM,files)
%save clustered images based on posterior assignments
%%%%%%%%%%%%%%
%input: [GMM,files]
%GMM: GMM model with specified SWITCH, Priors, Mu, Mu_s, Sigma,
%Sigma_s
%files: files array, each contains .mat file names and 
%r, which specified which images are used in clustering 
%%%%%%%%%%%%%%
%output: copy images to cluster directories
%%%%%%%%%%%%%%
global Data Label;
[N K] = size(Data);
dir_name = ['72458';'72460';'72463';'72467';'72468'...
	    ;'72469';'72470';'72471';'72472';'72473'...
	    ;'72474';'72475';'72476';'1';'2';'72412'];

%make directory and save directory names, one for each component
write_dir = [];
for i=1:GMM.num_comp
    temp_dir = strcat('CLUSTER_',int2str(i));
    mkdir(temp);
    write_dir = [write_dir;temp];
end
[Pzi] = LikelihoodEst(GMM.Priors,GMM.Mu,GMM.Mu_s,...
		      GMM.Sigma,GMM.Sigma_s,GMM.SWITCH,...
		      N,K,length(GMM.Priors));
%1xN matrix specifying which component it belongs to
[Pzi_max Pzi_hard] = max(Pzi');

clear Pzi; %not used any more?

%loop over all file names
count = 0;
for i=1:size(dir_name,1)
    sample_length = length(files(i).r);
    for j=1:sample_length	
	temp_name = files(i).files(files(i).r(j)).name;
	%get the png file name
	temp_png = strcat(temp_name(15:end-4),'.png');
	fprintf(strcat('reading dir',dir_name(i,:),...
		       ' png file ',temp_png,'.\n'));
	temp_img = get_image(temp_png,dir_name(i,:));
	%save the image into the cluster dir
	imwrite(temp_img,...
		strcat(write_dir(Pzi_hard(count+j),:),...
		       '/',int2str(i),'_',temp_png));
    end
    count = count+sample_length;
end
