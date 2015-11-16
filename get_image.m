function [img] = get_image(file_png,load_dir)
%input a .png image file name and where to load
%output: an image array

if strcmp(load_dir(end),'/')
  img = imread(strcat(load_dir,file_png));
else
  img = imread(strcat(load_dir,'/',file_png));
end


