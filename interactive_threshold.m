function [img, bin_image, ext_blob, int_blob] = interactive_threshold (img, img_original, initial_center, segmentation_threshold, canal_flag)
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu


    %case where colorbar is selected: throshold nad isolate blob
    [img, bin_image] = threshold_image(img_original, segmentation_threshold);
  
    %Create blob
    bin_image_original = im2bw(img,0); 
    %Fill holes in blobs (this is pretty slow)
    bin_image = imfill(bin_image_original,'holes');
    %Select bone blob from all others (also a slow command)
    bin_image = bwselect(bin_image,initial_center(1), initial_center(2));
    boundary = bwperim (bin_image);
    ext_blob = imfill(boundary,'holes');
    
    bin_image= ext_blob;
    
    if canal_flag==1
        %Invert image to detect intramedulary cavity and throw away peripheral pixels
        inv_bin_image = abs(1-bin_image_original).*bin_image;
        inv_bin_image = imfill(inv_bin_image,'holes');
        inv_bin_image = bwselect(inv_bin_image,initial_center(1), initial_center(2));
        inv_blob_area = bwarea(inv_bin_image);
        inv_blob_radius = sqrt(inv_blob_area/pi);
        %now inv_bin_image is the intramedulary canal mask
        int_blob = inv_bin_image;

        %subtract the intramedulary canal from bone blob
        bin_image = abs(ext_blob-int_blob);
    end

    %Threshold image
	img=img_original.*bin_image;

    if sum(sum(img))==0
      disp('threshold selected is too high (or too low), please try lower value')
      bin_image=[];
      center = [];
      ext_blob=[];
      int_blob=[];
     return
    end

    [size_x,size_y]=size(img);


if canal_flag==0   %only pass int_blob if there is an intramedulary canal
    int_blob=[];
end
