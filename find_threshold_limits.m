function limit = find_threshold_limits(A, canal_flag, scanning_steps);

temp_center = [size(A)/2];
img_size   = size(A);
img_area   = img_size(1)*img_size(2);
min_value = 0;
max_value = 0.95*max(max(A));

scanning_step  = (max_value-min_value)/scanning_steps;

threshold_level=0; limit=zeros(2,1);

h = waitbar(0,'Determining allowable thresholding limits')

for t=1:scanning_steps
    waitbar(t/scanning_steps,h)
    threshold_level=threshold_level+scanning_step;
   
    [img, bin_image] = threshold_image(A, threshold_level);
    bin_image = imfill(bin_image,'holes');
    bin_image = bwselect(bin_image,temp_center(1), temp_center(2));
    boundary = bwperim (bin_image);
    blob = imfill(boundary,'holes');

    if and(sum(sum(blob))>10,sum(sum(blob))<0.75*img_area)
        if limit(1)==0
            limit(1)=threshold_level;
        end
        if limit(1)>0
            limit(2)=threshold_level;
        end
    end
    
end
 close(h)
        