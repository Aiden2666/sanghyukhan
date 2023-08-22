% Program Main.m
% This calls all the routines needed to calculate stresses for a bone
% x-section

%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu

%PLEASE INSPECT THE FOLLOWING INPUT SECTION
close('all');clear all;clc


%% 0. INPUT
         
% Flexural Loads
Mx             = 0;             %in N.m
Mx=Mx*1000;
My             = 0;             %in N.mm
% OR (3pt bending)
Qx             = 707;           %in N
Qy             = 707;           %in N
z_pos          = 150;           %in mm
% Axial Loads
P              = 500;           %in N
T              = 1;
T=T*1000;
% Torsion displacement


    

%DONT forget to edit image2density.m to match your image to density calibration relation. See manual for examples
  

%% Model parameters

nodes_r  = 20;                  % number of nodes in the radial direction
nodes_th = 96;                  % number of nodes along the circumference



zone_width=11;                   %select search zone width for boundary search refinement
                                 
                                               
%%
nodes_th = 4*round(nodes_th /4); %Adjust number of nodes along circumference 
                                 %so that it is divided by 4 (for meshing
                                 %purposes)

%% 1. Import the image. 

[A, resolution, slice_z, filename]=read_section_image;
if A==0;return;end
A=double(A); A_original=A;

plot_image(A, 1, ['Subject: ' filename '  Slice Position: ' num2str(slice_z) '  - Original Scan [original units]'], [1,2]);

%Interactively select bone by cropping the original image.
max_A_value = max(max(A));
A = max_A_value*imcrop(A/max_A_value);
A_original = A;
plot_image(A, 1, ['Subject: ' filename '  Slice Position: ' num2str(slice_z) '  - Original Scan [original units]'], [1,2]);
bar_handle=colorbar;


%% 2. perform image processing to "clean" image 
% A is the original image in image units
% A_bin is the logical map used to mask the image. 
plot_image(A_original, 1, ['Subject: ' filename '  Slice Position: ' num2str(slice_z) '  - Original Scan [original units]'], [1,2]);
bar_handle=colorbar; 

canal_flag = determine_mesh_type;

temp_center = [size(A)/2];

limit = find_threshold_limits(A_original, canal_flag, 50);

segmentation_threshold=mean(limit);

[A, A_bin]  =   threshold_image(A_original, segmentation_threshold); %threshold image. segmentation_threshold is the level in HUplot_image(A, 2, 'Cropped and Thresholded Image - [image units]', [2,2]); hold on

set(bar_handle,'YColor','w','FontSize',8,'Box','off','YLim',[limit(1), limit(2)]);

%% 3. Threshold image interactively and Calculate center of xsection 
%Repeat thresholding and cropping.(left clicking) 
%By selecting a point on the colorbar, it automatically uses it as 
%a new threshold level and reiterates.
%By selecting a point on the image, it recalculates the center, and
%re-crops the image. Then the countour is found.
%Right clicking to confirm, exit this loop and proceed
but=1;

while (but==1)

figure(1); [selected_point(1),selected_point(2),but]=ginput(1);
    if but==1


        if selected_point(1)<1; 
            segmentation_threshold=selected_point(2);
        else
            temp_center=selected_point;
        end

        [A, A_bin, ext_blob, int_blob] = interactive_threshold (A, A_original, temp_center, segmentation_threshold, canal_flag); 
        center_A        = find_center (A);         %find center of image
        
        plot_image(A, 2, 'Cropped and Thresholded Image - [image units]', [2,2]); hold on
        figure(2);plot(center_A(1),center_A(2),'r+');

    % 3A. Find the internal & external nodes
        [node_ext_polar, node_ext_cart, node_int_polar, node_int_cart] = detect_boundary (center_A, nodes_th, canal_flag, ext_blob, int_blob);

    % Plot countour(s)
        ex_contour_handle=line(node_ext_cart(1:end-1,1), node_ext_cart(1:end-1,2), 'Marker','.', 'MarkerSize',8,'LineStyle','-','ButtonDownFcn','pickandmove(ex_contour_handle,[])');
        if canal_flag==1     %plot inner contour iff intramedulary canal is present
            in_contour_handle=line(node_int_cart(1:end-1,1), node_int_cart(1:end-1,2), 'Marker','.', 'MarkerSize',8,'LineStyle','-','ButtonDownFcn','pickandmove(in_contour_handle,[])');
        end
        drawnow

% 3B. Refine the boundary detection by creating a trust region around the
%previously detected boundary (zone), and then apply image partial voluming
%correction algorithms. Then re-detect boundary.
        [x_dimension, y_dimension ] = find_bounding_box(node_ext_cart, resolution);
        
       
 %define the width of the search zone based on whether there is thick enough
 %cortical bone on the periosteal surface. For distal/proximal bone,
 %cortical bone is thin, for midshaft cortical bone is thick therefore
 %search zone width can be wider
        if canal_flag==false
            zone_width = 7;
        end

        set(ex_contour_handle,'XData',node_ext_cart(1:end-1,1),'YData',node_ext_cart(1:end-1,2), 'Color','r','LineStyle','-','ButtonDownFcn','pickandmove(ex_contour_handle,[])');
        drawnow;
        
        [node_ext_polar] = refine_boundary(node_ext_polar,  A, center_A, zone_width,ex_contour_handle);
   
        node_ext_cart = polar2cart(node_ext_polar,center_A);
        set(ex_contour_handle,'XData',node_ext_cart(1:end-1,1),'YData',node_ext_cart(1:end-1,2), 'LineStyle','-','ButtonDownFcn','pickandmove(ex_contour_handle,[])');
        drawnow
    
    end % end if mouse click is left button : right click selected exit

end % end while loop



%Read back the vertex positions and Dismiss the events for the contour plots
node_ext_cart = [ get(ex_contour_handle,'XData'); get(ex_contour_handle,'YData') ]' ;
node_ext_cart = [node_ext_cart ; node_ext_cart(1,:)];

if canal_flag==1
    set(in_contour_handle,'ButtonDownFcn','')
    node_int_cart = [ get(in_contour_handle,'XData'); get(in_contour_handle,'YData') ]' ;
    node_int_cart = [node_int_cart ; node_int_cart(1,:)];
end


%% 4. Convert to cartesian coords and Mesh cross section.

%Convert seed node positions from cartesian to polar coordinates
  [node_ext_polar] = cart2polar (node_ext_cart, center_A);
%  node_ext_polar(:,1)=pi/2-node_ext_polar(:,1);       %the pi/2-  is there cause matlab is treating x,y in an inverse way

  if canal_flag==1
  [node_int_polar] = cart2polar (node_int_cart, center_A);
%   node_int_polar(:,1)=pi/2-node_int_polar(:,1);      %the pi/2-  is there cause matlab is treating x,y in an inverse way
  end
% 


%% Create the mesh based on the external and internal nodes
[node_polar, element] = create_mesh(node_ext_polar(1:end-1,:), node_int_polar(1:end-1,:), nodes_r, canal_flag, center_A);

% Convert to cartesian coordinates  
node = polar2cart (node_polar, center_A);


% Calculate bounding box for bone. This can be compared to direct caliper
% measurements or x-ray images in order to determine segmentation accuracy
[x_dimension, y_dimension ] = find_bounding_box(node_ext_cart, resolution)    
OD_mean = mean(node_ext_polar(:,2));


%% 5. Convert to material properties and find centers. 
% A_d is the density image
% A_YM is the Young's Modulus image
% A_v is the Poisson's ratio image

center_A_bin          = find_center (A_bin);
[A_d]                 = image2density(A,A_bin);          %convert HU to apparent density [gr/cm3]
[A_YM, A_G, A_v]      = density2elastic_props(A_d);       %convert apparent density [g/cm3] to YM [MPa] G[MPa]

% Calculate center for each image
center_d        = find_center (A_d);
center_YM       = find_center (A_YM);
center_G        = find_center (A_G);
   
 total_BM = sum(sum(A_d));
BMD = total_BM/(x_dimension*2.3); % calculate BMD by taking the sum of mineral and divide by the total area of the slice


%% 6. Convert nodes to real units: mm and assign material properties to each element. 
% Then, determine the cross sectional properties

node_real(:,1)=node(:,1)*resolution(1); %[node_real] are the coordinates in real units (mm)
node_real(:,2)=node(:,2)*resolution(2); %[node]      are the coordinates in image units (pixels)
node_ext_cart(:,1)=node_ext_cart(:,1)*resolution(1);
node_ext_cart(:,2)=node_ext_cart(:,2)*resolution(2);
 if canal_flag==1
	node_int_cart(:,1)=node_int_cart(:,1)*resolution(1);
	node_int_cart(:,2)=node_int_cart(:,2)*resolution(2);
 end
% This will give the YM, G in Pa and centroid,area (in square mm) for each element

[YM, v, G, rho, area, centroid, center_mesh, center_mesh_E] = assign_material (node, node_real, element, A_YM, A_G, A_v, A_d, resolution, nodes_th);

% %  For Benchmarking only
% %  YM=2*(1+0.33)*ones(length(YM),1);
% %     v=0.33*ones(length(v),1);
% %     G=1*ones(length(G),1);

% Plotting
figure(2); p_mesh=patch('Vertices', node, 'Faces', element,'FaceAlpha',0.3,'EdgeColor',[0.3 0.3 0.3]);
 set(p_mesh,'FaceVertexCData',repmat(YM,1,3)./max(YM), 'FaceColor', 'flat','FaceAlpha',1); drawnow;
 colormap gray; bar_handle=colorbar; set(bar_handle,'YColor','w','FontSize',8,'Box','off'); drawnow


%% 6. Calculate the sectional properties.

[I, J, I_principal, V,   I_E, J_E, I_E_principal, V_E,   I_G, J_G, I_G_principal, V_G] = section_properties (YM, G, rho, centroid, area, center_mesh, center_mesh_E);

equivalent_diameter = 2*sqrt(sum(area)/pi);

%Set 0,0 at the middle
 node_real = node_real - repmat(center_mesh_E,length(node_real),1);
 centroid  = centroid  - repmat(center_mesh_E,length(element),  1);
 
                   node_ext_cart=node_ext_cart-repmat(center_mesh_E,length(node_ext_cart),1);
 if canal_flag==1; node_int_cart=node_int_cart-repmat(center_mesh_E,length(node_ext_cart),1); end

%% 7. Solve Shear problem

for bb=1:2
    if bb==1
        theta=1;
    else
        theta=T/K_mod
    end


[shear_stress, resultant_shear_stress, K_mod] = solve_FE_problem (node_real,element, YM, G, v, I_E, Qx, Qy, theta);

torque = theta*K_mod

resultant_shear_stress_e = sqrt(shear_stress(:,1).^2+shear_stress(:,2).^2);
%resultant_shear_stress_e = resultant_shear_stress_e./((rho./1.92).^1.333);
resultant_shear_stress = elemental2nodal (resultant_shear_stress_e, node, element);
plot_contour(node_real , element, node_ext_cart,node_int_cart, resultant_shear_stress, 4, 'Shear Stresses' , [1 1], 'tau')
end
%quiver(centroid(:,1),centroid(:,2),shear_stress(:,1),shear_stress(:,2),0.65,'k'); drawnow

%% 8. Bending Stresses: use Beam theory to calculate bending stresses

[bending_stress_e] = calculate_bending (Mx, My, Qx, Qy, z_pos, P, I_E, node_real, area,centroid, YM,  [0 0]);
%bending_stress = bending_stress./((rho./1.92).^2);

[bending_stress] = elemental2nodal (bending_stress_e, node, element);

plot_contour(node_real, element, node_ext_cart, node_int_cart, bending_stress , 3, 'Out-of-Plane Normal Stresses' , [2 1], 'sigma'); drawnow

sigma_yield = 68 * rho.^2;
tau_yield   = 21.6 * rho.^1.65;

%Equivalent Stress (Tsai - Wu)
equivalent_stress = ( abs(bending_stress_e) ./ sigma_yield).^2 + (resultant_shear_stress_e ./ tau_yield).^2;

plot_contour(node_real, element, node_ext_cart, node_int_cart, equivalent_stress , 5, 'Equivalent Stress (Elliptical Criterion)' , [1 2], 'sigma'); drawnow

yield=find(equivalent_stress>1);
yield_area = sum(area(yield)) / sum(area);


yield_plot = ones(length(equivalent_stress),3);
yield_plot(yield,1) = 0.8; yield_plot(yield,2)=0.2; yield_plot(yield,3)=0.2;

plot_image(A, 6, 'Failed Elements (red)', [2,2]); hold on
        figure(2);plot(center_A(1),center_A(2),'r+');
 p_mesh=patch('Vertices', node, 'Faces', element,'FaceAlpha',0.3,'EdgeColor',[0.3 0.3 0.3],'EdgeColor','None');
  set(p_mesh,'FaceVertexCData',yield_plot, 'FaceColor', 'flat','FaceAlpha',1); drawnow;
  bar_handle=colorbar; set(bar_handle,'YColor','w','FontSize',8,'Box','off'); drawnow
 
 

break




%% 9. Trace cortical and trabecular bone
trabecular = find(YM<8000); 
cortical   = find(YM>8000);
marrow     = find(YM<2);

YM_plot = ones(length(YM),3);

YM_plot(trabecular,1) = 1; YM_plot(trabecular,2)=0.0; YM_plot(trabecular,3)=0.0;
YM_plot(cortical,2) = 1; YM_plot(cortical,1)=0.0; YM_plot(cortical,1)=0.0;
YM_plot(marrow,3) = 0; YM_plot(marrow,1)=0; YM_plot(marrow,2)=0;

plot_image(A, 7, 'Trabecular and Cortical Bone', [2,2]); hold on
        figure(2);plot(center_A(1),center_A(2),'r+');

p_mesh=patch('Vertices', node, 'Faces', element,'FaceAlpha',0.3,'EdgeColor',[0.3 0.3 0.3]);
 set(p_mesh,'FaceVertexCData',YM_plot, 'FaceColor', 'flat','FaceAlpha',1); drawnow;
 bar_handle=colorbar; set(bar_handle,'YColor','w','FontSize',8,'Box','off'); drawnow


 

trabecular_area = sum(YM_plot(:,1).*area);
cortical_area = sum(YM_plot(:,2).*area);
total_area =cortical_area+ trabecular_area
cortical_ratio = cortical_area/total_area
trabecular_ratio = 1-cortical_ratio





%% 10. ABAQUS export
extrusion_depth=100; extrusion_nodes=20; inhomogeneity=1;

abaqus_export2(node_real, element, YM, G, v, nodes_th, extrusion_depth, extrusion_nodes, 'test1.inp', inhomogeneity)