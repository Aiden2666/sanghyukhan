function plot_contour(node, element, contour, contour_in, magnitude , fig, name, location, greek)

%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu

% Find sizes of input data
num_nodes    = length(node);
num_elements = length(element);
num_data     = length(magnitude);

% Now, magnitude can be either nodal or element values 
% If the magnitude data size is equal to the node number, then
% magnitude refers to nodal values. If magnitude data size is equal to
% element number, magnitude refers to element values.
if num_elements==num_data                                   % elemental values
    mag       = zeros(num_nodes,1);
    mag_count = zeros(num_nodes,1);
    for t=1:num_elements                                    % convert elemental values to nodal  
        mag(element(t,:))=mag(element(t,:))+magnitude(t);   % values by averaging   
        mag_count(element(t,:))=mag_count(element(t,:))+1;  % element values over each node
    end
    mag=mag./mag_count;
magnitude = mag;
end

[max_value, max_node_index] = max(magnitude);               % find maximum nodal value 
[min_value, min_node_index] = min(magnitude);               % find minimum nodal value 
        
% Setting up the plot window
screen_size = get(0,'ScreenSize');
screen_size(1:2)=[];
horizontal_step =floor(0.89*screen_size(1)/2);              % acount for different 
vertical_step = floor(0.95*screen_size(2)/2);               % screen resolutions
cur_fig=figure(fig) ; clf; hold on
set(fig,'NumberTitle','off')
set(fig,'Name',name)
set(cur_fig,'position',[5+(location(1)-1)*(8+horizontal_step) -30+(location(2)-1)*vertical_step horizontal_step vertical_step],'Color',[ 0.05 0.05 0.05],'MenuBar','none')
set(get(gca,'Title'),'Color','w')
set(gca,'Color',[ 0.05 0.05 0.05])
main_axes_handle = gca;

% Plot fringe
p_mesh=patch('Vertices', node, 'Faces', element,'CData',magnitude, 'FaceColor', 'interp', 'EdgeColor','None');
% Plot Contour Line
contour_handle=line(contour(1:end,1), contour(1:end,2), 'LineStyle','-', 'LineWidth',2, 'Color','w');
% if not(isempty(contour_in))
% contour_handle2=line(contour_in(1:end,1), contour_in(1:end,2), 'LineStyle','-', 'LineWidth',1, 'Color','w');
% end

% Create colorbar. Different colormaps available


if min(magnitude)<0                             % Case where magnitude spans R
    load('positivenegativemap','mycmap')
    set(fig,'Colormap',mycmap) 
    set(main_axes_handle,'CLim',[-max(abs([max_value,min_value])) max(abs([max_value,min_value]))]);
else
    load('positivemap','mycmap')                % Case where magnitude spans R+
    set(fig,'Colormap',mycmap) 
end

bar_handle=colorbar;

% More Plotting functions. 
axis tight; axis image ; 
%title(name)

set(bar_handle,'YColor','w','FontSize',8,'Box','off')

% Find extends of countour
left_side = get(gca,'XLim'); left_side = left_side(1)*1.5;

% Make sure there is no overlaying of legends
if abs(node(max_node_index,2) - node(min_node_index,2))/(max(node(:,2))-min(node(:,2)))<0.1
    offset = 0.15*(max(node(:,2))-min(node(:,2)));
else
    offset = 0;
end

text(left_side, node(max_node_index,2)*1.15,[' \it\' greek '_{max}: ' num2str(max_value, '%6.2e' ) ], 'Color', [0.85 0.85 0.85], 'FontSize', 12, 'Clipping','off') ;
l_max=line ([left_side node(max_node_index,1)], [node(max_node_index,2) node(max_node_index,2)],  'Marker','.', 'MarkerSize',15,'LineStyle','-', 'LineWidth', 2, 'Color', [0.4 0.9 0.5], 'Clipping','off', 'EraseMode','none');
 

if min(magnitude)<0 
text(left_side, node(min_node_index,2)*1.15+offset,[' \it\' greek '_{min}: ' num2str(min_value, '%6.2e') ], 'Color', [0.85 0.85 0.85], 'FontSize', 12) ;
l_min=line ([left_side node(min_node_index,1)], [node(min_node_index,2)+offset node(min_node_index,2)],  'Marker','.', 'MarkerSize',15,'LineStyle','-',  'LineWidth', 2, 'Color', [0.4 0.9 0.5], 'Clipping','off', 'EraseMode','none');
end

view(0,-90)
zoom(0.8)
