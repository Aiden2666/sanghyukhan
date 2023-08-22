%% Register the event
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu
function pickandmove(src,eventdata)
h=get(0,'currentfigure');

% Hook the first mousedown event
set(h,'WindowButtonDownFcn',{@button_down,src});
% Call the main button down function 
button_down(src,eventdata,src);

%% Identify the selected point
function button_down(src,eventdata,lobj)
h=get(0,'currentfigure');
button_type=get(gcf,'selectiontype');

switch button_type
    case 'normal'                           % case where left button selection
        curpos=get(gca,'CurrentPoint');
        curpos=curpos(1,1:2)';              % get mouse current position
        
        X=get(lobj,'Xdata');                % get ploted points from line
        Y=get(lobj,'Ydata');                % object of the current graph
        
        % Find the point closest to the curent mouse click position
        [v,ind] = min(sqrt((curpos(1)-X).^2+(curpos(2)-Y).^2));
        single_pt_selection = false;
        
    case 'alt'                              % case where right button selection
        curpos=get(gca,'CurrentPoint');
        curpos=curpos(1,1:2)';              % get mouse current position
        
        X=get(lobj,'Xdata');                % get ploted points from line
        Y=get(lobj,'Ydata');                % object of the current graph
        
        % Find the point closest to the curent mouse click position
        [v,ind] = min(sqrt((curpos(1)-X).^2+(curpos(2)-Y).^2));
        single_pt_selection = true;
end

% Hook the rest of the mouse events 'up' and 'move'
set(h,'WindowButtonUpFcn',{@button_up},...
    'WindowButtonMotionFcn',{@mouse_move,lobj,ind, single_pt_selection});

%% On button up, release all events
function button_up(src,eventdata)
h=get(0,'currentfigure');
%Unhook all events
set(h,'WindowButtonDownFcn','','WindowButtonUpFcn','','WindowButtonMotionFcn','');


%% Update position for points based on the new position of the mouse
function mouse_move(src,eventdata,lobj,index, single_pt_selection)

curpos=get(gca,'CurrentPoint');
curpos=curpos(1,1:2)';                      % get mouse current position

X=get(lobj,'Xdata');                        % get ploted points from line
Y=get(lobj,'Ydata');                        % object of the current graph

num_pts = length(X);

% Calculate the change in position [DX,DY]
DX=X(index)-curpos(1);
DY=Y(index)-curpos(2);

% Calculate curvature of each point based on an examination region (ROE)
% Examination region is here defined in terms of number of points
% to the left and to the right. For a ROE of 2, for example, there will
% be 5 points to evaluate the curvature: the selected point and 2 from the
% left and 2 from the right. This can also employ an eucledian distance 
% based criterion 

ROE=3;

index_start=index-ROE;                      %region of Examination (region_e)
index_end  =index+ROE;                      %[index_start : index_end]
region_e = [index_start : index_end];

% Here are 2 tests to account for the closed loop numbering:

% 1. Case where the selected point is in the beginning of the numbering
% while the tail of the region is at the ending of the numbering
if index_start<1
    region_e = [index_start+num_pts:num_pts, 1:index_end];
end
% 2. Case where the selected point is in the ending of the numbering
% while the head of the region is at the beginning of the numbering
if index_end>num_pts
    region_e = [index_start:num_pts, 1:(index_end-num_pts)];
end

% Calculate the local curvature by fitting a circle and finding its radius
[R, center_local,error] = Calculate_Curvature (X(region_e)',Y(region_e)');

% Calculate Range of influence (ROI)
% The range is porportional to the square of the local curvature R, but it
% has to be relative to the dimensions of the object. So:
% First, we find the "equivalent radius" of the entire object. This is 
% done by fitting a circle to ALL the points of the contour.
[R_object,center_object,dummy] = Calculate_Curvature (X',Y');

if R>R_object                               % in case the local curvature
    R=R_object;                             % is larger than the "equivalent"
end                                         % curvature, its value is limited

%Override and apply penalty on the points that have a convex curvature
centers_distance = sqrt(sum((center_object-center_local).^2)); %eucledian distance
if centers_distance>R_object
    R=0.4*R;                                % here we apply the penalty to the convex point
end

R=R*((1-error)^16);

% Now R can be compared to R_object in order to find how much influence the
% particular selected point should have on the neighbors. 


% The ROI is a vector that contains the indices of the points that are to
% be influenced by the particular mouse move.
% It is defined in a way so that if the local curvature is maxed, half of
% the points of the object are going to follow the transformation (num_pts/2)
ROI = ceil ( (R^2/R_object^2) * (num_pts/2) );
ROI=  ceil (ROI/2);                         % get the one-sided range

t=-ROI:ROI;
DDX=cos (t*(pi/2)/ROI) * DX;                % define a cosine deformation 
DDY=cos (t*(pi/2)/ROI) * DY;                % function with a peak at the selected point 

index_start=index-ROI;                      % region of Influence (region_i)
index_end  =index+ROI;                      % [index_start : index_end]
region_i = [index_start : index_end];

% Again, Here are 2 tests to account for the closed loop numbering:

% 1. Case where the selected point is in the beginning of the numbering
% while the tail of the region is at the ending of the numbering
if index_start<1
    region_i = [index_start+num_pts:num_pts, 1:index_end];
end
% 2. Case where the selected point is in the ending of the numbering
% while the head of the region is at the beginning of the numbering
if index_end>num_pts
    region_i = [index_start:num_pts, 1:(index_end-num_pts)];
end


% Update positions for points within the ROI     
if single_pt_selection==false
    X(region_i)=X(region_i)-DDX;
    Y(region_i)=Y(region_i)-DDY;
else
    X(index)=X(index)-DX;
    Y(index)=Y(index)-DY;
end
% Set ploted points from line object of the current graph
set(lobj,'Xdata',X);
set(lobj,'Ydata',Y);

