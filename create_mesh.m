function [node, element] = create_mesh(node_ext, node_int, nodes_r, canal_flag, center_A)
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu


nodes_th = length(node_ext);        %number of circumferential nodes

a=0;

if canal_flag==0
    node_int=node_ext;
    radius_int = 0.40 * min(node_ext(:,2));
    node_int(:,2)=radius_int;
    nodes_side = nodes_th/4+1;
end


% Temporary shift angles to make things easier (will be reversed later)
node_ext(:,1)=node_ext(:,1)+pi/2;
node_int(:,1)=node_int(:,1)+pi/2;

%% Seed donut mesh nodes. 
for r=nodes_r+1:-1:1       %radial seed
    for t=1:nodes_th    %circumferential seed
        a=a+1; 
        %determine angle difference between the exterior and the interior node 
        angle = node_ext(t,1) - node_int(t,1);
        if angle>pi
            angle = -(2*pi - angle);
        end
        if angle<-pi
            angle =  (2*pi + angle);
        end
        node(a,1)=node_int(t,1) + angle * (r-1)/nodes_r;
        node(a,2)=node_int(t,2) + (node_ext(t,2) - node_int(t,2)) * (r-1)/nodes_r;      
    end
end

node(:,1)=node(:,1)-pi/2;

num_nodes=a;

%% In case of no intramedulary canal, seed nodes for a butterfly mesh
if canal_flag==false
    %create a rectangular quad mesh
    x=linspace(-radius_int*0.65,radius_int*0.65,nodes_side)';
    [xi,yi] = meshgrid(x,x);
    rec_node=[reshape(xi,[],1) reshape(yi,[],1)];
    rec_node=[rec_node zeros(size(rec_node,1),1)]; % add 3d dimension (will be removed later)
    % find perimeter nodes (counterclockwise)
    perim = 1:nodes_side-1; perim = [perim nodes_side:nodes_side:(nodes_side)^2-nodes_side]; perim = [perim nodes_side^2:-1:nodes_side^2-nodes_side+2]; perim = [perim nodes_side^2-nodes_side+1:-nodes_side:1+nodes_side]';
    perim=circshift (perim, round(-.5*nodes_side))';
    perim=fliplr(perim);
   
    % isolate nodes from the internal loop of the donut mesh
    circle_perim = [nodes_th*(nodes_r)+1:nodes_th*(nodes_r+1)]';

    % map (deform) rectangular mesh points to match circle
    [node_cart(:,1), node_cart(:,2)] = pol2cart(node(:,1), node(:,2));
    [xn yn zn]=scodef(rec_node(:,1),rec_node(:,2),rec_node(:,3),[node_cart(circle_perim,:) zeros(nodes_th,1)],[node_cart(circle_perim,:) zeros(nodes_th,1)]-rec_node(perim,:),4,radius_int*2.5,5);
%     xn =  rec_node(:,1);
%     yn =  rec_node(:,2);
    [th2,r2]=cart2pol(xn,yn);
   % th2(perim)=[]; r2(perim)=[];
    
    r2=r2*0.975;
    node(circle_perim,:) = [th2(perim) r2(perim)];
    node=[node ;th2, r2];
    
end



%% Define element connectivity for donut mesh
a=0;

for r=1:nodes_r
    for t=1:nodes_th-1
        a=a+1;
        element(a,4)=(r-1)* nodes_th + t ;
        element(a,3)=(r-1)* nodes_th + t + 1 ;
        element(a,2)= r   * nodes_th + t + 1;
        element(a,1)= r   * nodes_th + t ;
    end
        %last element for every loop has to link back to the loop's number one node
        a=a+1; 
        element(a,4)=(r-1)* nodes_th + t + 1;
        element(a,3)=(r-1)* nodes_th +     1;
        element(a,2)= r   * nodes_th +     1;
        element(a,1)= r   * nodes_th + t + 1;
end

%% create elements for the rectangular block
if canal_flag==false
    for x=1:nodes_side-1
        for y=1:nodes_side-1
            a=a+1;

            f_node_1 = num_nodes+y+nodes_side*(x-1);

            f_node_2 = f_node_1+nodes_side;

            f_node_3 = f_node_1+nodes_side+1;

            f_node_4 = f_node_1+1;

            element(a,1:4)=[f_node_4 f_node_3 f_node_2 f_node_1];
        end
    end
    
    
  % equivalence nodes   
    perim=perim+num_nodes;
  for t=1:a
      for g=1:4
          ii=find(perim==element(t,g));
        if isempty(ii)==0
            element(t,g) = circle_perim(ii);
        end
      end
  end
  
  % get rid of empty nodes and renumber connectivity
for t=length(node):-1:1
    ii=find(perim==t); 
    if isempty(ii)==0
     ix=find(element>t);
     element_vector = reshape(element,[],1);
     element_vector(ix)=element_vector(ix)-1;
     element = reshape(element_vector,[],4);
     end
end
node(perim,:)=[];

end        
