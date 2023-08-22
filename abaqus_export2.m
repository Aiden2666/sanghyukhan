function abaqus_export (node,element,YM, G, v, nodes_th, depth, z_nodes ,filename, inhomogeneity)
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu
fid = fopen(filename,'w');

z_nodes=2*round(z_nodes/2)+1;  %make sure that the number of node is always an odd number,so that element number is even 
nodes=length(node);
nodes_r = nodes/nodes_th;
elements=length(element);

if inhomogeneity==1 %no more than 999 properties allowed in abaqus. either limit number of elements or create binning routine
    if elements>999
        return
    end
end


if depth>0
%% Create node definitions
    fprintf(fid,'**VA_TWIST&BEND GENERATED FILE\r\n');
    fprintf(fid,'*Node, Nset=All_Nodes\r\n');
    a=0;
    for z=1:z_nodes
        for t=1:nodes
        a=a+1;
        fprintf(fid, '%i,%g,%g,%g\r\n', a,node(t,1), -node(t,2), -(z-1)*depth/(z_nodes-1));
        end
    end
    
    fprintf(fid,'100000 , 0 ,0 ,0 \r\n');   
    fprintf(fid,'100001 , 0 ,0 ,%g \r\n',-depth);   
    
    
%% Create element connectivity matrix
    fprintf(fid,'*Element, Type=C3D8, Elset=All_Elem\r\n');
    a=0; 
    for z=1:z_nodes-1
        z
        for t=1:elements
            a=a+1;
            row1=((z-1)*nodes);
            row2=(z*nodes);
            % % %old syntax   %fprintf(fid, '%i,%i,%i,%i,%i,%i,%i,%i,%i\r\n', a, element(t,1)+row1, element(t,2)+row1,element(t,3)+row1,element(t,4)+row1,element(t,1)+row2, element(t,2)+row2,element(t,3)+row2,element(t,4)+row2 ); 
            fprintf(fid, '%i,%i,%i,%i,%i,%i,%i,%i,%i\r\n', a, element(t,1)+row1,element(t,4)+row1,element(t,3)+row1,element(t,2)+row1,element(t,1)+row2,element(t,4)+row2,element(t,3)+row2,element(t,2)+row2); 
            
        end
    end
%% Create element and node sets.


%% Create section and material properties  
    
if inhomogeneity==1
    for t=1:elements
        fprintf(fid,'*Elset, Elset=set_%i, Generate \r\n',t);
        fprintf(fid,'%i,%i,%i\r\n',t,t+elements*(z_nodes-2),elements);
        fprintf(fid,'*Solid Section, elset=set_%i, material=mat_%i\r\n',t,t);
        fprintf(fid,'1.,\r\n');
        fprintf(fid,'*Material, name=mat_%i\r\n',t);
        fprintf(fid,'*Elastic\r\n');
        fprintf(fid,'%g, %g\r\n',YM(t)+1,v(t));
    end
else
    fprintf(fid,'*Solid Section, elset=All_Elem, material=Bone\r\n');
    fprintf(fid,'1.,\r\n');
    fprintf(fid,'*Material, name=Bone\r\n');
    fprintf(fid,'*Elastic\r\n');
    fprintf(fid,'1000000., 0.32\r\n');
end

%% find nodes for boundary conditions and Create sets
	
    %Impactor nodes 
    mid_node_1st  = (floor(z_nodes/2))*nodes;
    fprintf(fid,'*Nset, Nset=Force_Nodes, Generate\r\n');
    fprintf(fid,'%i,%i,%i\r\n',mid_node_1st+1,mid_node_1st+nodes,1);   

    fprintf(fid,'*Elset, Elset=Supp_prox_Elems, Generate\r\n');
    fprintf(fid,'%i,%i,%i\r\n',1, elements, 1);
    
    fprintf(fid,'*Elset, Elset=Supp_dist_Elems, Generate\r\n');
    fprintf(fid,'%i,%i,%i\r\n',elements*(z_nodes-2)+1, elements*(z_nodes-1), 1);
    
    fprintf(fid,'*Nset, Nset=Ref_Node_Prox\r\n');
    fprintf(fid,'100000\r\n');   
    
    fprintf(fid,'*Nset, Nset=Ref_Node_Dist\r\n');
    fprintf(fid,'100001\r\n');   
        
    fprintf(fid,'*Nset, Nset=Prox_Node, Generate\r\n');
    fprintf(fid,'%i,%i,%i\r\n',1,nodes,1);   
    fprintf(fid,'*Nset, Nset=Dist_Node, Generate\r\n');
    fprintf(fid,'%i,%i,%i\r\n',nodes*(z_nodes-1)+1,nodes*z_nodes,1); 
    
    fprintf(fid,'*Rigid Body, Ref Node=Ref_Node_Prox, Elset=Supp_prox_Elems\r\n');
    
    fprintf(fid,'*Rigid Body, Ref Node=Ref_Node_Dist, Elset=Supp_dist_Elems\r\n');
    
    quarter_element_1st  = (floor((z_nodes-1)/4))*elements;
    fprintf(fid,'*Elset, Elset=sample, Generate\r\n');
    fprintf(fid,'%i,%i,%i\r\n',quarter_element_1st+1,quarter_element_1st+elements,1); 
    
    
    %Step information and output requests
    fprintf(fid,'*Step, name=Step-1\r\n');
    fprintf(fid,'*Static\r\n');
    fprintf(fid,'1., 1., 1e-05, 1.\r\n');
    fprintf(fid,'*Boundary\r\n');
    fprintf(fid,' Ref_Node_Prox, 1,3\r\n');
    fprintf(fid,' Ref_Node_Prox, 5,6\r\n');
    fprintf(fid,' Ref_Node_Dist, 1,2\r\n');
    fprintf(fid,'***\r\n');
    fprintf(fid,' Ref_Node_Dist, 6,6,0.001\r\n');

    fprintf(fid,'*Cload\r\n');
    fprintf(fid,' Force_Nodes, 2, %g\r\n', -(1/nodes) ) ;
    
    fprintf(fid,'*Restart, write, frequency=0\r\n');
    fprintf(fid,'*Output, field, variable=PRESELECT\r\n');
    fprintf(fid,'*Output, history, variable=PRESELECT\r\n');
    fprintf(fid,'*El Print, elset=sample\r\n');
    fprintf(fid,' S13, S23, S33\r\n');
    fprintf(fid,'*End Step\r\n');
 
end

fclose(fid);
