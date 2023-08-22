function [shear_stress_e, resultant_shear_stress, K_mod] = solve_FE_problem (node, element, YM_e, G_e, v_e, I, Fx, Fy, theta)

%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu

nu=v_e(1);  %Poissons ration taken as constant

% initialize variables
n_nodes    = length(node)   ; 
n_elements = length(element);
icount = zeros(n_nodes,1); K_mod=0;

stress_x=zeros(4,1); 
stress_y=zeros(4,1);
shear_stress_n = zeros(n_nodes   ,2);
shear_stress_e = zeros(n_elements,2);
K       = zeros(n_nodes, n_nodes);    % Stiffness Matrix
F       = zeros(n_nodes,1);           % Force Vector

% Begin Caculations
Ix  = I(1,1);Iy  = I(2,2);Ixy = I(2,1);

% Find A and B as functions of F and I.  see Barber JR, Elasticity: eq. (17.6)
A=(Fx*Ix - Fy*Ixy) / (Ix*Iy - Ixy^2);
B=(Fy*Iy - Fx*Ixy) / (Ix*Iy - Ixy^2);

% Calculate normal vectors for each boundary node
boundary_node_map = find_boundary (node, element);

%YM_n = elemental2nodal (YM_e, node, element);
G_n  = elemental2nodal (G_e,  node, element);
%v_n =  elemental2nodal (v_e,  node, element);

%% Calculate shape functions and derivatives evaluated at gauss points. every line is one of the 4 shape functions

% Create isoparametric coords for isoparamatric nodes
ksi_0 = [-1  1  1  -1]; eta_0 = [-1  -1  1  1]; %changed it 8/2010
element=fliplr(element);%added 8/2010

% Gauss Integration Points
ksi = ksi_0/sqrt(3);  eta = eta_0/sqrt(3);

N       = ( 1 + ksi_0' * ksi ) .* ( 1 + eta_0' * eta )  / 4;% shape functions            [4x4]
N_ksi   = repmat(ksi_0',1,4)   .* ( 1 + eta_0' * eta )  / 4;% derivative wrt to x (parametric x)
N_eta   = repmat(eta_0',1,4)   .* ( 1 + ksi_0' * ksi )  / 4;% derivative wrt to y (parametric y)

% Gauss Integration Points for the edge (line) segments
e=[-1/sqrt(3); 1/sqrt(3)];
N_g=[0.5*(1-e),  0.5*(1+e)]; %Boundary shape functions

%% Form element stiffness and force vector
for t=1:n_elements                  % for every element
    
    k_e = zeros(4,4); f_f = zeros(4,1); f_g = zeros(4,1);

   %% Element Stiffness matrix    
    for j=1:4                       % for every gauss integration point
        
         %isoparametric element, therefore:   
         x      = N(:,j)'*node(element(t,:),1);        % x = sum (N * x_a)               [1x1]
         y      = N(:,j)'*node(element(t,:),2);        % y = sum (N * y_a)               [1x1]
         x_ksi  = N_ksi(:,j)'*node(element(t,:),1);    % x_ksi = sum(N_ksi * x_a)        [1x1]
         y_ksi  = N_ksi(:,j)'*node(element(t,:),2);    % y_ksi = sum(N_ksi * y_a)        [1x1]
         x_eta  = N_eta(:,j)'*node(element(t,:),1);    % x_eta = sum(N_eta * x_a)        [1x1]
         y_eta  = N_eta(:,j)'*node(element(t,:),2);    % y_eta = sum(N_eta * y_a)        [1x1]
         
         a_jacobian = x_ksi*y_eta-x_eta*y_ksi ;        % j=det[u_1,1 u_1,2]Hughes 3.9.10[1x1]

         X_ksi  = [x_ksi x_eta ; y_ksi y_eta];         %                                [2x2]
         X      = inv(X_ksi);                          %                                [2x2]
         
         N_     = [N_ksi(:,j) N_eta(:,j)] * X;         %                                [4x2]                                  
         N_x    = N_(:,1);                             %                                [4x1]
         N_y    = N_(:,2);                             %                                [4x1]   
         
         G   = N(:,j)'*G_n(element(t,:));               % value of G at the              [1x1]
                                                       % j integration point
         
         G_x = N_x' * G_n(element(t,:));               % derivative of G wrt x          [1x1]
         G_y = N_y' * G_n(element(t,:));               % derivative of G wrt y          [1x1]         
                   
         % based on the matrix form of the weak form, enter the components
         k_e = k_e + G*(N_x * N_x'  +   N_y * N_y'    ) * a_jacobian; %                 [4x4]                                                     [4x4]  
         
	     % 8/2010 changed the last component from (y -x) to [y -x]

         f_f = f_f +  N(:,j) * ( 2 * nu * x * y *  (B * G_x + A * G_y ) - 2 * G * (A * x + B * y)) * a_jacobian +  N_ * theta * G * ([ y; -x]) * a_jacobian;  %  [4x1]           
      %  f_f = f_f +  N(:,j) * ( 2 * nu * x * y *  (B * G_x + A * G_y ) - 2 * G * (A * x + B * y) + theta * (G_x * y - G_y *x) ) * a_jacobian;  % [4x1]           
        
          
    end

    %% Boundary condition applied only to the boundary nodes
for i=1:4
     if(boundary_node_map(element(t,i))==1)
         j=i+1; if j==5; j=1; end
        
         if(boundary_node_map(element(t,j))==1)
            e1=element(t,i); e2=element(t,j);           % This is the global node number of the boundary nodes
            l_jacobian=0.5*norm(node(e1,:)-node(e2,:)); % Line segment jacobian (=the distance between nodes)
            x_e=[node(e1,1);node(e2,1)];                % Pick up the node locations
            y_e=[node(e1,2);node(e2,2)];
            
            g_e=[G_n(e1);G_n(e2)];                      % Get shear modulus values for the two boundary nodes
            g=N_g'*g_e;                                 % Calculate shear mod at int points
            x=N_g'*x_e;                                 % Find int. point locations
            y=N_g'*y_e;
            
           % tt=[(theta*y+2*B*nu*x.*y)';(-theta*x+2*A*nu*x.*y)'];
            tt=[(2*B*nu*x.*y)';(2*A*nu*x.*y)'];
            ne = [ -(y_e(1)-y_e(2)) ;  (x_e(1)-x_e(2)) ];
            ne = -ne/norm(ne);
%            figure(12); quiver(mean(x_e), mean(y_e), ne(1), ne(2));hold on
             dummy=N_g'*g.*(tt'*ne)*l_jacobian;
             f_g(i)=dummy(1);
             f_g(j)=dummy(2);
        end
     end
end
    

    f_e =  f_f + f_g; %NORMALS ARE FLIPPED
    

%% Assemble in global stiffness matrix
    LM=element(t,:);
    K(LM,LM) = K(LM,LM) + k_e;
    F(LM)    = F(LM)    + f_e;

    
end %end of loop for each element
%%

% Solve the linear system of equations to find the warping function value
% for each int. point : K d = F

d=K\F;


%%

d_x = zeros(4,1);
d_y = zeros(4,1);

 %% PostProcessing
for t=1:n_elements

    for j=1:4
         %isoparametric element, therefore:   
         x      = N(:,j)'*node(element(t,:),1);        % x = sum (N * x_a)               [1x1]
         y      = N(:,j)'*node(element(t,:),2);        % y = sum (N * y_a)               [1x1]
         x_ksi  = N_ksi(:,j)'*node(element(t,:),1);    % x_ksi = sum(N_ksi * x_a)        [1x1]
         y_ksi  = N_ksi(:,j)'*node(element(t,:),2);    % y_ksi = sum(N_ksi * y_a)        [1x1]
         x_eta  = N_eta(:,j)'*node(element(t,:),1);    % x_eta = sum(N_eta * x_a)        [1x1]
         y_eta  = N_eta(:,j)'*node(element(t,:),2);    % y_eta = sum(N_eta * y_a)        [1x1]
         
         jacobian = x_ksi*y_eta-x_eta*y_ksi;           % j=det[u_1,1 u_1,2]Hughes 3.9.10[1x1]
         
         X_ksi  = [x_ksi x_eta ; y_ksi y_eta];         %                                [2x2]
         X      = inv(X_ksi);                          %                                [2x2]
         
         N_     = [N_ksi(:,j) N_eta(:,j)] * X;         %                                [4x2]                                  
         N_x    = N_(:,1);                             %                                [4x1]
         N_y    = N_(:,2);                             %                                [4x1]    
         
         d_x(j) = N_x' * d(element(t,:));                % derivative of d wrt x          [4x1]
         d_y(j) = N_y' * d(element(t,:));                % Sum(Sum(N_ii*xi)*u_i)          [4x1]        
         
         
         stress_x(j) = d_x(j) - theta*y - 2*B*nu*x*y;
         stress_y(j) = d_y(j) + theta*x - 2*A*nu*x*y;
     
         if theta~=0
             K_mod=K_mod +  (x*stress_y(j)-y*stress_x(j)) /theta * jacobian *  G_e(t);   %  4/2010 moved it before the end statement and changed it
         end
         
    end

    %extrapolate d_x and d_y values from gauss integration points to nodes
    stress_x_element=linsolve(N,stress_x).*G_n(element(t,:));
    stress_y_element=linsolve(N,stress_y).*G_n(element(t,:)); 
    
    %average over the nodes by summing contibution from each element and
    %then divide by the number of elements associated to this node
    shear_stress_n(element(t,:),1) = shear_stress_n(element(t,:),1) + stress_x_element;
    shear_stress_n(element(t,:),2) = shear_stress_n(element(t,:),2) + stress_y_element;  
    icount(element(t,:)) = icount(element(t,:)) + 1;
    
    shear_stress_e(t,1) = G_e(t)*mean(stress_x);
    shear_stress_e(t,2) = G_e(t)*mean(stress_y);
    
    
end


 shear_stress_n(:,1)=(shear_stress_n(:,1)./icount);
 shear_stress_n(:,2)=(shear_stress_n(:,2)./icount);

 resultant_shear_stress    = sqrt(shear_stress_n(:,1).^2+shear_stress_n(:,2).^2);

