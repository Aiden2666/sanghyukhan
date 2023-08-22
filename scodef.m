function [XI YI ZI]=scodef(X,Y,Z,C,dC,degree,radius,knots)

%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu

% Pack coordinates
Q=[X Y Z];

% Ci are the contraints
% for each constraint we are going to evaluate the fi function and then
% pack them all together in a matrix form.
% Ri is the of influence of each contraint
% knots is the knot vector
% nk the number of knots
% p is the order
% dCi is the displacment of the constrain
% For decreasing the execution time we are going first to evaluate the
% distance function (||Q-Ci||/Ri) for all Q's, then sort it, evaluate at
% once all fi's and then reverse the sorting. SPCOL funtion can evaluate
% recursively all tau's of the b-spline basis function.

[n_C,tmp]=size(C);
R=ones(n_C,1).*radius; % Radius of influence

nk=ones(1,n_C).*knots'; % Number of knots should nk>=p+1
% Calculate uniform knot vectors
for i=1:length(nk)
    allknots{i,:}=linspace(-1,1,nk(i)); % Knot vectors
end

p=ones(1,n_C).*degree'; % B-Spline degree degree should be p<=nk-1

fC=[];
fQ=[];

for i=1:n_C
    knots=allknots{i};
    
    % Calculate Ci(Q) where Q are the contraints
    CiQ=sqrt((C(:,1)-repmat(C(i,1),size(C,1),1)).^2+...
             (C(:,2)-repmat(C(i,2),size(C,1),1)).^2+...
             (C(:,3)-repmat(C(i,3),size(C,1),1)).^2)/R(i);

    % First sort constriants and keep the indices
    [CiQ,I]=sortrows(CiQ);
    
    % Kill dupes and keep their indices
    [CiQ,K,L]=unique(CiQ);
    
    % Estimate a basis function values for the ith contraint
    fiC=spcol(knots,p(i),CiQ);
    
    % Expand fiC by adding dupes that were earlier rejected with UNIQUE
    fiC=fiC(L,:);
    
    % resort fiQ with indices I
    fiC(I,:)=fiC;

    % Create the full matrix for the constraints
    fC=[fC fiC];

    % Calculate Ui(Q) where Q are the all the points
    UiQ=sqrt((Q(:,1)-repmat(C(i,1),size(Q,1),1)).^2+...
             (Q(:,2)-repmat(C(i,2),size(Q,1),1)).^2+...
             (Q(:,3)-repmat(C(i,3),size(Q,1),1)).^2)/R(i);
    
        
    % First sort UiQ and keep the indices
    [UiQ,J]=sortrows(UiQ);
    
    % Kill dupes and keep their indices
    [UiQ,K,L]=unique(UiQ);
   
    % Find which UiQ have a value less or equal to 1, as the knot vector
    % spans from -1 to 1, and outside all values are equal to zero
    ROI_UiQ=abs(UiQ)<=1;
    I=find(ROI_UiQ==1);
    
    % Initiate fiQ matrix
    fiQ=zeros(length(Q),length(knots)-p(i));
    
    
    % Estimate a unigue value for each contraint
    fiQ(I,:)=spcol(knots,p(i),UiQ(I));
    
    % expand fiQ by adding dupes that were earlier rejected with UNIQUE
    fiQ=fiQ(L,:);
    
    % re-sort fiQ with indices J
    fiQ(J,:)=fiQ;
    fQ=[fQ fiQ];
end
% Find deformation matrix M by inverse functions fc
M=pinv(fC)*dC;

% Multiply with weighted functions fQ and calculate the displacments for
% all points
dQ=fQ*M;

% Add displacments to initial points
Q=Q+dQ;

% number of vertices
n_v=length(X);
% reshape and draw
XI=reshape(Q(:,1),[],n_v)';
YI=reshape(Q(:,2),[],n_v)';
ZI=reshape(Q(:,3),[],n_v)';


