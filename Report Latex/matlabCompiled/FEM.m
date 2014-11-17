function [ u_sol ] = FEM( p, data, E, v, loadVector, newtonBoundary)
%FEM solve the linear elasticity problem using the finite element method
% Returns the solution as a 3*N vector (N in each spatial direction)
% If you want to split to 3*N matrix with displacement in each direction
% use U = [u(1:3:end), U(

N = length(p);

% Building the C matrix:

C1 =  E*v/((1+v)*(1-2*v))*ones(3,3) + E/(1+v)*eye(3);
C2 = E/(2*(1+v))*eye(3);
C = [ C1        , zeros(3,3);
    zeros(3,3), C2        ]; % pretty sure this should be inverted
% found at www.rpi.edu/~des/3DElasticity.ppt, slide 24 (inverse function)

A = zeros(3*N,3*N);
b = zeros(3*N,1);

for i = 1:length(data)
    nodes = data(i,1:4);
    P = p(nodes,:);       % Active points in tetraeder
    
    % Calculating area
    Q = [[1;1;1;1], P];
    vol = abs(det(Q))/6;
     
    %% Getting stiffness matetrx
    % Finding constants in phi (basis function = [1, x, y, z] * c)
%     c1 = Q\[1; 0; 0; 0];
%     c2 = Q\[0; 1; 0; 0];
%     c3 = Q\[0; 0; 1; 0];
%     c4 = Q\[0; 0; 0; 1];
%     c = [c1, c2, c3, c4];
    c = inv(Q);
      
    phiGrad = c(2:end,:);
    B = zeros(6,12);
    B([1,4,5],1:3:10) = phiGrad;
    B([4,2,6],2:3:11) = phiGrad;
    B([5,6,3],3:3:12) = phiGrad;
    
    
    
    % u_e: displacement field:
    map(1:3:3*length(nodes)) = 3*nodes-2;
    map(2:3:3*length(nodes)) = 3*nodes-1;
    map(3:3:3*length(nodes)) = 3*nodes;
    
    % Adding element stiffness matrix to main matrix
    A(map,map) = A(map,map) + B'*C*B*vol ;
    
    %% Getting b vector 
    % doing the gauss quadrature of order 1
    gaussQuad = loadVector(data(i,5))*vol;
    
    % Putting in right place:
    % only want to add gravity compononent to z-dir
    force = [0; 0; -9.81]*vol*gaussQuad;
    b(map) = b(map) + repmat(force,4,1);
    
    
end

%% Get A without boundary points
map2(1:3:3*length(newtonBoundary)) = 3*newtonBoundary-2;
map2(2:3:3*length(newtonBoundary)) = 3*newtonBoundary-1;
map2(3:3:3*length(newtonBoundary)) = 3*newtonBoundary;

% Setting cols of boundaryPoints equal to 0
A(map2, :) = 0;
b(map2) = 0;
A(map2, map2) = A(map2, map2) + speye(length(map2), length(map2));

% Making A sparse so linear system will be solved fast
Asp = sparse(A);
% Solving the linear system
u_sol = Asp\b;

end

