close all
clear all

addpath include

N = 30000;

[p tetr edge] = getSphere(N);



trep = TriRep(tetr, p);
[tr, Xb] = freeBoundary(trep);
trisurf(tr, Xb(:,1), Xb(:,2), Xb(:,3), 'FaceColor', 'red','FaceAlpha', 0.8); 
 

% Declaring functions
f = @(x,y,z) -12*pi*cos(2*pi*(x^2+y^2+z^2)) +16*pi^2*(x^2+y^2+z^2)*sin(2*pi*(x^2+y^2+z^2));
u = @(x) sin(2*pi*(x(1)^2 + x(2)^2 + x(3)^2));

A = sparse(N,N);
b = zeros(N,1);

% divpsi1 = [-1;-1;-1];
% divpsi2 = [1; 0; 0];
% divpsi3 = [0; 1; 0];
% divpsi4 = [0; 0; 1];
% divpsi = [divpsi1, divpsi2, divpsi3, divpsi4];

for i = 1:length(tetr)
    nodes = tetr(i,:);
    P = p(nodes,:);       % Active points in triangle
    
    % Calculating area
    Q = [[1;1;1;1], P];
    vol = abs(det(Q))/6;
    
    %% Getting stiffness matrix 
    % Finding constants in phi (basis function = [1, x, y] * c)
    c1 = Q\[1; 0; 0; 0];
    c2 = Q\[0; 1; 0; 0];
    c3 = Q\[0; 0; 1; 0];
    c4 = Q\[0; 0; 0; 1];
    
    divphi1 = c1(2:4);
    divphi2 = c2(2:4);
    divphi3 = c3(2:4);
    divphi4 = c4(2:4);
    divphi = [divphi1, divphi2, divphi3, divphi4];
    
    % Ak is element matrix for this element 
    Ak = divphi'*divphi*vol;
      
    % Put element matrix in right place
    A(nodes,nodes) = A(nodes,nodes) + Ak;
    
    %% Getting b vector 
    % Finding functions phi*f:
    f1 = @(x) ([1, x(1), x(2), x(3)]*c1)*f(x(1),x(2), x(3));
    f2 = @(x) ([1, x(1), x(2), x(3)]*c2)*f(x(1),x(2), x(3));
    f3 = @(x) ([1, x(1), x(2), x(3)]*c3)*f(x(1),x(2), x(3));
    f4 = @(x) ([1, x(1), x(2), x(3)]*c4)*f(x(1),x(2), x(3));
    
    % Getting values:
    val1 = quadrature3d(P(1,:), P(2,:), P(3,:), P(4,:), 5, f1);
    val2 = quadrature3d(P(1,:), P(2,:), P(3,:), P(4,:), 5, f2);
    val3 = quadrature3d(P(1,:), P(2,:), P(3,:), P(4,:), 5, f3);
    val4 = quadrature3d(P(1,:), P(2,:), P(3,:), P(4,:), 5, f4);
    
    % Putting in right place:
    b(nodes) = b(nodes) + [val1; val2; val3; val4];      
end

%% Get A without boundary points
AllboundaryPoints = [edge(:,1); edge(:,2); edge(:,3)];
boundaryPoints = unique(AllboundaryPoints);

% Setting rows and cols of boundaryPoints equal to 0
A(boundaryPoints, :) = 0;
% A(:, boundaryPoints) = 0;
b(boundaryPoints) = 0;
A(boundaryPoints, boundaryPoints) = speye(length(boundaryPoints), length(boundaryPoints));

% Solving the linear system
u_sol = A\b;


% Test of 3d viewing
writeVTF(p, tetr, u_sol, 'test.vtf')

% Finding reference solution in points:
u_ref = zeros(N,1);
for i = 1:length(p)
    point = p(i,:);
    u_ref(i) = u(point); 
end

figure
plot(u_ref, '*-black')
hold on
plot(u_sol, '*-b')
title('Values in points')


