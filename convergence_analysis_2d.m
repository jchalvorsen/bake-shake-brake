
clear all
close all

addpath include

% Declaring functions
f = @(x,y) -8*pi*cos(2*pi*(x^2+y^2)) +16*pi^2*(x^2+y^2)*sin(2*pi*(x^2+y^2));
u = @(x) sin(2*pi*(x(1)^2+x(2)^2));

% 
% N = 300;

e_n = 8;

error = zeros(e_n,1);
N_vec = zeros(e_n,1);

for n = 1:e_n
    
    N = 20*2^(n+1);

[p tri edge] = getDisk(N);

A = sparse(N,N);
b = zeros(N,1);

for i = 1:length(tri)
    nodes = tri(i,:);
    P = p(nodes,:);       % Active points in triangle
    
    % Calculating area
    Q = [[1;1;1], P];
    area = 0.5*abs(det(Q));
    
    %% Getting stiffness matrix 
    % Finding constants in phi (basis function = [1, x, y] * c)
    c1 = Q\[1; 0; 0];
    c2 = Q\[0; 1; 0];
    c3 = Q\[0; 0; 1];
    
    divphi1 = c1(2:3);
    divphi2 = c2(2:3);
    divphi3 = c3(2:3);
    divphi = [divphi1, divphi2, divphi3];
    
    % Ak is element matrix for this element 
    Ak = divphi'*divphi*area;
      
    % Put element matrix in right place
    A(nodes,nodes) = A(nodes,nodes) + Ak;
    
    %% Getting b vector 
    % Finding functions phi*f:
    f1 = @(x) ([1, x(1), x(2)]*c1)*f(x(1),x(2));
    f2 = @(x) ([1, x(1), x(2)]*c2)*f(x(1),x(2));
    f3 = @(x) ([1, x(1), x(2)]*c3)*f(x(1),x(2));
    
    % Getting values:
    val1 = quadrature2d(P(1,:), P(2,:), P(3,:), 4, f1);
    val2 = quadrature2d(P(1,:), P(2,:), P(3,:), 4, f2);
    val3 = quadrature2d(P(1,:), P(2,:), P(3,:), 4, f3);
    
    % Putting in right place:
    b(nodes) = b(nodes) + [val1; val2; val3];  
end
%% Get A without boundary points
boundaryPoints = edge(:,1);

% Setting rows and cols of boundaryPoints equal to 0
A(boundaryPoints, :) = 0;
% A(:, boundaryPoints) = 0;
b(boundaryPoints) = 0;
A(boundaryPoints, boundaryPoints) = speye(length(boundaryPoints), length(boundaryPoints));

% Solving the linear system
u_sol = A\b;

%Finding reference solution in points:
u_ref = zeros(N,1);
for i = 1:length(p)
    point = p(i,:);
    u_ref(i) = u(point); 
end





    error(n) = norm(u_ref-u_sol,2);
    N_vec(n) = 1/N;
    
end
%% Plotting:

figure
loglog(N_vec, error, '*-r');
title('Loglogplot of error');
hold on
loglog(N_vec, N_vec);

% %% Plotting error vector
% figure; loglog(error, '-*')
% title('Loglog-plot of error')
% 
% figure; plot(error, '-*r')
% title('Plot of error')
%% Finding solution on M*M grid and plotting


% figure
% plot(u_ref, '*-black')
% hold on
% plot(u_sol, '*-b')
% title('Values in points')
% 
% 
% figure;
% % Error between u_ref and u_sol
% plot(u_ref-u_sol, '*r')
% title('Error for values in points')



