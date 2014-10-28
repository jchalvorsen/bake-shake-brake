close all
clear all
% We want to solve the linear elasticity problem
% grad(o(u)) = -f
% on a 3D surface thingy
% where f is the weight of each element (force that acts upon it)
% and dirichlet homogenous BCs where structure is attached to ground

addpath include/queen

% loading blocked elements
edges = load('coarse_bndry.m');     % tetraeders of the surface
tetr = load('coarse_element.m');
p = load('coarse_nodes.m');




% Write structure to vtf file
%writeVTF(p, tetr, 0, 'queen.vtf')

% Can plot, but takes a lot of resources
%figure
%tetramesh(tetr, p)


N = length(p);

% Declaring material constants
E = 1;
v = 0.1;

% Building the C matrix:

C1 =  E*v/((1+v)*(1-2*v))*ones(3,3) + E/(1+v)*eye(3);
C2 = E/(2*(1+v))*eye(3);
C = inv([ C1        , zeros(3,3);
      zeros(3,3), C2        ]); % pretty sure this should be inverted
  % found at www.rpi.edu/~des/3DElasticity.ppt, slide 24 (inverse function)

density = @(x,y,z) 1;

A = sparse(3*N,3*N);
b = zeros(3*N,1);
h = waitbar(0, 'In progress');
for i = 1:length(tetr)
    if mod(i,1000) == 0
        waitbar(i/length(tetr),h, 'In progress')
    end
    nodes = tetr(i,:);
    P = p(nodes,:);       % Active points in tetrangle
    
    % Calculating area
    Q = [[1;1;1;1], P];
    vol = abs(det(Q))/6;
    
    
    %% Getting stiffness matetrx
    % Finding constants in phi (basis function = [1, x, y, z] * c)
    c1 = Q\[1; 0; 0; 0];
    c2 = Q\[0; 1; 0; 0];
    c3 = Q\[0; 0; 1; 0];
    c4 = Q\[0; 0; 0; 1];
    %c = inv(Q);
    c = [c1, c2, c3, c4];
    
    % looping over x and y indexing
    Ak = zeros(12);
    for q = 1:4
        for w = 1:4
            % Basis function x = [c1 + c2*x + c3*y, c4*z; 0; 0]
            % y = [0; c1 + c2*x + c3*y + c4*z; 0]
            % z = [0; 0; c1 + c2*x + c3*y + c4*z]
            %e = [c2; 0; 0.5*c3];
            
            %{
             The eps have form of [-; -; -; -; -; -] (1x6)
             and have the form:
             eps = Du where, u (1x3),
             D = [ delx,    0,    0;
                   0,    dely,    0;
                   0,       0, delz;
                   dely, delx,    0;
                   0,    delz, dely;
                   delz,    0, delx ]      
            %}
            
            ex = @(i) [c(2,i); 0; 0; c(3,i); 0; c(4,i)];
            ey = @(i) [0; c(3,i); 0; c(2,i); c(4,i); 0];
            ez = @(i) [0; 0; c(4,i); 0; c(3,i); c(2,i)];
            
            
            % q represents basis function #1, while w is #2         
            f11 = ex(q)'*C*ex(w);
            f12 = ex(q)'*C*ey(w);
            f13 = ex(q)'*C*ez(w);
            f22 = ey(q)'*C*ey(w);
            f23 = ey(q)'*C*ez(w);
            f33 = ez(q)'*C*ez(w);
            
            fk = [f11, f12, f13;
                  f12, f22, f23;
                  f13, f23, f33];
              
            Ak(3*q-2:3*q, 3*w-2:3*w) = fk*vol;
        end
    end
    %Ak
    % Put element matetrx in right place
    % Map to right place:
    map(1:3:3*length(nodes)) = 3*nodes-2;
    map(2:3:3*length(nodes)) = 3*nodes-1;
    map(3:3:3*length(nodes)) = 3*nodes;
    A(map,map) = A(map,map) + Ak;

    %% Getting b vector
    % Finding functions phi*f:
    
    f1 = @(x) ([1, x(1), x(2), x(3)]*c1)*density(x);
    f2 = @(x) ([1, x(1), x(2), x(3)]*c2)*density(x);
    f3 = @(x) ([1, x(1), x(2), x(3)]*c3)*density(x);
    f4 = @(x) ([1, x(1), x(2), x(3)]*c4)*density(x);
    
    % Getting values:
    val1 = quadrature3d(P(1,:), P(2,:), P(3,:), P(4,:), 5, f1);
    val2 = quadrature3d(P(1,:), P(2,:), P(3,:), P(4,:), 5, f2);
    val3 = quadrature3d(P(1,:), P(2,:), P(3,:), P(4,:), 5, f3);
    val4 = quadrature3d(P(1,:), P(2,:), P(3,:), P(4,:), 5, f4);
    
    % Putting in right place:
    for j = 0:2
        b(3*nodes-2+j) = b(3*nodes-2+j) + [val1; val2; val3; val4]*vol*-9.81;
    end

end
%% Get A without boundary points
boundaryPoints = find((p(:,3) == 0)); % Dirichlet homogenous BC: f(boundaryPoints) = 0

% Setting cols of boundaryPoints equal to 0
A(boundaryPoints, :) = 0;
%A(:, boundaryPoints) = 0;
b(boundaryPoints) = 0;
A(boundaryPoints, boundaryPoints) = A(boundaryPoints, boundaryPoints) + speye(length(boundaryPoints), length(boundaryPoints));

% Solving the linear system
u_sol = A\b;

%% Plotting:
% figure
% hold on
% % Scatterplot of our result points
% for i = 1:length(p)
%     point = p(i,:);
%     plot3(point(1), point(2), u_sol(i),'*')
% end
% title('Scatterplot of results')


% Finding reference solution in points:
% u_ref = zeros(2*N,1);
% for i = 1:length(p)
%     point = p(i,:);
%     u_ref(2*i) = u(point);
%     u_ref(2*i-1) =u(point);
% end

figure
%plot(u_ref, '*-black')
hold on
plot(u_sol, '*-b')
title('Values in points')


figure
%subplot(2,1,1)
trisurf(tetr,p(:,1),p(:,2),p(:,3),u_sol(1:3:end));
view(2),axis equal,colorbar,title('FEM solution')


U = [u_sol(1:3:end), u_sol(2:3:end), u_sol(3:3:end)];



