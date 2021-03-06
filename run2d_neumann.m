close all
clear all

addpath include

N = 700;

[p tri edge] = getDisk(N);

% Declaring functions
f = @(x,y) -8*pi*cos(2*pi*(x^2+y^2)) +16*pi^2*(x^2+y^2)*sin(2*pi*(x^2+y^2));
g = @(x,y) 4*pi*cos(2*pi*(x^2+y^2));
% Still same solution
u = @(x) sin(2*pi*(x(1)^2+x(2)^2));

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

% Splitting edge in neumann and dirichlet boundaries
boundaryPoints = edge(:,1);
actualBP = p(boundaryPoints,:);

% Might need refining
edgeN = edge(actualBP(:,2) >= 0,:);
edgeD = edge(actualBP(:,2) < 0,:);

% Neumann boundary conditions:
for i = 1:length(edge)
    nodes = edge(i,:);
    P = p(nodes,:); 
    Q = [[1;1], P];
    
    phi1 = @(x) ([1, x(1), x(2)]*(Q\[1; 0]));
    phi2 = @(x) ([1, x(1), x(2)]*(Q\[0; 1]));
    
    f1 = @(x) phi1(x)*g(x(1),x(2));
    f2 = @(x) phi2(x)*g(x(1),x(2));
    
    val1 = quadrature1D_modified(P(1,:), P(2,:), 4, f1);
    val2 = quadrature1D_modified(P(1,:), P(2,:), 4, f2);
    
    % Adding to b if we are in neumann area
    b(nodes(1)) = b(nodes(1)) + val1*(p(nodes(1),2) >= 0);
    b(nodes(2)) = b(nodes(2)) + val2*(p(nodes(2),2) >= 0);
end



% Dirichlet BC: Setting rows and cols of boundaryPoints equal to 0
boundaryPointsD = edgeD(:,1);
A(boundaryPointsD, :) = 0;
% A(:, boundaryPointsD) = 0;
b(boundaryPointsD) = 0;
A(boundaryPointsD, boundaryPointsD) = speye(length(boundaryPointsD), length(boundaryPointsD));

% Solving the linear system
u_sol = A\b;

%% Finding solution on M*M grid and plotting
M = 100;
U_sol = zeros(M,M);
z = linspace(-1,1,M);
for i = 1:length(tri)
    nodes = tri(i,:);
    P = p(nodes,:);       % Active points in triangle
   
    max_x = max(P(:,1));
    min_x = min(P(:,1));
    max_y = max(P(:,2));
    min_y = min(P(:,2));
    
    x = z(max_x >= z & z >= min_x);
    y = z(max_y >= z & z >= min_y);
    
    % Finding basis function:
    Q = [[1;1;1], P];
    phi1 = @(x) [1, x(1), x(2)]*(Q\[1; 0; 0]);
    phi2 = @(x) [1, x(1), x(2)]*(Q\[0; 1; 0]);
    phi3 = @(x) [1, x(1), x(2)]*(Q\[0; 0; 1]);
       
    uu = zeros(length(x), length(y));
    for j = 1:length(x)
        for k = 1:length(y)
            point = [x(j), y(k)];
            phi1v = phi1(point);
            phi2v = phi2(point);
            phi3v = phi3(point);
            q = [phi1v, phi2v, phi3v];
            % Check if point is inside the triangle
            if (phi1v <= 1 && phi1v >= 0) && (phi2v <= 1 && phi2v >= 0) && (phi3v <= 1 && phi3v >= 0)
                % Add if inside triangle
                uu(j,k) = q*u_sol(nodes);
            end
     
        end
    end   
    
    % merging uu into U_sol
    xstart = find(z ==x (1));
    xend = xstart + length(x) - 1;
    ystart = find(z == y(1));
    yend = ystart + length(y) - 1;

    U_sol(xstart:xend, ystart:yend) = uu + U_sol(xstart:xend, ystart:yend);
end
 
figure
subplot(1,2,1)
surf(z, z, U_sol)
title('FEM solution')


% Plotting reference solution:
U = zeros(100);
z = linspace(-1,1);
for i = 1:100
    for j = 1:100
        pz = [z(i), z(j)];
        if norm(pz,2) <= 1
            U(i,j) = u(pz);
        end
    end
end
subplot(1,2,2)
surf(z,z,U)
title('Reference solution')

figure
hold on
% Scatterplot of our result points
for i = 1:length(p)
    point = p(i,:);
    plot3(point(1), point(2), u_sol(i),'*')
end
title('Scatterplot of results')


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


figure
subplot(2,1,1)
trisurf(tri,p(:,1),p(:,2),0*p(:,1),u_sol,'edgecolor','k','facecolor','interp');
view(2),axis equal,colorbar,title('FEM solution')

subplot(2,1,2)
trisurf(tri,p(:,1),p(:,2),0*p(:,1),u_ref,'edgecolor','k','facecolor','interp');
view(2),axis equal,colorbar, title('Analytical solution')

