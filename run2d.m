close all
clear all

addpath include

N = 100;

[p tri edge] = getDisk(N);
%triplot(tri,p(:,1), p(:,2))

f = @(x,y) -8*pi*cos(2*pi*(x^2+y^2)) +16*pi^2*(x^2+y^2)*sin(2*pi*(x^2+y^2));
u = @(x) sin(2*pi*(x(1)^2+x(2)^2));

A = sparse(N,N);

divpsi1 = [-1;-1];
divpsi2 = [1; 0];
divpsi3 = [0; 1];
divpsi = [divpsi1, divpsi2, divpsi3];

for i = 1:length(tri)
    thisTri = tri(i,:);
    P = p(thisTri,:);       % Active points in triangle
    
    % Calculating area
    areaMatrix = [[1;1;1], P];
    area = 0.5*abs(det(areaMatrix));
    
    % Constructing jacobian and solving shit
    Jac = [P(2,:) - P(1,:); P(3,:) - P(1,:)];
    
    % Ak is element matrix for this element
    Ak = transpose(Jac\divpsi)*(Jac\divpsi)*area;
    
    % Put elemenet matrix in right place
    A(thisTri,thisTri) = A(thisTri,thisTri) + Ak;
    
end

%% Trying to get b vector:

b = zeros(N,1);
for i = 1:length(tri)
    thisTri = tri(i,:);
    P = p(thisTri,:);       % Active points in triangle
    p1 = P(1,:);
    p2 = P(2,:);
    p3 = P(3,:);
    
    
    % Finding basis function:
    Q = [[1;1;1], P];
    phi1 = @(x) ([1, x(1), x(2)]*(Q\[1; 0; 0]))*f(x(1),x(2));
    phi2 = @(x) ([1, x(1), x(2)]*(Q\[0; 1; 0]))*f(x(1),x(2));
    phi3 = @(x) ([1, x(1), x(2)]*(Q\[0; 0; 1]))*f(x(1),x(2));
    
    % Getting values:
    val1 = quadrature2d(p1, p2, p3, 4, phi1);
    val2 = quadrature2d(p1, p2, p3, 4, phi2);
    val3 = quadrature2d(p1, p2, p3, 4, phi3);
    
    % Putting in right place:
    b(thisTri) = b(thisTri) + [val1; val2; val3];
    
end
% Get A without boundary points (internal A):

allPoints = 1:N;
boundaryPoints = edge(:,1);
internalPoints = setdiff(allPoints,boundaryPoints);

% Setting rows and cols of boundaryPoints equal to 0
A(boundaryPoints, :) = 0;
A(:, boundaryPoints) = 0;
b(boundaryPoints) = 0;
A(boundaryPoints, boundaryPoints) = speye(length(boundaryPoints), length(boundaryPoints));


% Getting internal nodes of b:
b_int = b(internalPoints);
det(A)
u_int = A\b;

% Finding solution on M*M grid
M = 100;
U_sol = zeros(M,M);
z = linspace(-1,1,M);
for i = 1:length(tri)
    thisTri = tri(i,:);
    P = p(thisTri,:);       % Active points in triangle
    p1 = P(1,:);
    p2 = P(2,:);
    p3 = P(3,:);
   
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
            p1v = phi1(point);
            p2v = phi2(point);
            p3v = phi3(point);
            q = [p1v, p2v, p3v];
            % Check if point is inside the triangle
            if (p1v <= 1 && p1v >= 0) && (p2v <= 1 && p2v >= 0) && (p3v <= 1 && p3v >= 0)
                % Add if inside triangle
                uu(j,k) = q*u_int(thisTri);
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
    plot3(point(1), point(2), u_int(i),'*')
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
plot(u_int, '*-b')
title('Values in points')


figure
subplot(2,1,1)
trisurf(tri,p(:,1),p(:,2),0*p(:,1),u_int,'edgecolor','k','facecolor','interp');
view(2),axis equal,colorbar,title('FEM solution')

subplot(2,1,2)
trisurf(tri,p(:,1),p(:,2),0*p(:,1),u_ref,'edgecolor','k','facecolor','interp');
view(2),axis equal,colorbar, title('Analytical solution')

