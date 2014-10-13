close all
clear all

addpath include

N = 50;

[p tri edge] = getDisk(N);
triplot(tri,p(:,1), p(:,2))

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
    Jac = [P(2,:)' - P(1,:)', P(3,:)' - P(1,:)'];
    
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

A_int = A(internalPoints,internalPoints);

% Getting internal nodes of b:
b_int = b(internalPoints);

u_int = A_int\b_int;

U_sol = zeros(100,100);
z = linspace(-1,1);
for i = 1:1
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
    
    
    
    
end




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
figure
surf(z,z,U)

figure
hold on
% plot our result function
for i = 1:length(internalPoints)
    point = p(internalPoints(i),:);
    plot3(point(1), point(2), u_int(i))
end
    