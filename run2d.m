close all
clear all

addpath include

N = 50;

[p tri edge] = getDisk(N);
triplot(tri,p(:,1), p(:,2))

f = @(x,y) -8*pi*cos(2*pi*(x^2+y^2)) +16*pi^2*(x^2+y^2)*sin(2*pi*(x^2+y^2));
u = @(x) sin(2*pi*(x(1)^2+x(2)^2));

Z = tri(1,:);
P = p(Z,:);

M = [[1;1;1], P];
coef = inv(M);
psi1 = @(x,y) coef(1,1) + coef(2,1)*x + coef(3,1)*y;

QQ = zeros(100);
z = linspace(-1,1);
for i = 1:100
    for j = 1:100
        QQ(i,j) = psi1(z(i),z(j));
    end
end
%surf(QQ)


%figure
hold on
for i = 1:length(edge)
    left = p(edge(i,1),:);
    right = p(edge(i,2),:);
    x = linspace(left(1),right(1));
    y = linspace(left(2), right(2));
    %plot(x,y)
end


%% Trying to make all in a loop:
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

u_sol = zeros(length(internalPoints),1);
for i = internalPoints
    point = p(i,:);
    %norm(point,2)
    u_sol(i) = u(point);
end
% Solving system:
u_is = A_int\b_int;

figure
plot(u_is, '*-r')
hold on
plot(u_sol, 'black')
