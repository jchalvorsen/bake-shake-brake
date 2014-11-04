close all
clear all

addpath include

% Declaring functions
f = @(x,y,z) -12*pi*cos(2*pi*(x^2+y^2+z^2)) +16*pi^2*(x^2+y^2+z^2)*sin(2*pi*(x^2+y^2+z^2));
u = @(x) sin(2*pi*(x(1)^2 + x(2)^2 + x(3)^2));


e_n = 4;
error = zeros(e_n, 1);

for n = 1:e_n
    
    N = 1000*2^(n);
    
    [p tetr edge] = getSphere(N);
    
    trep = TriRep(tetr, p);
    [tr, Xb] = freeBoundary(trep);
    trisurf(tr, Xb(:,1), Xb(:,2), Xb(:,3), 'FaceColor', 'red','FaceAlpha', 0.8);
    
    
    A = sparse(N,N);
    b = zeros(N,1);
    
    divpsi1 = [-1;-1;-1];
    divpsi2 = [1; 0; 0];
    divpsi3 = [0; 1; 0];
    divpsi4 = [0; 0; 1];
    divpsi = [divpsi1, divpsi2, divpsi3, divpsi4];
    
    for i = 1:length(tetr)
        nodes = tetr(i,:);
        P = p(nodes,:);       % Active points in triangle
        
        %% Getting stiffness matrix
        % Calculating area
        Q = [[1;1;1;1], P];
        vol = abs(det(Q))/6;
        
        % Constructing jacobian and solving shit
        Jac = [P(2,:) - P(1,:); P(3,:) - P(1,:); P(4,:) - P(1,:)];
        
        % Ak is element matrix for this element
        Ak = transpose(Jac\divpsi)*(Jac\divpsi)*vol;
        
        % Put element matrix in right place
        A(nodes,nodes) = A(nodes,nodes) + Ak;
        
        %% Getting b vector
        % Finding basis function:
        
        phi1 = @(x) ([1, x(1), x(2), x(3)]*(Q\[1; 0; 0; 0]))*f(x(1), x(2), x(3));
        phi2 = @(x) ([1, x(1), x(2), x(3)]*(Q\[0; 1; 0; 0]))*f(x(1), x(2), x(3));
        phi3 = @(x) ([1, x(1), x(2), x(3)]*(Q\[0; 0; 1; 0]))*f(x(1), x(2), x(3));
        phi4 = @(x) ([1, x(1), x(2), x(3)]*(Q\[0; 0; 0; 1]))*f(x(1), x(2), x(3));
        
        
        % Getting values:
        val1 = quadrature3d(P(1,:), P(2,:), P(3,:), P(4,:), 5, phi1);
        val2 = quadrature3d(P(1,:), P(2,:), P(3,:), P(4,:), 5, phi2);
        val3 = quadrature3d(P(1,:), P(2,:), P(3,:), P(4,:), 5, phi3);
        val4 = quadrature3d(P(1,:), P(2,:), P(3,:), P(4,:), 5, phi4);
        
        % Putting in right place:
        b(nodes) = b(nodes) + [val1; val2; val3; val4];
        
    end
    
    %% Get A without boundary points
    AllboundaryPoints = [edge(:,1); edge(:,2); edge(:,3)];
    boundaryPoints = unique(AllboundaryPoints);
    
    % Setting rows and cols of boundaryPoints equal to 0
    A(boundaryPoints, :) = 0;
    b(boundaryPoints) = 0;
    A(boundaryPoints, boundaryPoints) = speye(length(boundaryPoints), length(boundaryPoints));
    
    % Solving the linear system
    u_sol = A\b;
    
    
    
    % Finding reference solution in points:
    u_ref = zeros(N,1);
    for i = 1:length(p)
        point = p(i,:);
        u_ref(i) = u(point);
    end
    
    error(n) = norm(u_ref-u_sol,2);
end


figure
loglog(error);

% figure
% plot(u_ref, '*-black')
% hold on
% plot(u_sol, '*-b')
% title('Values in points')


