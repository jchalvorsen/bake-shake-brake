close all
clear all

% We want to solve the linear elasticity problem
% grad(o(u)) = -f
% with bondaries on x = -1, 1 or y = -1,1 and dirichlet homogenous BCs

addpath include

n_it = 3;
error = zeros(n_it,1);
n_vec = zeros(n_it,1);

for k = 1:n_it
    
    n = 5*2^(k)
    
    [p, tri, edge] = getPlate(n); % n: number of nodes in each spatial direction
    N = length(p);
    
    % Declaring material constants
    E = 1;
    v = 0.1;
    
    % Building the C matrix:
    C = E/(1-v^2)*[1, v, 0; v, 1, 0; 0, 0, (1-v)/2];
    
    % Declaring functions
    fx = @(x,y) E/(1-v^2) * (-2*y^2 - x^2 + v*x^2 - 2*v*x*y -2*x*y + 3 - v);
    fy = @(x,y) E/(1-v^2) * (-2*x^2 - y^2 + v*y^2 - 2*v*x*y -2*x*y + 3 - v);
    u = @(x) (x(1)^2-1)*(x(2)^2-1); % The same in both directions
    
    
    A = sparse(2*N,2*N);
    b = zeros(2*N,1);
    
    for i = 1:length(tri)
        nodes = tri(i,:);
        P = p(nodes,:);       % Active points in triangle
        
        % Calculating area
        Q = [[1;1;1], P];
        area = 0.5*abs(det(Q));
        
        
        %% Getting stiffness matrix
        % Finding constants in phi (basis function = [1, x, y] * c)
        % c1 = Q\[1; 0; 0];
        % c2 = Q\[0; 1; 0];
        % c3 = Q\[0; 0; 1];
        % c = [c1, c2, c3];
        c = inv(Q);
        
        % looping over x and y indexing
        Ak = zeros(6);
        for q = 1:3
            for w = 1:3
                % Basis function x = [c1 + c2*x + c3*y; 0] and y = [0; c1 + c2*x + c3*y]
                %e = [c2; 0; 0.5*c3];
                e1x = [c(2,q); 0; c(3,q)];
                e2x = [c(2,w); 0; c(3,w)];
                
                %ey = [0; c(3); 0.5*c(2)];
                e1y = [0; c(3,q); c(2,q)];
                e2y = [0; c(3,w); c(2,w)];
                
                
                f11 = e1x'*C*e2x;
                f12 = e1x'*C*e2y;
                f21 = e1y'*C*e2x;
                f22 = e1y'*C*e2y;
                
                Ak(2*q-1,2*w-1) = f11*area;
                Ak(2*q-1,2*w) = f12*area;
                Ak(2*q  ,2*w-1) = f21*area;
                Ak(2*q  ,2*w) = f22*area;
            end
        end
        %Ak
        % Put element matrix in right place
        % Map to right place:
        map(1:2:2*length(nodes)) = 2*nodes-1;
        map(2:2:2*length(nodes)) = 2*nodes;
        A(map,map) = A(map,map) + Ak;
        
        %% Getting b vector
        % Finding functions phi*f:
        for j = 0:1
            if j == 0
                f1 = @(x) ([1, x(1), x(2)]*c(:,1))*fx(x(1),x(2));
                f2 = @(x) ([1, x(1), x(2)]*c(:,2))*fx(x(1),x(2));
                f3 = @(x) ([1, x(1), x(2)]*c(:,3))*fx(x(1),x(2));
            else
                f1 = @(x) ([1, x(1), x(2)]*c(:,3))*fy(x(1),x(2));
                f2 = @(x) ([1, x(1), x(2)]*c(:,2))*fy(x(1),x(2));
                f3 = @(x) ([1, x(1), x(2)]*c(:,1))*fy(x(1),x(2));
            end
            % Getting values:
            val1 = quadrature2d(P(1,:), P(2,:), P(3,:), 1, f1);
            val2 = quadrature2d(P(1,:), P(2,:), P(3,:), 1, f2);
            val3 = quadrature2d(P(1,:), P(2,:), P(3,:), 1, f3);
            
            % Putting in right place:
            b(2*nodes-1+j) = b(2*nodes-1+j) + [val1; val2; val3];
        end
    end
    %% Get A without boundary points
    boundaryPoints = [2*edge'-1; 2*edge'];
    
    % Setting cols of boundaryPoints equal to 0
    A(boundaryPoints, :) = 0;
    b(boundaryPoints) = 0;
    A(boundaryPoints, boundaryPoints) = A(boundaryPoints, boundaryPoints) + speye(length(boundaryPoints), length(boundaryPoints));
    
    % Solving the linear system
    u_sol = A\b;
    
    
    % Finding reference solution in points:
    u_ref = zeros(2*N,1);
    for i = 1:length(p)
        point = p(i,:);
        u_ref(2*i) = u(point);
        u_ref(2*i-1) = u(point);
    end
    
    error(k) = norm(u_ref-u_sol,Inf);
    n_vec(k) = 1/(n);
    
end

%% Convergence plot:
figure
loglog(n_vec, error, '*-r');
title('Loglogplot of error in 2D');
hold on
loglog(n_vec, n_vec, 'b');
loglog(n_vec, n_vec.^2, 'g');
legend('Error', '1', '2');

U_sol = sqrt(u_sol(1:2:end).^2 + u_sol(2:2:end).^2);
U_ref = sqrt(u_ref(1:2:end).^2 + u_ref(2:2:end).^2);

figure
trisurf(tri, p(:,1), p(:,2), U_sol)
figure
trisurf(tri, p(:,1), p(:,2), U_ref)
figure
trisurf(tri, p(:,1), p(:,2), U_sol - U_ref)

