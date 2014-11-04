function [ u_sol ] = FEM( p, data, E, v, loadVector)
%FEM solve stuff vie the finite element method
%   Detailed explanation goes here

N = length(p);

% Building the C matrix:

C1 =  E*v/((1+v)*(1-2*v))*ones(3,3) + E/(1+v)*eye(3);
C2 = E/(2*(1+v))*eye(3);
C = [ C1        , zeros(3,3);
    zeros(3,3), C2        ]; % pretty sure this should be inverted
% found at www.rpi.edu/~des/3DElasticity.ppt, slide 24 (inverse function)

%density = @(x) 7750 + (x(3) > 8)*500000; % kg/ m^3

A = zeros(3*N,3*N);
b = zeros(3*N,1);
h = waitbar(0, 'In progress');
for i = 1:length(data)
    if mod(i,1000) == 0
        waitbar(i/length(data),h, 'In progress')
    end
    nodes = data(i,1:4);
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
            
            ex1 = [c(2,q); 0; 0; c(3,q); 0; c(4,q)];
            ex2 = [c(2,w); 0; 0; c(3,w); 0; c(4,w)];
            ey1 = [0; c(3,q); 0; c(2,q); c(4,q); 0];
            ey2 = [0; c(3,w); 0; c(2,w); c(4,w); 0];
            ez1 = [0; 0; c(4,q); 0; c(3,q); c(2,q)];
            ez2 = [0; 0; c(4,w); 0; c(3,w); c(2,w)];
            
            e1 = [ex1, ey1, ez1]';
            e2 = [ex2, ey2, ez2];
            
            fk = e1*C*e2;
            
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
    
    
    % new try without quadratures and function handles:
    midpoint = 1/4*ones(1,4)*P;
    
    val = [1, midpoint]*c*loadVector(data(i,5))*vol;
    
    
    % Putting in right place:
    % only want to add gravity compononent to z-dir
    b(3*nodes) = b(3*nodes) + val'*vol*-9.81; % not sure about constants
    
    
end
close(h)

%% Get A without boundary points
boundaryPoints = find((p(:,3) == 0)); % Dirichlet homogenous BC: f(boundaryPoints) = 0

map2(1:3:3*length(boundaryPoints)) = 3*boundaryPoints-2;
map2(2:3:3*length(boundaryPoints)) = 3*boundaryPoints-1;
map2(3:3:3*length(boundaryPoints)) = 3*boundaryPoints;

% Setting cols of boundaryPoints equal to 0
A(map2, :) = 0;
%A(:, boundaryPoints) = 0;
b(map2) = 0;
A(map2, map2) = A(map2, map2) + speye(length(map2), length(map2));


% Removing extra elements because of double nodes:
% finding rows equal to zero:
z = zeros(1,length(A));
for i = 1:length(A)
    if isequal(A(i,:),z)
        A(i,i) = 1;
    end    
end

% extraNodes = find(mapping ~=0) + length(pts);
% map3(1:3:3*length(extraNodes)) = 3*extraNodes-2;
% map3(2:3:3*length(extraNodes)) = 3*extraNodes-1;
% map3(3:3:3*length(extraNodes)) = 3*extraNodes;
% A(map3, map3) = eye(length(map3));


% Making A sparse so linear system will be solved fast
Asp = sparse(A);
% Solving the linear system
u_sol = Asp\b;

%     % Plotting inside loop:
%     U1 = [u_sol(1:3:end,1), u_sol(2:3:end,1), u_sol(3:3:end,1)];
%     clf
%     trisurf(tetr,p(:,1)+U1(:,1),p(:,2)+U1(:,2),p(:,3)+U1(:,3),U1(:,3));
%     axis equal,colorbar,title(['FEM solution at time = ' num2str(ts(t))])
%     view(-65 - 2*ts(t), 24);
%     saveas(fig, ['figures/time' num2str(ts(t)) '.png']);
%
%



end

