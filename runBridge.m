close all
clear all
addpath minecraft/bridge
addpath minecraft

% load bridge
load elements_bridge.m
load nodes_bridge.m
[pts el] = getLargestConnectedDomain(nodes_bridge, elements_bridge);
data = hex2tetr(el);
tetr = data(:,1:4);
p = [pts(:,1), pts(:,3), pts(:,2)];

%normalizing p:
m = min(p);
p = [p(:,1) - m(1), p(:,2) - m(2), p(:,3) - m(3)];


%load car:
addpath minecraft/car
load elements_car.m
load nodes_car.m
[pts_car el_car] = getLargestConnectedDomain(nodes_car, elements_car);
data_car = hex2tetr(el_car);
tetr_car = data_car(:,1:4);
p_car = [pts_car(:,1), pts_car(:,3), pts_car(:,2)];

%normalizing p_car:
p_car = [p_car(:,1) - m(1), p_car(:,2) - m(2), p_car(:,3) - m(3)];



% adding car to bridge:
%p = [p; p_car];
%tetr = [tetr; tetr_car + length(p)];

writeVTF(p, tetr, data(:,5), 'bridge.vtf')

% We then have tetr and pts and are ready to do the analysis


N = length(p);

% Declaring material constants
E = 29e6; % 52 gPa
v = 0.2;

% Building the C matrix:

C1 =  E*v/((1+v)*(1-2*v))*ones(3,3) + E/(1+v)*eye(3);
C2 = E/(2*(1+v))*eye(3);
C = [ C1        , zeros(3,3);
      zeros(3,3), C2        ]; % pretty sure this should be inverted
  % found at www.rpi.edu/~des/3DElasticity.ppt, slide 24 (inverse function)

density = @(x,y,z) 7750; % kg/ m^2

A = zeros(3*N,3*N);
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

    val = [1, midpoint]*c*density(midpoint)*vol;
    
    
    % Putting in right place:
    % only want to add gravity compononent to z-dir
    b(3*nodes) = b(3*nodes) + val'*vol*-9.81/4; % not sure about constants
    

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

% Making A sparse so linear system will be solved fast 
Asp = sparse(A);
% Solving the linear system
u_sol = Asp\b;

figure
%plot(u_ref, '*-black')
hold on
plot(u_sol, '*-b')
title('Values in points')

figure
subplot(1,3,1)
trisurf(tetr,p(:,1),p(:,2),p(:,3),u_sol(1:3:end));
view(2),axis equal,colorbar,title('FEM solution in x dir')
subplot(1,3,2)
trisurf(tetr,p(:,1),p(:,2),p(:,3),u_sol(2:3:end));
view(2),axis equal,colorbar,title('FEM solution in y dir')
subplot(1,3,3)
trisurf(tetr,p(:,1),p(:,2),p(:,3),u_sol(3:3:end));
view(2),axis equal,colorbar,title('FEM solution in z dir')

U = [u_sol(1:3:end), u_sol(2:3:end), u_sol(3:3:end)];


figure
%subplot(2,1,1)
trisurf(tetr,p(:,1)+U(:,1),p(:,2)+U(:,2),p(:,3)+U(:,3),U(:,3));
view(2),axis equal,colorbar,title('FEM solution')

% Export to glview
writeVTF2(p, tetr, 'Displacement', U, 'FileName', 'bridge_displacement.vtf');



% wanting to find stresses:




