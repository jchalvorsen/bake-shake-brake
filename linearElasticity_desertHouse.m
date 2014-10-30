close all
clear all
% We want to solve the linear elasticity problem
% grad(o(u)) = -f
% on a 3D surface thingy
% where f is the weight of each element (force that acts upon it)
% and dirichlet homogenous BCs where structure is attached to ground

addpath include/desertHouse

% loading blocked elements
data = load('tetr.m');
tetr = data(:,1:4);
points = load('pts.m');
p = [points(:,1), points(:,3), points(:,2)];
block_id = data(:,5);


% Write structure to vtf file
%writeVTF(p, tetr, 0, 'queen.vtf')

% Can plot, but takes a lot of resources
%figure
%tetramesh(tetr, p)


N = length(p);
id = unique(block_id);

% {grass, plank, gravel, dispenser, moss stone, stone brick}
% keySet =   {2, 5, 13, 23, 48, 99}; 
% density = [300, 700, 2400, 1600, 2740, 1920];
% Young = 10e9 * [10e-9, 11, 10e-9, 17, 17, 20];
% rho_map = containers.Map(keySet,density);
% E_map = containers.Map(keySet,weight);

    
    % Declaring material constants
    E = 10000;
    v = 0.3;
    density = 300;
    

 % pretty sure this should be inverted
  % found at www.rpi.edu/~des/3DElasticity.ppt, slide 24 (inverse function)


A = spalloc(3*N,3*N, 12*12*length(tetr));
b = zeros(3*N,1);
h = waitbar(0, 'In progress');
for i = 1:length(tetr)
    if mod(i,1000) == 0
        waitbar(i/length(tetr),h, 'In progress')
    end
    

    % Building the C matrix:
    
    C1 =  E*v/((1+v)*(1-2*v))*ones(3,3) + E/(1+v)*eye(3);
    C2 = E/(2*(1+v))*eye(3);
    C = inv([ C1        , zeros(3,3);
        zeros(3,3), C2        ]);
    
    nodes = tetr(i,:);
    P = p(nodes,:);       % Active points in tetrahedron
    
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
            % e = [c2; 0; 0.5*c3];
            
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
            
            
            ex1 = [c(2,q), 0, 0, c(3,q), 0, c(4,q)];
            ex2 = [c(2,w); 0; 0; c(3,w); 0; c(4,w)];
            ey1 = [0, c(3,q), 0, c(2,q), c(4,q), 0];
            ey2 = [0; c(3,w); 0; c(2,w); c(4,w); 0];
            ez1 = [0, 0, c(4,q), 0, c(3,q), c(2,q)];
            ez2 = [0; 0; c(4,w); 0; c(3,w); c(2,w)];
            
            e1 = [ex1; ey1; ez1];
            e2 = [ex2, ey2, ez2];
            
            fk = e1*C*e2;
%             
%           % q represents basis function #1, while w is #2         

              
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
    
    midpoint = 1/4*ones(1,4)*P;
    val = [1, midpoint]*c*density*vol;
  
    % Putting in right place:
    % only want to add gravity compononen to z-dir
    b(3*nodes) = b(3*nodes) + val'*vol*-9.81;
    

end
close(h)
%% Get A without boundary points
min_z = min(p(:,3));
boundaryPoints = find((p(:,3) == min_z)); % Dirichlet homogenous BC: f(boundaryPoints) = 0

map2(1:3:3*length(boundaryPoints)) = 3*boundaryPoints-2;
map2(2:3:3*length(boundaryPoints)) = 3*boundaryPoints-1;
map2(3:3:3*length(boundaryPoints)) = 3*boundaryPoints;

% Setting cols of boundaryPoints equal to 0
A(map2,:)  = 0;
b(map2) = 0;
A(map2, map2) = A(map2, map2) + speye(length(map2), length(map2));

% Solving the linear system
u_sol = A\b;

%% Plotting:

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


figure
%subplot(2,1,1)
trisurf(tetr,p(:,1)+U(:,1),p(:,2)+U(:,2),p(:,3)+U(:,3),u_sol(1:3:end));
view(2),axis equal,colorbar,title('FEM solution')

% Export to glview
writeVTF(p, tetr, u_sol(1:3:end), 'desertHouse.vtf')
writeVTF2(p, tetr, 'Displacement', U, 'FileName', 'desertHouse_displacement.vtf');
disp('hei')