function [ normedNodeStress ] = stressRecovery( p, tetr, E, v, u_sol)
%STRESSRECOVERY Recovers Von Mises stress from a set of displacements
%   Returns Von Mises stress

N = length(p);
C1 =  E*v/((1+v)*(1-2*v))*ones(3,3) + E/(1+v)*eye(3);
C2 = E/(2*(1+v))*eye(3);
C = [ C1        , zeros(3,3);
    zeros(3,3), C2        ];

normedElementStress = zeros(length(tetr),1);
for i = 1:length(tetr)
    nodes = tetr(i,:);
    P = p(nodes,:);       % Active points in tetrangle
    
    % Calculating area
    Q = [[1;1;1;1], P]; 
    
    % Getting stiffness matetrx
    % Finding constants in phi (basis function = [1, x, y, z] * c)
    % c1 = Q\[1; 0; 0; 0];
    % c2 = Q\[0; 1; 0; 0];
    % c3 = Q\[0; 0; 1; 0];
    % c4 = Q\[0; 0; 0; 1];
    % c = [c1, c2, c3, c4];
    c = inv(Q);
    
    % B: strain displacement matrix
    phiGrad = c(2:end,:);
    B = zeros(6,12);
    B([1,4,5],1:3:10) = phiGrad;
    B([4,2,6],2:3:11) = phiGrad;
    B([5,6,3],3:3:12) = phiGrad;
    
    % u_e: displacement field:
    map(1:3:3*length(nodes)) = 3*nodes-2;
    map(2:3:3*length(nodes)) = 3*nodes-1;
    map(3:3:3*length(nodes)) = 3*nodes;
    
    u_e = u_sol(map);
    
    % s: stress vector for element:  
    s = C*B*u_e;
    
    
    % calculating the von mises stress:
    stress_sqaured = 0.5*((s(1)-s(2))^2 + (s(2)-s(3))^2 + (s(3)-s(1))^2) ...
                   + 3*(s(4)^2 + s(5)^2 + s(6)^2);
    
    normedElementStress(i) = normedElementStress(i) + sqrt(stress_sqaured);   
end

%% Averaging out stresses for each node:
normedNodeStress = zeros(N,1);
for i = 1:N
    [lines, ~] = find(tetr == i);
    normedNodeStress(i) = mean(normedElementStress(lines));
end
end

