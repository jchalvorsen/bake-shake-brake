close all
clear all

addpath include

N = 900;

[p tri edge] = getSphere(N);

trep = TriRep(tri, p);
[tr, Xb] = freeBoundary(trep);
trisurf(tr, Xb(:,1), Xb(:,2), Xb(:,3), 'FaceColor', 'red','FaceAlpha', 0.8); 
 

% Declaring functions
f = @(x,y,z) -12*pi*cos(2*pi*(x^2+y^2+z^2)) +16*pi^2*(x^2+y^2+z^2)*sin(2*pi*(x^2+y^2+z^2));
u = @(x) sin(2*pi*(x(1)^2+x(2)^2));

A = sparse(N,N);
b = zeros(N,1);

divpsi1 = [-1;-1];
divpsi2 = [1; 0];
divpsi3 = [0; 1];
divpsi = [divpsi1, divpsi2, divpsi3];
