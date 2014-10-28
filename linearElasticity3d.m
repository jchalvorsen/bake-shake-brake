close all
clear all
% We want to solve the linear elasticity problem
% grad(o(u)) = -f
% on a 3D surface thingy
% where f is the weight of each element (force that acts upon it)
% and dirichlet homogenous BCs where structure is attached to ground

addpath include/queen

% loading blocked elements
edges = load('coarse_bndry.m');     % tetraeders of the surface
tetr = load('coarse_element.m');
p = load('coarse_nodes.m');


zeronodes = find((p(:,3) == 0)); % Dirichlet homogenous BC: f(zeronodes) = 0

% Write structure to vtf file
%writeVTF(p, tetr, 0, 'queen.vtf')

% Can plot, but takes a lot of resources
%figure
%tetramesh(tetr, p)






