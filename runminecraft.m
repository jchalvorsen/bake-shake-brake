close all
clear all
addpath minecraft/house

%load pts.m
%load tetr.m
%writeVTF(pts, tetr(:,1:4), tetr(:,5), 'testHouse.vtf')

addpath minecraft/random
addpath minecraft
load elements.m
load nodes.m
[pts el] = getLargestConnectedDomain(nodes, elements);
tetr = hex2tetr(el);

writeVTF(pts, tetr(:,1:4), tetr(:,5), 'testRandom.vtf')

% We then have tetr and pts and are ready to do the analysis