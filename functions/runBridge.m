clear all
close all


addpath ../minecraft

% load bridge
addpath ../minecraft/bridge
load elements_bridge.m
load nodes_bridge.m
[pts_bridge, el_bridge] = getLargestConnectedDomain(nodes_bridge, elements_bridge);

%load car:
addpath ../minecraft/car
load elements_car.m
load nodes_car.m
[pts_car, el_car] = getLargestConnectedDomain(nodes_car, elements_car);

%% constants
E = 29e6;
v = 0.2;
loadVector = @(x) 7750 + (x(3) > 8)*500000;

ts = 0:5:25;

fig = figure;
% loop over t:
N = length(pts_car) + length(pts_bridge);
M = length(ts);
u_sol = zeros(3*N,M);
stress = zeros(N,M);
for i = 1:M
    %t = 25;
    [ pts, newelements ] = mergeBridgeCar( pts_bridge, el_bridge, pts_car, el_car, ts(i) );
    data = hex2tetr(newelements);
    
    u_sol(:,i) = FEM( pts, data, E, v, loadVector);
    %%
    stress(:,i) = stressRecovery( pts, data(:,1:4), E, v, u_sol(:,i) );
    
    % Plotting inside loop:
    U = [u_sol(1:3:end,i), u_sol(2:3:end,i), u_sol(3:3:end,i)];
    clf
    trisurf(data(:,1:4),pts(:,1)+U(:,1),pts(:,2)+U(:,2),pts(:,3)+U(:,3),stress(:,i));
    axis equal, colorbar,title(['FEM solution at time = ' num2str(ts(i))])
    view(-65 - 2*ts(i), 24);
    saveas(fig, ['../figures/time' num2str(ts(i)) '.png']);
    
    
end

