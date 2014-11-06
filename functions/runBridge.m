clear all
close all


addpath ../minecraft

%% loading
% load bridge
addpath ../minecraft/bridge
load elements_bridge.m
load nodes_bridge.m
%[pts_bridge, el_bridge] = getLargestConnectedDomain(nodes_bridge, elements_bridge);
% ditching this because slow
pts_bridge = nodes_bridge;
el_bridge = [elements_bridge(:,1:8) + 1, elements_bridge(:,9:end)]; % fix 0-indexing


%load car:
addpath ../minecraft/car
load elements_car.m
load nodes_car.m
%[pts_car, el_car] = getLargestConnectedDomain(nodes_car, elements_car);
pts_car = nodes_car;
el_car = [elements_car(:,1:8) + 1, elements_car(:,9:end)]; % fix 0-indexing


%% constants
E = 29e6;
v = 0.2;
loadVector = @(x) 7750*(x == 1) + (x ~= 1)*100000;  % density of steel if blockid is 1, heavy else (car)

% min stepsize is 0.5 because of the placement of the points
ts = 25;

% loop over t:
N = length(pts_car) + length(pts_bridge);
M = length(ts);

%% matrixes to save all data to from paralell computing
u_sol2 = zeros(3*N,M);
stress2 = zeros(N,M);
pts2 = zeros(N,3,M);
data2 = zeros(26592,6,M);
%% running parallell loop
for i = 1:M
    i
    %t = 25;
    
    
    
    [ pts, newelements ] = mergeBridgeCar( pts_bridge, el_bridge, pts_car, el_car, ts(i) );
    data = hex2tetr(newelements);

    % removing unused elements
    [pts, data] = removeUnused(pts, data);
    
    % Finding boundary points:
    newtonBoundary = find((pts(:,3) == 0)); % Dirichlet homogenous BC at z = 0;
    u_sol = FEM( pts, data, E, v, loadVector, newtonBoundary);
    stress = stressRecovery( pts, data(:,1:4), E, v, u_sol );
    
    n = length(stress);
    % Saving to use later
    stress2(1:n,i) = stress;
    u_sol2(1:3*n,i) = u_sol;
    pts2(1:n,:,i) = pts;
    data2(1:length(data),:,i) = data;
    
end

%% plotting loop:
fig = figure;
%set(fig, 'Visible','off')
colormap(jet)
for i = 1:M
    U = [u_sol2(1:3:end,i), u_sol2(2:3:end,i), u_sol2(3:3:end,i)];
    clf
    
    trisurf(data2(:,1:4,i),pts2(:,1,i)+U(:,1),pts2(:,2,i)+U(:,2),pts2(:,3,i)+U(:,3),stress2(:,i));
    axis equal, colorbar%title(['FEM solution at time = ' num2str(ts(i))])
    view(-65 - 2*ts(i), 24)
    shading interp
    set(gca,'clim',[0 5e5])  % setting color gradient
    axis off
    saveas(fig, sprintf('../figures/time%03d.png',10*ts(i)));
    
    % to save as gif, go to terminal in figure folder and:
    % convert -delay 10 -loop 0 *.png myimage.gif
    
end