clear all
close all


addpath ../minecraft

%% loading
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
loadVector = @(x) 7750*(x == 1) + (x~=1)*100000;  % density of steel if blockid is 1, heavy else (car)

% min stepsize is 0.5 because of the placement of the points
ts = 0:10:32;

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
    
    u_sol = FEM( pts, data, E, v, loadVector);
    stress = stressRecovery( pts, data(:,1:4), E, v, u_sol );
    
    
    % Saving to use later
    stress2(:,i) = stress;
    u_sol2(:,i) = u_sol;
    pts2(:,:,i) = pts;
    data2(:,:,i) = data;
    
end

%% plotting loop:
fig = figure;
%set(fig, 'Visible','off')
colormap(jet)
for i = 1:M
    U = [u_sol2(1:3:end,i), u_sol2(2:3:end,i), u_sol2(3:3:end,i)];
    clf
    
    trisurf(data2(:,1:4,i),pts2(:,1,i)+U(:,1),pts2(:,2,i)+U(:,2),pts2(:,3,i)+U(:,3),stress2(:,i));
    axis equal, %colorbar%title(['FEM solution at time = ' num2str(ts(i))])
    view(-65 - 2*ts(i), 24)
    shading interp
    set(gca,'clim',[0 5e5])  % setting color gradient
    axis off
    saveas(fig, sprintf('../figures/time%03d.png',10*ts(i)));
    
    % to save as gif, go to terminal and:
    % convert -delay 20 -loop 0 *.png myimage.gif
    
end