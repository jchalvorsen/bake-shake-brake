function [ newpoints, newelements ] = mergeBridgeCar( pts_bridge, el_bridge, pts_car, el_car, z )
%MERGEBRIDGECAR merges a brigde and car pts and elements
% Returns the new points and elements. Runs extremely fast


% skewing bridge to start at origo:
m = min(pts_bridge);
pts_b = [pts_bridge(:,1) - m(1), pts_bridge(:,2) - m(2), pts_bridge(:,3) - m(3)];

% Setting position of car:
z_dir = min(pts_car(:,3)) - z;
pts_car_new = [pts_car(:,1) - m(1), pts_car(:,2) - m(2), pts_car(:,3) - z_dir];

% intersect returns indexes IA from pts_b and IB from pts_car_new of
% intersecting lines (common lines)
[~,IA,IB] = intersect(pts_b, pts_car_new, 'rows');

% Adding a scalar to indexes of el_car to be able to merge el_car and
% el_bridge later. Also preserves block metadata and id
el_car2 = [el_car(:,1:8) + length(pts_b), el_car(:,9:end)];

% mapping back everything:
for i = 1:length(IB)
    % swap all points in el_car2 with intersecting points in el_bridge
    el_car2(el_car(:,1:8) == IB(i)) = IA(i);
end

% making new data with both car and bridge
newpoints = [pts_b; pts_car_new];
newelements = [el_bridge; el_car2];

 % fixing minecraft y-z axis
newpoints = [newpoints(:,1), newpoints(:,3), newpoints(:,2)];
end

