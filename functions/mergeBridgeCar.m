function [ newpoints, newelements ] = mergeBridgeCar( pts_bridge, el_bridge, pts_car, el_car, y )
%MERGEBRIDGECAR merges a brigde and car pts and elements
%   This is a mess, but it works at least.

% skewing bridge to start at origo:
m = min(pts_bridge);
pts_b = [pts_bridge(:,1) - m(1), pts_bridge(:,2) - m(2), pts_bridge(:,3) - m(3)];

% Setting position of car:
z_dir = min(pts_car(:,3)) - y;
pts_car_new = [pts_car(:,1) - m(1), pts_car(:,2) - m(2), pts_car(:,3) - z_dir];


% combining car and bridge mappings:
mapping = zeros(length(pts_car_new),1);
for i = 1:length(pts_car_new)
    for j = 1:length(pts_b)
        if sum(pts_car_new(i,:) == pts_b(j,:)) == 3
            mapping(i) = j;
        end
    end
end

el_car2 = [el_car(:,1:8) + length(pts_b), el_car(:,9:end)];  % fixing indexing for for the new points, preserving block metadata and id

% mapping back everything:
for i = 1:length(mapping)
    if mapping(i) > 0
        el_car2(el_car(:,1:8) == i) = mapping(i);
    end
end

% making new data with both car and bridge
newpoints = [pts_b; pts_car_new];
newpoints = [newpoints(:,1), newpoints(:,3), newpoints(:,2)]; % fixing minecraft y-z axis
newelements = [el_bridge; el_car2];

end

