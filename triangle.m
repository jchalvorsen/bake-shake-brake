% Check if point (x,y) is in triangle defined by P

% All combinations of gridpoints
z = linspace(-1,1);
x = repmat(z,100,1);
y = repmat(z',1,100);

% Ptri defines the closed triangle (start to end point)
Ptri = [P;P(1,:)];

% Checks if (x(i),y(i)) is inside triangle for i = 1:100
[in on] = inpolygon(x,y,Ptri(:,1),Ptri(:,2));

% Only interior points:
% in = in+on includes both interior and boundary points

% Plot triangle and interior points
hold on
plot(Ptri(:,1),Ptri(:,2));
plot(x(in),y(in),'rx')
plot(x(on),y(on),'gx')