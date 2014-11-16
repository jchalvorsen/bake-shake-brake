function I = quadrature3d(p1 ,p2 ,p3, p4, Nq ,g )
%QUADRATURE2D Summary of this function goes here
%   Detailed explanation goes here

assert(Nq ~= 2)
assert(Nq ~= 3)

a = 1/4 + 3*sqrt(5)/20;
b = 1/4 - sqrt(5)/20;
chi = 0;
rho = 0;
switch (Nq)
    case 1
        chi = 1/4*ones(1,4);
        rho = 1;

    case 4
        chi = [ a, b, b, b;
                b, a, b, b;
                b, b, a, b;
                b, b, b, a];                
        rho = 1/4*ones(4,1);
        
    case 5
        chi = [1/4, 1/4, 1/4, 1/4;
               1/2, 1/6, 1/6, 1/6;
               1/6, 1/2, 1/6, 1/6;
               1/6, 1/6, 1/2, 1/6;
               1/6, 1/6, 1/6, 1/2];
        rho = [-4/5, + 9/20*[1 1 1 1]]';
end
% Time in volume to answer
vol = abs(1/6*det([p1'-p4', p2'-p4', p3'-p4']));
x = chi*[p1; p2; p3; p4];
y = zeros(3, Nq);
for i = 1:Nq   
    y(:,i) = rho(i)*g(x(i,:));
end
I = vol*sum(y,2);
