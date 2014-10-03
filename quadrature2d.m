function I = quadrature2D(p1 ,p2 ,p3 ,Nq ,g )
%QUADRATURE2D Summary of this function goes here
%   Detailed explanation goes here

assert(Nq ~= 2)

C1 = 0;
C2 = 0;
C3 = 0;
rho = 0;
switch (Nq)
    case 1
        C1 = 1/3;
        C2 = 1/3;
        C3 = 1/3;
        rho = [1];

    case 3
        C1 = [0.5, 0.5, 0  ];
        C2 = [0.5, 0  , 0.5];
        C3 = [0  , 0.5, 0.5];
        rho = [1/3, 1/3, 1/3 ];
        
    case 4
        C1 = [1/3, 3/5, 1/5, 1/5];
        C2 = [1/3, 1/5, 3/5, 1/5];
        C3 = [1/3, 1/5, 1/5, 3/5];
        rho = [-9/16, 25/48*[1, 1, 1 ]];
end
x = zeros(Nq,length(p1));
partial_sums = zeros(Nq,1);
for i = 1:Nq
    x(i,:) = C1(i).*p1 + C2(i).*p2 + C3(i).*p3;
    partial_sums(i) = rho(i)*g(x(i,:));
end
I = sum(partial_sums);
