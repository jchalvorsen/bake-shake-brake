function I = quadrature1D_modified(a ,b ,Nq ,g )
%QUADRATURE1D Summary of this function goes here
%   Detailed explanation goes here


z = 0;
rho = 0;
switch (Nq)
    case 1
        z = [0.5];
        rho = [1];
    case 2
        z = [0.5 - 1/sqrt(3),
            0.5 + 1/sqrt(3) ];
        rho = [0.5,
            0.5];
    case 3
        z = [0.5 - sqrt(15)/10,
            0.5,
            0.5 + sqrt(15)/10];
        rho = [5/18,
            4/9,
            5/18];
    case 4
        z = [0.5 - sqrt(525 + 70*sqrt(30))/70,
            0.5 - sqrt(525 - 70*sqrt(30))/70,
            0.5 + sqrt(525 - 70*sqrt(30))/70,
            0.5 + sqrt(525 + 70*sqrt(30))/70];
        rho = [ (18 - sqrt(30))/72,
            (18 + sqrt(30))/72,
            (18 + sqrt(30))/72,
            (18 - sqrt(30))/72];
end
% This mapping came to me in a dream

% Write line segment (a,b) as vector v(t), t = [0,1]
% Replace x by v(1) and y by v(y)
v = @(t) [a(1) + t*(b(1)-a(1)), a(2) + t*(b(2)-a(2))];

I = norm(b-a,2) * rho.*g(v(z));

end

