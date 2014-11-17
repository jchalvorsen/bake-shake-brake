close all
clear all

% using material constants for steel
E = 29e6;
v = 0.2;

% declaring load and analytical function
Q = @(x) x^2-1;
u = @(x, y, z) Q(x)*Q(y)*Q(z);
loadfunction = @(x, y, z)  E*v/((1+v)*(1-2*v)) ...
    *[ -2*Q(y)*Q(z) + (2*v-1)*Q(x)*(Q(y)+Q(z)) - 2*x*(y*Q(z)+z*Q(y));
    -2*Q(x)*Q(z) + (2*v-1)*Q(y)*(Q(x)+Q(z)) - 2*y*(x*Q(z)+z*Q(x));
    -2*Q(x)*Q(y) + (2*v-1)*Q(z)*(Q(x)+Q(y)) - 2*z*(x*Q(y)+y*Q(x))];
l = @(x) loadfunction(x(1), x(2), x(3));

Ns = 4:2:16;
error = zeros(length(Ns),1);
for ii = 1:length(Ns);
    N = Ns(ii)
    p = zeros(N,3);
    
    % getting the triangulation
    edge = [];
    k = 1;
    for j=1:N,
        for i=1:N
            for z=1:N
                p(k,:) = [(z-1)/(N-1), (i-1)/(N-1), (j-1)/(N-1)] * 2 - 1;
                if max((p(k,:) == -1) +  (p(k,:) == +1)) == 1
                    edge = [edge, k];
                end
                k = k+1;                
            end
        end
    end   
    quad  = delaunay(p(:,1), p(:,2), p(:,3));   
    
    %% solving   
    u_sol = FEM( p, quad, E, v, l, edge);
    
    %% fixing analytical solution and calculating error
    u_e = zeros(3*length(p),1);
    for i = 1:length(p)
        u_e(3*i-2:3*i) = u(p(i,1), p(i,2), p(i,3));
    end
    error(ii) = norm(u_sol-u_e,'inf'); 
end
%%

figure % plot solution
subplot(1,2,1)
trisurf(quad,p(:,1) + u_sol(1:3:end), p(:,2) + u_sol(2:3:end), p(:,3) + u_sol(3:3:end));
view(-0.5, 34)
title('FEM solution')
subplot(1,2,2)
trisurf(quad,p(:,1) + u_e(1:3:end), p(:,2) + u_e(2:3:end), p(:,3) + u_e(3:3:end));
view(-0.5, 34)
title('Analytical solution')

% plot error
figure
loglog(1./Ns, error, '*-r');
title('Loglogplot of error');
hold on
loglog(1./Ns, 1./Ns);
grid on
legend('Error', 'Order 1')
xlabel('step size')
