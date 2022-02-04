function [A, b] = getSysMat(u_past, kappa, tau, h, s_j, s_jm1, ...
    ds_jdt, ds_jm1dt, alpha, g_jm1, g_j)
%GETSYSTEMMATRIX Вычисление матрицы системы для решения задачи
%Стефана методом конечных разностей с применением метода выпрямления
%фронта. Граничные условия учитываются.
%   u_past - распределение температур на предыдущем временном шаге

h_sq = h*h;
inv_h_sq = 1/h_sq;
C1 = kappa*(s_j - s_jm1)^2/tau;
Np = length(u_past);
C2_1 = kappa*(s_j - s_jm1)/(2*h);
C2_2 = (ds_jdt - ds_jm1dt);

e = ones(Np, 1);
x = (0:Np-1)'*h;
C2 = C2_1*( ds_jm1dt + x*C2_2 );
A = sparse(Np);
A = spdiags([ inv_h_sq - [C2(2:end); 0; 0],...
             -[e;0]*(2/h_sq + C1),...
              inv_h_sq + [0;C2]], [-1 0 1], Np, Np);
b = zeros(Np, 1);
b(2:Np-1) = -C1*u_past(2:Np-1);

C3 = alpha(1, 2)/(2*h);
A(1, 1) = alpha(1, 1)*(s_j - s_jm1) - 3*C3;
A(1, 2) = 4*C3;
A(1, 3) = -C3;
b(1) = g_jm1*(s_j - s_jm1);

C4 = alpha(2, 2)/(2*h);
A(Np, Np - 2) = C4;
A(Np, Np - 1) = -4*C4;
A(Np, Np) = 3*C4 + alpha(2, 1)*(s_j - s_jm1);
b(Np) = g_j*(s_j - s_jm1);
end

