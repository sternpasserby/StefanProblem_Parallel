function [ld, md, ud, b] = getSysMat(u_past, kappa, tau, ksi, s_j, s_jm1, ...
    ds_jdt, ds_jm1dt, alpha, g_jm1, g_j)
%GETSYSTEMMATRIX Вычисление матрицы системы для решения задачи
%Стефана методом конечных разностей с применением метода выпрямления
%фронта. Граничные условия учитываются.
%   u_past - распределение температур на предыдущем временном шаге

% Константы
h = ksi(2:end)-ksi(1:end-1);
Np = length(u_past);
C1 = kappa*(s_j - s_jm1)^2;
Ksi = kappa*(s_j - s_jm1)*( (1-ksi)*ds_jm1dt + ksi*ds_jdt );

% Формирование основной части матрицы А
ld = [(Ksi(2:Np-1).*h(1:Np-2) - 2)./( h(1:Np-2).*(h(1:Np-2) + h(2:Np-1)) ) 0]';
md = [ 0 C1/tau + 2./( h(1:Np-2).*h(2:Np-1) ) 0]';
ud = [ 0 -( Ksi(2:Np-1).*h(2:Np-1) + 2)./( h(2:Np-1).*(h(2:Np-1) + h(1:Np-2)) ) ]';
%A = spdiags([ ld md ud ], [-1 0 1], Np, Np);
b = C1/tau*u_past;

% Левое г.у.
C3 = alpha(1, 2)/(s_j - s_jm1);
A_11 = alpha(1, 1) - C3*(2*h(1) + h(2)) / ( h(1)*( h(1) + h(2) ) );
A_12 = C3*( h(1)+h(2) ) / ( h(1)*h(2) );
A_13 = -C3*h(1) / ( h(2)*( h(1) + h(2) ) );
b_1 = g_jm1;
temp = -A_13/ud(2);
md(1) = A_11 + ld(1)*temp;
ud(1) = A_12 + md(2)*temp;
b(1) = b_1 + b(2)*temp;

% Правое г.у.
C5 = alpha(2, 2) / (s_j - s_jm1);
A_Np_Npm2 = C5*h(Np-1) / ( h(Np-2)*( h(Np-2)+h(Np-1) ) );
A_Np_Npm1 = -C5*( h(Np-2) + h(Np-1) ) / ( h(Np-1)*h(Np-2) );
A_Np_Np = alpha(2, 1) + C5*( h(Np-2) + 2*h(Np-1) ) / ( h(Np-1)*( h(Np-2) + h(Np-1) ) );
b_Np = g_j;
temp = -A_Np_Npm2/ld(Np-2);
ld(Np-1) = A_Np_Npm1 + md(Np-1)*temp;
md(Np) = A_Np_Np + ud(Np-1)*temp;
b(Np) = b_Np + b(Np-1)*temp;

% A = sparse(Np);
% b = zeros(Np, 1);
% 
% C2 = alpha(1, 1)*(s_j-s_jm1);
% C3 = alpha(1, 2);
% A(1, 1) = C2 - C3*( 2*h(1)+h(2) )/( h(1)*(h(1)+h(2)) );
% A(1, 2) = C3*( h(1)+h(2) )/( h(1)*h(2) );
% A(1, 3) = -C3*h(1)/( h(2)*(h(1)+h(2)) );
% b(1) = g_jm1*(s_j-s_jm1);
% for i = 2:Np-1
%     C1 = kappa*(s_j - s_jm1)^2;
%     Ksi_i = kappa*(s_j - s_jm1)*( (1-ksi(i))*ds_jm1dt + ksi(i)*ds_jdt );
%     A(i, i - 1) = (Ksi_i*h(i-1) - 2) / ( h(i-1)*(h(i)+h(i-1)) );
%     A(i, i) = C1/tau + 2 / (h(i)*h(i-1));
%     A(i, i + 1) = - ( Ksi_i*h(i)+2 ) / ( h(i)*( h(i) + h(i-1) ) );
%     b(i) = C1*u_past(i)/tau;
% end
% C4 = alpha(2,1)*( s_j-s_jm1 );
% C5 = alpha(2, 2);
% A(Np, Np - 2) = C5*h(Np-1)/( h(Np-2)*(h(Np-2)+h(Np-1)) );
% A(Np, Np - 1) = -C5*( h(Np-2) + h(Np-1) ) / ( h(Np-1)*h(Np-2) );
% A(Np, Np) = C4 + C5*( h(Np-2) + 2*h(Np-1) ) / ( h(Np-1)*(h(Np-2)+h(Np-1)) );
% b(Np) = g_j*( s_j-s_jm1 );
end

