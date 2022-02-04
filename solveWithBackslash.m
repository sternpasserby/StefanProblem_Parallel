function x = solveWithBackslash(A, b)
%SOLVEWITHTHOMAS Summary of this function goes here
%   Detailed explanation goes here

N = length(b);
x = zeros(N, 1);

% Получения нуля на месте A(2, 1) и A(N-1, N)
C1 = -A(2, 1)/A(1, 1);
A(2, :) = A(2, :) + A(1, :)*C1;
b(2) = b(2) + b(1)*C1;
C2 = -A(N - 1, N) / A(N, N);
A(N-1, :) = A(N-1, :) + A(N, :)*C2;
b(N-1) = b(N-1) + b(N)*C2;

% Решение системы функцией матлаба и явным вычислением x(1) и x(N)
x(2:N-1) = A(2:N-1, 2:N-1) \ b(2:N-1);
x(1) = (b(1) - A(1, 2)*x(2) - A(1, 3)*x(3))/A(1, 1);
x(N) = ( b(N) - A(N, N-1)*x(N-1) - A(N, N-2)*x(N-2))/A(N, N);

end

