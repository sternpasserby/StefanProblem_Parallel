function [s, t, U, X, T] = StefanProblemSolver(pc, bc, ic, h, tau, tMax, Np_save, tau_save)
%STEFANPROBLEMSOLVER Решатель трехфазной задачи Стефана
%   На вход подаются 3 структуры, число узлов сетки для каждой фазы Np,
%   размер временнОго шага tau и время моделирования tMax. Возвращает 3
%   матрицы U, X, T и вектор s.

%%% Распаковка структур
lambda1 = pc.lambda1;
c1 = pc.c1;
rho1 = pc.rho1;
a1_sq = pc.a1_sq;
lambda2 = pc.lambda2;
c2 = pc.c2;
rho2 = pc.rho2;
a2_sq = pc.a2_sq;
Uf = pc.Uf;
qf = pc.qf;
rho = (rho1 + rho2)/2;
alpha = bc.alpha;
accumRate = ic.accumRate;

% Задание числа узлов на каждую фазу
Np_min = 100;
Np1 = max(Np_min, ceil((ic.s1 - ic.s0)/h) );
Np2 = max(Np_min, ceil((ic.s2 - ic.s1)/h) );
Np3 = max(Np_min, ceil((ic.s3 - ic.s2)/h) );
h1 = 1/(Np1 - 1);
h2 = 1/(Np2 - 1);
h3 = 1/(Np3 - 1);

% Задание характерных параметров для обезразмеривания
x0 = ic.s3 - ic.s0;
%x0 = 1;
t0 = rho1*c1*x0*x0/lambda1;
U0 = Uf;
beta = qf*rho / (rho1*c1*U0);
kappa = c2*rho2/lambda2 * lambda1 / (c1*rho1); 

% Обезразмеривание исходных данных
g0 = @(t)(bc.g0(t*t0)/U0);
g1 = @(t)(bc.g1(t*t0)/U0);
g2 = @(t)(bc.g2(t*t0)/U0);
g3 = @(t)(bc.g3(t*t0)/U0);
g4 = @(t)(bc.g4(t*t0)/U0);
g5 = @(t)(bc.g5(t*t0)/U0);
u1 = ic.u1/U0;
u2 = ic.u2/U0;
u3 = ic.u3/U0;
tau = tau/t0;
tau_save = tau_save/t0;
tMax = tMax/t0;
Uf = Uf/U0;
alpha(:, 2) = bc.alpha(:, 2)/x0;

M = ceil(tMax/tau); % Число шагов по времени
s0 = zeros(1, M + 1);
s1 = zeros(1, M + 1);
s2 = zeros(1, M + 1);
s3 = zeros(1, M + 1);
 t = zeros(1, M + 1);
s0(1) = ic.s0/x0;
s1(1) = ic.s1/x0;
s2(1) = ic.s2/x0;
s3(1) = ic.s3/x0;

nRows = 3*min(min(min(Np1, Np2), Np3), Np_save);
nCols = min(M, round(tMax/tau_save) );
X = zeros(nRows, nCols);
U = zeros(nRows, nCols);
T = zeros(nRows, nCols);

ksi1 = linspace(0, 1, Np1)';
ksi2 = linspace(0, 1, Np2)';
ksi3 = linspace(0, 1, Np3)';
ksi_save = linspace(0, 1, nRows/3)';
A1 = sparse(Np1);
A2 = sparse(Np2);
A3 = sparse(Np3);
b1 = zeros(Np1, 1);
b2 = zeros(Np2, 1);
b3 = zeros(Np3, 1);

isUpperPhase = true;
isLowerPhase = true;
if s0(1) == s1(1)
    isLowerPhase = false;
end
if s3(1) == s2(1)
    isUpperPhase = false;
end

if isLowerPhase
    x1 = linspace(s0(1), s1(1), length(u1));
    x1q = s0(1) + ksi1.*(s1(1) - s0(1));
    u1 = interp1(x1, u1, x1q);
else
    u1 = ksi1.*NaN;
end
x2 = linspace(s1(1), s2(1), length(u2));
x2q = s1(1) + ksi2.*(s2(1) - s1(1));
u2 = interp1(x2, u2, x2q);
if isUpperPhase
    x3 = linspace(s2(1), s3(1), length(u3));
    x3q = s2(1) + ksi3.*(s3(1) - s2(1));
    u3 = interp1(x3, u3, x3q);
else
    u3 = ksi3.*NaN;
end

% Запись начальных условий в выходные массивы
if isLowerPhase
    x1 = s0(1) + ksi1.*(s1(1) - s0(1));
    x1q = s0(1) + ksi_save.*(s1(1) - s0(1));
    u1q = interp1(x1, u1, x1q);
else
    x1q = ksi_save.*NaN;
    u1q = ksi_save.*NaN;
end
x2 = s1(1) + ksi2.*(s2(1) - s1(1));
x2q = s1(1) + ksi_save.*(s2(1) - s1(1));
u2q = interp1(x2, u2, x2q);
if isUpperPhase
    x3 = s2(1) + ksi3.*(s3(1) - s2(1));
    x3q = s2(1) + ksi_save.*(s3(1) - s2(1));
    u3q = interp1(x3, u3, x3q);
else
    x3q = ksi_save.*NaN;
    u3q = ksi_save.*NaN;
end
U(:, 1) = [u1q; u2q; u3q];
X(:, 1) = [x1q; x2q; x3q];

saveTime = tau_save;
saveId = 1;

time = 0;
tau0 = tau;     % Шаг по времени, заданный пользователем
n = 1;
dsMin = 0.01/x0;
dlMin = 1*1e-3/x0; % Нижняя граница толщины новой фазы

tic;
while time <= tMax
    u1_past = u1;
    u2_past = u2;
    u3_past = u3;
    
    % Вычисление скоростей движения границ
    C1 = 1/beta/( s1(n) - s0(n) );
    C2 = lambda2/(lambda1*beta*( s2(n) - s1(n) ));
    C3 = 1/beta/( s3(n) - s2(n) );
    ds0dt = 0;
    ds1dt = C2/(2*h2)*(-3*u2_past(1) + 4*u2_past(2) - u2_past(3)) - ...
        C1/(2*h1)*(3*u1_past(end) - 4*u1_past(end-1) + u1_past(end-2));
    ds2dt = C2/(2*h2)*(-u2_past(end) + 4*u2_past(end-1) - 3*u2_past(end-2)) - ...
        C3/(2*h3)*(-3*u3_past(1) + 4*u3_past(2) - u3_past(3));
    ds3dt = 0;
    
    % Интегрирование уравнений движения границ
    s0(n+1) = s0(n) + tau*ds0dt;
    s1(n+1) = s1(n) + tau*ds1dt;
    s2(n+1) = s2(n) + tau*ds2dt;
    s3(n+1) = s3(n) + tau*ds3dt;
    
    if abs(s2(n+1) - s2(n)) > dsMin && isUpperPhase
        tau = tau/2;
        continue;
    end
    if abs(s1(n+1)-s1(n)) > dsMin && isLowerPhase
        tau = tau/2;
        continue;
    end
    
    if ~isUpperPhase
        s2(n+1) = s3(n+1);
        ds2dt = ds3dt;
    end
    if ~isLowerPhase
        s1(n+1) = s0(n+1);
        ds1dt = ds0dt;
    end
    
    % Вырождение верхней и нижней фаз
    if s2(n+1) >= s3(n+1) && isUpperPhase
        s2(n+1) = s3(n+1);
        ds2dt = ds3dt;
        isUpperPhase = false;
    end
    if s1(n+1) <= s0(n+1) && isLowerPhase
        s1(n+1) = s0(n+1);
        ds1dt = ds0dt;
        isLowerPhase = false;
    end
    time = time + tau;
    
    % Получение распределения тепла для первой фазы (если она есть)
    if isLowerPhase
        [A1, b1] = getSysMat(u1_past, 1, tau, h1, s1(n+1), s0(n+1), ds1dt, ds0dt, ...
           alpha(1:2, :), g0(time), g1(time));
        %u1 = A \ b;
        u1 = solveWithBackslash(A1, b1);
        %u1 = solveWithThomas(A, b);
    end
    
    % Получение распределения тепла для второй фазы
    if isUpperPhase
        alphaUpper = alpha(4, :);
        gUpper = g3;
    else
        alphaUpper = alpha(6, :);
        gUpper = g5;
    end
    if isLowerPhase
        alphaLower = alpha(3, :);
        gLower = g2;
    else
        alphaLower = alpha(1, :);
        gLower = g0;
    end
    [A2, b2] = getSysMat(u2_past, kappa, tau, h2, s2(n+1), s1(n+1), ds2dt, ds1dt, ...
        [alphaLower; alphaUpper], gLower(time), gUpper(time));
    %u2 = A \ b;
    u2 = solveWithBackslash(A2, b2);
    %u2 = solveWithThomas(A, b);
    
    % Получение распределения тепла для третьей фазы
    if isUpperPhase
        [A3, b3] = getSysMat(u3_past, 1, tau, h3, s3(n+1), s2(n+1), ds3dt, ds2dt, ...
           alpha(5:6, :), g4(time), g5(time));
        %u3 = A \ b;
        u3 = solveWithBackslash(A3, b3);
        %u3 = solveWithThomas(A, b);
    end
    
    % Зарождение верхней фазы
    if u2(end) > Uf && u2(end-1) > Uf && ~isUpperPhase
        % Поиск номера последнего узла, который должен быть водой
        id = length(u2);
        for i = length(u2):-1:1
            if u2(i) < Uf
                id = i + 1;
                break;
            end
        end
        
        % Вычисление толщины новой фазы
        x = s1(n+1) + ksi2.*(s2(n+1) - s1(n+1));
        dl = c2*rho2/qf/rho1*trapz(x(id:end), abs(u2(id:end) - Uf)*U0);
        
        if dl >= dlMin
            s2(n+1) = s3(n+1) - dl;
            u3 = ones(Np3, 1)*Uf;
            u2(id:end) = Uf;
            u2 = interp1(x, u2, s1(n+1) + ksi2.*(s2(n+1) - s1(n+1)), 'linear', 'extrap');
            isUpperPhase = true;
        end
       
    end
    
    % Зарождение нижней фазы
    if u2(1) > Uf && u2(2) > Uf && ~isLowerPhase
        % Поиск номера последнего узла, который должен быть водой
        id = 1;
        for i = 1:length(u2)
            if u2(i) < Uf
                id = i - 1;
                break;
            end
        end
        
        % Вычисление толщины новой фазы
        x = s1(n+1) + ksi2.*(s2(n+1) - s1(n+1));
        dl = c2*rho2/qf/rho1*trapz(x(1:id), abs(u2(1:id) - Uf)*U0);
        
        if dl >= dlMin
            s1(n+1) = s0(n+1) + dl;
            u1 = ones(Np1, 1)*Uf;
            u2(1:id) = Uf;
            u2 = interp1(x, u2, s1(n+1) + ksi2.*(s2(n+1) - s1(n+1)), 'linear', 'extrap');
            isLowerPhase = true;
        end
    end
    
    % Проседание льда из-за различной плотности льда и воды
    dL = (s1(n+1) - s1(n) - (s2(n+1) - s2(n)) )*( 1 - rho1/rho2 );
    s2(n+1) = s2(n+1) + dL;
    s3(n+1) = s3(n+1) + dL;
    
    % Аккумуляция
    dL = accumRate/(365.25*24*3600)*tau*t0/rho2/x0;
    s2(n+1) = s2(n+1) + dL;
    s3(n+1) = s3(n+1) + dL;
    
    % Запись результатов
    t(n + 1) = time;
    if (saveTime < time)
        fprintf("Progress: %4.2f%%\n", saveTime/tMax*100);
        saveTime = saveTime + tau_save;
        saveId = saveId + 1;
        
        x1 = s0(n+1) + ksi1.*(s1(n+1) - s0(n+1));
        x2 = s1(n+1) + ksi2.*(s2(n+1) - s1(n+1));
        x3 = s2(n+1) + ksi3.*(s3(n+1) - s2(n+1));
        if isLowerPhase
            x1q = s0(n+1) + ksi_save.*(s1(n+1) - s0(n+1));
            u1q = interp1(x1, u1, x1q, 'linear', 'extrap');
        else
            x1q = ksi_save.*NaN;
            u1q = ksi_save.*NaN;
        end
        x2q = s1(n+1) + ksi_save.*(s2(n+1) - s1(n+1));
        u2q = interp1(x2, u2, x2q, 'linear', 'extrap');
        if isUpperPhase
            x3q = s2(n+1) + ksi_save.*(s3(n+1) - s2(n+1));
            u3q = interp1(x3, u3, x3q, 'linear', 'extrap');
        else
            x3q = ksi_save.*NaN;
            u3q = ksi_save.*NaN;
        end
        U(:, saveId) = [u1q; u2q; u3q];
        X(:, saveId) = [x1q; x2q; x3q];
        T(:, saveId) = ones(nRows, 1)*time;
    end
    
    n = n+1;
    tau = tau0;
end
s = [s0;s1;s2;s3];

% Масшабирование к исходной размерности
X = X*x0;
T = T*t0;
U = U*U0;
s = s*x0;
t = t*t0;

elapsedTime = toc;
fprintf("Elapsed time for Stefan Problem Solver: %4.2f sec.\n", elapsedTime);
end

