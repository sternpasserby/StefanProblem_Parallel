function [s, t, U, X, T] = StefanProblemSolver(pc, bc, varargin)
%STEFANPROBLEMSOLVER Решатель трехфазной задачи Стефана
%   На вход подаются 3 структуры, число узлов сетки для каждой фазы Np,
%   размер временнОго шага tau и время моделирования tMax. Возвращает 3
%   матрицы U, X, T и вектор s.

%%% ПАРСИНГ ВХОДНЫХ ПАРАМЕТРОВ
% Задание значений входных параметров по умолчанию
defaultNp = [100 10000 100];
defaultIc = struct('s', [-130; -129; 1100; 1100], ...
                   'dsdt', zeros(4, 1), ...
                   'x1', linspace(-130, -129, defaultNp(1)), ...
                   'u1', 273.15 + zeros(defaultNp(1), 1), ...
                   'x2', linspace(-129, 1100, defaultNp(2)), ...
                   'u2', 273.15 - 2*ones(defaultNp(2), 1), ...
                   'x3', ones(defaultNp(3), 1)*1100, ...
                   'u3', 273.15 + ones(defaultNp(3), 1), ...
                   'tInit', 0);
defaultAccumRate = 85.5;
defaultTau = 3600*24*14;
defaultTMax = 20*365.25*24*3600; 
defaultChLength = defaultIc.s(3) - defaultIc.s(2);
defaultChTemperature = pc.Uf;
defaultGridType = 'Uniform';
defaultNpSave = [100 200 100];
defaultTauSave = defaultTau;
defaultDs = 0.01;
defaultMinNewPhaseThickness = 1*1e-3;

% Настройка объекта parserObj типа InputParser
parserObj = inputParser;
parserObj.StructExpand = false;
addRequired(parserObj, 'pc');
addRequired(parserObj, 'bc');
addOptional(parserObj, 'ic', defaultIc);
addParameter(parserObj, 'Np', defaultNp);
addParameter(parserObj, 'accumRate', defaultAccumRate);
addParameter(parserObj, 'tau', defaultTau);
addParameter(parserObj, 'tMax', defaultTMax);
addParameter(parserObj, 'chLength', defaultChLength);
addParameter(parserObj, 'chTemperature', defaultChTemperature);
addParameter(parserObj, 'gridType', defaultGridType);
addParameter(parserObj, 'NpSave', defaultNpSave);
addParameter(parserObj, 'tauSave', defaultTauSave);
addParameter(parserObj, 'minDs', defaultDs);
addParameter(parserObj, 'minNewPhaseThickness', defaultMinNewPhaseThickness);

% Распаковка структур
lambda1 = pc.lambda1;
     c1 = pc.c1;
   rho1 = pc.rho1;
%  a1_sq = pc.a1_sq;
lambda2 = pc.lambda2;
     c2 = pc.c2;
   rho2 = pc.rho2;
%  a2_sq = pc.a2_sq;
     Uf = pc.Uf;
     qf = pc.qf;
  alpha = bc.alpha;
rho = (rho1 + rho2)/2;

% Получение значений параметров от объекта parserObj типа InputParser
parse(parserObj, pc, bc, varargin{:});
           accumRate = parserObj.Results.accumRate;
                  ic = parserObj.Results.ic;
                  Np = parserObj.Results.Np;
                 tau = parserObj.Results.tau;
                tMax = parserObj.Results.tMax;
                  x0 = parserObj.Results.chLength;
                  U0 = parserObj.Results.chTemperature;
            gridType = parserObj.Results.gridType;
              NpSave = parserObj.Results.NpSave;
             tauSave = parserObj.Results.tauSave;
               minDs = parserObj.Results.minDs;
minNewPhaseThickness = parserObj.Results.minNewPhaseThickness;
         
% Вычисление безразмерных параметров
t0 = rho1*c1*x0*x0/lambda1;
beta = qf*rho / (rho1*c1*U0);
kappa = c2*rho2/lambda2 * lambda1 / (c1*rho1); 

% Обезразмеривание исходных данных
g0 = @(t)(bc.g0(t*t0)/U0);
g1 = @(t)(bc.g1(t*t0)/U0);
u1 = ic.u1/U0;
u2 = ic.u2/U0;
u3 = ic.u3/U0;
S = ic.s/x0;
tInit = ic.tInit/t0;
tau = tau/t0;
tauSave = tauSave/t0;
tMax = tMax/t0;
Uf = Uf/U0;
alpha(:, 2) = bc.alpha(:, 2)/x0;
minDs = minDs/x0;
minNewPhaseThickness = minNewPhaseThickness/x0;
s0 = ic.s(1)/x0;

% Настройка function handle на функцию, которая генерирует заданную сетку
switch gridType
    case 'Uniform'
        getGrid = @reconstructGrid_uniform;
    case 'SplineBased'
        getGrid = @reconstructGrid_spline;
    case 'SigmoidBased'
        getGrid = @(x)( reconstructGrid_sigmoid(x, 0.01*1e-3/(ic.s(3)-ic.s(2)) ) );
    otherwise
        error("Unknown grid type %s", gridType);
end

% Интерполяция начального распределения температуры для каждой фазы
ksiNew = getGrid(Np(1));
if ~(ic.s(1)==ic.s(2))
    ksi = ( ic.x1 - ic.s(1) ) ./ ( ic.s(2) - ic.s(1) );
    %u1 = csInterp(ksi, u1, ksiNew, [alpha(1, :); 1 0], g0(tInit), Uf);
    u1 = interp1(ksi, u1, ksiNew, 'linear', 'extrap')';
else
    ksi = ic.x1.*NaN;
    u1 = u1.*NaN;
end
ph1 = struct('ksi', ksiNew, ...
             'u', u1, ...
             'uPast', u1, ...
             's', S(2), ...
             'sPast', S(2), ...
             'dsdt', 0, ...
             'exists', ~(ic.s(1)==ic.s(2)) );
%ph1.s(1) = ic.s(2)/x0;
ksiNew = getGrid(Np(2));
ksi = ( ic.x2 - ic.s(2) ) ./ ( ic.s(3) - ic.s(2) );
u2 = interp1(ksi, u2, ksiNew, 'linear', 'extrap')';
ph2 = struct('ksi', ksiNew, ...
             'u', u2, ... %'u', csInterp(ksi, u2, ksiNew, [1 0; 1 0], Uf, Uf), ...
             'uPast', u2, ...
             's', S(3), ...
             'sPast', S(3), ...
             'dsdt', 0, ...
             'exists', true );
%ph2.s(1) = ic.s(3)/x0;
ksiNew = getGrid(Np(3));
if ~(ic.s(3)==ic.s(4))
    ksi = ( ic.x3 - ic.s(3) ) ./ ( ic.s(4) - ic.s(3) );
    %u3 = csInterp(ksi, u3, ksiNew, [1 0; alpha(2, :)], Uf, g1(tInit));
    u3 = interp1(ksi, u3, ksiNew, 'linear', 'extrap')';
else
    ksi = ic.x3.*NaN;
    u3 = u3.*NaN;
end
ph3 = struct('ksi', ksiNew, ...
             'u', u3, ...
             'uPast', u3, ...
             's', S(4), ...
             'sPast', S(4), ...
             'dsdt', 0, ...
             'exists', ~(ic.s(3)==ic.s(4)) );
%ph3.s(1) = ic.s(4)/x0;

numOfTimeSteps = ceil(tMax/tau);
% t = zeros(1, numOfTimeSteps + 1);
% s = zeros(4, numOfTimeSteps + 1);
t = zeros(1, 10*numOfTimeSteps);
s = zeros(4, 10*numOfTimeSteps);

nRows = sum(NpSave);
nCols = min(numOfTimeSteps, ceil(tMax/tauSave) ) + 1;
X = zeros(nRows, nCols);
U = zeros(nRows, nCols);
T = zeros(nRows, nCols);

% A1 = sparse(Np(1));
% A2 = sparse(Np(2));
% A3 = sparse(Np(3));
ld1 = zeros(Np(1) - 1, 1);
ld2 = zeros(Np(2) - 1, 1);
ld3 = zeros(Np(3) - 1, 1);
md1 = zeros(Np(1), 1);
md2 = zeros(Np(2), 1);
md3 = zeros(Np(3), 1);
ud1 = zeros(Np(1) - 1, 1);
ud2 = zeros(Np(2) - 1, 1);
ud3 = zeros(Np(3) - 1, 1);
b1 = zeros(Np(1), 1);
b2 = zeros(Np(2), 1);
b3 = zeros(Np(3), 1);

% Запись начальных условий в выходные массивы
ksiSave1 = getGrid(NpSave(1));
if ph1.exists
    %u1q = csInterp(ph1.ksi, ph1.u, ksiSave1, [alpha(1, :); 1 0], g0(tInit), Uf);
    u1q = interp1(ph1.ksi, ph1.u, ksiSave1, 'linear', 'extrap')';
    x1q = S(1) + ksiSave1*( ph1.s(1) - S(1) );
else
    x1q = ksiSave1.*NaN;
    u1q = ksiSave1'.*NaN;
end
ksiSave2 = getGrid(NpSave(2));
%u2q = csInterp(ph2.ksi, ph2.u, ksiSave2, [1 0; 1 0], Uf, Uf);
u2q = interp1(ph2.ksi, ph2.u, ksiSave2, 'linear', 'extrap')';
x2q = ph1.s(1) + ksiSave2*( ph2.s(1) - ph1.s(1) );
ksiSave3 = getGrid(NpSave(3));
if ph3.exists
    %u3q = csInterp(ph3.ksi, ph3.u, ksiSave3, [1 0; alpha(2, :)], Uf, g1(tInit));
    u3q = interp1(ph3.ksi, ph3.u, ksiSave3, 'linear', 'extrap')';
    x3q = ph2.s(1) + ksiSave3*( ph3.s(1) - ph2.s(1) );
else
    x3q = ksiSave3.*NaN;
    u3q = ksiSave3'.*NaN;
end
U(:, 1) = [u1q; u2q; u3q];
X(:, 1) = [x1q'; x2q'; x3q'];
s(:, 1) = S';
t(1) = tInit;

saveTime = tauSave;
saveId = 2;

time = tInit;
tau0 = tau;     % Шаг по времени, заданный пользователем
n = 2;
numOfActualTimeSteps = 0;

%tauAr = zeros(size(t));

clear printProgressBar;
tic;
while time <= tMax
%while n <= numOfTimeSteps + 1
    ph1.uPast = ph1.u;
    ph2.uPast = ph2.u;
    ph3.uPast = ph3.u;
    ph1.sPast = ph1.s;
    ph2.sPast = ph2.s;
    ph3.sPast = ph3.s;
    
    % Вычисление скоростей движения границ
    C1 = 1/beta/( ph1.s - s0 );
    C2 = lambda2/(lambda1*beta*( ph2.s - ph1.s ));
    C3 = 1/beta/( ph3.s - ph2.s );
    h0 = ph2.ksi(2) - ph2.ksi(1);
    h1 = ph2.ksi(3) - ph2.ksi(2);
    h_Npm1 = ph1.ksi( Np(1) ) - ph1.ksi( Np(1)-1 );
    h_Npm2 = ph1.ksi( Np(1)-1 ) - ph1.ksi( Np(1)-2 );
    ph1.dsdt = C2/( h0*h1*(h0+h1) )*( ( ph2.uPast(2) - ph2.uPast(1) )*(h0+h1)^2 - (ph2.uPast(3)-ph2.uPast(1))*h0^2 ) - ...
        C1/( h_Npm1*h_Npm2*(h_Npm1+h_Npm2) )*...
        ( ( ph1.uPast(Np(1)) - ph1.uPast(Np(1)-1) )*(h_Npm2+h_Npm1)^2 + (ph1.uPast(Np(1)-2) - ph1.uPast(Np(1)))*h_Npm1^2 );
    h_Npm1 = ph2.ksi( Np(2) ) - ph2.ksi( Np(2)-1 );
    h_Npm2 = ph2.ksi( Np(2)-1 ) - ph2.ksi( Np(2)-2 );
    h0 = ph3.ksi(2) - ph3.ksi(1);
    h1 = ph3.ksi(3) - ph3.ksi(2);
    ph2.dsdt = C2/( h_Npm1*h_Npm2*(h_Npm1+h_Npm2) )*...
        ( ( ph2.uPast(Np(2)) - ph2.uPast(Np(2)-1) )*(h_Npm2+h_Npm1)^2 + (ph2.uPast(Np(2)-2)-ph2.uPast(Np(2)))*h_Npm1^2 ) - ...
        C3/( h0*h1*(h0+h1) )*( ( ph3.uPast(2) - ph3.uPast(1) )*(h0+h1)^2 - (ph3.uPast(3)-ph3.uPast(1))*h0^2 );
    ph3.dsdt = 0;
    
    if ~ph3.exists
        ph2.s = ph3.s;
        ph2.dsdt = ph3.dsdt;
    end
    if ~ph1.exists
        ph1.s = s0;
        ph1.dsdt = 0;
    end
    
    %fprintf("%10.4e %10.4e %10.4e %10.4e\n", tau, ph1.dsdt, ph2.dsdt, ph3.dsdt);
    
    if (abs( tau*ph2.dsdt ) > minDs && ph3.exists) || ...
            (abs( tau*ph3.dsdt ) > minDs && ph1.exists)
        tau = tau/2;
        continue;
    end
    
    % Интегрирование уравнений движения границ
    ph1.s = ph1.sPast + tau*ph1.dsdt;
    ph2.s = ph2.sPast + tau*ph2.dsdt;
    ph3.s = ph3.sPast + tau*ph3.dsdt;
    
    % Вырождение верхней и нижней фаз
    if ph2.s >= ph3.s && ph3.exists
        ph2.s = ph3.s;
        ph2.dsdt = ph3.dsdt;
        ph3.exists = false;
    end
    if ph1.s <= s0 && ph1.exists
        ph1.s = s0;
        ph1.dsdt = 0;
        ph1.exists = false;
    end
    time = time + tau;
    
    % Температура фазового перехода с поправкой на давление толщи льда
    Uf_adj = (273.15 - 7.43*1e-8*rho2*9.81*( ph2.s - ph1.s )*x0)/U0;     % adj - adjusted
    
    % Получение распределения тепла для первой фазы (если она есть)
    if ph1.exists
        [ld1, md1, ud1, b1] = getSysMat(ph1.uPast, 1, tau, ph1.ksi, ph1.s, s0, ph1.dsdt, 0, ...
           [alpha(1, :); 1 0], g0(time), Uf_adj);
        ph1.u = mex_TDMA(md1, ud1, ld1, b1);
    end
    
    % Получение распределения тепла для второй фазы
    if ph3.exists
        alphaUpper = [1 0];
        gUpper = Uf;
    else
        alphaUpper = alpha(2, :);
        gUpper = g1(time);
    end
    if ph1.exists
        alphaLower = [1 0];
        gLower = Uf_adj;
    else
        alphaLower = alpha(1, :);
        gLower = g0(time);
    end
    [ld2, md2, ud2, b2] = getSysMat(ph2.uPast, kappa, tau, ph2.ksi, ph2.s, ph1.s, ph2.dsdt, ph1.dsdt, ...
        [alphaLower; alphaUpper], gLower, gUpper);
    ph2.u = mex_TDMA(md2, ud2, ld2, b2);
    
    % Получение распределения тепла для третьей фазы
    if ph3.exists
        [ld3, md3, ud3, b3]  = getSysMat(ph3.uPast, 1, tau, ph3.ksi, ph3.s, ph2.s, ph3.dsdt, ph2.dsdt, ...
           [1 0; alpha(2, :)], Uf, g1(time));
        ph3.u = mex_TDMA(md3, ud3, ld3, b3);
    end
    
    % Зарождение верхней фазы
    if ph2.u(end) > Uf && ph2.u(end-1) > Uf && ~ph3.exists
        % Поиск номера последнего узла, который должен быть водой
        id = Np(2);
        for i = Np(2):-1:1
            if ph2.u(i) < Uf
                id = i + 1;
                break;
            end
        end
        
        % Вычисление толщины новой фазы
        x = ph1.s + ph2.ksi.*(ph2.s - ph1.s);
        dl = c2*rho2/qf/rho1*trapz(x(id:end), abs(ph2.u(id:end) - Uf)*U0);
        
        if dl >= minNewPhaseThickness
            ph2.s = ph3.s - dl;
            ph3.u = ones(Np(3), 1)*Uf;
            ph2.u(id:end) = Uf;
            %ph2.u = csInterp(x, ph2.u, ph1.s(n+1) + ph2.ksi.*(ph2.s(n+1) - ph1.s(n+1)), [1 0; 1 0], Uf_adj, Uf);
            ph2.u = interp1(x, ph2.u, ph1.s + ph2.ksi.*(ph2.s - ph1.s), 'linear', 'extrap')';
            ph3.exists = true;
        end
       
    end
    
    % Зарождение нижней фазы
    if ph2.u(1) > Uf_adj && ph2.u(2) > Uf_adj && ~ph1.exists
        % Поиск номера последнего узла, который должен быть водой
        id = 1;
        for i = 1:Np(2)
            if ph2.u(i) < Uf_adj
                id = i - 1;
                break;
            end
        end
        
        % Вычисление толщины новой фазы
        x = ph1.s + ph2.ksi.*(ph2.s - ph1.s);
        dl = c2*rho2/qf/rho1*trapz(x(1:id), abs(ph2.u(1:id) - Uf_adj)*U0);
        
        if dl >= minNewPhaseThickness
            ph1.s = s0 + dl;
            ph1.u = ones(Np(1), 1)*Uf_adj;
            ph2.u(1:id) = Uf_adj;
            ph2.u = interp1(x, ph2.u, ph1.s + ph2.ksi.*(ph2.s - ph1.s), 'linear', 'extrap')';
            ph1.exists = true;
        end
    end
    
    % Проседание льда из-за различной плотности льда и воды
    dL = (ph1.s - ph1.sPast - (ph2.s - ph2.sPast) )*( 1 - rho1/rho2 );
    ph2.s = ph2.s + dL;
    ph3.s = ph3.s + dL;
    
    % Аккумуляция
    dL = accumRate/(365.25*24*3600)*tau*t0/rho2/x0;
    ph2.s = ph2.s + dL;
    ph3.s = ph3.s + dL;
    
%     x1 = s0 + ph1.ksi*( ph1.s - s0 );
%     x2 = ph1.s + ph2.ksi*( ph2.s - ph1.s );
%     x3 = ph2.s + ph3.ksi*( ph3.s - ph2.s );
%     plot(x1*x0, ( ph1.u - 1 )*U0)
%     hold on
%     plot(x2*x0, ( ph2.u - 1 )*U0)
%     plot(x3*x0, ( ph3.u - 1 )*U0)
%     hold off
    
%     if abs(time - tInit - (n-1)*tau0) < 1e-6*tau0
%         tau = min(tau0, tMax - time);
%         %t(n) = time;
%         %s(:, n) = [s0; ph1.s; ph2.s; ph3.s];
%         n = n+1;
%     else
%         tau = tInit + (n-1)*tau0 - time;
%     end
    
    numOfActualTimeSteps = numOfActualTimeSteps + 1;
    t(numOfActualTimeSteps + 1) = time;
    s(:, numOfActualTimeSteps + 1) = [s0; ph1.s; ph2.s; ph3.s];
    tauAr(numOfActualTimeSteps) = tau;
    
    % Запись результатов
    %t(n + 1) = time;
    %if (saveTime <= time)
    if saveTime - time <= 1e-6*tauSave || n == numOfTimeSteps + 2
        %fprintf("Progress: %4.2f%%\n", saveTime/tMax*100);
        printProgressBar(saveTime, tMax);
        
        %ksiSave1 = getGrid(NpSave(1));
        if ph1.exists
            %u1q = csInterp(ph1.ksi, ph1.u, ksiSave1, [alpha(1, :); 1 0], g0(time), Uf);
            u1q = interp1(ph1.ksi, ph1.u, ksiSave1, 'linear', 'extrap')';
            x1q = s0 + ksiSave1*( ph1.s - s0 );
        else
            x1q = ksiSave1.*NaN;
            u1q = ksiSave1'.*NaN;
        end
        ksiSave2 = getGrid(NpSave(2));
        %u2q = csInterp(ph2.ksi, ph2.u, ksiSave2, [1 0; 1 0], Uf, Uf);
        u2q = interp1(ph2.ksi, ph2.u, ksiSave2, 'linear', 'extrap')';
        x2q = ph1.s + ksiSave2*( ph2.s - ph1.s );
        ksiSave3 = getGrid(NpSave(3));
        if ph3.exists
            %u3q = csInterp(ph3.ksi, ph3.u, ksiSave3, [1 0; alpha(2, :)], Uf, g1(time));
            u3q = interp1(ph3.ksi, ph3.u, ksiSave3, 'linear', 'extrap')';
            x3q = ph2.s + ksiSave3*( ph3.s - ph2.s );
        else
            x3q = ksiSave3.*NaN;
            u3q = ksiSave3'.*NaN;
        end
        U(:, saveId) = [u1q; u2q; u3q];
        X(:, saveId) = [x1q'; x2q'; x3q'];
        T(:, saveId) = ones(nRows, 1)*time;
        
        saveTime = saveTime + tauSave;
        saveId = saveId + 1;
    end
    
    tau = tau0;
    n = n+1;
end
id = find(t == 0);
id(1) = [];
t(id) = [];
s(:, id) = [];
%tauAr(id) = [];

% Масшабирование к исходной размерности
X = X*x0;
T = T*t0;
U = U*U0;
s = s*x0;
t = t*t0;

elapsedTime = toc;
fprintf("Elapsed time for Stefan Problem Solver: %4.2f sec.\n", elapsedTime);
end

function xNew = reconstructGrid_uniform(Np)
    xNew = linspace(0, 1, Np);
end

function xNew = reconstructGrid_spline(Np)
    X = [0 1/4 3/4 1];
    H = 0.1;
    ppF = struct();
    ppF.form = 'pp';
    ppF.breaks = X;
    c1 = ( (1-H)/3 + H )*(X(2)-X(1));
    ppF.coefs = [ (1-H)/3/(X(2)-X(1))^2 (H-1)/(X(2)-X(1)) 1 0; ...
                   0 0 H c1; ...
                   (1-H)/3/(X(4)-X(3))^2 0 H c1 + H*(X(3)-X(2))];
    ppF.pieces = 3;
    ppF.order = 4;
    ppF.dim = 1;
    
    der_ppF = fnder(ppF, 1);
    
    y = linspace( ppval(ppF, 0), ppval(ppF, 1), Np );
    xNew = linspace(0, 1, Np);
    for i = 1:10
        xNew = xNew - (ppval(ppF, xNew)-y)./(ppval(der_ppF, xNew));
        %fprintf("Error: %6.2e\n", max(abs(xNew-temp)));
        %temp = xNew;
    end
    
    if ~issorted(xNew)
        error("xNew is not monotonically increasing!")
    end
    if ~all(xNew >= 0 & xNew <= 1)
        error("Some of xNew elements lie outside of [0, 1] domain!");
    end
end

function xNew = reconstructGrid_sigmoid(Np, dx1)
    x1 = -10;
    alpha = -1/x1*log( 1/dx1 - 1 );
    x = linspace(x1, -x1, Np-2);
    xNew = [0 1 ./ (1 + exp(-alpha*x)) 1];
    xNew(2) = (xNew(3)-xNew(1))/2;
    xNew(Np-1) = xNew(Np-2) + (xNew(Np)-xNew(Np-2))/2;
end

function yNew = csInterp(x, y, xq, alpha, g0, g1)
    conds = [1 1];
    if alpha(1, 2) == 0
        conds(1) = 0;
        C0 = g0/alpha(1, 1);
    else
        conds(1) = 1;
        C0 = (g0 - alpha(1, 1)*y(1)) / alpha(1, 2);
    end
    if alpha(2, 2) == 0
        conds(2) = 0;
        C1 = g1/alpha(2, 1);
    else
        conds(2) = 1;
        C1 = (g1 - alpha(2, 1)*y(end)) / alpha(2, 2);
    end
    pp = csape(x, [C0; y; C1]', conds);
    yNew = ppval(pp, xq)';
end

function printProgressBar(n, N)
%PRINTPROGRESSBAR Summary of this function goes here
%   Detailed explanation goes here

persistent barLength;
persistent leftSymbol;
persistent rightSymbol;
persistent completedSymbol;
persistent todoSymbol;
persistent reverseStr;

fraction = n/N;
isInit = (isempty(reverseStr) && fraction <= 1);
if isInit
          barLength = 25;
         leftSymbol = '[';
        rightSymbol = ']';
    completedSymbol = '=';
         todoSymbol = ' ';
         reverseStr = '';
end

numOfSegments = floor(fraction*barLength);
completedStr = repmat(completedSymbol, 1, numOfSegments);
todoStr = repmat(todoSymbol, 1, barLength - numOfSegments);
%msg = sprintf([leftSymbol completedStr todoStr  rightSymbol ' %5.2f%% (%d / %d)'], fraction*100, n, N);
msg = sprintf(...
    [leftSymbol ... 
    completedStr ... 
    todoStr ... 
    rightSymbol, ...
    ' %5.2f%%'], fraction*100);
    %' %5.2f%% (%d/%d)'], fraction*100, n, N);
disp([reverseStr msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg)+1);
end


