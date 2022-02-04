function times = measureT(pool, numOfTasks, numOfRuns)

%%% Физические константы
pc = struct;                  % pc - problem constants, константы задачи
pc.lambda1 = 0.6;             % Коэффициент теплопроводности воды, Вт / (м * K)
pc.c1 = 4180.6;               % Коэффициеент удельной теплоёмкости воды, Дж / (кг * К)
pc.rho1 = 1000;               % Плотность воды, кг/м^3
pc.a1_sq = pc.lambda1/...     % Коэффициент температуропроводности воды, м^2/с
    pc.c1/pc.rho1;            
pc.lambda2 = 2.33;            % Коэффициент теплопроводности льда, Вт / (м * K)
pc.c2 = 2110.0;               % Коэффициеент удельной теплоёмкости льда, Дж / (кг * К)
pc.rho2 = 916.7;              % Плотность льда, кг/м^3
pc.a2_sq = pc.lambda2/...     % Коэффициент температуропроводности льда, м^2/с
    pc.c2/pc.rho2;            
pc.qf = 330*1e3;              % Удельная теплота плавления льда, Дж / кг
%pc.rho = (rho1 + rho2)/2;     % Средняя плотность
%pc.L = 1;                     % Длина стержня, м
pc.Uf = 273.15;               % Температура фазового перехода, К

%%% Смешанные краевые условия
% Формат краевых условий:
% alpha00*u1(0, t) + alpha01*du1/dx(0, t) = g0(t) - для левого конца
% alpha10*u2(L, t) + alpha11*du2/dx(L, t) = g1(t) - для правого конца
bc = struct;                  % bc - boundary conditions
bc.alpha = zeros(6, 2);
bc.alpha = [0 -pc.lambda1; 
            1 0
            1 0
            1 0
            1 0
            1 0];
bc.g0 = @(t)(0.05);
bc.g1 = @(t)(pc.Uf);
bc.g2 = @(t)(pc.Uf);
bc.g3 = @(t)(pc.Uf);
bc.g4 = @(t)(pc.Uf);
bc.g5 = @(t)(- 4.3 + 8*sin(2*pi*t/31556952 + pi/2) + 273.15); %% Внимание, здесь сдвиг по фазе на pi/2
%bc.g5 = @(t)(1 + 273.15);

%%% Параметры численного решения
Np = 1000;            % Число узлов сетки для каждой фазы
tMax = 10*365.25*24*3600;        % Время, до которого необходимо моделировать, с
tau = 3600*24;     % Шаг по времени, с
tauSave = 3600*24*365.25/2;

%%% Начальные условия
ic = struct;                 % ic - initial conditions
ic.s0 = 0;            % Начальные положения границы раздела сред, м
ic.s1 = 0;
ic.s2 = 9;
ic.s3 = 10;
ic.u1 = zeros(1, Np) + 273.15 + 0;
ic.u2 = zeros(1, Np) + 273.15 - 1;
ic.u3 = zeros(1, Np) + 273.15 + 1;

% Сформировать массив индексов строк в таблице с данными
load('2021_03_30 AntarcticaBM2_parsed.mat', '-mat');
iAr = zeros(1, numOfTasks);
i = 1;
while i <= numOfTasks
%     if mod(i-1, 50) ~= 0
%         continue;
%     end
%     x(end + 1) = Data.X(i);
    bedrock = Data.Bedrock_m(i);
    iceSurf = Data.Surface_m(i);
    iceThickness = Data.IceThickness_m(i);
    GHF = Data.GHF_Martos_mWm2(i);
    accumRate = Data.AccumRate_kg1m2a1(i);
    if bedrock > 0 && (iceSurf - iceThickness ~= bedrock) % величины должны быть целыми, так что можно не исхищряться 
        continue;
    end
    if bedrock < 0 && iceSurf - iceThickness > 0
        continue;
    end
    iAr(i) = i;
    i = i + 1;
end

for k = 1:numOfRuns
    j = 0;
    time = tic;
    for i = iAr
        bedrock = Data.Bedrock_m(i);
        iceSurf = Data.Surface_m(i);
        iceThickness = Data.IceThickness_m(i);
        GHF = Data.GHF_Martos_mWm2(i);
        accumRate = Data.AccumRate_kg1m2a1(i);

        ic.s0 = bedrock;
        ic.s1 = iceSurf - iceThickness;
        ic.s2 = iceSurf;
        ic.s3 = iceSurf;
        ic.accumRate = accumRate;
        bc.g0 =  @(t)(GHF/1000);
        j = j + 1;
        F(j) = parfeval(pool, @StefanProblemSolver, 2, pc, bc, ic, Np, tau, tMax, 100, tauSave);
        %[s, t] = StefanProblemSolver(pc, bc, ic, Np, tau, tMax, 100, tauSave);
        %S{end + 1} = s;
    end

    S = cell(size(F));
    t = 0:tauSave:tMax;
    for i = 1:length(S)
        % fetchNext blocks until next results are available.
        [completedIdx,value] = fetchNext(F);
        S{completedIdx} = interp1(F(completedIdx).OutputArguments{2}, F(completedIdx).OutputArguments{1}', t)';
        %fprintf('Got result with index: %d.\n', completedIdx);
    end
    times(k) = toc(time);
end

end

