clear;
delete(gcp('nocreate'))

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
tMax = 200*365.25*24*3600;        % Время, до которого необходимо моделировать, с
tau = 3600*24*30;     % Шаг по времени, с
tauSave = 3600*24*365.25;

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
x = [];
i = 0;
k = 1;
pool = parpool();
for j = 1:length(Data.X)
    if Data.Y(j) ~= -3000
        continue;
    end
    bedrock = Data.Bedrock_m(j);
    iceSurf = Data.Surface_m(j);
    iceThickness = Data.IceThickness_m(j);
    GHF = Data.GHF_Martos_mWm2(j);
    accumRate = Data.AccumRate_kg1m2a1(j);
    if bedrock > 0 && (iceSurf - iceThickness ~= bedrock) % величины должны быть целыми, так что можно не исхищряться 
        continue;
    end
    if bedrock < 0 && iceSurf - iceThickness > 0
        continue;
    end
    if iceThickness == 0
        continue
    end
    i = i + 1;
    if mod(i-1, 50) ~= 0
        continue;
    end
    x(end + 1) = Data.X(j);
    
    ic.s0 = bedrock;
    ic.s1 = iceSurf - iceThickness;
    ic.s2 = iceSurf;
    ic.s3 = iceSurf;
    ic.accumRate = accumRate;
    bc.g0 =  @(t)(GHF/1000);
    F(k) = parfeval(pool, @StefanProblemSolver, 2, pc, bc, ic, 0.25, tau, tMax, 100, tauSave);
    k = k + 1;
end

S = cell(size(F));
t = 0:tauSave:tMax;
for i = 1:length(S)
    [completedIdx,value] = fetchNext(F);
    S{completedIdx} = interp1(F(completedIdx).OutputArguments{2}, F(completedIdx).OutputArguments{1}', t)';
    printProgressBar(i/length(S));
    %fprintf('Got result with index: %d.\n', completedIdx);
end

vid = VideoWriter("vid");
open(vid);
X = unique(Data.X);
s = zeros(length(S), 4);
for i = 1:length(t)
    for j = 1:length(S)
        s(j, 1) = S{j}(1, i);
        s(j, 2) = S{j}(2, i);
        s(j, 3) = S{j}(3, i);
        s(j, 4) = S{j}(4, i);
    end
    plot(x, s, '-o');
    %disp(y(1));
    %bar(X, s, 1, 'Stacked');
    %axis([-inf inf 0 1]);
    drawnow();
    title(sprintf("Time: %8.2f year", t(i)/3600/24/365.25));
    frame = getframe(gcf());
    writeVideo(vid, frame);
end
close(vid);
