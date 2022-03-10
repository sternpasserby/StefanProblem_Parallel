function runGlacierModelling(pool, resFolderPath, initDataFilename, points_id)
%RUNGLACIERMODELLING Summary of this function goes here
%   Detailed explanation goes here

dirName = resFolderPath + "\\Data.bin";
numOfPoints = length(points_id);
if isfolder(resFolderPath)
    fprintf("Folder ""%s"" already exists. Loading completed points indices.\n", resFolderPath);
    fid = fopen(dirName, "rb");
    M = fread(fid, 1, 'int');
    fseek(fid, M*4, 0);
    completedPoints_id = zeros(1, M);
    i = 1;
    while true
        id = fread(fid, 1, 'int');
        if isempty(id)
            break;
        else
            L = fread(fid, 1, 'int');
            completedPoints_id(i) = id;
            fseek(fid, 5*L*8, 0);
            i = i + 1;
        end
    end
    fclose(fid);
    completedPoints_id(i:end) = [];
    for i = 1:length(completedPoints_id)
        id = find( points_id == completedPoints_id(i), 1 );
        if ~isempty(id)
            points_id(id) = [];
        end
    end
    numOfPoints = length(points_id);
else
    fprintf("Folder ""%s"" does not exists. Creating folder.\n", resFolderPath);
    mkdir(resFolderPath);
    fid = fopen(dirName, "wb");
    fwrite(fid, numOfPoints, 'int');
    fwrite(fid, points_id, 'int');
    fclose(fid);
end
numOfPoints = length(points_id);

load(initDataFilename, 'Data');
batchSize = pool.NumWorkers;

pc = getPhysicalConstants();

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
tMax = 20*365.25*24*3600;        % Время, до которого необходимо моделировать, с
tau = 3600*24*30;     % Шаг по времени, с
tauSave = 3600*24*365.25;

%%% Начальные условия
ic = struct;                 % ic - initial conditions
ic.s0 = 0;            % Начальные положения границы раздела сред, м
ic.s1 = 0;
ic.s2 = 9;
ic.s3 = 10;
ic.u1 = zeros(1, Np) + 273.15 + 0;
ic.u2 = zeros(1, Np) + 273.15 - 3;
ic.u3 = zeros(1, Np) + 273.15 + 1;

taskInd = 0;
taskInd2pInd = zeros(batchSize, 1);
batchSize = min(batchSize, numOfPoints);
fprintf('Progress: ');
pb = ConsoleProgressBar();
pb.setProgress( 0, numOfPoints );
for i = 1:batchSize
    k = points_id(i);
    ic.s0 = Data.Bedrock_m(k);
    ic.s1 = Data.Surface_m(k) - Data.IceThickness_m(k);
    ic.s2 = Data.Surface_m(k);
    ic.s3 = Data.Surface_m(k);
    ic.accumRate = Data.AccumRate_kg1m2a1(k);
    bc.g0 =  @(t)(Data.GHF_Martos_mWm2(k)/1000);
    F(i) = parfeval(pool, @StefanProblemSolver, 2, pc, bc, ic, 0.25, tau, tMax, 100, tauSave);
    taskInd2pInd(i) = k;
end

numOfCompPoints = 0;
for i = batchSize+1:numOfPoints + batchSize
    
    % Получение результатов, запись их на диск
    taskInd = fetchNext(F);
    k = taskInd2pInd(taskInd);
    fid = fopen(dirName, "ab");
    L = length( F(taskInd).OutputArguments{2} );
    fwrite(fid, k, 'int');
    fwrite(fid, L, 'int');
%     fwrite(fid, L, 'int');
%     fwrite(fid, k*ones(5, L), 'double');
    fwrite(fid, F(taskInd).OutputArguments{2}, 'double');
    fwrite(fid, F(taskInd).OutputArguments{1}, 'double');
    fclose(fid);
    
    numOfCompPoints = numOfCompPoints + 1;
    pb.setProgress( numOfCompPoints, numOfPoints );
    
    % Загрузка новых точек
    if i <= length(points_id)
        k = points_id(i);
        ic.s0 = Data.Bedrock_m(k);
        ic.s1 = Data.Surface_m(k) - Data.IceThickness_m(k);
        ic.s2 = Data.Surface_m(k);
        ic.s3 = Data.Surface_m(k);
        ic.accumRate = Data.AccumRate_kg1m2a1(k);
        bc.g0 =  @(t)(Data.GHF_Martos_mWm2(k)/1000);
        F(taskInd) = parfeval(pool, @StefanProblemSolver, 2, pc, bc, ic, 0.25, tau, tMax, 100, tauSave);
        taskInd2pInd(taskInd) = k;
    end
    
end

end

function createResultsFile(filename, InitDataFilename, pointIndices)
    isCompleted = false;
    completedPoints = zeros(length(pointIndices), 1);
    Results = cell(length(pointIndices), 1);
    CreateDateTime = datetime();
    FinishDateTime = [];
    save(filename, ...
        'isCompleted', ...
        'pointIndices', ...
        'completedPoints', ...
        'Results', ...
        'CreateDateTime', ...
        'FinishDateTime', ...
        'InitDataFilename', ...
        '-v7.3', '-nocompression');
end

function pc = getPhysicalConstants()
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
end
