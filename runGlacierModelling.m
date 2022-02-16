function runGlacierModelling(pool, resFilename, initDataFilename, x, y)
%RUNGLACIERMODELLING Summary of this function goes here
%   Detailed explanation goes here

rfo = matfile('');  % Results File Object
isPointCompleted = [];
if isfile(resFilename)
    fprintf("File %s already exists. I am going to check it's contents.\n", resFilename);
    rfo = matfile(resFilename, 'Writable', true);
    if rfo.isCompleted
        error('The specified file already contains completed computations. I will not rewrite it.')
    end
    rfo = matfile(resFilename, 'Writable', true);
    if ~isequal(rfo.x, x)
        error('The vector x from the file does not match your vector x');
    end
    if ~isequal(rfo.y, y)
        error('The vector y from the file does not match your vector y');
    end
    if ~isequal(rfo.InitDataFilename, initDataFilename)
        error('The initDataFilename from the file does not match your initDataFilename');
    end
    isPointCompleted = rfo.isPointCompleted;
else
    fprintf("File %s does not exists. I am going to create it.\n", resFilename);
    createResultsFile(resFilename, initDataFilename, x, y);
    rfo = matfile(resFilename, 'Writable', true);
    isPointCompleted = false(length(y), length(x));
end

load(initDataFilename, 'Data');
batchSize = pool.NumWorkers;

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

taskNumber = 0;
ijArray = zeros(batchSize, 2);
fprintf('Progress: ');
pb = ConsoleProgressBar();
for i = 1:length(y)
    for j = 1:length(x)
        if isPointCompleted(i, j)
            continue;
        end
        
        k = find((Data.X == x(j)) & (Data.Y == y(i)));
        if isempty(k)
            isPointCompleted(i, j) = true;
            rfo.isPointCompleted(i, j) = true;
            rfo.Results(i, j) = {'There is no initial data for this point'};
            continue;
        end
        
        bedrock = Data.Bedrock_m(k);
        iceSurf = Data.Surface_m(k);
        iceThickness = Data.IceThickness_m(k);
        GHF = Data.GHF_Martos_mWm2(k);
        accumRate = Data.AccumRate_kg1m2a1(k);
        if bedrock > 0 && (iceSurf - iceThickness ~= bedrock) % величины должны быть целыми, так что можно не исхищряться 
            isPointCompleted(i, j) = true;
            rfo.isPointCompleted(i, j) = true;
            rfo.Results(i, j) = {'There is an empty space between ice and bedrock'};
            continue;
        end
        if bedrock < 0 && iceSurf - iceThickness > 0
            isPointCompleted(i, j) = true;
            rfo.isPointCompleted(i, j) = true;
            rfo.Results(i, j) = {'There is an empty space between ice and water at z = 0'};
            continue;
        end
        if iceThickness == 0
            isPointCompleted(i, j) = true;
            rfo.isPointCompleted(i, j) = true;
            rfo.Results(i, j) = {'Ice thickness is zero for this point'};
            continue
        end
        
        if taskNumber ~= batchSize
            taskNumber = taskNumber + 1;
            ic.s0 = bedrock;
            ic.s1 = iceSurf - iceThickness;
            ic.s2 = iceSurf;
            ic.s3 = iceSurf;
            ic.accumRate = accumRate;
            bc.g0 =  @(t)(GHF/1000);
            F(taskNumber) = parfeval(pool, @StefanProblemSolver, 2, pc, bc, ic, 0.25, tau, tMax, 100, tauSave);
            ijArray(taskNumber, :) = [i j];
        end
        if taskNumber == batchSize || (i == length(y) && j == length(x))
            for k = 1:taskNumber
                completedId = fetchNext(F);
                iTemp = ijArray(completedId, 1);
                jTemp = ijArray(completedId, 2);
                rfo.Results(iTemp, jTemp) = {F(completedId).OutputArguments};
                isPointCompleted(iTemp, jTemp) = true;
                rfo.isPointCompleted(iTemp, jTemp) = true;
            end
            pb.setProgress(((i-1)*length(x) + j), (length(x)*length(y)));
            %printProgressBar( ((i-1)*length(x) + j)/(length(x)*length(y)) );
            taskNumber = 0;
        end
    end
end

rfo.FinishDateTime = datetime();
rfo.isCompleted = true;

end

function createResultsFile(filename, InitDataFilename, x, y)
    isCompleted = false;
    isPointCompleted = false(length(y), length(x));
    Results = cell(length(y), length(x));
    CreateDateTime = datetime();
    FinishDateTime = [];
    save(filename, ...
        'isCompleted', ...
        'x', ...
        'y', ...
        'isPointCompleted', ...
        'Results', ...
        'CreateDateTime', ...
        'FinishDateTime', ...
        'InitDataFilename', ...
        '-v7.3', '-nocompression');
end

