function runGlacierModelling(pool, resFolderPath, initDataFilename, points_id, varargin)
%RUNGLACIERMODELLING Summary of this function goes here
%   Detailed explanation goes here

%%% ПАРСИНГ ВХОДНЫХ ПАРАМЕТРОВ
% Задание значений входных параметров по умолчанию
defaultNp = [100 1000 100];
defaultTau = 3600*24*14;
defaultTMax = 20*365.25*24*3600;
defaultGridType = 'Uniform';
defaultNpSave = [100 200 100];
defaultTauSave = defaultTau;
defaultmaxPartSize = 1024*1024*1024;     % Размер тома в байтах
defaultShowInfo = true;

% Настройка объекта parserObj типа InputParser
parserObj = inputParser;
parserObj.StructExpand = false;
addRequired(parserObj, 'pool');
addRequired(parserObj, 'resFolderPath');
addRequired(parserObj, 'initDataFilename');
addRequired(parserObj, 'points_id');
addParameter(parserObj, 'Np', defaultNp);
addParameter(parserObj, 'tau', defaultTau);
addParameter(parserObj, 'tMax', defaultTMax);
addParameter(parserObj, 'gridType', defaultGridType);
addParameter(parserObj, 'NpSave', defaultNpSave);
addParameter(parserObj, 'tauSave', defaultTauSave);
addParameter(parserObj, 'maxPartSize', defaultmaxPartSize);
addParameter(parserObj, 'showInfo', defaultShowInfo);

parse(parserObj, pool, resFolderPath, initDataFilename, points_id, varargin{:});
         Np = parserObj.Results.Np;
        tau = parserObj.Results.tau;
       tMax = parserObj.Results.tMax;
   gridType = parserObj.Results.gridType;
     NpSave = parserObj.Results.NpSave;
    tauSave = parserObj.Results.tauSave;
maxPartSize = parserObj.Results.maxPartSize;
   showInfo = parserObj.Results.showInfo;

%%% Создание папки и файла для результатов (или корректировка плана по частично выполненным расчётам)
partBaseName = "Data";
numOfPoints = length(points_id);
if isfolder(resFolderPath)
    if showInfo
        fprintf("Folder ""%s"" already exists. Loading completed points indices.\n", resFolderPath);
    end
    
    dirInfo = dir(resFolderPath + "/*.bin");
    numOfParts = length(dirInfo);
    partBaseName = string( dirInfo(1).name );
    partBaseName = extractBetween(partBaseName, 1, ...
        strlength(partBaseName) - 4 - strlength(regexp(partBaseName,'\d*','Match')) );
    
    dirName = resFolderPath + "/" + partBaseName + 1 + ".bin";
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
    
    if numOfParts > 1
        for j = 2:numOfParts
            dirName = resFolderPath + "/" + partBaseName + j + ".bin";
            fid = fopen(dirName, "rb");
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
        end
    end
    
    completedPoints_id(i:end) = [];
    for i = 1:length(completedPoints_id)
        id = find( points_id == completedPoints_id(i), 1 );
        if ~isempty(id)
            points_id(id) = [];
        end
    end
    numOfPoints = length(points_id);
else
    if showInfo
        fprintf("Folder ""%s"" does not exist. Creating folder.\n", resFolderPath);
    end
    mkdir(resFolderPath);
    numOfParts = 1;
    dirName = resFolderPath + "/" + partBaseName + numOfParts + ".bin";
    fid = fopen(dirName, "wb");
    fwrite(fid, numOfPoints, 'int');
    fwrite(fid, points_id, 'int');
    fclose(fid);
end
numOfPoints = length(points_id);
if numOfPoints == 0
    if showInfo
        fprintf("No more grid points to calculate! Aborting...\n");
    end
    return;
end

load(initDataFilename, 'Data');
batchSize = pool.NumWorkers;

pc = getPhysicalConstants();

%%% Смешанные краевые условия
% Формат краевых условий:
% alpha00*u1(0, t) + alpha01*du1/dx(0, t) = g0(t) - для левого конца
% alpha10*u2(L, t) + alpha11*du2/dx(L, t) = g1(t) - для правого конца
bc = struct;                  % bc - boundary conditions
bc.alpha = [0 -pc.lambda1; 1 0];

taskInd = 0;
taskInd2pInd = zeros(batchSize, 1);
batchSize = min(batchSize, numOfPoints);
dateStart = datetime(now,'ConvertFrom','datenum');
if showInfo
    fprintf('Progress: ');
    pb = ConsoleProgressBar();
    pb.setProgress( 0, numOfPoints );
end
for i = 1:batchSize
    k = points_id(i);

    % Граничные условия
    bc.g0 = @(t)(Data.GHF_Martos_mWm2(k)/1000);
    bc.g1 = @(t)( Data.T_Average_C(k) + 0.082/(10*365.25*24*3600)*t + Data.dT_Average_C(k)*sin(2*pi*t/31556952) + 273.15);
    %F(i) = parfeval(pool, @StefanProblemSolver, 2, pc, bc, ic, 0.25, tau, tMax, 100, tauSave);
    
    % Начальные условия
    s = [ Data.Bedrock_m(k); Data.Surface_m(k) - Data.IceThickness_m(k); Data.Surface_m(k); Data.Surface_m(k) ];
    x2 = linspace(s(2), s(3), Np(2));
    Uf_adj = (273.15 - 7.43*1e-8*pc.rho2*9.81*( s(3) - s(2) ));
    u2 = linspace(Uf_adj, bc.g1(0), Np(2));
    ic = struct('s', s, ...
            'dsdt', zeros(4, 1), ...
            'x1', linspace(s(1), s(2), Np(1)), ...
            'u1', 273.15 + zeros(Np(1), 1), ...
            'x2', x2, ...
            'u2', u2, ...
            'x3', linspace(s(3), s(4), Np(3)), ...
            'u3', 273.15 + ones(Np(3), 1), ...
            'tInit', 0);
    
    F(i) = parfeval(pool, @StefanProblemSolver, 2, pc, bc, ic, ... 
                                              'tau', tau, ...
                                              'tauSave', tauSave, ...
                                              'tMax', tMax, ...
                                              'Np', Np,...
                                              'gridType', gridType, ...
                                              'NpSave', NpSave, ...
                                              'accumRate', Data.AccumRate_kg1m2a1(k));
    taskInd2pInd(i) = k;
end

numOfCompPoints = 0;

fid = fopen(dirName, "ab");
for i = batchSize+1:numOfPoints + batchSize
    dirInfo = dir(dirName);
    if dirInfo.bytes >= maxPartSize
        fclose(fid);
        dirName = resFolderPath + "/" + partBaseName + numOfParts + ".bin";
        numOfParts = numOfParts + 1;
        fid = fopen(dirName, "ab");
    end
    
    % Получение результатов, запись их на диск
    taskInd = fetchNext(F);
    
    %fprintf("Elapsed time for point %d: %6.2f sec\n%", taskInd2pInd(taskInd), ...
    %    seconds(F(taskInd).FinishDateTime - F(taskInd).StartDateTime));
    
    k = taskInd2pInd(taskInd);
    L = length( F(taskInd).OutputArguments{2} );
    fwrite(fid, k, 'int');
    fwrite(fid, L, 'int');
%     fwrite(fid, L, 'int');
%     fwrite(fid, k*ones(5, L), 'double');
    fwrite(fid, F(taskInd).OutputArguments{2}, 'double');
    fwrite(fid, F(taskInd).OutputArguments{1}, 'double');
    
    numOfCompPoints = numOfCompPoints + 1;
    if showInfo
        pb.setProgress( numOfCompPoints, numOfPoints );
    end
    
    % Загрузка новых точек
    if i <= length(points_id)
        k = points_id(i);
        
        % Граничные условия
        bc.g0 = @(t)(Data.GHF_Martos_mWm2(k)/1000);
        bc.g1 = @(t)( Data.T_Average_C(k) + 0.082/(10*365.25*24*3600)*t + Data.dT_Average_C(k)*sin(2*pi*t/31556952) + 273.15);
        %F(i) = parfeval(pool, @StefanProblemSolver, 2, pc, bc, ic, 0.25, tau, tMax, 100, tauSave);
        
        % Начальные условия
        s = [ Data.Bedrock_m(k); Data.Surface_m(k) - Data.IceThickness_m(k); Data.Surface_m(k); Data.Surface_m(k) ];
        x2 = linspace(s(2), s(3), Np(2));
        Uf_adj = (273.15 - 7.43*1e-8*pc.rho2*9.81*( s(3) - s(2) ));
        u2 = linspace(Uf_adj, bc.g1(0), Np(2));
        ic = struct('s', s, ...
                'dsdt', zeros(4, 1), ...
                'x1', linspace(s(1), s(2), Np(1)), ...
                'u1', 273.15 + zeros(Np(1), 1), ...
                'x2', x2, ...
                'u2', u2, ...
                'x3', linspace(s(3), s(4), Np(3)), ...
                'u3', 273.15 + ones(Np(3), 1), ...
                'tInit', 0);
        
        F(taskInd) = parfeval(pool, @StefanProblemSolver, 2, pc, bc, ic, ...
                                                  'tau', tau, ...
                                                  'tauSave', tauSave, ...
                                                  'tMax', tMax, ...
                                                  'Np', Np,...
                                                  'gridType', gridType, ...
                                                  'NpSave', NpSave, ...
                                                  'accumRate', Data.AccumRate_kg1m2a1(k));
        taskInd2pInd(taskInd) = k;
    end
    
end
fclose(fid);

dateEnd = datetime(now,'ConvertFrom','datenum');
if showInfo
    fprintf("Elapsed time for glacier modelling: "); 
    disp(dateEnd - dateStart);
end

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
