% Скрипт для запуска параллельных расчётов для Антарктиды с использованием
% трёхфазной задачи Стефана.

clear;

initDataFilename = '../AntarcticData/2021_03_30 AntarcticaBM2_parsed.mat';
Np = [500 5000 500];
NpSave = [100 1000 100];
tMax = 1000*365.25*24*3600;
tau = 3600*24*365.25/3;
tauSave = 3600*24*365.25*10;
NpBoundsSave = 100;

pool = gcp('nocreate');
if isempty(pool)
    pool = parpool(3);
end

% Если нет скомпилированного mex-файла, скомпилировать
if ~(isfile("mex_TDMA.mexw64") || isfile("mex_TDMA.mexa64"))
    mex -largeArrayDims mex_TDMA.cpp
end

% Присоединить mex файл к пулу воркеров, на всякий случай
if isfile("mex_TDMA.mexw64")
    addAttachedFiles(pool, "mex_TDMA.mexw64");
elseif isfile("mex_TDMA.mexa64")
    addAttachedFiles(pool, "mex_TDMA.mexa64");
else
    error("Can't find compiled mex file!");
end
    
Data = load(initDataFilename);

% Выбор точек грида данных, для которых будет проводится расчёт
points_id = [];
for i = 1:length(Data.X)
    
    %%% Первичное отсеивание, нужно для дебага. При настоящем расчёте эти
    %%% два оператора if надо закомментировать
    % Отсеять все точки кроме тех, у которых координата Y = -3000
    if Data.Y(i) ~= -3000
        continue
    end
    
    % Брать каждую тридцатую точку из отсеянных до этого момента
    if mod(i, 30) ~= 0
        continue
    end
    
    %%% Отсеивание всех точек, у которых между нижней кромкой льда и горной
    %%% породой есть пространство, и точек с нулевой толщиной льда
    bedrock = Data.Bedrock(i);
    iceSurf = Data.Surface(i);
    iceThickness = Data.IceThickness(i);
    if (iceSurf - iceThickness - bedrock ~= 0)
        continue
    end
    if iceThickness == 0
        continue
    end
    points_id(end+1) = i;
end

parentDir = "Results/";     % Папка, куда будут складываться папки с результатами расчётов
mkdir(parentDir);
resFolderName = "Five111";    % Имя папки для результатов расчёта
runGlacierModelling(pool, parentDir + resFolderName, initDataFilename, points_id, ...
    'tau', tau, ...
    'tauSave', tauSave, ...
    'tMax', tMax, ...
    'Np', Np,...
    'gridType', 'SigmoidBased', ...
    'NpSave', NpSave, ...
    'showInfo', true, ...
    'NpBoundsSave', NpBoundsSave);

