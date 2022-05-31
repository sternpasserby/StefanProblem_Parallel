% Скрипт для исследование масштабируемости параллельного алгоритма
% моделирований для всей Антарктиды, реалированного в функции runGlacierModelling.m

clear;

wAr = [1 3];         % Для правильного графика ускорения обязательно надо, чтобы присутствовало значение 1
numOfRuns = 2;       % Сколько раз измерить время для каждого числа воркеров
numOfPoints = 20;    % Число точек грида исходных данных, для которых будет запускаться параллельный расчёт

resultsFolder = "Speedup Data/" +...
    string( datestr(now, 'yy_mm_dd-HHMMSS') ) + "/";         % Папка для результатов исследования масштабирования
tempfolder = "temp";                                         % Папка, куда будут складываться временные результаты моделирования
initDataFilename = '../AntarcticData/2021_03_30 AntarcticaBM2_parsed.mat';    % Имя файла с исходными данными для моделирования

if isfolder(tempfolder)
    rmdir(tempfolder, 's');
end

mkdir(resultsFolder);

% Параметры моделирования для каждой точки
Np = [500 2000 500];
NpSave = [100 1000 100];
tMax = 500*365.25*24*3600;
tau = 3600*24*365.25/3;
tauSave = 3600*24*365.25*10;
NpBoundsSave = 100;

points_id = getPoints_id(numOfPoints, initDataFilename);

% Если нет скомпилированного mex-файла, скомпилировать
if ~(isfile("mex_TDMA.mexw64") || isfile("mex_TDMA.mexa64"))
    mex -largeArrayDims mex_TDMA.cpp
end

times = zeros(numOfRuns, length(wAr));
evalc( "delete(gcp('nocreate'))" );     % Здесь и далее некоторые функции используются через evalc, чтобы подавить их вывод в консоль
fprintf("\n%12s    %15s    %7s\n","NumOfWorkers", "mean(Time), sec", "Speedup")
for i = 1:length(wAr)
    
    % Запустить пул воркеров
    [~, pool] = evalc( sprintf("parpool(%d)", wAr(i)) );
    if isfile("mex_TDMA.mexw64")
        addAttachedFiles(pool, "mex_TDMA.mexw64");
    elseif isfile("mex_TDMA.mexa64")
        addAttachedFiles(pool, "mex_TDMA.mexa64");
    else
        error("Can't find compiled mex file!");
    end
    
    % Измерение времени моделирования numOfRuns раз
    for j = 1:numOfRuns
        time = tic();
        runGlacierModelling(pool, tempfolder, initDataFilename, points_id, ...
            'tau', tau, ...
            'tauSave', tauSave, ...
            'tMax', tMax, ...
            'Np', Np,...
            'gridType', 'SigmoidBased', ...
            'NpSave', NpSave, ...
            'showInfo', false, ...
            'NpBoundsSave', NpBoundsSave);
        times(j, i) = toc(time);
        rmdir(tempfolder, 's');
    end
    evalc( "delete(pool)" );

    fprintf("%12d    %15.4f    %7.4f\n", wAr(i), mean(times(:, i)), mean(times(:, 1))/mean(times(:, i)));
end

%%% Построение графиков и сохранение их в папку с результатами исследования
%%% масштабирования вместе с замерами времени
figure
plot(wAr, mean(times), '-s');
xlabel('Number of Workers')
ylabel('Time, sec')
savePlot(gca, resultsFolder + "Time");

figure
plot(wAr, mean(times(:, 1))./mean(times), '-s');
xlabel('Number of Workers')
ylabel('Speedup')
savePlot(gca, resultsFolder + "Speedup");

save(resultsFolder + 'data.mat', 'times', 'wAr', 'numOfRuns', 'numOfPoints', 'Np',...
    'NpSave', 'tMax', 'tau', 'tauSave', 'NpBoundsSave', '-mat');

function savePlot(h, filename)
    savefig(h.Parent, filename);
    print(h.Parent, filename, '-dpng', '-r300');
    print(h.Parent, filename, '-depsc');
end

% Функция для выбора точек грида исходных данных, на которых будет
% запускаться параллельный расчёт
function points_id = getPoints_id(numOfPoints, initDataFilename)
    Data = load(initDataFilename);
    points_id = [];
    for i = 1:length(Data.X)
        
        % Отсеять все точки, где толщина льда 0 или где между нижней
        % кромкой льда и горными породами есть пространство
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
        if length(points_id) == numOfPoints
            break;
        end
    end
end