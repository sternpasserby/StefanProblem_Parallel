clear;

wAr = [1 3];     % Для правильного графика ускорения обязательно надо, чтобы присутствовало значение 1
numOfRuns = 2;
numOfPoints = 20;
resultsFolder = "Speedup Data/" +...
    string( datestr(now, 'yy_mm_dd-HHMMSS') ) + "/";         % Папка для результатов исследования масштабирования
tempfolder = "temp";                                         % Папка, куда будут складываться временные результаты моделирования
initDataFilename = '2021_03_30 AntarcticaBM2_parsed.mat';    % Имя файла с исходными данными для моделирования
if isfolder(tempfolder)
    rmdir(tempfolder, 's');
end

mkdir(resultsFolder);

Np = [500 2000 500];
NpSave = [100 1000 100];
tMax = 500*365.25*24*3600;
tau = 3600*24*365.25/3;
tauSave = 3600*24*365.25*10;
NpBoundsSave = 100;

points_id = getPoints_id(numOfPoints, initDataFilename);

mex -largeArrayDims mex_TDMA.cpp

times = zeros(numOfRuns, length(wAr));
evalc( "delete(gcp('nocreate'))" );
%pool = parpool(max(wAr));
fprintf("\n%12s    %15s    %7s\n","NumOfWorkers", "mean(Time), sec", "Speedup")
for i = 1:length(wAr)
%     for j = max(wAr)-wAr(i):-1:1
%         f(j) = parfeval(pool, @pause, 0, inf);
%     end
    [~, pool] = evalc( sprintf("parpool(%d)", wAr(i)) );
    if isfile("mex_TDMA.mexw64")
        addAttachedFiles(pool, "mex_TDMA.mexw64");
    elseif isfile("mex_TDMA.mexa64")
        addAttachedFiles(pool, "mex_TDMA.mexa64");
    else
        error("Can't find compiled mex file!");
    end
    
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
%     for j = max(wAr)-wAr(i):-1:1
%         cancel(f(j));
%     end

    fprintf("%12d    %15.4f    %7.4f\n", wAr(i), mean(times(:, i)), mean(times(:, 1))/mean(times(:, i)));
end

figure
h = plot(wAr, mean(times), '-s');
xlabel('Number of Workers')
ylabel('Time, sec')
savePlot(h, resultsFolder + "Time");

figure
h = plot(wAr, mean(times(:, 1))./mean(times), '-s');
xlabel('Number of Workers')
ylabel('Speedup')
savePlot(h, resultsFolder + "Speedup");

save(resultsFolder + 'data.mat', 'times', 'wAr', 'numOfRuns', 'numOfPoints', 'Np',...
    'NpSave', 'tMax', 'tau', 'tauSave', 'NpBoundsSave', '-mat');

function savePlot(h, filename)
    savefig(filename);
    print(filename, '-dpng', '-r300');
    print(filename, '-depsc');
end

function points_id = getPoints_id(numOfPoints, initDataFilename)
    load(initDataFilename, 'Data');
    points_id = [];
    for i = 1:length(Data.X)
        bedrock = Data.Bedrock_m(i);
        iceSurf = Data.Surface_m(i);
        iceThickness = Data.IceThickness_m(i);
        %GHF = Data.GHF_Martos_mWm2(i);
        %accumRate = Data.AccumRate_kg1m2a1(i);

    %     if bedrock > 0 && (iceSurf - iceThickness ~= bedrock) % величины должны быть целыми, так что можно не исхищряться 
    %         continue;
    %     end
    %     if bedrock < 0 && iceSurf - iceThickness > 0
    %         continue;
    %     end
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