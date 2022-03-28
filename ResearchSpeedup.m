clear;

wAr = [1 3 6];     % Для правильного графика ускорения обязательно надо, чтобы присутствовало значение 1
numOfRuns = 2;
numOfPoints = 20;
resFolderPath = "temp";                                      % Папка, куда будут складываться результаты моделирования
initDataFilename = '2021_03_30 AntarcticaBM2_parsed.mat';    % Имя файла с исходными данными для моделирования
if isfolder(resFolderPath)
    rmdir(resFolderPath, 's');
end

Np = [500 2000 500];
NpSave = [100 1000 100];
tMax = 500*365.25*24*3600;
tau = 3600*24*365.25/3;
tauSave = 3600*24*365.25*10;

points_id = getPoints_id(numOfPoints, initDataFilename);

times = zeros(numOfRuns, length(wAr));
delete(gcp('nocreate'));
%pool = parpool(max(wAr));
fprintf("%20s%20s\n","NumOfWorkers", "mean(Time), sec")
for i = 1:length(wAr)
%     for j = max(wAr)-wAr(i):-1:1
%         f(j) = parfeval(pool, @pause, 0, inf);
%     end
    pool = parpool(wAr(i));
    
    for j = 1:numOfRuns
        time = tic();
        runGlacierModelling(pool, resFolderPath, initDataFilename, points_id, ...
            'tau', tau, ...
            'tauSave', tauSave, ...
            'tMax', tMax, ...
            'Np', Np,...
            'gridType', 'SigmoidBased', ...
            'NpSave', NpSave, ...
            'showInfo', false);
        times(j, i) = toc(time);
        rmdir(resFolderPath, 's');
    end
    delete(pool);
%     for j = max(wAr)-wAr(i):-1:1
%         cancel(f(j));
%     end
    fprintf("%20d%20.4f\n", wAr(i), mean(times(:, i)));
end

figure
plot(wAr, mean(times), '-s')
xlabel('Number of Workers')
ylabel('Time, sec')

figure
plot(wAr, mean(times(:, 1))./mean(times), '-s')
xlabel('Number of Workers')
ylabel('Speedup')

save([datestr(now, 'yy_mm_dd-HHMMSS') '.mat'], 'times', 'wAr', 'numOfRuns',...
    'numOfPoints', '-mat');

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