clear; close all;
addpath( genpath('..\ForCoastlinePlotting') ); 

folderName = "Results/Four/";
initDataFilename = '2021_03_30 AntarcticaBM2_parsed.mat';

sec2years = 1/(3600*24*365.25);

fprintf("\n=============== Read glacier modelling results ===============\n");

if ~exist(folderName + "Pics", 'dir')
    mkdir(folderName + "Pics");
end

fprintf("Finding last common time moment... ");
tEnd = getLastCommonTimeMoment(folderName);
fprintf("Done. tEnd = %.6e sec (%.2f years)\n", tEnd, tEnd * sec2years);

fprintf("Getting phase boundaries at t = 0... ");
[S0, S1, S2, S3, x, y] = getPhaseCoordinates(folderName, initDataFilename, 0);
fprintf("Done.\n")

% Получение скорости аккумуляции (в м/с)
AccumRate = getAccumulationSpeed(initDataFilename, folderName, size(S0));

fprintf("Getting phase boundaries at t = %.2f years... ", tEnd * sec2years);
[~, S1End, S2End, S3End, ~, ~] = getPhaseCoordinates(folderName, initDataFilename, tEnd);
fprintf("Done.\n")
dS1 = S1End - S1;
fprintf("Volume of melted water after %.2f years: %.4e km^3\n", tEnd*sec2years, ...
    100*sum( dS1( ~isnan(dS1) )/1000, 'all') );      % шаг грида данных 10 км, отсюда площадь на точку данных 100 км2
fprintf("Number of glacier points: %d\n", numel( S0( ~isnan(S0) ) ) );

% dH = (S1End - S0)*(1000/rho2-1);    % Вклад проседания

% % Скорость донного таяния
% figure
% h = imagesc([x(1) x(end)], [y(1) y(end)], 1000*(S1End - S0)/(tEnd*sec2years) );
% setupPlot(h, "Basal Melting Rate", "mm/year")
% plotAntarcticCoastLines(h.Parent, "Color", "k", "LineWidth", 0.5);
% caxis([0 15])
% savePlot(h, folderName + "Pics/" + "BasalMeltingRate");
% 
% % Скорость приповерхностного таяния
% figure
% h = imagesc([x(1) x(end)], [y(1) y(end)], 1000*( S3End - S2End - (S3 - S2) )/sec2years );
% setupPlot(h, "Average Upper Melting Rate", "mm/year")
% plotAntarcticCoastLines(h.Parent, "Color", "k", "LineWidth", 0.5);
% %caxis([0 1500])
% savePlot(h, folderName + "Pics/" + "AverageUpperMeltingRate");
% 
% % Конечная толщина льда
% figure
% h = imagesc([x(1) x(end)], [y(1) y(end)], S2End-S1End );
% setupPlot(h, "End Ice Thickness", "m")
% plotAntarcticCoastLines(h.Parent, "Color", "k", "LineWidth", 0.5);
% caxis([0 inf])
% savePlot(h, folderName + "Pics/" + "EndIceThickness");
% 
% % Изменение толщины льда с вычетом вклада аккумуляции
% figure
% h = imagesc([x(1) x(end)], [y(1) y(end)], (S2End - AccumRate*tEnd - S1End) - (S2 - S1) );
% setupPlot(h, sprintf("Change in ice thickness\n(no accumulation)"), "m")
% plotAntarcticCoastLines(h.Parent, "Color", "k", "LineWidth", 0.5);
% %caxis([0 500])
% savePlot(h, folderName + "Pics/" + "ChangeInIceThickness_noAccum");
% 
% % Изменения координаты верхней кромки ледника с вычетом вклада аккумуляции
% figure
% h = imagesc([x(1) x(end)], [y(1) y(end)], S2End - AccumRate*tEnd - S2 );
% setupPlot(h, sprintf("Change in surface elevation\n(no accumulation)"), "m")
% plotAntarcticCoastLines(h.Parent, "Color", "k", "LineWidth", 0.5);
% %caxis([0 500])
% savePlot(h, folderName + "Pics/" + "ChangeInSurfElevation_noAccum");

%%% Генерация картинок в разные моменты времени
t = 0:250/sec2years:tEnd; t(1) = [];

if ~exist(folderName + "Pics/Dynamics", 'dir')
    mkdir(folderName + "Pics/Dynamics");
end
if ~exist(folderName + "Pics/Dynamics/BasalMelt", 'dir')
    mkdir(folderName + "Pics/Dynamics/BasalMelt");
end
if ~exist(folderName + "Pics/Dynamics/ChangeInIceThickness_noAccum", 'dir')
    mkdir(folderName + "Pics/Dynamics/ChangeInIceThickness_noAccum");
end
if ~exist(folderName + "Pics/Dynamics/ChangeInSurfElevation_noAccum", 'dir')
    mkdir(folderName + "Pics/Dynamics/ChangeInSurfElevation_noAccum");
end

figure; ax1 = gca;
figure; ax2 = gca;
figure; ax3 = gca;
figure; ax4 = gca;
% vMelted = zeros(size(t));

for i = 1:length(t)
    fprintf("Getting phase boundaries at t = %.2f years... ", t(i) * sec2years);
    [~, S1End, S2End, S3End, ~, ~] = getPhaseCoordinates(folderName, initDataFilename, t(i));
    fprintf("Done.\n")
    
    basalMelt = S1End - S0;
%     vMelted(i) = 100*sum( basalMelt( ~isnan(basalMelt) ), 'all' );
    
    % Донное таяние
    h1 = imagesc(ax1, [x(1) x(end)], [y(1) y(end)], basalMelt );
    setupPlot(h1, sprintf("Basal melting\nt = %.2f years", t(i)*sec2years), "m")
    caxis(ax1, [0 15])
    plotAntarcticCoastLines(h1.Parent, "Color", "k", "LineWidth", 0.5);
    savePlot(h1, folderName + "Pics/Dynamics/BasalMelt/" + "BasalMelt" + i);

    % Изменение толщины льда с вычетом вклада аккумуляции
    h2 = imagesc(ax2, [x(1) x(end)], [y(1) y(end)], (S2End - AccumRate*t(i) - S1End) - (S2 - S1) );
    setupPlot(h2, sprintf("Change in ice thickness, t = %.2f years\n(no accumulation)", t(i)*sec2years), "m")
    caxis(ax2, [-17 0])
    plotAntarcticCoastLines(h2.Parent, "Color", "k", "LineWidth", 0.5);
    savePlot(h2, folderName + "Pics/Dynamics/ChangeInIceThickness_noAccum/" + "ChangeInIceThickness_noAccum" + i);
    
%     % Толщина льда
%     h3 = imagesc(ax3, [x(1) x(end)], [y(1) y(end)], S2End-S1End );
%     setupPlot(h3, sprintf("End Ice Thickness, t = %.2f years", t(i)*sec2years), "m")
%     caxis(ax3, [0 5000])
%     plotAntarcticCoastLines(h3.Parent, "Color", "k", "LineWidth", 0.5);
%     savePlot(h3, folderName + "Pics/Dynamics/IceThickness/" + "IceThickness" + i);
    
    % Изменения координаты верхней кромки ледника с вычетом вклада аккумуляции
    h4 = imagesc(ax4, [x(1) x(end)], [y(1) y(end)], S2End - AccumRate*t(i) - S2 );
    setupPlot(h4, sprintf("Change in surface elevation, t = %.2f years\n(no accumulation)", t(i)*sec2years), "m")
    caxis(ax4, [-1.5 0])
    plotAntarcticCoastLines(h4.Parent, "Color", "k", "LineWidth", 0.5);
    savePlot(h4, folderName + "Pics/Dynamics/ChangeInSurfElevation_noAccum/" + "ChangeInSurfElevation_noAccum" + i);
end

% figure
% plot(t*sec2years, vMelted*1e-6, '-o')
% ax = gca;
% ax.XLabel.String = "Time, years";
% ax.YLabel.String = "V of melted water under ice, $km^3$";

function setupPlot(h, titleName, cbTitleName)
    ax = h.Parent;
    ax.Colormap = colormap(jet);
    %colormap jet
    %ax.YDir = 'normal';
    set(ax ,'YDir','normal') 
    set(h, 'AlphaData', ~isnan(h.CData))
    ax.XLabel.String = "X, m";
    ax.YLabel.String = "Y, m";
    ax.Title.String = titleName;
    cb=colorbar(ax); cb.TickLabelInterpreter='latex'; cb.Title.String = cbTitleName; 
        cb.Title.Interpreter = "Latex";
    axis(ax,'equal'); % axis equal
end

function savePlot(h, filename)
    ax = h.Parent;
    savefig(ax.Parent, filename)
    print(ax.Parent, filename, '-dpng', '-r300')
    print(ax.Parent, filename, '-depsc')
end

function [S0, S1, S2, S3, x, y] = getPhaseCoordinates(folderName, initDataFilename, time)
    load(initDataFilename, 'Data');
    X = Data.X;
    Y = Data.Y;
    clear Data;
    
    x = unique(X);
    y = unique(Y);
    
    numOfRows = length(y);
    numOfCols = length(x);
    S0 = zeros(numOfRows, numOfCols)*NaN;
    S1 = zeros(numOfRows, numOfCols)*NaN;
    S2 = zeros(numOfRows, numOfCols)*NaN;
    S3 = zeros(numOfRows, numOfCols)*NaN;
    
    % Получение информации о томах
    dirInfo = dir(folderName + "/*.bin");
    numOfParts = length(dirInfo);
    partBaseName = string( dirInfo(1).name );
    partBaseName = extractBetween(partBaseName, 1, ...
        strlength(partBaseName) - 4 - strlength(regexp(partBaseName,'\d*','Match')) );
    
    readPoints = 0;
    pb = ConsoleProgressBar();
    for i = 1:numOfParts
        filename = folderName + partBaseName + i + ".bin";
        if i == 1
            fid = fopen(filename, "rb");
            N = fread(fid, 1, 'int');
            points_id = fread(fid, [1, N], 'int');
        else
            fid = fopen(filename, "rb");
        end

        while true
            id = fread(fid, 1, 'int');
            if isempty(id)
                break
            end
            L = fread(fid, 1, 'int');
            t = fread(fid, [1, L], 'double');
            
            % idx - номер интервала на сетке t ( подразумевается нумерация, начинающаяся с единицы)
            [~, idx] = min( abs(t - time) );
            %idx = min(idx, L-1);
            if t(idx) > time
                idx = idx-1;
            end
            if time == t(end)
                idx = idx-1;
            end
            
            % Далее в файле идёт матрица 4 на L. Строки - координаты границ
            % s_i в моменты времени, которые хранятся в массиве t
            fseek(fid, 32*(idx-1), 0);  % 32, потому что 4 элемента в столбце, и каждый занимает 8 байт
            A = fread(fid, [4, 2], 'double');
            fseek(fid, 32*( L-idx-1 ), 0);
            
            %A = fread(fid, [4, L], 'double');

            i1 = find(y == Y(id));
            j1 = find(x == X(id));
            
            % Линейная интерполяция
            temp = (time - t(idx)) / ( t(idx+1)-t(idx) );
            S0(i1, j1) = ( A(1, 2) - A(1, 1) )*temp + A(1, 1);
            S1(i1, j1) = ( A(2, 2) - A(2, 1) )*temp + A(2, 1);
            S2(i1, j1) = ( A(3, 2) - A(3, 1) )*temp + A(3, 1);
            S3(i1, j1) = ( A(4, 2) - A(4, 1) )*temp + A(4, 1);
            
%             S0(i1, j1) = interp1(t(idx:idx+1), A(1, :), time);
%             S1(i1, j1) = interp1(t(idx:idx+1), A(2, :), time);
%             S2(i1, j1) = interp1(t(idx:idx+1), A(3, :), time);
%             S3(i1, j1) = interp1(t(idx:idx+1), A(4, :), time);
            
            readPoints = readPoints + 1;
            if mod(readPoints, 5000) == 0
                pb.setProgress( readPoints, N );
            end
        end
        fclose(fid);
    end
    delete(pb);
end

function res = getLastCommonTimeMoment(folderName)

    % Получение информации о томах
    dirInfo = dir(folderName + "/*.bin");
    numOfParts = length(dirInfo);
    partBaseName = string( dirInfo(1).name );
    partBaseName = extractBetween(partBaseName, 1, ...
        strlength(partBaseName) - 4 - strlength(regexp(partBaseName,'\d*','Match')) );
    
    readPoints = 0;
    pb = ConsoleProgressBar();
    res = inf;
    for i = 1:numOfParts
        filename = folderName + partBaseName + i + ".bin";
        if i == 1
            fid = fopen(filename, "rb");
            N = fread(fid, 1, 'int');
            points_id = fread(fid, [1, N], 'int');
        else
            fid = fopen(filename, "rb");
        end

        while true
            id = fread(fid, 1, 'int');
            if isempty(id)
                break
            end
            L = fread(fid, 1, 'int');
            t = fread(fid, [1, L], 'double');
            fseek(fid, 32*L, 0);
            %A = fread(fid, [4, L], 'double');*uint8
            %A = fread(fid, [4, L], 'double');
            
            if t(end) < res
                res = t(end);
            end

            readPoints = readPoints + 1;
            if mod(readPoints, 5000) == 0
                pb.setProgress( readPoints, N );
            end
        end
        fclose(fid);
    end
end

function AccumRate = getAccumulationSpeed(initDataFilename, folderName, size)
    AccumRate = NaN*zeros( size );
    
    dirInfo = dir(folderName + "/*.bin");
    filename = string( dirInfo(1).name );
    fid = fopen(folderName + filename, "rb");
    N = fread(fid, 1, 'int');
    points_id = fread(fid, [1, N], 'int');
    
    load(initDataFilename, 'Data');
    x = unique(Data.X);
    y = unique(Data.Y);
    rho2 = 916.7;
    temp = 1/rho2/(3600*24*365.25);
    for i = 1:N
        id = points_id(i);
        i1 = find( y == Data.Y(id) );
        j1 = find( x == Data.X(id) );

        AccumRate(i1, j1) = Data.AccumRate_kg1m2a1(id)*temp;
    end
    fclose(fid);
    clear Data 
end

function printArray(A)
    [numOfLines, ~] = size(A);
    for i = 1:numOfLines
        fprintf("%.4e ", A(i, :));
        fprintf("\n");
    end
end



% vw = VideoWriter('newfile.avi');
% vw.FrameRate = 1;
% open(vw)
% h = imagesc([x(1) x(end)], [y(1) y(end)], S3-S2 );
% setupPlot(h, sprintf("S3(t) - S2(t), t  = %d years", 0), "m")
% axis tight manual 
% set(gca,'nextplot','replacechildren'); 
% t = linspace(0, tEnd, 20);
% pb = ConsoleProgressBar();
% for i = 1:length(t)
%     [~, ~, S2, S3, ~, ~] = getPhaseCoordinates(folderName, initDataFilename, t(i));
%     imagesc([x(1) x(end)], [y(1) y(end)], S3-S2 );
%     title( sprintf( "S3(t) - S2(t), t  = %d years", floor(t(i)/(3600*24*365.25)) ) );
%     writeVideo(vw, getframe(gcf) );
%     pb.setProgress( i, length(t) );
% end
% close(vw)