clear; close all;

folderName = "Results/Four/";
initDataFilename = '2021_03_30 AntarcticaBM2_parsed.mat';
mkdir(folderName + "Pics");

tEnd = getLastCommonTimeMoment(folderName);
[S0, S1, S2, S3, x, y] = getPhaseCoordinates(folderName, initDataFilename, 0);
[~, S1End, S2End, S3End, ~, ~] = getPhaseCoordinates(folderName, initDataFilename, tEnd);

% Получение вклада аккумуляции
load(initDataFilename, 'Data');
dAccum = NaN*zeros(size(S0));
dirInfo = dir(folderName + "/*.bin");
filename = string( dirInfo(1).name );
fid = fopen(folderName + filename, "rb");
N = fread(fid, 1, 'int');
points_id = fread(fid, [1, N], 'int');
x = unique(Data.X);
y = unique(Data.Y);
rho2 = 916.7;
for i = 1:N
    id = points_id(i);
    i1 = find( y == Data.Y(id) );
    j1 = find( x == Data.X(id) );
    
    dAccum(i1, j1) = Data.AccumRate_kg1m2a1(id)/rho2/(3600*24*365.25)*tEnd;
end
fclose(fid);

dH = (S1End - S0)*(1000/rho2-1);    % Вклад проседания

figure
h = imagesc([x(1) x(end)], [y(1) y(end)], 1000*(S1End - S0)/tEnd*3600*24*365.25 );
setupPlot(h, "Basal Melting Rate", "mm/year")
caxis([0 15])
savePlot(h, folderName + "Pics/" + "BasalMeltingRate");

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

figure
h = imagesc([x(1) x(end)], [y(1) y(end)], 1000*( S3End - dAccum - S3 + dH )/tEnd*3600*24*365.25 );
setupPlot(h, "Average Upper Melting Rate", "mm/year")
%caxis([0 1500])
savePlot(h, folderName + "Pics/" + "AverageUpperMeltingRate");

figure
h = imagesc([x(1) x(end)], [y(1) y(end)], S2End-S1End );
setupPlot(h, "End Ice Thickness", "m")
caxis([0 inf])
savePlot(h, folderName + "Pics/" + "EndIceThickness");

figure
h = imagesc([x(1) x(end)], [y(1) y(end)], (S2End - dAccum - S1End) - (S2 - S1) );
setupPlot(h, "Change in ice thickness (no accumulation)", "m")
%caxis([0 500])
savePlot(h, folderName + "Pics/" + "ChangeInIceThickness_noAccum");

function setupPlot(h, titleName, cbTitleName)
    colormap jet
    set(gca ,'YDir','normal') 
    set(h, 'AlphaData', ~isnan(h.CData))
    xlabel("X, m");
    ylabel("Y, m")
    title(titleName);
    cb=colorbar; cb.TickLabelInterpreter='latex'; title(cb, cbTitleName)
    axis equal
end

function savePlot(h, filename)
    savefig(filename)
    print(filename, '-dpng', '-r300')
end

function [S0, S1, S2, S3, x, y] = getPhaseCoordinates(folderName, initDataFilename, time)
    load(initDataFilename, 'Data');
    x = unique(Data.X);
    y = unique(Data.Y);
    
%     isS0 = false;
%     isS1 = false;
%     isS2 = false;
%     isS3 = false;
%     for i = 1:length(chosenVarsArray)
%         if chosenVarsArray(i) == "S0"
%             isS0 = true;
%         end
%         if chosenVarsArray(i) == "S1"
%             isS1 = true;
%         end
%         if chosenVarsArray(i) == "S2"
%             isS2 = true;
%         end
%         if chosenVarsArray(i) == "S3"
%             isS3 = true;
%         end
%     end
    
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
            A = fread(fid, [4, L], 'double');

            i1 = find(y == Data.Y(id));
            j1 = find(x == Data.X(id));
            
            S0(i1, j1) = interp1(t, A(1, :), time);
            S1(i1, j1) = interp1(t, A(2, :), time);
            S2(i1, j1) = interp1(t, A(3, :), time);
            S3(i1, j1) = interp1(t, A(4, :), time);
            
            %s1End = interp1(t, A(2, :), tVid(k));
            %image(i1, j1) = s1End - A(2, 1);
            
            readPoints = readPoints + 1;
        end
        pb.setProgress( readPoints, N );
        fclose(fid);
    end
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
            A = fread(fid, [4, L], 'double');
            
            if t(end) < res
                res = t(end);
            end

            readPoints = readPoints + 1;
        end
        pb.setProgress( readPoints, N );
        fclose(fid);
    end
end

function printArray(A)
    [numOfLines, ~] = size(A);
    for i = 1:numOfLines
        fprintf("%.4e ", A(i, :));
        fprintf("\n");
    end
end