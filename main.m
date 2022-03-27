clear;

pool = gcp();
load('2021_03_30 AntarcticaBM2_parsed.mat', 'Data');
points_id = [];
for i = 1:length(Data.X)
    if Data.Y(i) ~= -3000
        continue
    end
    if mod(i, 30) ~= 0
        continue
    end
    
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
end

%parentDir = "Results\\" + datestr(now, 'yy_mm_dd-HHMMSS') + "\\";
parentDir = "Results/";
mkdir(parentDir);
resFolderName = "One";
runGlacierModelling(pool, parentDir + resFolderName, '2021_03_30 AntarcticaBM2_parsed.mat', points_id)
% 
% tMax = 200*365.25*24*3600;
% tauSave = 3600*24*365.25;
% load(resFilename, 'Results', 'pointIndices');
% t = 0:tauSave:tMax;
% S = cell(size(Results));
% for i = 1:length(S)
%     S{i} = interp1(Results{i}{2}, Results{i}{1}', t);
% end

% vid = VideoWriter("vid");
% open(vid);
% x = Data.X(pointIndices);
% s = zeros(length(S), 1);
% f = figure('WindowState','maximized');
% for i = 1:length(t)
%     for j = 1:length(S)
%         s(j) = S{j}(i, 2) - S{j}(i, 1);
%     end
%     plot(x, s, '-o', 'LineWidth', 2);
%     legend("Relative subglacial water level")
%     xlabel("X, m")
%     ylabel("Water level, m")
%     %bar(x, s, 1);
%     axis([-inf inf 0 1.7]);
%     set(gca, 'FontSize', 20)
%     drawnow();
%     title(sprintf("Time: %12.2f year", t(i)/3600/24/365.25));
%     frame = getframe(gcf());
%     writeVideo(vid, frame);
% end
% close(vid);