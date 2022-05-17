clear;

initDataFilename = '2021_03_30 AntarcticaBM2_parsed.mat';
Np = [500 5000 500];
NpSave = [100 1000 100];
tMax = 1000*365.25*24*3600;
tau = 3600*24*365.25/3;
tauSave = 3600*24*365.25*10;

pool = gcp();
addAttachedFiles(pool, "mex_TDMA.mexa64");
load(initDataFilename, 'Data');
points_id = [];
for i = 1:length(Data.X)
%    if Data.Y(i) ~= -3000
%        continue
%    end
%     if mod(i, 30) ~= 0
%         continue
%     end
    
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
resFolderName = "Three";
runGlacierModelling(pool, parentDir + resFolderName, initDataFilename, points_id, ...
    'tau', tau, ...
    'tauSave', tauSave, ...
    'tMax', tMax, ...
    'Np', Np,...
    'gridType', 'SigmoidBased', ...
    'NpSave', NpSave, ...
    'showInfo', true);

