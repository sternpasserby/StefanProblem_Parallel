clear;

folderName = "Results\One\";

dirInfo = dir(folderName + "\\*.bin");
numOfParts = length(dirInfo);
partBaseName = string( dirInfo(1).name );
partBaseName = extractBetween(partBaseName, 1, ...
    strlength(partBaseName) - 4 - strlength(regexp(partBaseName,'\d*','Match')) );

readPoints = 0;
for i = 1:numOfParts
    filename = folderName + partBaseName + i + ".bin";
    if i == 1
        fid = fopen(filename, "rb");
        N = fread(fid, 1, 'int'); fprintf("%d\n", N);
        points_id = fread(fid, [1, N], 'int'); fprintf("%d ", points_id); fprintf("\n");
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

        fprintf("%d\n%d\n", id, L);
        printArray(t);
        printArray(A)
        readPoints = readPoints + 1;
    end
    
    fclose(fid);
end
fprintf("Number of points read: %d\n", readPoints);


function printArray(A)
    [numOfLines, ~] = size(A);
    for i = 1:numOfLines
        fprintf("%.4e ", A(i, :));
        fprintf("\n");
    end
end