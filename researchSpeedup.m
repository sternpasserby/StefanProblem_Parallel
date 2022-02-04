clear;

wAr = 1:6;
numOfRuns = 5;
numOfTasks = 10;
times = zeros(numOfRuns, length(wAr));
%fprintf("%20s%20s\n","NumOfWorkers", "Time, sec")
for i = 1:length(wAr)
    pool = parpool(wAr(i));
    times(:, i) = measureT(pool, numOfTasks, numOfRuns)';
    %fprintf("%20d%20.4f\n", wAr(i), times(i));
    delete(pool);
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
    'numOfTasks', '-mat');