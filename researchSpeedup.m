clear;

wAr = 1:6;
numOfRuns = 5;
numOfTasks = 20;
times = zeros(numOfRuns, length(wAr));
delete(gcp('nocreate'));
pool = parpool(max(wAr));
fprintf("%20s%20s%20s\n","NumOfWorkers", "mean(Time), sec")
for i = 1:length(wAr)
    for j = i+1:max(wAr)
        f(j) = parfeval(pool, @pause, 0, inf);
    end
%     pool = parpool(wAr(i));
    times(:, i) = measureT(pool, numOfTasks, numOfRuns)';
%     delete(pool);
    for j = i+1:max(wAr)
        cancel(f(j));
    end
    temp = memory;
    fprintf("%20d%20.4f%20.4f\n", wAr(i), times(i));
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