classdef ConsoleProgressBar < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        barLength = 25;
        
        leftSymbol = '[';
        rightSymbol = ']';
        completedSymbol = '=';
        todoSymbol = ' ';
        reverseStr = '';
        
        setProgressExecutedFirstTime = true;
        startTick = 0;
        elapsedTime = 0;
    end
    
    methods
        function obj = ConsoleProgressBar()
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            %obj.Property1 = inputArg1 + inputArg2;
        end
        
        function setProgress(obj, n, N)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if obj.setProgressExecutedFirstTime == true
                obj.setProgressExecutedFirstTime = false;
                obj.startTick = tic;
            end
            obj.elapsedTime = toc(obj.startTick);
            
            progressFraction = n/N;
            
            numOfSegments = floor(progressFraction*obj.barLength);
            completedStr = repmat(obj.completedSymbol, 1, numOfSegments);
            todoStr = repmat(obj.todoSymbol, 1, obj.barLength - numOfSegments);
            msg = sprintf(...
                [obj.leftSymbol... 
                completedStr... 
                todoStr... 
                obj.rightSymbol... 
                ' %5.2f%% (%d/%d)\n'...
                'Elapsed time: %.2f sec\n'...
                '         ETA: %.2f sec'],... 
                progressFraction*100, n, N, obj.elapsedTime, obj.elapsedTime*(N/n-1) );
            
            disp([obj.reverseStr, msg]);
            obj.reverseStr = repmat(sprintf('\b'), 1, length(msg)+1 );
        end
        
%         function delete(obj)
%             fprintf([obj.reverseStr ]);
%         end
    end
end

