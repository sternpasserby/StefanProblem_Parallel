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
            progressFraction = n/N;
            numOfSegments = floor(progressFraction*obj.barLength);
            completedStr = repmat(obj.completedSymbol, 1, numOfSegments);
            todoStr = repmat(obj.todoSymbol, 1, obj.barLength - numOfSegments);
            msg = sprintf(...
                [obj.leftSymbol... 
                completedStr... 
                todoStr... 
                obj.rightSymbol... 
                ' %5.2f%% (%d/%d)'],... 
                progressFraction*100, n, N);
            disp([obj.reverseStr, msg]);
            obj.reverseStr = repmat(sprintf('\b'), 1, length(msg)+1);
        end
        
%         function delete(obj)
%             disp(obj.reverseStr);
%         end
    end
end

