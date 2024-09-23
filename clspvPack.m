%%
% class for a pack of PVs. The class is aimed to simplify the get  data 
% from archieved server. The orginal was in Python.
% FYWang, FCCA, 9/20/2024
%%

classdef clspvPack
    properties         
        pvNames                         % list of PVs
        alignSetup                      % alignment settings
        endTime                         % end time, in the format '12/11/2022 06:40:00'
        duration_hour                   % duration in hours
        sysDeltaTime_sec = 7 * 3600;    % EPICS time delay in seconds
        rawData                         % raw data from archiver, rawData{i} first column is time in seconds
        synData                         % synchronized data
    end
    
    methods
        % 0. Initialization
        function obj = clspvPack(pvNames, endTime, duration_hour, alignSetup)
            if nargin < 4
                alignSetup = struct('base_id', 1, 'val_range', [0, 5], 'disTimeAddBack_sec', 1, 'Trim', true);
                if nargin < 3
                    duration_hour = 1;
                    if nargin < 2
                        endTime = '05/19/2023 20:00:00';
                        if nargin < 1      
                            % default pv list
                            pvNames = [];
                        end
                    end
                end
            end
                    
            obj.pvNames = pvNames;
            obj.endTime = endTime;
            obj.duration_hour = duration_hour;
            obj.alignSetup = alignSetup;
            obj.rawData = {};
            obj.synData = {};
        end
        
        % 1. Extract history data of PVs at a given time window
        function obj = getHistory(obj)            
            startTime = datestr(datetime(obj.endTime) - obj.duration_hour/24, 'mm/dd/yyyy HH:MM:SS');
            timeRange = {startTime, obj.endTime};
            obj.rawData = {};            
            
            % get data by pvNames        
            % setup waitbar
            progress = 0;
            hWaitbar = waitbar(progress, ' Getting data ...'); % initialize the waitbar            
            for i = 1:length(obj.pvNames)
                [obj.rawData{end+1}, ~] = get_history(obj.pvNames{i}, timeRange);   
    
                progress = i/length(obj.pvNames);
                waitbar(progress, hWaitbar, sprintf('Getting data: %d%d', round(progress*100)));
            end            
            close(hWaitbar);
        end
        
        % 2. Set base channel to a given PV
        function obj = setBase(obj, pv, val_range, disTimeAddBack, Trim)
            if nargin > 1 && ~isempty(pv)
                idx = find(strcmp(obj.pvNames, pv));
                if ~isempty(idx)
                    obj.alignSetup.base_id = idx ; 
                else
                    disp(['The PV - ', pv, ' does not exist!']);
                    return;
                end
            end
            if nargin > 2 && ~isempty(val_range)    
                % or use default
                % multiple val_range could be set, like [[0, 1]; [-3, -1]; ...]
                obj.alignSetup.val_range = val_range;
            end
            if nargin > 3 && ~isempty(disTimeAddBack)   
                % or use default, 1 second
                obj.alignSetup.disTimeAddBack_sec = disTimeAddBack;
            end
            if nargin > 4   % or use default
                obj.alignSetup.Trim = Trim;
            end
        end
        
        % 3. Align history data in time by assigning a base channel
        function obj = alignHistory(obj, getHist)
            if nargin < 2
                getHist = false;
            end
            if isempty(obj.rawData)
                disp('No history data ...');
                if ~getHist
                    return;
                else
                    disp('Getting history data ...');
                    obj = obj.getHistory();
                end
            end
            
            synData = struct('absTime', [], 'relTime', [], 'vals', []);
            atime = obj.rawData{obj.alignSetup.base_id}(:,1); % Time reference
            
            if ~obj.alignSetup.Trim
                % no Trim
                time_cum = atime - atime(1); % Cumulated time
                idt = 1:length(atime);
            else
                % cut off out of range data
                idt = [];
                for k = 1:size(obj.alignSetup.val_range, 1)
                    idt_k = obj.rawData{obj.alignSetup.base_id}(:,2) >= obj.alignSetup.val_range(k, 1) & ...
                            obj.rawData{obj.alignSetup.base_id}(:,2) <= obj.alignSetup.val_range(k, 2);
                    idt = [idt, find(idt_k)];
                end
                % sort and remove duplicates
                idt = unique(sort(idt));
            end
            
            if isempty(idt)
                disp('No valid data!');
                return;
            end
            
            time_inv = [0; diff(atime(idt))];
            time_inv(time_inv > obj.alignSetup.disTimeAddBack_sec) = obj.alignSetup.disTimeAddBack_sec;
            time_cum = cumsum(time_inv);
            
            synData.absTime = atime(idt);
            synData.relTime = time_cum;
            
            % Align data
            for i = 1:length(obj.rawData)
                if i == obj.alignSetup.base_id
                    synData.vals(:,i) = obj.rawData{i}(idt, 2);
                else
                    synData.vals(:,i) = interp1(obj.rawData{i}(:,1), obj.rawData{i}(:,2), atime(idt));
                end
            end
            obj.synData = synData;
        end
        
        % 4. Plot data
        function pltHistory(obj, pvNames, raw_not_syn)
            if nargin < 2 || isempty(pvNames)
                pv_ids = 1:length(obj.pvNames);
                pvNames = obj.pvNames;
            else
                pv_ids = find(ismember(obj.pvNames, pvNames));
            end
            if nargin < 3
                raw_not_syn = true;
            end
            
            
            figure();
            hold on;
            
            if raw_not_syn
                if isempty(obj.rawData)
                    disp('No raw data available!');
                    return;
                end
                for i = pv_ids
                    plot(obj.rawData{i}(:,1), obj.rawData{i}(:,2));
                end
                xlabel('time [sec]');
                legend(pvNames);
                title(['raw data ',...
                     datestr(datenum(obj.endTime) - obj.duration_hour/24, ...
                    'mm/dd/yyyy HH:MM:SS' ), ...
                    '-', obj.endTime]);
            else
                if isempty(obj.synData)
                    disp('No synchronized data available!');
                    return;
                end
                plot(obj.synData.relTime, obj.synData.vals);
                xlabel('time [sec]');
                legend(pvNames);
                title(['synchronized data with disTimeAddBack = ', num2str(obj.alignSetup.disTimeAddBack_sec), ' s']);
            end
            grid on;
            hold off;
        end
         
        % 5. test function
        

    end

    methods(Static)
        function runTest()
            pvGun = {'GUN:GUNB:100:DFACT', 'GUN:GUNB:100:FWD:PWR', 'GUN:GUNB:100:CATHENDCAPCENTEMP', 'GUN:GUNB:100:REV1:PWR'};
            endTime = '09/28/2023 18:00:00';
            duration_hour = 10;        
    
            alignSet=struct('base_id', 2, 'val_range', [50, 100e3], 'disTimeAddBack_sec', 1, 'Trim',true);
    
            obj = clspvPack(pvGun, endTime, duration_hour, alignSet);
    
            % reset align channel to pvGun[1]
            % mpvPack.setBase(pvGun{1});    
    
            obj = obj.getHistory();
            obj.pltHistory();
            obj = obj.alignHistory();
            obj.pltHistory([],false); %plot aligned data

        end
    end
end

