% 
% Class for handling a pack of PVs (Process Variables). 
% The class is designed to simplify data retrieval from an archived server. 
% Originally implemented in Python, this version is adapted for MATLAB.
%
% properties:
%   pvNames         : Cell array of PV names (Process Variables)
%   webServer       : String, Web server for fetching data
%   alignSetup      : Structure containing settings for data alignment
%   endTime         : End time for data collection ('MM/DD/YYYY HH:MM:SS')
%   duration_hour   : Duration of data collection window (in hours)
%   rawData = {}    : Cell array to store raw data from the server (time, values)
%   synData = {}    : Cell array to store synchronized/aligned data
%
% Methods: 
% 1. getHistory   : Fetch history data of PVs from the archive.
% 2. setBase      : Set a base channel for data alignment.
% 3. alignHistory : Align history data in time to the assigned base channel.
% 4. pltHistory   : Plot history data, both raw and synchronized.
% 5. getInds      : Get the index of PVs in the list of PV names.
% 
% Author: FYWang, FCCA, 9/20/2024
%%

classdef clsPVpack
    properties
        pvNames         % Cell array of PV names (Process Variables)
        webServer       % String, Web server for fetching data
        alignSetup      % Structure containing settings for data alignment
        endTime         % End time for data collection ('MM/DD/YYYY HH:MM:SS')
        duration_hour   % Duration of data collection window (in hours)
        rawData = {}    % Cell array to store raw data from the server (time, values)
        synData = {}    % Cell array to store synchronized/aligned data
    end

    methods
        % 0. Initialization: Constructor for the class
        function obj = clsPVpack(varargin)            
            % Constructor initializes the class with optional input parameters
            %
            % Parameters:
            % - pvNames       : List of PV names
            % - webServer     : Server for fetching data
            % - endTime       : End time for data collection
            % - duration_hour : Duration of the data window in hours
            % - alignSetup    : Structure for alignment settings
            %
            % Example:
            % obj = clsPVpack('pvNames', {PV_list}, 'endTime', '12/11/2022 06:40:00', 'duration_hour', 6)

            
            p = inputParser;

            % default values LCLSII gun
            pvNames = {'GUN:GUNB:100:FWD:PWR', ...
                       'GUN:GUNB:100:DFACT', ...
                       'GUN:GUNB:100:REV1:PWR'};
            endTime = '06/05/2023 08:08:08';
            duration_hour = 6;
            alignSetup = struct('base_id', 1, ...               % index of PV for the base
                                'base_pv', pvNames{1}, ...
                                'val_range', [1e3, 1e5], ...    % valid range of the base [[a, b]; [c, d]; ...]
                                'disTimeAddBack_sec', 1, ...    % add back time interval (s) for discontinuety in base 
                                'dtResample_sec', 1, ...        % resample time interval in seconds
                                'Trim', true);                  % cut off out of range data or not. true for CUT, faluse for NO cut
  
            % Parameter validation
            addParameter(p, 'pvNames', pvNames);
            addParameter(p, 'webServer', 'LCLS');
            addParameter(p, 'endTime', endTime);
            addParameter(p, 'duration_hour', duration_hour);         
            addParameter(p, 'alignSetup', alignSetup);
            
            % Parse input arguments
            parse(p, varargin{:});

            % Assign parsed values to the object
            obj.pvNames = p.Results.pvNames;
            obj.webServer = p.Results.webServer;
            obj.endTime = p.Results.endTime;
            obj.duration_hour = p.Results.duration_hour;
            obj.alignSetup = p.Results.alignSetup;
            obj.rawData = {};
            obj.synData = {};
        end
        
        % 1. Fetch historical data for PVs over a given time window
        function obj = getHistory(obj)
            % Fetch historical data for the specified PVs within a time range
            %
            % Example:
            % obj = obj.getHistory(); 

            startTime = datestr(datetime(obj.endTime, ...
                                'InputFormat','MM/dd/yyyy HH:mm:ss') ...
                              - obj.duration_hour/24,  ...
                                'mm/dd/yyyy HH:MM:SS');
            timeRange = {startTime, obj.endTime};
            obj.rawData = {};         
                        
            % get data by pvNames        
            % setup waitbar
            progress = 0;
            % initialize the waitbar 
            hWaitbar = waitbar(progress, ' Getting data ...', ...
                               'Name', 'Get History Data');            
            for i = 1:length(obj.pvNames)
                try
                    [obj.rawData{end+1}, ~] = get_history(...
                        obj.pvNames{i}, timeRange, obj.webServer);   
                catch
                    warning(['The PV {' obj.pvNames{i}  ...
                             '} is not valided!']);
                    obj.rawData{end+1} = [];
                end
                %update progress bar
                progress = i/length(obj.pvNames);
                waitbar(progress, hWaitbar, ...
                    sprintf('[%s]: %d %%', obj.pvNames{i}, ...
                                                    round(progress*100)));
            end            
            close(hWaitbar);
        end
               
        % 2. Set base channel to a given PV, or id according to pvNames
        function obj = setBase(obj, varargin)
            %pv, val_range, disTimeAddBack, Trim)

            p = inputParser;
            % define expected parameters
            addParameter(p, 'base_pv', [], @isstr);                 % base_pv: should be a string
            addParameter(p, 'base_id', [], @isnumeric);             % base_id should be integer
            addParameter(p, 'val_range', [], @isnumeric);           % val_range: shoulde be numeric 
            addParameter(p, 'disTimeAddBack_sec', [], @isnumeric);  % disTimeAddBack should be numeric
            addParameter(p, 'dtResample_sec', [], @isnumeric);      % dtResample_sec should be numeric
            addParameter(p, 'Trim', [], @islogical);                % Trim should be logical

            % Parse the input
            parse(p, varargin{:});

            % Process 'base_pv' parameter
            if ~isempty(p.Results.base_pv)
                idx = find(strcmp(obj.pvNames, p.Results.base_pv));
                if ~isempty(idx)
                    obj.alignSetup.base_id = idx;
                    obj.alignSetup.base_pv = p.Results.base_pv;
                else
                    disp(['The PV - ', pv, ' does not exist!']);
                    return;
                end
            end

            % Process 'base_id' if base_id is used 
            if ~isempty(p.Results.base_id)
                if p.Results.base_id >=1 && p.Results.base_id <= length(obj.pvNames)
                    obj.alignSetup.base_id = p.Results.base_id;
                    obj.alignSetup.base_pv = obj.pvNames{p.Results.base_id};
                else
                    error('Error. base_id is invalid. ');
                end
            end

            % Process 'val_range' parameter
            if ~isempty(p.Results.val_range)
                if size(p.Results.val_range, 2) == 2
                    % multiple val_range could be set, 
                    % like [[0, 1]; [-3, -1]; ...]
                    obj.alignSetup.val_range = p.Results.val_range;
                else
                    error('Error. val_range: dimension err.')
                end
            end
    
            % Process 'disTimeAddBack_sec' parameter
            if ~isempty(p.Results.disTimeAddBack_sec)
                obj.alignSetup.disTimeAddBack_sec = ...
                    p.Results.disTimeAddBack_sec;
            end

            % Process 'dtResample_sec' parameter
             if ~isempty(p.Results.dtResample_sec)
                obj.alignSetup.dtResample_sec = ...
                    p.Results.drResample_sec;
            end
    
            % Process 'Trim' parameter
            if ~isempty(p.Results.Trim)
                obj.alignSetup.Trim = p.Results.Trim;
            end

        end
        
        % 3. Align history data in time to the assigned base channel
        function obj = alignHistory(obj, varargin)
            % Syntax:
            % obj = alignHistory('getHistory', true)

            p = inputParser;
            % Parameter validation
            addParameter(p, 'getHistory', true, @islogical); %getHistort should be logical

            % Parse input arguments
            parse(p, varargin{:});

            % check if there is history data or not
            
            if isempty(obj.rawData)
                disp('No history data ...');
                if ~p.Results.getHistory
                    return;
                else
                    disp('Getting history data ...');
                    obj = obj.getHistory();
                end
            end
            
            % align data

            obj.synData = struct('startTime', [], 'relTime', [], 'vals', []);
            timeRange = obj.rawData{obj.alignSetup.base_id}(:,1); % Time reference

            % converte to relative time
            timeRange = (timeRange - timeRange(1))*24*3600; % in seconds
           
            if ~obj.alignSetup.Trim
                % no Trim
                time_cum = timeRange - timeRange(1); % Cumulated time
                idt = 1:length(timeRange);
            else
                % cut off out of range data
                idt = [];       % in range data points
                for k = 1:size(obj.alignSetup.val_range, 1)
                    idt_k = obj.rawData{obj.alignSetup.base_id}(:,2) >= obj.alignSetup.val_range(k, 1) & ...
                            obj.rawData{obj.alignSetup.base_id}(:,2) <= obj.alignSetup.val_range(k, 2);
                    idt = [idt, find(idt_k)];
                end
                % sort and remove duplicates
                idt = unique(sort(idt));
            end
            
            if isempty(idt)
                disp('All data are out of range!');
                return;
            end
            
            % ajust time for cut off
            diff_idt = 1 + find(diff(idt) > 1);   % find where it was cut
            time_inv = [0; diff(timeRange(idt))]; % reltime in seconds
            time_inv(diff_idt) = obj.alignSetup.disTimeAddBack_sec;
            time_cum = cumsum(time_inv);
            
            if isempty(diff_idt)
                obj.synData.startTime = datestr(obj.rawData{obj.alignSetup.base_id}(1,1));
            else
                obj.synData.startTime = datestr(obj.rawData{obj.alignSetup.base_id}(diff_idt,1));
            end

            obj.synData.relTime = time_cum;
            
            % Align data to the baseID with nearest data point
            for i = 1:length(obj.rawData)
                if i == obj.alignSetup.base_id
                    obj.synData.vals(:,i) = obj.rawData{i}(idt, 2);
                else
                    if isscalar(obj.rawData{i}(:,1))
                        % only one data point -- bad 
                        obj.synData.vals(:, i) = obj.rawData{i}(1,2);
                    else
                        obj.synData.vals(:,i) = interp1((obj.rawData{i}(:,1) - obj.rawData{i}(1,1))*24*3600, ...
                        obj.rawData{i}(:,2), timeRange(idt), 'nearest');
                    end
                end
            end

            % resample the data with time interval = reSample_sec
            reSample = obj.synData.relTime(1):obj.alignSetup.dtResample_sec:obj.synData.relTime(end);
            
            obj.synData.vals = interp1(obj.synData.relTime, obj.synData.vals, reSample, 'nearest');
            obj.synData.relTime = reSample;

            disp(['Total time of data in hours: ' num2str(time_cum(end)/3600, '%.1f')]);

        end
        
        % 4. Plot data
        function pltHistory(obj, varargin)
            % pltHisory('pvNames', pvNames, 'raw_not_syn', true_or_false,'figNum', figNum, 'norm', true_or_false)
            p = inputParser;

            % default
            addParameter(p, 'pvNames', obj.pvNames);
            addParameter(p, 'raw_not_syn', true, @islogical);       % true for raw data
            addParameter(p, 'figNum', 10, @isnumeric);              % figure number
            addParameter(p, 'norm', false);                         % normlize the data for plot

            parse(p, varargin{:});

            pv_ids = find(ismember(obj.pvNames, p.Results.pvNames));

            if isempty(pv_ids)
                disp(['The PVs are invalid. +++ ' p.Results.pvNames]);
                return;
            end
           
            if p.Results.raw_not_syn
                if isempty(obj.rawData)
                    disp('No raw data available!');
                    return;
                end
                figure(p.Results.figNum);       
                %hFig = findobj('Type', 'figure', 'Number', p.Results.figNum);
                
                for i = pv_ids
                    if p.Results.norm
                        % normalization
                         stairs((obj.rawData{i}(:,1)-obj.rawData{i}(1,1))*24, ...
                               (obj.rawData{i}(:,2)-min(obj.rawData{i}(:,2)))/...
                               (1e-9+max(obj.rawData{i}(:,2))-min(obj.rawData{i}(:,2))));
                    else
                        % not 
                        stairs((obj.rawData{i}(:,1)-obj.rawData{i}(1,1))*24, obj.rawData{i}(:,2));
                    end
                    hold on;
                end
                grid on;
                hold off;
                xlabel('time [hour]');
                legend(p.Results.pvNames, 'location','best');
                title(['raw data starts @ ',...
                     datestr(datenum(obj.endTime) - obj.duration_hour/24, ...
                    'mm/dd/yyyy HH:MM:SS' ), ...
                    ' --- for ', num2str(obj.duration_hour), ' hours.']);
            else
                if isempty(obj.synData)
                    error('No synchronized data available!');
                else
                    if isempty(obj.synData.vals)
                        error('No synchronized data available!');
                    end
                end
                figure(p.Results.figNum);
                
                if p.Results.norm
                    % norm
                    plot(obj.synData.relTime , ...
                        (obj.synData.vals(:, pv_ids) -  ...
                        min(obj.synData.vals(:,pv_ids)))./...
                        (1e-9+max(obj.synData.vals(:,pv_ids)) - ...
                        min(obj.synData.vals(:,pv_ids))));
                else
                    % not
                    plot(obj.synData.relTime/3600, obj.synData.vals(:, pv_ids));
                end
                grid on;
                xlabel('time [hour]');
                legend(p.Results.pvNames);
                title(['synchronized data with disTimeAddBack = ', num2str(obj.alignSetup.disTimeAddBack_sec), ' s']);
            end
            
        end
         
        % 5. get index of pvs
        function pv_ids = getInds(obj, pvNames)
             pv_ids = find(ismember(obj.pvNames, pvNames));
        end
        % 6. moving window scatter plot 
        function mwScatter(obj, varargin)
            % Scatter plots of pvs_x, pvs_y at moving time window mw_hour
            % from begining of the data till the end of data. The plots are
            % group in figue (figNum) with subplots defined by num_row and
            % num_col.
            % data.mwScatter('pvs_x', data.pvNames(1), 'pvs_y', data.pvNames(2:3), 'num_row', 1, 'mw_hour', 1, 'strLegend', {'ref', 'df'})

            % default values
            mw_hour = 8;    
            pvs_x = {};     % Cell array to store pv list for x-axis
            scx = 1;        % scale factor for pvs_x
            scy = 1;        % scale factor for pvs_y
            pvs_y = {};     % Cell array to store pv list for y-axis
            num_row = 1;    % row num. for subplot
            num_col = [];   % col num. for subplot
            figNum = 19;    % figure num.
            sctOpt = {};    % cell for scatter plot options
                        
            p = inputParser;

            % load the default parameters
            addParameter(p, 'mw_hour', mw_hour, @(x) isnumeric(x) && x>0); 
            addParameter(p, 'pvs_x', pvs_x);
            addParameter(p, 'pvs_y', pvs_y);
            addParameter(p, 'scx', scx, @isnumeric);
            addParameter(p, 'scy', scy, @isnumeric);
            addParameter(p, 'num_row', num_row, @(x) isnumeric(x) && x>0 && x==round(x));
            addParameter(p, 'num_col', num_col, @(x) isnumeric(x) && x>0 && x==round(x));
            addParameter(p, 'figNum', figNum, @(x) isnumeric(x) && x>=0 && x==round(x));
            addParameter(p, 'sctOpt', sctOpt, @iscell);
            addParameter(p, 'xlabel', []);
            addParameter(p, 'ylabel', []);
            addParameter(p, 'legend', {});

            % Parse the input
            parse(p, varargin{:});
             
            % process parameters
            % pv_x
            if isempty(p.Results.pvs_x) 
                error('pvs_x has no pvs!');                
            else
                pvs_x_id = getInds(obj, p.Results.pvs_x);
                if length(pvs_x_id) < length(pvs_x)
                    error('pvs_x has invalid pvs!');
                end
            end
            % pv_y
            if isempty(p.Results.pvs_y) 
                error('pvs_y has no pvs!');
            else
                pvs_y_id = getInds(obj, p.Results.pvs_y);
                if length(pvs_y_id) < length(pvs_y)
                    error('pvs_y has invalid pvs!');
                 end
            end
            %
            if length(p.Results.scx) > 1 && length(p.Results.scx) ~= length(p.Results.pvs_x)
                error('scaling factor x do not match with pvs_x!');                
            end

            if length(p.Results.scy) > 1 && length(p.Results.scy) ~= length(p.Results.pvs_y)
                error('scaling factor x do not match with pvs_x!');                
            end

            % synData
            if isempty(obj.synData)
                error('No synchronized data!');
            else
                if isempty(obj.synData.vals)
                    error('No synchronized data!');
                end
            end

            % Np: how many plots
            Np = ceil(obj.synData.relTime(end)/3600/p.Results.mw_hour);

            if isempty(p.Results.num_col)
                num_col = ceil(Np/p.Results.num_row);
            else
                num_col = p.Results.num_col;
            end

            tstart = 0;
            tend = p.Results.mw_hour*3600;

            figure(p.Results.figNum);

            for plt = 1:Np
                if plt > p.Results.num_row*p.Results.num_col
                    warning('Number of plots are not enough for all data!')
                    break;
                end
                % do plots
                ind = obj.synData.relTime >= tstart & obj.synData.relTime < tend;
                subplot(p.Results.num_row, num_col, plt);
                scatter(obj.synData.vals(ind, pvs_x_id).*p.Results.scx, ...
                        obj.synData.vals(ind, pvs_y_id).*p.Results.scy, p.Results.sctOpt{:});
                grid on;
                title([num2str(tstart/3600) ' -- ' num2str(tend/3600) 'hrs']);

                if ~isempty(p.Results.xlabel)
                    xlabel(p.Results.xlabel);
                end

                if ~isempty(p.Results.ylabel)
                    ylabel(p.Results.ylabel);
                end
                
                if ~isempty(p.Results.legend)
                    legend(p.Results.legend{:});
                end
                
                tstart = tend;
                tend = tstart + p.Results.mw_hour*3600;
            end
            
            %
            
        end

    end

    methods(Static)
        function testRun()
            pvGun = {'GUN:GUNB:100:FWD:PWR', ...
                     'GUN:GUNB:100:DFACT', ...
                     'GUN:GUNB:100:REV1:PWR'};
            webServer = 'LCLS';
            endTime = '09/28/2023 18:00:00';
            duration_hour = 24;        
    
            alignSet=struct('base_id', 1, 'val_range', [50, 100e3], ...
                            'dtResample_sec', 1, ...
                            'disTimeAddBack_sec', 1, 'Trim',true);
    
            obj = clsPVpack('pvNames', pvGun, 'webServer', webServer, ...
                            'endTIme', endTime, ...
                            'duration_hour', ...
                            duration_hour, 'alignSet', alignSet);
    
            % reset align channel to pvGun[1]
            % mpvPack.setBase(pvGun{1});    
    
            obj = obj.getHistory();          
            obj = obj.alignHistory();

            %plot raw data,
            obj.pltHistory('raw_not_syn', true, 'figNum', 101);   

            %plot aligned data, and normalized
            obj.pltHistory('raw_not_syn', false, 'norm', true, 'figNum', 102);

            % moving window scatter plot
            obj.mwScatter('pvs_x', pvGun(1), 'pvs_y', pvGun(2), ...
                          'scx', 1e-3, 'scy', 1e-3, 'mw_hour', 2, ...
                          'num_row', 2, 'figNum', 103, ...
                          'sctOpt',{40, '+'}, ...
                          'xlabel','FwdPwr [kW]', 'ylabel', 'df [kHz]');
            
        end
    end
end

