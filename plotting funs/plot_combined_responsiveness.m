function plot_combined_responsiveness(baseFolder)
    % get all the files starting with days
    dayFolders = dir(fullfile(baseFolder, 'day*'));
    validDays = [];
    allResults = {};
    
    % go through every day file
    for i = 1:length(dayFolders)
        dayPath = fullfile(baseFolder, dayFolders(i).name);
        
        % look for _distance.mat files
        matFiles = dir(fullfile(dayPath, '*_distance.mat'));
        
        % need at least three mat files
        if length(matFiles) >= 2
            dayResults = struct();
            dayData = [];
            
            % combine all the data within the day file
            for j = 1:length(matFiles)
                filePath = fullfile(dayPath, matFiles(j).name);
                try
                    N1 = load(filePath);
                    N = N1.analysisResults;
                    loadedData.analysisResults.centre = N.centre;
                    loadedData.analysisResults.distances = N.distances;
                    loadedData.analysisResults.bin_edges = N.bin_edges;
                    loadedData.analysisResults.bin_centers = N.bin_centers;
                    loadedData.analysisResults.response_percentage = N.response_percentage;
                    loadedData.analysisResults.bin_stats = N.bin_stats;
                    loadedData.analysisResults.neuronData = N.neuronData;
                    dayData = [dayData; loadedData.analysisResults];
                catch ME
                    warning('unable to load file %s: %s', matFiles(j).name, ME.message);
                end
            end
            
            % calculate stats for the day
            if ~isempty(dayData)
                dayNum = str2double(regexp(dayFolders(i).name, '\d+', 'match'));
                
                % normalise the day bin range
                all_bin_edges = arrayfun(@(x) x.bin_edges, dayData, 'UniformOutput', false);
                common_bins = get_common_bins(all_bin_edges);
                
                % calculate means and stds for every bin
                [mean_resp, std_resp, bin_centers] = calculate_day_stats(dayData, common_bins);
                
                % save results
                dayResults.day = dayNum;
                dayResults.mean_resp = mean_resp;
                dayResults.std_resp = std_resp;
                dayResults.bin_centers = bin_centers;
                dayResults.n = length(dayData);
                dayResults.color = rand(1,3); % assign colour for each day
                
                allResults{end+1} = dayResults;
                validDays = [validDays, dayNum];
            end
        end
    end
    
    % sort by days' values
    [~, sortIdx] = sort(validDays);
    allResults = allResults(sortIdx);
    
    % plotting results
    if ~isempty(allResults)
        plot_responsiveness_curves(allResults, baseFolder);
    else
        warning('unable to find enough data number（every day requires ≥3 distance.mat files）');
    end
end

%% functions to normalise the bin range
function common_bins = get_common_bins(all_bin_edges)
    % find all the bin range
    min_edge = min(cellfun(@min, all_bin_edges));
    max_edge = max(cellfun(@max, all_bin_edges));
    
    % use the finest bin size
    bin_widths = cellfun(@(x) x(2)-x(1), all_bin_edges);
    finest_width = min(bin_widths);
    
    common_bins = min_edge:finest_width:max_edge;
    if common_bins(end) < max_edge
        common_bins = [common_bins, common_bins(end)+finest_width];
    end
end

%% calculate the means and stds for every day
function [mean_resp, std_resp, bin_centers] = calculate_day_stats(dayData, common_bins)
    bin_centers = common_bins(1:end-1) + diff(common_bins)/2;
    nBins = length(bin_centers);
    nFiles = length(dayData);
    
    % initization
    all_resp = nan(nFiles, nBins);
    
    % align every day's data into the common 
    for i = 1:nFiles
        [~, bin_idx] = histc(dayData(i).bin_centers, common_bins);
        valid_bins = bin_idx > 0 & bin_idx <= nBins;
        all_resp(i, bin_idx(valid_bins)) = dayData(i).response_percentage(valid_bins);
    end
    
    % stats
    mean_resp = nanmean(all_resp, 1);
    std_resp = nanstd(all_resp, 0, 1);
end

%% make response curve
function plot_responsiveness_curves(allResults, savePath)
    figure('Position', [100 100 1000 600], 'Color', 'w');
    hold on;
    
    % pre-define colours
    colors = lines(length(allResults));
    
    % plot every line
    legendEntries = cell(1, length(allResults));
    for i = 1:length(allResults)
        dayData = allResults{i};
        valid_bins = ~isnan(dayData.mean_resp);
        
        % make error bar for std and mean
        errorbar(dayData.bin_centers(valid_bins), ...
                dayData.mean_resp(valid_bins), ...
                dayData.std_resp(valid_bins), ...
                'o-', 'Color', colors(i,:), ...
                'LineWidth', 1.5, 'MarkerSize', 6, ...
                'MarkerFaceColor', colors(i,:), ...
                'CapSize', 5);
        
        legendEntries{i} = sprintf('Day %d (n=%d)', dayData.day, dayData.n);
    end
    
   
    xlabel('Distance from center (μm)', 'FontSize', 12);
    ylabel('Responsive Percentage (%)', 'FontSize', 12);
    title('Responsiveness vs Distance Across Days', 'FontSize', 14);
    legend(legendEntries, 'Location', 'bestoutside', 'FontSize', 10);
    grid on;
    box on;
    
    % set axies
   lims = [min(cellfun(@(x) min(x.bin_centers), allResults)), ...
         max(cellfun(@(x) max(x.bin_centers), allResults))];
    xlim(lims);
    y1 = max(cellfun(@(x) max(x.mean_resp), allResults));
    y2 = max(cellfun(@(x) max(x.std_resp), allResults));
    ylim([0 ceil(y1+y2)]);
    
    % save image
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end
    saveName = fullfile(savePath, 'multiday_responsiveness.png');
    saveas(gcf, saveName);
    fprintf('image has been saved to: %s\n', saveName);
    close(gcf);
end