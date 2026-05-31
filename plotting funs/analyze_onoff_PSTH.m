function analyze_onoff_PSTH(folderPath)
% This function summarize and analyze the ON/OFF PSTH responses of retinal neurons of mice in the same experimental group/day of age
%% Usage example:
% analyze_onoff_PSTH('D:/rd10_baseline/day N')
%% Input description:
% folderPath - the folder path containing *_totalneuronsV2.mat (cell classification) and the corresponding *_psth.mat (PSTH data)
%% Flow:
% Combine all neurons into four categories: On-Transient, On-Sustained, Off-Transient, and Off-Sustained.
% Calculate the population median and 25%/75% interquartile range, calculate the peak index at the retinal level and save it.
    if nargin < 1 || isempty(folderPath)
        folderPath = pwd;
    end

    % --- Peak / baseline windows in seconds ---
    win = struct();
    win.OnTrans.baseline  = [0, 2];   win.OnTrans.response  = [-2, 0];
    win.OnSus.baseline    = [0, 2];   win.OnSus.response    = [-2, 0];
    win.OffTrans.baseline = [-2, 0];   win.OffTrans.response = [ 0, 2];
    win.OffSus.baseline   = [-2, 0];   win.OffSus.response   = [ 0, 2];


    %======================================================================
    % 1. find all *_totalneuronsV2.mat files（each file corresponds to one
    % retina's data）
    %======================================================================
    v2Files = dir(fullfile(folderPath, '*_totalneuronsV2.mat'));
    if isempty(v2Files)
        error('In %s No *_totalneuronsV2.mat file is found', folderPath);
    end

    % prepare four response groups and add PSTHs from all the same-aged
    % mice retinas
    groupNamesList = {'OnTrans', 'OnSus', 'OffTrans', 'OffSus'};
    allGroupTraces = struct();  % Every retina：allGroupTraces.(gName) = [nNeurons x nBins]
    for gi = 1:numel(groupNamesList)
        allGroupTraces.(groupNamesList{gi}) = [];
    end
    time = [];
    binSize_ms = [];
    retinaMetrics = struct();  % one entry per retina
    retinaCount = 0;
    retinaCellStats = struct();


    %======================================================================
    % 2. go through every *_totalneuronsV2 + *_psth in current day folder 
    %======================================================================
    for f = 1:numel(v2Files)
        v2File = v2Files(f).name;
        [~, baseName] = fileparts(v2File);

        % "_totalneuronsV2" file
        prefix = erase(baseName, '_totalneuronsV2');
        
        retinaCount = retinaCount + 1;
        retinaMetrics(retinaCount).prefix = prefix;
        retinaCellStats(retinaCount).prefix = prefix;


        psthFileName = [prefix '_psth.mat'];
        psthFilePath = fullfile(folderPath, psthFileName);
        v2FilePath   = fullfile(folderPath, v2File);

        if ~exist(psthFilePath, 'file')
            warning('Unable to find matched PSTH file：%s，Skip this retina', psthFileName);
            continue;
        end

        fprintf('  Mouse prefix "%s":\n    v2:   %s\n    psth: %s\n', ...
                prefix, v2FilePath, psthFilePath);

        %------------------------------------------------------------------
        % 2.1  read totalneurons structure
        %------------------------------------------------------------------
        v2Data = load(v2FilePath, 'totalneurons');
        if ~isfield(v2Data, 'totalneurons')
            warning('In %s no valid totalneurons structure，skip.', v2FilePath);
            continue;
        end
        tn = v2Data.totalneurons;

        try
            on_trans_0  = tn.num_OnNeurons.trans;
            on_sus_0    = tn.num_OnNeurons.sus;
            off_trans_0 = tn.num_OffNeurons.trans;
            off_sus_0   = tn.num_OffNeurons.sus;
        catch
            warning(' In %s totalneurons structure unecpected, skip. ', v2FilePath);
            continue;
        end

        % 0-based -> 1-based MATLAB index
        on_trans_idx  = on_trans_0  + 1;
        on_sus_idx    = on_sus_0    + 1;
        off_trans_idx = off_trans_0 + 1;
        off_sus_idx   = off_sus_0   + 1;

        nOnTrans  = numel(on_trans_0);
        nOnSus    = numel(on_sus_0);
        nOffTrans = numel(off_trans_0);
        nOffSus   = numel(off_sus_0);
        
        nOn  = nOnTrans + nOnSus;
        nOff = nOffTrans + nOffSus;
        nAll_classified = nOn + nOff;

         all_idx_0 = tn.allneurons;  
         nAll_units = numel(all_idx_0);
         fracOn_inAll         = nOn  / nAll_units;
         fracOff_inAll        = nOff / nAll_units;       % OFF perc

        retinaCellStats(retinaCount).nOnTrans = nOnTrans;
        retinaCellStats(retinaCount).nOnSus   = nOnSus;
        retinaCellStats(retinaCount).nOffTrans= nOffTrans;
        retinaCellStats(retinaCount).nOffSus  = nOffSus;
        retinaCellStats(retinaCount).nOn  = nOn;
        retinaCellStats(retinaCount).nOff = nOff;
        retinaCellStats(retinaCount).nAll_classified = nAll_classified;
        retinaCellStats(retinaCount).nAll_units      = nAll_units;
        retinaCellStats(retinaCount).fracOn_inAll = fracOn_inAll;
        retinaCellStats(retinaCount).fracOff_inAll        = fracOff_inAll;



        %------------------------------------------------------------------
        % 2.2 read pltcurve in psth file  
        %------------------------------------------------------------------
        psthData = load(psthFilePath, 'pltcurve');  
        if ~isfield(psthData, 'pltcurve')
            warning('In %s no pltcurve structure，skip. ', psthFilePath);
            continue;
        end

        trialPSTHs = psthData.pltcurve.trialPSTHs;
        edges      = psthData.pltcurve.PSTH_binEdges;
        edges      = unwrap_edges(edges);

        if numel(edges) < 2
            warning('The length of PSTH_binEdges is not enough，skip %s。', psthFilePath);
            continue;
        end

        % bin size（s）and convert to ms
        this_binSize_sec = edges(2) - edges(1);
        this_binSize_ms  = this_binSize_sec * 1000;
        this_time        = edges(1:end-1) + this_binSize_sec/2;

        % set up global time / binSize
        if isempty(time)
            time       = this_time;
            binSize_ms = this_binSize_ms;
        else
            % sanity check：bin number consistent
            if numel(this_time) ~= numel(time)
                warning('In file %s, PSTH bin number inconsistent，skip the retina', psthFilePath);
                continue;
            end
        end

        nNeuronsTotal = numel(trialPSTHs);

        % define index in four categories (in current retina data)
        groups = struct();
        groups.OnTrans.name  = 'On-Transient';
        groups.OnTrans.idx   = on_trans_idx;
        groups.OnSus.name    = 'On-Sustained';
        groups.OnSus.idx     = on_sus_idx;
        groups.OffTrans.name = 'Off-Transient';
        groups.OffTrans.idx  = off_trans_idx;
        groups.OffSus.name   = 'Off-Sustained';
        groups.OffSus.idx    = off_sus_idx;

        %------------------------------------------------------------------
        % 2.3 for the retina， trial-average PSTH 
        %------------------------------------------------------------------
        factor = 1000 / binSize_ms;  % counts/bin → spikes/s

        for gi = 1:numel(groupNamesList)
            gName = groupNamesList{gi};
            gInfo = groups.(gName);

            idx = gInfo.idx(:)';  % vectorize
            % prevent index exceeding limit
            idx = idx(idx >= 1 & idx <= nNeuronsTotal);

            if isempty(idx)
                continue;
            end

            % preallication: neuron numhers in one retina 
            % every neuron corresponds to 1 × nBins (spikes/s)
            nLocalNeurons = numel(idx);
            nBins = numel(time);
            localTraces = nan(nLocalNeurons, nBins);

            for ii = 1:nLocalNeurons
                nid = idx(ii);
                psth_i = trialPSTHs{nid}{2};  % nTrials × nBins

                if size(psth_i, 2) ~= nBins
                    warning('Neuron %d in file %s has inconsistent bin number，skip this neuron.', ...
                             nid, psthFilePath);
                    continue;
                end

                trialAvg = mean(psth_i, 1, 'omitnan');  % counts/bin
                localTraces(ii, :) = trialAvg * factor; % spikes/s
            end

            % Add traces into retina group
            allGroupTraces.(gName) = [allGroupTraces.(gName); localTraces];
            % ---- NEW: retina-level metric for statistics ----
            validRows = all(~isnan(localTraces), 2);
            localTraces2 = localTraces(validRows, :);
            if size(localTraces2,1) >= 10
                medTrace = median(localTraces2, 1, 'omitnan');
                mm = compute_peak_metric(time, medTrace, win.(gName));

                retinaMetrics(retinaCount).(gName).ok        = true;
                retinaMetrics(retinaCount).(gName).nNeurons  = size(localTraces2,1);
                retinaMetrics(retinaCount).(gName).peakDelta = mm.peakDelta; % <-- use this for age comparisons
                retinaMetrics(retinaCount).(gName).peakAbs   = mm.peakAbs;
                retinaMetrics(retinaCount).(gName).baseMean  = mm.baseMean;
            else
                retinaMetrics(retinaCount).(gName).ok = false;
            end
        end

        clear psthData trialPSTHs;
    end

    %======================================================================
    % 3. After processing all the files -> Calculate median + 25/75%
    % and plot images based on the four reponse type
    %======================================================================
    for gi = 1:numel(groupNamesList)
        gName = groupNamesList{gi};
        traces = allGroupTraces.(gName);  % nNeuronsAll × nBins
        nNeuronsAll = size(traces, 1);
        if nNeuronsAll > 10
            traces_used = traces;
        else
            traces_used = [];
        end

        if isempty(traces_used)
            fprintf('Group %s in %s has np valid/enough neuron。\n', gName, folderPath);
            continue;
        end

        medianTrace = median(traces_used, 1, 'omitnan');
        p25         = prctile(traces_used, 25, 1);
        p75         = prctile(traces_used, 75, 1);

        % average neuron's mean firing rate
        meanRateEach = mean(traces_used, 2, 'omitnan');

        % plotting
        fig = figure('Name', sprintf('%s - %s', folderPath, gName), ...
                     'Color', 'w');
        hold on;

        x = time(:)';
        y1 = p25;
        y2 = p75;
        % shadow
        fill([x fliplr(x)], [y1 fliplr(y2)], [0.8 0.8 0.8], ...
             'EdgeColor', 'none', 'FaceAlpha', 0.4);
        % median trace
        plot(time, medianTrace, 'k', 'LineWidth', 2);

        yL = [0 45];
        %plot([0 0], yL, ':k', 'LineWidth', 1);
        ylim(yL);
        xlim([-2 0]);

        %xlabel('Time (s)');
        %ylabel('Firing rate (spikes/s)');
        nNeuronsAll = size(traces_used,1);
        title(sprintf('%s (n=%d)', gName, nNeuronsAll));
        box off;

        % save as either PNG or TIF
        figBaseName = sprintf('%s_%s_PSTH_median', ...
                              get_folder_shortname(folderPath), gName);
        %exportgraphics(fig,fullfile(folderPath, [figBaseName '.pdf']),'ContentType','vector');
        close(fig);

        %{
        % save stats in meanRateEach
        stats.meanRateEach_spkPerSec = meanRateEach;
        stats.medianTrace_spkPerSec  = medianTrace;
        stats.p25_spkPerSec          = p25;
        stats.p75_spkPerSec          = p75;
        %save(fullfile(folderPath, [figBaseName '_stats.mat']), 'stats');
        %}
        save(fullfile(folderPath, 'retina_level_peak_metrics.mat'), ...
     'retinaMetrics', 'win', 'time', 'binSize_ms');
        save(fullfile(folderPath,'retina_level_cell_fraction.mat'), 'retinaCellStats');

    end
end
function edges = unwrap_edges(edgesIn)
    edges = edgesIn;
    if isstruct(edges)
        fns = fieldnames(edges);
        if ~isempty(fns)
            edges = edges.(fns{1});
        else
            error('PSTH_binEdges is an empty struct.');
        end
    end
    if iscell(edges)
        edges = edges{1};
    end
    if ~isvector(edges)
        error('PSTH_binEdges is not vectorized。');
    end
    edges = edges(:)';   % converted to vector
end

% convert D:/rd10_baseline/day N to dayN
function name = get_folder_shortname(folderPath)
    [~, name] = fileparts(folderPath);
end
function mm = compute_peak_metric(time, trace, win)
    baseMask = time >= win.baseline(1) & time <= win.baseline(2);
    respMask = time >= win.response(1)  & time <= win.response(2);

    baseMean = mean(trace(baseMask), 'omitnan');
    peakAbs  = max(trace(respMask), [], 'omitnan');

    mm.baseMean  = baseMean;
    mm.peakAbs   = peakAbs;
    mm.peakDelta = peakAbs - baseMean;
end
