function T = build_retina_table(ageFolders, ageLabels)
% This function aggregate retina-level metrics across multiple age groups into a tidy MATLAB Table.
%
% Inputs:
%   ageFolders - Cell array of strings/chars containing full directory paths for each age group.
%   ageLabels  - Numeric vector corresponding to the chronological age of each group.
%
% Example Input Format:
%   ageFolders = {
%       'D:\Data\MEA_recording\rd10_baselines\dayN'
%       'D:\Data\MEA_recording\rd10_baselines\dayN+1'
%   };
%   ageLabels = [23 45 60 90 120 150 200];
%
% Output:
%   T - A MATLAB table with columns: {'Retina', 'Age', 'Group', 'DeltaPeak'}
%       optimized for downstream statistical analysis (e.g., ANOVA, linear mixed models).
rows = {};

for ai = 1:numel(ageFolders)
    S = load(fullfile(ageFolders{ai}, 'retina_level_peak_metrics.mat'), 'retinaMetrics');
    R = S.retinaMetrics;

    for r = 1:numel(R)
        retinaID = R(r).prefix;

        groupNames = {'OnTrans','OnSus','OffTrans','OffSus'};
        for gi = 1:numel(groupNames)
            g = groupNames{gi};

            if isfield(R(r), g) && isfield(R(r).(g), 'ok') && R(r).(g).ok
                rows(end+1,:) = {
                    retinaID, ...                    % Retina ID
                    ageLabels(ai), ...               % Age
                    g, ...                           % Cell type group
                    R(r).(g).peakDelta               % Δpeak
                };
            end
        end
    end
end

T = cell2table(rows, ...
    'VariableNames', {'Retina','Age','Group','DeltaPeak'});

% Convert to categorical where appropriate
T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);

T.Retina = categorical(T.Retina);
T.Group  = categorical(T.Group);
T.Age    = categorical(T.Age);

disp(T.Properties.VariableNames)
end
