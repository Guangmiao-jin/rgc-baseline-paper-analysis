function Tf = collect_fraction_all(ageFolders, ageLabels, denominatorType)
% This function aggregates raw RGC subtype counts and assigns the total 
% trial size for subsequent Binomial GLMM regression.
%
% Input Arguments:
%   ageFolders      - Cell array of directory paths for each age group.
%   ageLabels       - Numeric vector of chronological ages.
%   denominatorType - String/Char specifying the 'BinomialSize' baseline:
%                     'all_units'        -> Total recorded MEA units (nAllCells = nAll_units)
%                     'classified_only' -> Total light-responsive cells (nAllCells = Sum of 4 subtypes)

    % Default to total MEA units if not specified
    if nargin < 3 || isempty(denominatorType)
        denominatorType = 'all_units'; 
    end

    rows = {};

    for ai = 1:numel(ageFolders)
        statsFile = fullfile(ageFolders{ai}, 'retina_level_cell_fraction.mat');
        if ~exist(statsFile, 'file'), continue; end
        
        S = load(statsFile, 'retinaCellStats');
        R = S.retinaCellStats;
        
        for r = 1:numel(R)
            % Extract raw baseline counts
            nOnTrans  = double(R(r).nOnTrans);
            nOnSus    = double(R(r).nOnSus);
            nOffTrans = double(R(r).nOffTrans);
            nOffSus   = double(R(r).nOffSus);
            nAllUnits = double(R(r).nAll_units);
            
            % Determine what the total trials (BinomialSize) should be
            switch lower(denominatorType)
                case 'all_units'
                    totalTrials = nAllUnits;
                case 'classified_only'
                    totalTrials = nOnTrans + nOnSus + nOffTrans + nOffSus;
                otherwise
                    error('Invalid denominatorType! Use ''all_units'' or ''classified_only''.');
            end

            % Append row only if there are valid cells to model
            if totalTrials > 0
                % Prevent random-effect naming cross-contamination across ages
                uniqueRetinaID = sprintf('%d_%s', ageLabels(ai), R(r).prefix);
                
                rows(end+1,:) = { ...
                    uniqueRetinaID, ...          % Retina ID (Categorical)
                    ageLabels(ai), ...           % Age Group
                    nOnTrans, ...                % Counts for OnTrans GLMM
                    nOnSus, ...                  % Counts for OnSus GLMM
                    nOffTrans, ...               % Counts for OffTrans GLMM
                    nOffSus, ...                 % Counts for OffSus GLMM
                    totalTrials ...              % Plugs into 'BinomialSize'
                };
            end
        end
    end

    % Construct tidy table
    Tf = cell2table(rows, 'VariableNames', ...
        {'Retina','Age','OnTrans','OnSus','OffTrans','OffSus','nAllCells'});

    % Ensure continuous numerical data types for GLMM stability
    Tf.OnTrans   = double(Tf.OnTrans);
    Tf.OnSus     = double(Tf.OnSus);
    Tf.OffTrans  = double(Tf.OffTrans);
    Tf.OffSus    = double(Tf.OffSus);
    Tf.nAllCells = double(Tf.nAllCells);

    % Convert factors to categorical tracking groups
    Tf.Retina = categorical(Tf.Retina);
    Tf.Age    = categorical(Tf.Age);
    Tf.Age    = removecats(Tf.Age);
end