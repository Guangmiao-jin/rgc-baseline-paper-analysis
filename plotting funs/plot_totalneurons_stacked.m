function plot_totalneurons_stacked(filedir)
    % Prepare containers
    matFiles = dir(fullfile(filedir,'*_totalneuronsV2.mat'));
    basename = matFiles(1).folder;
    for k = 1:numel(matFiles)
        matfiles{k} = fullfile(matFiles(k).folder, matFiles(k).name);
    end
    nFiles = numel(matfiles);
    results = zeros(nFiles, 5); % Each row: [ON OFF ONOFF unconventional unresponsive]

    % Go through each file
    for i = 1:nFiles
        S = load(matfiles{i});
        fn = fieldnames(S);
        d = S.(fn{1});
        
        % ON: sum of subcategories
        ON_count = length(d.num_OnNeurons.trans) + length(d.num_OnNeurons.sus);
        % OFF: sum of subcategories
        OFF_count = length(d.num_OffNeurons.trans) + length(d.num_OffNeurons.sus);
        % ONOFF
        ONOFF_count = length(d.num_OnOffNeurons);
        % unconventional
        unconventional_count = length(d.num_unconventional);
        % unresponsive
        unresponsive_count = length(d.num_notRespond);

        results(i,:) = [ON_count, OFF_count, ONOFF_count, unconventional_count, unresponsive_count];
    end

    % Calculate percentages
    perc = results ./ sum(results,2) * 100; % each row sums to 100

    % Bar colors: red, blue, green, orange, gray
    colors = [1 0 0; 0 0 1; 0 1 0; 1 0.6 0; 0.5 0.5 0.5];

    % --------- Plot stacked bar (per file) ----------
    if nFiles ~=1
    figure('Color','w');
    b = bar(perc, 'stacked', 'BarWidth',0.7);   
        for k = 1:5
            b(k).FaceColor = colors(k,:);
        end
    
    ylabel('Percentage (%)');
    set(gca,'XTickLabel',{}); % no label underneath
    legend({'ON','OFF','ON-OFF','Unconventional','Unresponsive'}, 'Location','eastoutside');
    ylim([0 100])
    title('Stacked Bar: Neuron Classification per File (percentage)');

    % Add percentage labels on bars
    xt = get(gca, 'XTick');
    for i = 1:nFiles
        y_bottom = 0;
        for k = 1:5
            pct = perc(i,k);
            if pct > 0
                y = y_bottom + pct/2;
                text(i, y, sprintf('%.1f%%', pct),...
                     'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',9,'FontWeight','bold');
            end
            y_bottom = y_bottom + pct;
        end
    end

    % Save per-file figure
    saveas(gcf, [basename '\stackedbar_perfile.png']);
    close(gcf)
    end
    % --------- Summary plot: pool all values ----------
    total_counts = sum(results,1);
    total_perc = total_counts / sum(total_counts) * 100;

    figure('Color','w');
    b2 = bar(1, total_perc, 'stacked', 'BarWidth',0.7);
    for k = 1:5
        b2(k).FaceColor = colors(k,:);
    end
    ylabel('Percentage (%)');
    set(gca,'XTickLabel',{'Summary'}); % just one bar
    ylim([0 100])
    legend({'ON','OFF','ON-OFF','Unconventional','Unresponsive'}, 'Location','eastoutside');
    title('Stacked Bar: Combined Summary (percentage)');

    % Add percentage labels
    y_bottom = 0;
    for k = 1:5
        pct = total_perc(k);
        if pct > 0
            y = y_bottom + pct/2;
            text(1, y, sprintf('%.1f%%', pct),...
                'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10,'FontWeight','bold');
        end
        y_bottom = y_bottom + pct;
    end

    % Save summary figure
    saveas(gcf, [basename '\stackedbar_summary.png']);
    close(gcf)
end
