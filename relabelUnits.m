function [res_firings, res_metrics, res_isol_metrics, res_isol_pair_metrics] = relabelUnits(firings, metrics, isol_metrics, isol_pair_metrics)

    res_firings = firings;
    res_metrics = metrics;
    res_isol_metrics = isol_metrics;
    res_isol_pair_metrics = isol_pair_metrics;

    orig_unique_units = sort(unique(firings(3,:)));
    
    num_unique_units = length(orig_unique_units);
    relabel_unique_units = 1:num_unique_units;
    
    %relabel firings
    for i = 1:length(orig_unique_units)
       
        unit_spikes = (firings(3,:) == orig_unique_units(i));
        res_firings(3,unit_spikes) = relabel_unique_units(i);
        
    end
    
    %metrics
    for i = 1:length(metrics.clusters)
       
        res_metrics.clusters(i).label = relabel_unique_units( metrics.clusters(i).label == orig_unique_units);
        
    end
    
    %isol_metrics
    for i = 1:length(isol_metrics.clusters)
       
        res_isol_metrics.clusters(i).label = relabel_unique_units( isol_metrics.clusters(i).label == orig_unique_units);
        
    end
    
    %isol_pair_metrics
    for i = 1:length(isol_pair_metrics.cluster_pairs)
       
        pair = isol_pair_metrics.cluster_pairs(i).label;
        pair_split = strsplit(pair, ',');
        left_pair = str2num(pair_split{1});
        right_pair = str2num(pair_split{2});
        
        res_left_pair = relabel_unique_units( left_pair == orig_unique_units );
        res_right_pair = relabel_unique_units( right_pair == orig_unique_units );

        res_isol_pair_metrics.cluster_pairs(i).label = sprintf('%d,%d', res_left_pair, res_right_pair);
        
    end
    
    
end