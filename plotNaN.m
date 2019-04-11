function plotNaN(current_chan_names, current_processed, samplingFreq, refset_savedir)

    fprintf('plotting NaNs\n');

    figure();
    set(gcf,'color','w');
    set(gcf,'visible','off');
    
    output_fpath = [refset_savedir 'nan.png'];
    
    nan_vec = isnan(current_processed(:));
    not_nan_vec = ~isnan(current_processed(:));
    
    current_processed(nan_vec) = 0;
    current_processed(not_nan_vec) = 1;
    
    imagesc(current_processed, [0 1]);
    CMap=[1,1,1; 0,0,0];
    colormap(CMap);
    
    title('white squares indicate NaNs');
    
    yticks(1:length(current_chan_names));
    yticklabels(current_chan_names);
    
    min_marks = floor(xticks()/(60*samplingFreq));
    
    xticklabels(strsplit(num2str(min_marks), ' '));
    
    xlabel('time (min)');
    
    set(gca, 'FontSize', 5);
    
    tightfig();
    set(gcf, 'InvertHardcopy', 'off','PaperUnits','inches','PaperPosition',[0,0,25,19],'PaperPositionMode','auto');

    print(gcf,output_fpath,'-dpng','-r100');


end