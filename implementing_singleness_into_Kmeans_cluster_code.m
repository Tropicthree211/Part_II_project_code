


%% Bring in info from Singleness measure 
load /home/orie3770/Modelling/Singlenessmeasure2individual.mat %make sure to use with NoCUTs data, has only 5579 values in. doesn't include Cuts
singleness_filtered = singlenessbygenemeasure2(list_remaining_genes,1);
% might get issues with a bunch of NaNs in the data but we'll see how it
% comes out.

%put this just before the big cluster 


%% Boxplot of singleness within each of the clusters 
        filenamex = append(filename, 'singleness_boxplot');
        figure
        boxplot(singleness_filtered, classes) 
        %Think i can replace classes with any other division into classes. 
        title(filenamet)
        ylabel('singleness score')
        xlabel('cluster')
        
    %% save boxplot
        saveas(gcf, filenamex, 'fig')
        saveas(gcf, filenamex, 'jpeg')
        
        
    % this within the for loop in which clustering and repeats occur. 