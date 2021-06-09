%% 
klist=2:10;%the number of clusters to try
myfunc = @(X ,K)(kmeans(X, K));
%% evaluate clusters calinskiHarabasz
eva = evalclusters(Normalised_NETseq_Bins_Matrix ,myfunc,'CalinskiHarabasz','Klist',klist) 
%% evaluate clusters Silhouette
eva = evalclusters(Normalised_NETseq_Bins_Matrix ,myfunc,'silhouette','Klist',klist) 

%% perform the final kmeans cluster
mutant = 'nogene500_WTforwardstrand_TIFseqANNOTATIONS';
for K= 2:6
    for z = 1:5
        classes=kmeans(Normalised_NETseq_Bins_Matrix ,K);
        
        %% Pull out the classes from the clustering
        Class_sorted_NET_seq_matrix = zeros(size(Normalised_NETseq_Bins_Matrix, 1), size(Normalised_NETseq_Bins_Matrix, 2));
        Classmetagenes = zeros(size(Normalised_NETseq_Bins_Matrix, 2),K) ;
        classize_matrix = zeros(K, 1);
        
        
        for n =1:K
            %sorted matrix for each class
            Classn = Normalised_NETseq_Bins_Matrix(classes==n, :) ;
            
            prevClasses = Normalised_NETseq_Bins_Matrix(classes<n, :);
            
            classize = size(Classn, 1);
            classize_matrix((n+1), 1) =  classize
            prevclassizes = size(prevClasses,1);
            
            Class_sorted_NET_seq_matrix(prevclassizes+1:classize+prevclassizes, :) = Classn ;
            
            
            
            
            
            %metagene for each class
            Classmetagenes(:, n) = sum(Class_sorted_NET_seq_matrix(prevclassizes+1:classize+prevclassizes, :))/classize ;
            
            %plot each metagene one at a time
            %     figure
            %     plot(Classmetagenes(:,n))
            
        end
        filename = append(mutant, string(ybefore), '_before_', string(xafter), '_after_', string(K), '_clusters_', string(z));
        figure
        plot(Classmetagenes)
        title(filename) ;
        legend
        
        %% save as fig and jpeg
        saveas(gcf, filename, 'fig')
        
        saveas(gcf, filename, 'jpeg')
        disp('saved 3end clustering data')
        
        %% Alternate sorting method Class_sorted_NET_seq_matrix = Class_sorted_NET_seq_matrix2 showing that ix contains original position info of genes.
        [sorted_classes, ix] = sort(classes);
        Class_sorted_NET_seq_matrix2 = Normalised_NETseq_Bins_Matrix(ix,:) ;
        
        Check_sorting_works = Normalised_NETseq_Bins_Matrix(ix, :) == Class_sorted_NET_seq_matrix ;
        
        clear original_order
        original_order(ix, :) = Class_sorted_NET_seq_matrix ; % unsorts the above, so orginal_order == Normalised_NETseq_Bins_Matrix
        
        % %% Getting metagenes from sorted matrix
        % [p, index] = unique(sorted_classes)
        % Classmetagenes2 = zeros(size(Class_sorted_NET_seq_matrix2, 2),K);
        % for n=K:-1:1
        %
        %     Classmetagenes2(:, n) = sum(Class_sorted_NET_seq_matrix(size(sorted_classes,1):-1:index(n,1), :))/classize ;
        % end
        
        
        %% Applying 3' end cluster sorting to 5' end data
        %Apply same selection to 5' end
        NETseq_3primeexprnSelected_Normalised_Bins_Matrix_start = Normalised_NETseq_Bins_Matrix_start(list_remaining_genes, :);
        %Apply same sorting to 5' end as above using ix
        NETseq_matrix_3prime_Class_sorted_fiveprime_end = NETseq_3primeexprnSelected_Normalised_Bins_Matrix_start(ix,:);
        
        
        
        
        
        %% Find metagenes at 5' end
        Classmetagenes5prime_end = zeros(size(Normalised_NETseq_Bins_Matrix_start, 2), K);
        Cumulative_Classize = cumsum(classize_matrix) ;
        for q = 2:(K+1)
            % Classmetagenes5prime_end(:, (q-1)) = sum(NETseq_matrix_3prime_Class_sorted_fiveprime_end( (classize_matrix((q-1),1)+1):(classize_matrix((q-1),1)+classize_matrix(q,1)), :)))/classize_matrix(q, 1) ;
            Classmetagenes5prime_end(:, (q-1)) = nansum(NETseq_matrix_3prime_Class_sorted_fiveprime_end((Cumulative_Classize((q-1), 1) + 1): Cumulative_Classize(q , 1) , :))/classize_matrix(q, 1) ;
            %had to use nansum due to NaNs, need to find source see if they're
            %impacting results
        end
        filename2 = append(mutant, '5_end_sorted_based_on_3_end_clustering', string(ybefore), '_before', string(xafter), '_after_', string(K), '_clusters_', string(z));
        figure
        plot(Classmetagenes5prime_end)
        title(filename2);
        legend
        %% save as fig and jpeg
        saveas(gcf, filename2, 'fig')
        
        saveas(gcf, filename2, 'jpeg')
        
        disp('saved 5end clusters')
    end
    
end
unix('mv nogene500* Current_Cluster_Job2')
%% heatmap class metagenes
imagesc(Class_sorted_NET_seq_matrix)
title(append(filename, 'heatmap'))
truesize
%% Look for patterns in heatmap
figure
imagesc(NETseq_matrix_3prime_Class_sorted_fiveprime_end)
truesize