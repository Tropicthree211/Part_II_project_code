%% 
klist=2:10;%the number of clusters to try
myfunc = @(X ,K)(kmeans(X, K));
%% evaluate clusters calinskiHarabasz
eva = evalclusters(Normalised_NETseq_Bins_Matrix ,myfunc,'CalinskiHarabasz','Klist',klist) 
%% evaluate clusters Silhouette
eva = evalclusters(Normalised_NETseq_Bins_Matrix ,myfunc,'silhouette','Klist',klist) 

%% perform the final kmeans cluster
mutant = 'HARRYWTTIFseqANNOTATIONS';
Cluster_Job = 14 ;
JOB = append( 'Current_Cluster_Job', string(Cluster_Job));
unix(append('mkdir ', JOB))


 for K= 2:3
    for repeat = 1:5
        classes=kmeans(Normalised_NETseq_Bins_Matrix ,K);
        
        %% Pull out the classes from the clustering
        Class_sorted_NET_seq_matrix = zeros(size(Normalised_NETseq_Bins_Matrix, 1), size(Normalised_NETseq_Bins_Matrix, 2));
        Classmetagenes = zeros(size(Normalised_NETseq_Bins_Matrix, 2),K) ;
        classize_matrix = zeros(K, 1);
        classize_matrix0 = zeros(K+1, 1) ;
        classize_matrix_string = strings(K, 1);
        
        for n =1:K
            %sorted matrix for each class
            Classn = Normalised_NETseq_Bins_Matrix(classes==n, :) ;
            
            prevClasses = Normalised_NETseq_Bins_Matrix(classes<n, :);
            
            classize = size(Classn, 1);
            classize_matrix((n), 1) =  classize 
            classize_matrix0((n+1), 1) = classize; %this allows the 5' end clustering to work outside of the loop 
            
            prevclassizes = size(prevClasses,1);
            [sorted_classes, ix] = sort(classes);
            Class_sorted_NET_seq_matrix = Normalised_NETseq_Bins_Matrix(ix,:);
%             Class_sorted_NET_seq_matrix(prevclassizes+1:classize+prevclassizes,
%             :) = Classn ; old way of sorting use sort function allows
%             index to be pulled out that can sort 
            
            classize_matrix_string(n, 1) = string(classize_matrix(n, 1) ) ;
            
            
            
            %metagene for each class
            Classmetagenes(:, n) = sum(Class_sorted_NET_seq_matrix(prevclassizes+1:classize+prevclassizes, :))/classize ;
            
            %plot each metagene one at a time
            %     figure
            %     plot(Classmetagenes(:,n))
            
        end
        filename = append(mutant, string(ybefore), 'before', string(xafter), 'after', string(K), 'clusters', string(repeat));
        figure
        plot(Classmetagenes)
        title(filename) ;
        legend(classize_matrix_string)
        
        %% save as fig and jpeg
        saveas(gcf, filename, 'fig')
        
        saveas(gcf, filename, 'jpeg')
        disp('saved 3end clustering data')
        
        %% Alternate sorting method Class_sorted_NET_seq_matrix = Class_sorted_NET_seq_matrix2 showing that ix contains original position info of genes.
%         [sorted_classes, ix] = sort(classes);
%         Class_sorted_NET_seq_matrix2 = Normalised_NETseq_Bins_Matrix(ix,:) ;
%         
%         Check_sorting_works = Normalised_NETseq_Bins_Matrix(ix, :) == Class_sorted_NET_seq_matrix ;
%         
%         clear original_order
%         original_order(ix, :) = Class_sorted_NET_seq_matrix ; % unsorts the above, so orginal_order == Normalised_NETseq_Bins_Matrix
        
        % %% Getting metagenes from sorted matrix
        % [p, index] = unique(sorted_classes)
        % Classmetagenes2 = zeros(size(Class_sorted_NET_seq_matrix2, 2),K);
        % for n=K:-1:1
        %
        %     Classmetagenes2(:, n) = sum(Class_sorted_NET_seq_matrix(size(sorted_classes,1):-1:index(n,1), :))/classize ;
        % end
        
        
        %% Applying 3' end cluster sorting to 5' end data
        %load in 5' end NETseqMatrix
        %load('Unfiltered_5_end_genes.mat') %need to either have in the current path or save it somewhere specific and specify the path here.Saves me having to run the 5' end every time.
        %make sure to check the annotations files used match up. 
        %Apply same selection to 5' end
        NETseq_3primeexprnSelected_Normalised_Bins_Matrix_start = Normalised_NETseq_Bins_Matrix_start(list_remaining_genes, :);
        %Apply same sorting to 5' end as above using ix
        NETseq_matrix_3prime_Class_sorted_fiveprime_end = NETseq_3primeexprnSelected_Normalised_Bins_Matrix_start(ix,:);
        
        
        
        
        
        %% Find metagenes at 5' end
        Classmetagenes5prime_end = zeros(size(Normalised_NETseq_Bins_Matrix_start, 2), K);
        Cumulative_Classize = cumsum(classize_matrix0) ;
        for q = 2:K+1
            % Classmetagenes5prime_end(:, (q-1)) = sum(NETseq_matrix_3prime_Class_sorted_fiveprime_end( (classize_matrix((q-1),1)+1):(classize_matrix((q-1),1)+classize_matrix(q,1)), :)))/classize_matrix(q, 1) ;
            Classmetagenes5prime_end(:, (q-1)) = nansum(NETseq_matrix_3prime_Class_sorted_fiveprime_end((Cumulative_Classize((q-1), 1) + 1): Cumulative_Classize(q , 1) , :))/classize_matrix0(q, 1) ;
            %had to use nansum due to NaNs, need to find source see if they're
            %impacting results
        end
        filename2 = append(mutant, '5endsortedon3endclustering', string(ybefore), 'before', string(xafter), 'after', string(K), 'clusters', string(repeat));
        figure
        plot(Classmetagenes5prime_end)
        title(filename2);
        legend(classize_matrix_string)
        %% save as fig and jpeg
        saveas(gcf, filename2, 'fig')
        
        saveas(gcf, filename2, 'jpeg')
        
        disp('saved 5end clusters')
        %% heatmap class metagenes
        figure
        imagesc(Class_sorted_NET_seq_matrix)
        title(append(filename, 'heatmap'))
        truesize
        filename3 = append(filename, 'heatmap');
        saveas(gcf, filename3, 'jpeg')
        saveas(gcf, filename3, 'fig')
    end
    
end
unix(append('mv HARRY* ',JOB))

%% Look for patterns in heatmap
figure
imagesc(NETseq_matrix_3prime_Class_sorted_fiveprime_end)
truesize

%% Attach gene names 
load('WT_actively_expressing_genes_list.mat') %In current path or need to specify here 
%take care to use correct gene annotations file. 

%construct list of gene names in the same way as NET-seq matrix to apply
%filter to. 
%%
gene_names_list = strings(Total_gene_number, 1) ;%empty presized list to add all the gene names into
%% forward names 
w = 0;
for cctr = 1:1:NoChromosomes
    x = size(gene_names{cctr,1} , 2);
    for gctr = 1:1:x
        w = w+1;
        gene_names_list(w,1) = gene_names{cctr,1}(1, gctr);
    end
end

%% reverse names
% got error code in the middle and ended up duplicating 100 names. fixed
% now.
for cctr = 1:1:NoChromosomes
    x = size(gene_names{cctr,2} , 2);
    for gctr = 1:1:x
        w = w+1;
        gene_names_list(w,1) = gene_names{cctr,2}(1, gctr);
        if w>Total_gene_number
            disp('too many gene names')
        end
    end
end
%%
filtered_gene_names = gene_names_list(list_remaining_genes);

cluster_sorted_gene_names = filtered_gene_names(ix, :);
    



