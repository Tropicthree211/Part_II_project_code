%% 

% Before running: change cluster job number 
%change mutant loaded in by pulling out genes code 
%change unfiltered data to be sorted by pulled out genes automated now

%Should try to automate to minimise the number of manual changes needed to save
%time 
klist=2:10;%the number of clusters to try
myfunc = @(X ,K)(kmeans(X, K));
%% evaluate clusters calinskiHarabasz
eva = evalclusters(Normalised_NETseq_Bins_Matrix ,myfunc,'CalinskiHarabasz','Klist',klist) 
figure 
plot(eva.CriterionValues)
%% evaluate clusters Silhouette
eva = evalclusters(Normalised_NETseq_Bins_Matrix ,myfunc,'silhouette','Klist',klist) 
figure
plot(eva.CriterionValues)
%% elbow method implenented in code outside of MATLAB package 
dim=10;
% default number of test to get minimun under differnent random centroids
test_num=10;
distortion=zeros(dim,1);
for k_temp=1:dim(1)
    [~,~,sumd]=kmeans(Normalised_NETseq_Bins_Matrix,k_temp,'emptyaction','drop');
    distortion_temp=sum(sumd);
    % try differnet tests to find minimun disortion under k_temp clusters
    for test_count=2:test_num
        [~,~,sumd]=kmeans(Normalised_NETseq_Bins_Matrix,k_temp,'emptyaction','drop');
        distortion_temp=min(distortion_temp,sum(sumd));
    end
    distortion(k_temp,1)=distortion_temp;
end
variance=distortion(1:end-1)-distortion(2:end); %good way to pull out 
distortion_percent=cumsum(variance)/(distortion(1)-distortion(end));
figure
plot(distortion_percent,'b*--')
xlabel('Number of clusters')
ylabel('Percentage Distortion explained')

%% Eigenvalues

figure, sgtitle('set2')

subplot(1,2,1)

semilogy(diag(S),'k-o'), grid on

xlabel('r')

ylabel('\sigma_r')

set(gca)

subplot(1,2,2)

plot(cumsum(diag(S))./sum(diag(S)),'k-o') 

grid on

xlabel('r')

ylabel('Cumulative Energy')

set(gca)

%%
% variance_plot = zeros(1, 9);
% for K = 2:10 
%     clustering = kmeans(Normalised_NETseq_Bins_Matrix, K);
%     %Pull out each cluster 
%     clustern = Normalised_NETseq_Bins_Matrix(clustering == K, :) ;%might be wrong way round check
%     %calculate the variance of clustern, not sure if I need to do this for
%     %every single cluster and take an average or not.
%     
%     clusternvar = var(clustern, 1, [1 2]) ; % [1 2] here hopefully specifies the 1st and 2nd dimension of the cluster which will give the variance between 
%     %individual genes in the cluster rather than columns 
%     variance_plot(1,(K-1)) = clusternvar ;
% end   
% figure
% plot(variance_plot) %might need to transpose can't remember which dimension plot function likes. 
% title('Elbow method to determine optimal cluster number') 
% xlabel('Cluster number')
% ylabel('Variance')

%% perform the final kmeans cluster
Cluster_Job = 23 ;
JOB = append( 'Current_Cluster_Job', string(Cluster_Job));
unix(append('mkdir ', JOB))

clusters = append(datatype, mutant1, annotations)
% ufdata = 'NETseqspt4delxTIF';
%ufdata = 'NETseqdst1delxTIF';
%ufdata = 'NETseqxrn1delxTIF' ;
%ufdata = 'NETseqspt4AAxTIF';

%ufdata = 'NETseqspt5AA';
%ufdata = 'NETseqspt45DMSO';
%ufdata = 'NETseqset1del5endxTIF';

%ufdata = 'NETseqeaf3del5endxTIF';
%ufdata = 'PROseqWT'
%ufdata = 'PROseqspt4del';
% generic ufdata below: 
%ufdata = append(datatype1, mutant, annotations) %need to redo

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
        Cumulative_Classize = cumsum(classize_matrix0)
        filename = append(clusters, string(ybefore), 'before', string(xafter), 'after', string(K), 'clusters', string(repeat));
        figure
        plot(Classmetagenes)
        title(filename) ;
        legend(classize_matrix_string)
        
        %% save as fig and jpeg
        saveas(gcf, filename, 'fig')
        
        saveas(gcf, filename, 'jpeg')
        disp('saved 3end clustering data') 
%Don't need to visualise WT
%         data again 
        
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
        
        
%         
%         
%         %% Applying 3' end WT clustering to 
%         load(append(mutant1,'_actively_expressing_genes_list', annotations, '.mat'))
%         ufdataloaddir = '/home/orie3770/Unfiltered_data/';
%         %replace loadfilex with mutant to load 
%         %loadfilex ='Unfiltered_PROseq_TIF_nogene500pc100b200a.mat';
%         %loadfilex ='Unfiltered_NETseqspt4del_TIF_nogene500pc100b200a.mat';  
%         %loadfilex = 'Unfiltered_NETseqdst1del_TIF_nogene500pc100b200a.mat';
%         %loadfilex = 'Unfiltered_NETseqxrn1del_TIF_nogene500pc100b200a.mat';
%         %loadfilex ='Unfiltered_NETseqspt4AA_TIF_nogene500pc100b200a.mat';
%         %loadfilex = 'Unfiltered_PROseqspt4del_TIF_nogene500pc100b200a.mat';
%         %loadfilex = 'Unfiltered_NETseqspt5AA_TIF_nogene500pc100b200a.mat';
%         %loadfilex = 'Unfiltered_NETseqspt45DMSO_TIF_nogene500pc100b200a.mat';
%         %loadfilex = 'Unfiltered_5_endset1delnogene500pcTIFseq.mat'; 
%         %loadfilex = 'Unfiltered_5_endset2delnogene500pcTIFseq.mat';% need
%         loadfilex = append('Unfiltered', datatype1, mutant, annotations, '.mat');
%         %to redo clustering with this loaded in for set2 forgot to change 
%         %loadfilex = 'Unfiltered_5_endeaf3delnogene500pcTIFseq.mat';
%         %Can just load in different unfiltered files here, make sure
%         %annotation stage is the same.
%         load([ufdataloaddir loadfilex])
%         %Apply same expression selection to PROseq data
%         NETseq_3exprn_selected_ufdata = Normalised_ufdata_Bins_Matrix(list_remaining_genes, :);
%         %Apply same sorting to PROseq 
%         NETseq_3class_sorted_ufdata = NETseq_3exprn_selected_ufdata(ix, :);
%         
%         %% Find ufdata metagenes with NETseq clusters 
%         
%         Classmetagenesufdata = zeros(size(Normalised_ufdata_Bins_Matrix, 2), K);
%         Cumulative_Classize = cumsum(classize_matrix0) ;
%         for q = 2:K+1
%             % Classmetagenes5prime_end(:, (q-1)) = sum(NETseq_matrix_3prime_Class_sorted_fiveprime_end( (classize_matrix((q-1),1)+1):(classize_matrix((q-1),1)+classize_matrix(q,1)), :)))/classize_matrix(q, 1) ;
%             Classmetagenesufdata(:, (q-1)) = nansum(NETseq_3class_sorted_ufdata((Cumulative_Classize((q-1), 1) + 1): Cumulative_Classize(q , 1) , :))/classize_matrix0(q, 1) ;
%             %had to use nansum due to NaNs, need to find source see if they're
%             %impacting results
%         end
%         filename4 = append(clusters, ufdata, string(ybefore), 'before', string(xafter), 'after', string(K), 'clusters', string(repeat));
%         figure
%         plot(Classmetagenesufdata)
%         title(filename4);
%         legend(classize_matrix_string)
%         %% save as fig and jpeg
%         saveas(gcf, filename4, 'fig')
%         
%         saveas(gcf, filename4, 'jpeg')
%         
%         disp('saved ufdata clusters')
%         %% heatmap ufdata metagenes
%         figure
%         imagesc(NETseq_3class_sorted_ufdata)
%         title(append(filename4, 'heatmap'))
%         truesize
%         filename5 = append(filename4, 'heatmap');
%         saveas(gcf, filename5, 'jpeg')
%         saveas(gcf, filename5, 'fig')
%         disp('saved ufdata clustered heatmap')
        %% Applying 3' end cluster sorting to 5' end data
        
       % OLD now incorporate 5end data filtering above. 
       %can copy and paste above and edit to ufdata1 if ever want to filter
       %two or more things at once.
       
%         %load in 5' end NETseqMatrix
%         load('Unfiltered_5_end_genes.mat') %need to either have in the current path or save it somewhere specific and specify the path here.Saves me having to run the 5' end every time.
%         %make sure to check the annotations files used match up. 
%         %Apply same selection to 5' end
%         NETseq_3primeexprnSelected_Normalised_Bins_Matrix_start = Normalised_NETseq_Bins_Matrix_start(list_remaining_genes, :);
%         %Apply same sorting to 5' end as above using ix
%         NETseq_matrix_3prime_Class_sorted_fiveprime_end = NETseq_3primeexprnSelected_Normalised_Bins_Matrix_start(ix,:);
        
%         %% Find metagenes at 5' end
%         Classmetagenes5prime_end = zeros(size(Normalised_NETseq_Bins_Matrix_start, 2), K);
%         Cumulative_Classize = cumsum(classize_matrix0) ;
%         for q = 2:K+1
%             % Classmetagenes5prime_end(:, (q-1)) = sum(NETseq_matrix_3prime_Class_sorted_fiveprime_end( (classize_matrix((q-1),1)+1):(classize_matrix((q-1),1)+classize_matrix(q,1)), :)))/classize_matrix(q, 1) ;
%             Classmetagenes5prime_end(:, (q-1)) = nansum(NETseq_matrix_3prime_Class_sorted_fiveprime_end((Cumulative_Classize((q-1), 1) + 1): Cumulative_Classize(q , 1) , :))/classize_matrix0(q, 1) ;
%             %had to use nansum due to NaNs, need to find source see if they're
%             %impacting results
%         end
%         filename2 = append(mutant, '5endsortedon3endclustering', string(ybefore), 'before', string(xafter), 'after', string(K), 'clusters', string(repeat));
%         figure
%         plot(Classmetagenes5prime_end)
%         title(filename2);
%         legend(classize_matrix_string)
%         %% save as fig and jpeg
%         saveas(gcf, filename2, 'fig')
%         
%         saveas(gcf, filename2, 'jpeg')
%         
%         disp('saved 5end clusters')
%% Attach gene names
%load(append(mutant,'_actively_expressing_genes_list', annotations, '.mat')) %In current path or need to specify here
%moved earlier to avoid issues.
%take care to use correct gene annotations file.

%construct list of gene names in the same way as NET-seq matrix to apply
%filter to.
%%
gene_names_clustered = strings(Total_gene_number, K);
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

for clctr = 2:1:(K+1)
    w=0;
    for pctr = (Cumulative_Classize((clctr-1),1)+1):Cumulative_Classize(clctr,1)
        w=w+1;
    gene_names_clustered(w,(clctr-1)) = ...
    cluster_sorted_gene_names(pctr,1);
    end
end
%% save gene names

save(append(annotations, datatype, 'clustersortedgenenames', string(K), 'clusters', string(repeat)), 'cluster_sorted_gene_names') % removed ufdata as just WT clusters
save(append(annotations, datatype, 'clusterednamesforGO', string(K), 'clusters', string(repeat)), 'gene_names_clustered')%ditto
disp('saved clustered gene names')
%% heatmap class metagenes
figure
imagesc(Class_sorted_NET_seq_matrix)
title(append(filename, 'heatmap'))
truesize
filename3 = append(filename, 'heatmap');
saveas(gcf, filename3, 'jpeg')
saveas(gcf, filename3, 'fig')
disp('saved NETseq clustered heatmap')

%% Gene size info into box plot of clusters
% cluster_sorted_gene_sizes = gene_sizes_matrix_filtered(ix,:);
% %prespecify size of gene_sizes_clustered
% gene_sizes_clustered = cell(K,1);
% 
% for clctr = 2:1:(K+1)
%     w=0;
%     for pctr = (Cumulative_Classize((clctr-1),1)+1):Cumulative_Classize(clctr,1)
%         w=w+1;
%     gene_sizes_clustered{(clctr-1),1}(w,1)  = ...
%     cluster_sorted_gene_sizes(pctr,1);
%     end
% end
%Trim down gene size matrix to remove zeros.
% for clctr = 2:1:K
%     for gctr = 1:1:Total_gene_number
%         if gene_sizes_clustered(gctr, clctr) == 0 
%             gene_sizes_clustered(gctr,clctr) = [];
%         end
%     end
% end
    %earlier code didn't realise how simple it was see below:        
%% make boxplot of gene sizes within the clusters 
filename4 = append(filename, 'boxplot');
figure
boxplot(gene_sizes_matrix_filtered, classes)
title(filename4)
ylabel('Gene Size')
xlabel('cluster')


% boxplot(gene_sizes_clustered{K,1})
% hold(ax1) 
% if K-1>=1
%     
% boxplot(gene_sizes_clustered{(K-1),1})
% hold on %will work for just 2 and 3 clusters 
% end
% if K-2>=1
%     figure
%     boxplot(gene_sizes_clustered{(K-2),1})
%     hold on
% end
%% save boxplot
saveas(gcf, filename4, 'fig')
saveas(gcf, filename4, 'jpeg')

    end
    
 end
unix(append('mv nogene500* ', JOB))
unix(append('mv NETseq', mutant, '* ', JOB))

%% Look for patterns in heatmap
% figure
% imagesc(NETseq_matrix_3prime_Class_sorted_fiveprime_end)
% truesize


