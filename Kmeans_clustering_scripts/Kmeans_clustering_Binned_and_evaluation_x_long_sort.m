%% 
klist=2:10;%the number of clusters to try
myfunc = @(X ,K)(kmeans(X, K));
%% evaluate clusters calinskiHarabasz
eva = evalclusters(Normalised_NETseq_Bins_Matrix_start ,myfunc,'CalinskiHarabasz','Klist',klist) 
%% evaluate clusters Silhouette
eva = evalclusters(Normalised_NETseq_Bins_Matrix_start ,myfunc,'silhouette','Klist',klist) 

%% perform the final kmeans cluster 
K=3 ;
classes=kmeans(Normalised_NETseq_Bins_Matrix_start ,K);

%% Pull out the classes from the clustering
% Class_sorted_NET_seq_matrix = zeros(size(Normalised_NETseq_Bins_Matrix_start, 1), size(Normalised_NETseq_Bins_Matrix_start, 2));
Classmetagenes = zeros(size(Normalised_NETseq_Bins_Matrix_start, 2),K) ;
classize_matrix = zeros(K, 1);


for n =1:K
    %sorted matrix for each class
    Classn = Normalised_NETseq_Bins_Matrix_start(classes==n, :) ;
    
    prevClasses = Normalised_NETseq_Bins_Matrix_start(classes<n, :);
    
    classize = size(Classn, 1)
    classize_matrix(n, 1) = classize
    prevclassizes = size(prevClasses,1);
    
    Class_sorted_NET_seq_matrix(prevclassizes+1:classize+prevclassizes, :) = Classn ;
    
    
    


    %metagene for each class 
    Classmetagenes(:, n) = sum(Class_sorted_NET_seq_matrix(prevclassizes+1:classize+prevclassizes, :))/classize ;
    
    %plot each metagene one at a time 
%     figure 
%     plot(Classmetagenes(:,n))
    
end
figure
plot(Classmetagenes)


%% Alternate sorting method Class_sorted_NET_seq_matrix = Class_sorted_NET_seq_matrix2 showing that ix contains original position info of genes.
[sorted_classes, ix] = sort(classes);
Class_sorted_NET_seq_matrix2 = Normalised_NETseq_Bins_Matrix_start(ix,:) ;

%Normalised_NETseq_Bins_Matrix_start(ix, :) == Class_sorted_NET_seq_matrix
original_order(ix, :) = Class_sorted_NET_seq_matrix ; % unsorts the above, so orginal_order == Normalised_NETseq_Bins_Matrix_start

% %% Getting metagenes from sorted matrix 
% [p, index] = unique(sorted_classes)
% Classmetagenes2 = zeros(size(Class_sorted_NET_seq_matrix2, 2),K);
% for n=K:-1:1
%     
%     Classmetagenes2(:, n) = sum(Class_sorted_NET_seq_matrix(size(sorted_classes,1):-1:index(n,1), :))/classize ;
% end
