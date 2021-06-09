%% 
klist=2:10;%the number of clusters to try
myfunc = @(X ,K)(kmeans(X, K));
%% evaluate clusters calinskiHarabasz
eva = evalclusters(Normalised_NETseq_Bins_Matrix ,myfunc,'CalinskiHarabasz','Klist',klist) 
%% evaluate clusters Silhouette
eva = evalclusters(Normalised_NETseq_Bins_Matrix ,myfunc,'silhouette','Klist',klist) 

%% perform the final kmeans cluster 
K=5 ;
classes=kmeans(Normalised_NETseq_Bins_Matrix ,K);

%% Pull out the classes from the clustering
Class_sorted_NET_seq_matrix = zeros(size(Normalised_NETseq_Bins_Matrix, 1), size(Normalised_NETseq_Bins_Matrix, 2));
Classmetagenes = zeros(size(Normalised_NETseq_Bins_Matrix, 2),K) ;

for n = 1:K
    %sorted matrix for each class 
    Classn = Normalised_NETseq_Bins_Matrix(classes==n, :) ;
    
    prevClasses = Normalised_NETseq_Bins_Matrix(classes<n, :);
    
    classize = size(Classn, 1)  
    prevclassizes = size(prevClasses,1);
    
    Class_sorted_NET_seq_matrix(prevclassizes+1:classize+prevclassizes, :) = Classn ;
    
    %variance within classes
    
    
    %metagene for each class 
    Classmetagenes(:, n) = sum(Class_sorted_NET_seq_matrix(prevclassizes+1:classize+prevclassizes, :))/classize ;
    
    %plot each metagene one at a time 
%     figure 
%     plot(Classmetagenes(:,n))
    
end
figure
plot(Classmetagenes)


    