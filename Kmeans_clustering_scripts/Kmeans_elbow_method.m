%%  Kmeans Clustering elbow method 
%% perform the kmeans cluster 
K=4;
classes=kmeans(Normalised_NETseq_Matrix ,K);

%% Pull out the classes from the clustering
Class_sorted_NET_seq_matrix = zeros(Total_gene_number, 1000);
Classmetagenes = zeros(1000,K) ;
variance_matrix = zeros(1, K) ;
%remove NaNs
% Normalised_NETseq_Matrix_NoNaN = Normalised_NETseq_Matrix ;
% Normalised_NET_Seq_Matrix_NoNaN(isnan(Normalised_NETseq_Matrix_NoNaN)) = [];
for n = 1:K
    %sorted matrix for each class 
    Classn = Normalised_NETseq_Matrix(classes==n, :) ;
    
    prevClasses = Normalised_NETseq_Matrix(classes<n, :);
    
    classize = size(Classn, 1)  
    prevclassizes = size(prevClasses,1);
    
    Class_sorted_NET_seq_matrix(prevclassizes+1:classize+prevclassizes, :) = Classn ;
    
    %variance within classes
    variance_matrix(1, n) = var(Class_sorted_NET_seq_matrix(prevclassizes+1:classize+prevclassizes, :),0, 2) ;
    
    %metagene for each class 
    Classmetagenes(:, n) = sum(Class_sorted_NET_seq_matrix(prevclassizes+1:classize+prevclassizes, :))/classize ;
    
%     %plot each metagene one at a time 
%     figure 
%     plot(Classmetagenes(:,n))
    
end
figure
plot(Classmetagenes)
    
    