%% load in
load ulku_WT_NETseq_1.mat
load gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat
all_genes = zeros(1250,1) ;
%% calculating total gene number

gene_number_forward = 0 ; %start with matrix of zeros
gene_number_reverse = 0 ;
for d = 1:16
    
    gene_number_forward = gene_number_forward + sum(length(gene_positions{d,1}(:,1)))  ;
    gene_number_reverse = gene_number_reverse + sum(length(gene_positions{d,2}(:,1))) ;
    %iteratively add up all the genes present
    
end
Total_gene_number = gene_number_forward + gene_number_reverse ;
%% forward
NET_Seq_matrix = zeros(gene_number_forward,1250 ) ;
o = 0 ;
for i = 1:16
         
    %Isolating one chromosome%
    chromosomei_NETseq = Ulku_WT_NETseq_1(i,1) ;
    chromosomei_genepos = gene_positions(i,1) ;
    
    
    %pulling out the start point of genes%
    chromosomei_startpoints = chromosomei_genepos{1,1}(:,1) -250 ;
    %start pos + 999 gives first 1000 bases%
    chromosomei_endpoints = chromosomei_startpoints + 1249 ;
    
    %pulling out the individual genes%
    z = size(gene_positions{i,1}(:,1), 1) ;
   
    for k = 1:z 
        
        o = o + 1 ; 
        %add all genes together in single matrix for averaging later
        all_genes = all_genes + chromosomei_NETseq{1,1}(chromosomei_startpoints(k,1):chromosomei_endpoints(k,1),1) ;
        
        %construct a matrix with all gene profiles in
        NET_Seq_matrix(o, :) = chromosomei_NETseq{1,1}(chromosomei_startpoints(k,1):chromosomei_endpoints(k,1),1) ;
    
    end
end
    
%% reverse
NET_Seq_matrix2 = zeros( gene_number_reverse,1250 ) ;
w = 0 ;
for j = 1:16
   
    %Isolating one chromosome%
    chromosomej_NETseq = Ulku_WT_NETseq_1(j,2) ;
    chromosomej_genepos = gene_positions(j,2) ;
    
    
    %pulling out the start point of genes%
    chromosomej_startpoints = chromosomej_genepos{1,1}(:,2) + 250 ;
    %start points - 1249 gives first 1250 bases of reversed gene%
    chromosomej_endpoints = chromosomej_startpoints - 1249 ;
    
    %pulling out the individual genes%
    x = size(gene_positions{j,2}(:,1), 1) ;
    
    for k = 1:x
        
        
        w = w + 1  ;
       
        %add to total signal matrix for all genes
        all_genes = all_genes + chromosomej_NETseq{1,1}(chromosomej_startpoints(k,1):-1:chromosomej_endpoints(k,1),1) ;
    
     
        %construct a matrix with all gene profiles in
        NET_Seq_matrix2(w, :) = chromosomej_NETseq{1,1}(chromosomej_startpoints(k,1):-1:chromosomej_endpoints(k,1),1);
    end
end


%% combine reverse and forward matrices
Total_NET_Seq_Matrix = vertcat(NET_Seq_matrix,NET_Seq_matrix2) ;

%% average gene calculation
meta_gene = all_genes / Total_gene_number ;


