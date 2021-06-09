%% load in
load ulku_WT_NETseq_1.mat
load gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat
all_genes = zeros(1000,1) ;
%% calculating total gene number

gene_number = 0 ; %start with matrix of zeros
for d = 1:16
    
    gene_number = gene_number + sum(length(gene_positions{d,1}(:,1))) + sum(length(gene_positions{d,2}(:,1))) ;
    %iteratively add up all the genes present
end
%% forward
NET_Seq_matrix = zeros(1000, 1) ;
for i = 1:16
    %Isolating one chromosome%
    chromosomei_NETseq = Ulku_WT_NETseq_1(i,1) ;
    chromosomei_genepos = gene_positions(i,1) ;
    
    
    %pulling out the start point of genes%
    chromosomei_startpoints = chromosomei_genepos{1,1}(:,1) ;
    %start pos + 999 gives first 1000 bases%
    chromosomei_endpoints = chromosomei_startpoints + 999 ;
    
    %pulling out the individual genes%
    z = size(gene_positions{i,1}(:,1), 1) ;
    for k = 1:z
        %construct a matrix with all gene profiles in 
        NET_Seq_matrix = [NET_Seq_matrix chromosomei_NETseq{1,1}(chromosomei_startpoints(k,1):chromosomei_endpoints(k,1),1)];
        %add all genes together in single matrix for averaging later 
        all_genes = all_genes + chromosomei_NETseq{1,1}(chromosomei_startpoints(k,1):chromosomei_endpoints(k,1),1) ;
    end
end
NET_Seq_matrix(:, 1) = [] ;
%% reverse
NET_Seq_matrix2 = zeros(1000,1) ;
for j = 1:16
    %Isolating one chromosome%
    chromosomej_NETseq = Ulku_WT_NETseq_1(j,2) ;
    chromosomej_genepos = gene_positions(j,2) ;
    
    
    
    %pulling out the start point of genes%
    chromosomej_startpoints = chromosomej_genepos{1,1}(:,2) ;
    %start points - 999 gives first 1000 bases of reversed gene%
    chromosomej_endpoints = chromosomej_startpoints - 999 ;
    
    %pulling out the individual genes%
    x = size(gene_positions{j,2}(:,1), 1) ;
    
    for k = 1:x
        %construct a matrix with all gene profiles in 
        NET_Seq_matrix2 = [NET_Seq_matrix2, chromosomej_NETseq{1,1}(chromosomej_startpoints(k,1):-1:chromosomej_endpoints(k,1),1)] ;
        %add to total signal matrix for all genes
        all_genes = all_genes + chromosomej_NETseq{1,1}(chromosomej_startpoints(k,1):-1:chromosomej_endpoints(k,1),1) ;
    end
end
NET_Seq_matrix2(:, 1) = [] ;

%% combine reverse and forward matrices 
Total_NET_Seq_Matrix = [NET_Seq_matrix NET_Seq_matrix2] ;

%% average gene calculation
meta_gene = all_genes / gene_number


