%% load in
load ulku_WT_NETseq_1.mat
load gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat
all_genes = zeros(1000,1) ;
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
NET_Seq_matrix = zeros(gene_number_forward,1000 ) ;
o = 0 ;
for i = 1:16

    %Isolating one chromosome%
    chromosomei_NETseq = Ulku_WT_NETseq_1(i,1) ;
    chromosomei_genepos = gene_positions(i,1) ;
    
    
    %pulling out 500bp before the end of gene%
    chromosomei_startpoints = chromosomei_genepos{1,1}(:,2) - 500 ;
    
    %start pos + 999 gives last 500 bases + 500 after%
    chromosomei_endpoints = chromosomei_startpoints + 999 ;
    
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
NET_Seq_matrix2 = zeros( gene_number_reverse,1000 ) ;
w = 0 ;
for j = 1:16
   
    %Isolating one chromosome%
    chromosomej_NETseq = Ulku_WT_NETseq_1(j,2) ;
    chromosomej_genepos = gene_positions(j,2) ;
    
    
    %pulling out the 500bp before end of genes%
    chromosomej_startpoints = chromosomej_genepos{1,1}(:,1) + 499 ;
    %start points - 999 gives last 500 bases of reversed gene + 500 after%
    chromosomej_endpoints = chromosomej_startpoints - 999 ;
    
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

%% Remove negative values 
NETseq_matrix_no_negative = Total_NET_Seq_Matrix ; 
NETseq_matrix_no_negative(Total_NET_Seq_Matrix<0) = 0;


%g=size(Total_NET_Seq_Matrix(Logic_matrix_remove_negative));
% for h=1:g
%     Total_NET
%     
% end
% NETseq_matrix_no_negative = Total_NET_Seq_Matrix(Logic_matrix_remove_negative)==0 
% %transpose to allow matrix multiplication
% no_negativity = Logic_matrix_remove_negative';
%generate 
% NETseq_matrix_no_negative = Total_NET_Seq_Matrix*no_negativity ;

%% Normalise no negative NETseq Matrix
Normalised_NETseq_Matrix = Shape_normalisation_function(NETseq_matrix_no_negative);
% remove NaNs
Normalised_NETseq_Matrix(isnan(Normalised_NETseq_Matrix)) = 0;

%% Bin the data prior to normalisation 
%% averaging every 10bps Bins
%make the matrix 
NETseq_Bins10bp_avg = zeros(Total_gene_number, 100);
for k = 1:Total_gene_number
    for x = 1:100
        Bin = mean(NETseq_matrix_no_negative(k,10*x-9:10*x));
        NETseq_Bins10bp_avg(k, x) = Bin  ;
    end
end


%failed code 
%n = 10 ;
% reshape(NETseq_matrix_no_negative,3129,100,10)

% for k= 1:size(NETseq_matrix_no_negative, 1)
%     gene_noneg = NETseq_matrix_no_negative(k, :)
%     n = 10
%     arrayfun(@i
%% Normalise Binned NETseq
Normalised_NETseq_Bins_Matrix = Shape_normalisation_function(NETseq_Bins10bp_avg);
%remove NaNs 
Normalised_NETseq_Bins_Matrix(isnan(Normalised_NETseq_Bins_Matrix)) = 0;
%% average gene calculation
meta_gene = all_genes / Total_gene_number ;


