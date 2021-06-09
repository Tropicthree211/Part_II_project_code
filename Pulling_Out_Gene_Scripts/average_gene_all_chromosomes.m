avg_gene = zeros(1000,1) ;
%Isolating the chromosome%
 
 for   i = 1:16
%chromosomei_NETseq  
Ulku_WT_NETseq_1(i,1) 
%chromosomei_genepos = 
gene_positions(i,1) ;

%looking at the array data in numeric form%
% numericVector = cell2mat(chromosomei_genepos) ;

%pulling out the start point of genes%
gene_positions(i,1){1,1}(:,1) 
%start pos + 999 gives first 1000 bases%
%chromosomei_endpoints 
gene_positions(i,1){1,1}(:,1) + 999 ; 
 end
%pulling out the individual genes%

for j = 1:last
avg_gene = avg_gene + Ulku_WT_NETseq_1(i,1){1,1}(gene_positions(i,1){1,1}(:,1)(j,1):(gene_positions(i,1){1,1}(:,1)+999)(j,1),1) ;
end

%
end_product = avg_gene / sum(
%

