avg_gene = zeros(1000,1) ;
%Isolating just one of the first chromosome pair%
chromosome1_NETseq = Ulku_WT_NETseq_1(1,1) ;
chromosome1_genepos = gene_positions(1,1) ;
%looking at the array data in numeric form%
numericVector = cell2mat(chromosome1_genepos) ;
%pulling out the start point of genes%
chromosome1_startpoints = chromosome1_genepos{1,1}(:,1) ;
%start pos + 999 gives first 1000 bases%
metagene_endpoints = chromosome1_startpoints + 999 ; 
%pulling out the individual genes%
for i = 1:22
avg_gene = avg_gene + chromosome1_NETseq{1,1}(chromosome1_startpoints(i,1):metagene_endpoints(i,1),1) ;
end

%
chromosome1_end_product = avg_gene / 22 ;
%

