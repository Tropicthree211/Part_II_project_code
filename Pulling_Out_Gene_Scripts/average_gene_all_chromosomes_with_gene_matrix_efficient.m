%% load in
load ulku_WT_NETseq_1.mat
load gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat
all_genes = zeros(1000,1) ;
%% Start with WT 

dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/notag_normalised/WT_nnt/';

NETseq_sample1 = importdata([dataloaddir 'ulku_WT_NETseq_1.mat']);
NETseq_sample2 = importdata([dataloaddir 'ulku_WT_NETseq_2.mat']);

disp('Loaded in UU WT NETseq data')

%NOTES:
%The correlation between the two samples was very good
%The average level at each gene fell onto the diagonal reasonably well

%leave things generic such that other data can be loaded into the same
%variables and processed in the same way!

%Currently saving the output with a generic name to the load directory 
% - this is fine when the original data is separated into different
% directories (for data sets where they are all in one - consider splitting
% up the data -- or at the very least, put the processed data into clearly
% labellled subdirectories)

%mutants to load can be found in andrew's code 
%% Load in the processed annotation data

%NB when using alternative annotations, gene positions should always be 
%leftmost to rightmost irrespective of strand
%If this is not the case, need to reprocess the names here to conform

annotationloaddir = ...
    '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ENSEMBL/S_CEREVISIAE_ANNOTATIONS/';

extrasavepiece = '';


loadfile1 = 'gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat';
loadfile2 = 'gene_ids_Saccharomyces_cerevisiae__R64_1_1__101.mat';

load([annotationloaddir loadfile1])
load([annotationloaddir loadfile2])

disp('Loaded annotation files.')

 load gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat
 NoChromosomes = 16 ;
 % extra annotations can be loaded found in andrew's code 
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
    
    
    %pulling out the start point of genes%
    chromosomei_startpoints = chromosomei_genepos{1,1}(:,1) ;

    %start pos + 999 gives first 1000 bases%
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
    
    
    %pulling out the start point of genes%
    chromosomej_startpoints = chromosomej_genepos{1,1}(:,2) ;
    %start points - 999 gives first 1000 bases of reversed gene%
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

%% Normalise no negative NETseq Matrix
Normalised_NETseq_Matrix = Shape_normalisation_function(NETseq_matrix_no_negative);
% remove NaNs
Normalised_NETseq_Matrix(isnan(Normalised_NETseq_Matrix)) = 0;
%% averaging every 10bps Bins
%make the matrix 
NETseq_Bins10bp_avg = zeros(Total_gene_number, 100);
for k = 1:Total_gene_number
    for x = 1:100
        Bin = mean(NETseq_matrix_no_negative(k,10*x-9:10*x));
        NETseq_Bins10bp_avg(k, x) = Bin  ;
    end
end
%% Normalise Binned NETseq
Normalised_NETseq_Bins_Matrix = Shape_normalisation_function(NETseq_Bins10bp_avg);
%remove NaNs 
Normalised_NETseq_Bins_Matrix(isnan(Normalised_NETseq_Bins_Matrix)) = 0;
%% average gene calculation
meta_gene = all_genes / Total_gene_number ;



