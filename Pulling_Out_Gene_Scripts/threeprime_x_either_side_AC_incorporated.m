%% load in annotations
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


length_pulled_out = 1000;
all_genes = zeros(length_pulled_out,1) ;
%% calculating total gene number

gene_number_forward = 0 ; %start with matrix of zeros
gene_number_reverse = 0 ;
for d = 1:NoChromosomes
    
    gene_number_forward = gene_number_forward + sum(length(gene_positions{d,1}(:,1)))  ;
    gene_number_reverse = gene_number_reverse + sum(length(gene_positions{d,2}(:,1))) ;
    %iteratively add up all the genes present
    
end
Total_gene_number = gene_number_forward + gene_number_reverse ;
%% forward
NET_Seq_matrix_forward = zeros(gene_number_forward,length_pulled_out ) ;
o = 0 ;
for i = 1:NoChromosomes

    %Isolating one chromosome%
    chromosomei_NETseq = NETseq_sample1(i,1) ;
    chromosomei_genepos = gene_positions(i,1) ;
    
    
    %pulling out 500bp before the end of gene%
    chromosomei_startpoints = chromosomei_genepos{1,1}(:,2) - length_pulled_out/2 ;
    
    %start pos + 999 gives last 500 bases + 500 after%
    chromosomei_endpoints = chromosomei_startpoints + (length_pulled_out-1) ;
    
    %pulling out the individual genes%
    z = size(gene_positions{i,1}(:,1), 1) ;
   
    for k = 1:z 
        
        o = o + 1 ; 
        %add all genes together in single matrix for averaging later
        all_genes = all_genes + chromosomei_NETseq{1,1}(chromosomei_startpoints(k,1):chromosomei_endpoints(k,1),1) ;
        
        %construct a matrix with all gene profiles in
        NET_Seq_matrix_forward(o, :) = chromosomei_NETseq{1,1}(chromosomei_startpoints(k,1):chromosomei_endpoints(k,1),1) ;
    
    end
end
    
%% reverse
NET_Seq_matrix_reverse = zeros(gene_number_reverse,length_pulled_out ) ;
w = 0 ;
for j = 1:NoChromosomes
   
    %Isolating one chromosome%
    chromosomej_NETseq = NETseq_sample1(j,2) ;
    chromosomej_genepos = gene_positions(j,2) ;
    
    
    %pulling out the 250bp before end of genes%
    chromosomej_startpoints = chromosomej_genepos{1,1}(:,1) + (length_pulled_out/2 - 1) ;
    %start points - 499 gives last 250 bases of reversed gene + 250 after%
    chromosomej_endpoints = chromosomej_startpoints - (length_pulled_out - 1) ;
    
    %pulling out the individual genes%
    x = size(gene_positions{j,2}(:,1), 1) ;
    
    for k = 1:x
        
        
        w = w + 1  ;
       
        %add to total signal matrix for all genes
        all_genes = all_genes + chromosomej_NETseq{1,1}(chromosomej_startpoints(k,1):-1:chromosomej_endpoints(k,1),1) ;
    
     
        %construct a matrix with all gene profiles in
        NET_Seq_matrix_reverse(w, :) = chromosomej_NETseq{1,1}(chromosomej_startpoints(k,1):-1:chromosomej_endpoints(k,1),1);
    end
end


%% Look at metagene forward and reverse 

metagene_forward = zeros(1, length_pulled_out);
for k= 1:gene_number_forward
    metagene_forward = metagene_forward + NET_Seq_matrix_forward(k, :) ;
end
metagene_forward = metagene_forward ./ gene_number_forward ;

metagene_reverse = zeros(1, length_pulled_out);
for k= 1:gene_number_reverse
    metagene_reverse = metagene_reverse + NET_Seq_matrix_reverse(k, :) ;
end
metagene_reverse = metagene_reverse ./ gene_number_reverse ;

%% plot metagenes 
figure 
plot(metagene_forward')

figure
plot(metagene_reverse')


%% combine reverse and forward matrices
NET_Seq_Matrix1 = vertcat(NET_Seq_matrix_forward,NET_Seq_matrix_reverse) ;

%% forward  repeat for biological repeat then average to get overall NETseq matrix 
NET_Seq_matrix_forward2 = zeros(gene_number_forward,length_pulled_out ) ;
o = 0 ;
for i = 1:NoChromosomes

    %Isolating one chromosome%
    chromosomei2_NETseq = NETseq_sample2(i,1) ;
    chromosomei2_genepos = gene_positions(i,1) ;
    
    
    %pulling out 250bp before the end of gene%
    chromosomei2_startpoints = chromosomei2_genepos{1,1}(:,2) - length_pulled_out/2 ;
    
    %start pos + 499 gives last 500 bases + 500 after%
    chromosomei2_endpoints = chromosomei2_startpoints + (length_pulled_out-1) ;
    
    %pulling out the individual genes%
    z = size(gene_positions{i,1}(:,1), 1) ;
   
    for k = 1:z 
        
        o = o + 1 ; 
        %add all genes together in single matrix for averaging later
        all_genes = all_genes + chromosomei2_NETseq{1,1}(chromosomei2_startpoints(k,1):chromosomei2_endpoints(k,1),1) ;
        
        %construct a matrix with all gene profiles in
        NET_Seq_matrix_forward2(o, :) = chromosomei2_NETseq{1,1}(chromosomei2_startpoints(k,1):chromosomei2_endpoints(k,1),1) ;
    
    end
end
    
%% reverse
NET_Seq_matrix_reverse2 = zeros(gene_number_reverse,length_pulled_out ) ;
w = 0 ;
for j = 1:NoChromosomes
   
    %Isolating one chromosome%
    chromosomej2_NETseq = NETseq_sample2(j,2) ;
    chromosomej2_genepos = gene_positions(j,2) ;
    
    
    %pulling out the 250bp before end of genes%
    chromosomej2_startpoints = chromosomej2_genepos{1,1}(:,1) + (length_pulled_out/2 - 1) ;
    %start points - 499 gives last 250 bases of reversed gene + 250 after%
    chromosomej2_endpoints = chromosomej2_startpoints - (length_pulled_out - 1) ;
    
    %pulling out the individual genes%
    x = size(gene_positions{j,2}(:,1), 1) ;
    
    for k = 1:x
        
        
        w = w + 1  ;
       
        %add to total signal matrix for all genes
        all_genes = all_genes + chromosomej2_NETseq{1,1}(chromosomej2_startpoints(k,1):-1:chromosomej2_endpoints(k,1),1) ;
    
     
        %construct a matrix with all gene profiles in
        NET_Seq_matrix_reverse2(w, :) = chromosomej2_NETseq{1,1}(chromosomej2_startpoints(k,1):-1:chromosomej2_endpoints(k,1),1);
    end
end

%% Look at metagene forward and reverse 
metagene_forward2 = zeros(1, length_pulled_out);
for k= 1:gene_number_forward
    metagene_forward2 = metagene_forward2 + NET_Seq_matrix_forward(k, :) ;
end
metagene_forward2 = metagene_forward2 ./ gene_number_forward ;

metagene_reverse2 = zeros(1, length_pulled_out);
for k= 1:gene_number_reverse
    metagene_reverse2 = metagene_reverse2 + NET_Seq_matrix_reverse(k, :) ;
end
metagene_reverse2 = metagene_reverse2 ./ gene_number_reverse ;

 %% plot metagenes 
figure 
plot(metagene_forward2)

figure
plot(metagene_reverse2)


%% combine reverse and forward matrices
NET_Seq_Matrix2 = vertcat(NET_Seq_matrix_forward2,NET_Seq_matrix_reverse2) ;


%% Average the two samples 
NET_Seq_Matrix = (NET_Seq_Matrix1 + NET_Seq_Matrix2) ./ 2 ;

%% Remove negative values 
NETseq_matrix_no_negative = NET_Seq_Matrix ; 
NETseq_matrix_no_negative(NET_Seq_Matrix<0) = 0;


%g=size(NET_Seq_Matrix(Logic_matrix_remove_negative));
% for h=1:g
%     Total_NET
%     
% end
% NETseq_matrix_no_negative = NET_Seq_Matrix(Logic_matrix_remove_negative)==0 
% %transpose to allow matrix multiplication
% no_negativity = Logic_matrix_remove_negative';
%generate 
% NETseq_matrix_no_negative = Total_NET_Seq_Matrix*no_negativity ;

%% Remove low expression genes 

%% Calculate mean across genome as 'background' 
total_across_genome = 0;
genome_size = 0 ;
NETseq_composite = cell(size(NETseq_sample1));
for cctr = 1:1:NoChromosomes %chromosme counter 
    for sctr = 1:1:2 %strand counter 
        
        NETseq_composite{cctr,sctr} = ...
            nanmean( ...
            [NETseq_sample1{cctr,sctr},NETseq_sample2{cctr,sctr}],2 ...
            );
        
    end
end
%calculate mean from composite cell array 
for cctr = 1:NoChromosomes
    for sctr = 1:2
        total_across_genome = total_across_genome + sum(NETseq_composite{cctr, sctr}) ;
        
        genome_size = genome_size + sum(numel(NETseq_composite{cctr, sctr})) ;
        
        
    end
end

mean_across_genome = total_across_genome ./ genome_size ;
%Combine samples

%% 

for g = size(NETseq_matrix_no_negative, 1):-1:1 
    if nanmean(NETseq_matrix_no_negative(g, :), 2)<mean_across_genome %make this the
    
       NETseq_matrix_no_negative(g, :) = [] ;
    end      
end

 

%% Normalise no negative NETseq Matrix
Normalised_NETseq_Matrix = Shape_normalisation_function(NETseq_matrix_no_negative); %only works if shape normalisation function is in current folder, need to find out how to specify path to function but for now will leave


% remove NaNs
 for k= size(Normalised_NETseq_Matrix, 1):-1:1 
    
    if  sum(isnan(Normalised_NETseq_Matrix(k, :)), 2)>0 %gives logical showing which genes contain any NaNs 
        Normalised_NETseq_Matrix(k,:) = [] ;
%     else 
%         Normalised_NETseq_Matrix(k,:) = Normalised_NETseq_Matrix(k,:);
    end
    
 end
 
 %% Find Normalised metagene 

 metagene_Normalised = zeros(1, length_pulled_out);
for k= 1:size(Normalised_NETseq_Matrix, 1)
    metagene_Normalised = metagene_Normalised + Normalised_NETseq_Matrix(k, :) ;
end
metagene_Normalised = metagene_Normalised ./ size(Normalised_NETseq_Matrix, 1) ;

%% plot normalised metagene 
figure 
plot(metagene_Normalised)

%% NaN to zeros 
% Normalised_NETseq_Matrix(isnan(Normalised_NETseq_Matrix)) = 0; converts
% Nans to zeros rather than removing

%% Bin the data prior to normalisation 
%% averaging every 10bps Bins
%make the matrix 
NETseq_Bins10bp_avg = zeros(Total_gene_number, length_pulled_out/10);
for k = 1:Total_gene_number
    for x = 1:length_pulled_out/10
        Bin = mean(NETseq_matrix_no_negative(k,10*x-9:10*x));
        NETseq_Bins10bp_avg(k, x) = Bin  ;
    end
end


%% Normalise Binned NETseq
Normalised_NETseq_Bins_Matrix = Shape_normalisation_function(NETseq_Bins10bp_avg);
%remove NaNs
for k= size(Normalised_NETseq_Bins_Matrix, 1):-1:1 
    
    if  sum(isnan(Normalised_NETseq_Bins_Matrix(k, :)), 2)>0 %gives logical showing which genes contain any NaNs 
        Normalised_NETseq_Bins_Matrix(k,:) = [] ;
%     else not needed 
%         Normalised_NETseq_Matrix(k,:) = Normalised_NETseq_Matrix(k,:);
    end
    
end

%% average gene calculation
meta_gene = all_genes ./ (2*Total_gene_number) ;


