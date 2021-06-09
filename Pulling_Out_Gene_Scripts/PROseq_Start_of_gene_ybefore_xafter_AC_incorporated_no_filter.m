% Have chagned every variable to end _start so this code can be run
% alongside the threeprime code without deleting each other's stuff

% Will have to alter the kmeans clustering code to accept _start variables

%all the loading sections are irrelevant echoes of how the code was made by
%modifying the 3end script
%% Load in the processed annotation data

%NB when using alternative annotations, gene positions should always be 
%leftmost to rightmost irrespective of strand
%If this is not the case, need to reprocess the names here to conform

% annotationloaddir = ...
%     '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ENSEMBL/S_CEREVISIAE_ANNOTATIONS/';
% 
% extrasavepiece = '';
% 
% 
% loadfile1 = 'gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat';
% loadfile2 = 'gene_ids_Saccharomyces_cerevisiae__R64_1_1__101.mat';
% 
% load([annotationloaddir loadfile1])
% load([annotationloaddir loadfile2])
% 
% disp('Loaded annotation files.')
% 
%  load gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat
%  NoChromosomes = 16 ;
 % extra annotations can be loaded found in andrew's code 
%% Start with WT 

% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ULKU_PRO_SEQ/BigWiG_spikein_normalised2019_PROandTEF-seq/notag_normalised/WT_nnt/';
% 
% PROseq_sample = importdata([dataloaddir 'ulku_WT_PROseq_1.mat']);
% PROseq_sample2 = importdata([dataloaddir 'ulku_WT_PROseq_2.mat']);
% 
% disp('Loaded in UU WT PROseq data')

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


 
xafter_start = 1000;
ybefore_start = 0;
length_pulled_out_start = (xafter_start+ybefore_start);
all_genes_start = zeros(length_pulled_out_start,1) ;
%% spt4 delete

% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ULKU_PRO_SEQ/BigWiG_spikein_normalised2019_PROandTEF-seq/notag_normalised/spt4del_nnt/';
% 
% PROseq_sample = importdata([dataloaddir 'ulku_spt4del_PROseq_1.mat']);
% PROseq_sample2 = importdata([dataloaddir 'ulku_spt4del_PROseq_2.mat']);
% 
% disp('Loaded in UU spt4 delete PROseq data')

%NOTES:
%The correlation between the two samples was very good
%The average levels in the genes lie significantly off the diagonal - this
%might indicate a levelling issue and a more careful average might be
%needed in the longer run.
% xafter_start = 1000;
% ybefore_start = 250;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes_start = zeros((xafter_start+ybefore_start),1) ;


%% dst1 delete

% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ULKU_PRO_SEQ/BigWiG_spikein_normalised2019_PROandTEF-seq/notag_normalised/dst1del_nnt/';
% 
% PROseq_sample = importdata([dataloaddir 'ulku_dst1del_PROseq_1.mat']);
% PROseq_sample2 = importdata([dataloaddir 'ulku_dst1del_PROseq_2.mat']);
% 
% disp('Loaded in UU dst1 delete PROseq data')
% 
% %NOTES:
% %The correlation between the two samples was very good
% %The average levels in each were reasonably close to the diagonal
% xafter_start = 1000;
% ybefore_start = 250;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes_start = zeros((xafter_start+ybefore_start),1) ;

%% xrn1 delete

% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ULKU_PRO_SEQ/BigWiG_spikein_normalised2019_PROandTEF-seq/notag_normalised/xrn1del_nnt/';
% 
% PROseq_sample = importdata([dataloaddir 'ulku_xrn1del_PROseq_1.mat']);
% PROseq_sample2 = importdata([dataloaddir 'ulku_xrn1del_PROseq_2.mat']);
% 
% disp('Loaded in UU xrn1 delete PROseq data')
% 
% %NOTES:
% %The correlation between the two samples was very good
% %The average levels at genes were very close to the diagonal
% 
% 
%  
% xafter_start = 1000;
% ybefore_start = 250;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes_start = zeros((xafter_start+ybefore_start),1) ;
%% SPT4 Anchor Away
% 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ULKU_PRO_SEQ/BigWiG_spikein_normalised2019_PROandTEF-seq/AnchorAway/AA_PROCESSED_PROSEQ/SPT4AA/';
% 
% PROseq_sample = importdata([dataloaddir 'ulku_SPT4AA_PROseq_1.mat']);
% PROseq_sample2 = importdata([dataloaddir 'ulku_SPT4AA_PROseq_2.mat']);
% 
% disp('Loaded in UU SPT4 Anchor Away PROseq data')

%NOTES:
%The correlation between the two samples is very good
%The average levels ag genes were very close to the diagonal
% xafter_start = 1000;
% ybefore_start = 250;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes_start = zeros((xafter_start+ybefore_start),1) ;

%% SPT5 Anchor Away

% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ULKU_PRO_SEQ/BigWiG_spikein_normalised2019_PROandTEF-seq/AnchorAway/AA_PROCESSED_PROSEQ/SPT5AA/';
% 
% PROseq_sample = importdata([dataloaddir 'ulku_SPT5AA_PROseq_1.mat']);
% PROseq_sample2 = importdata([dataloaddir 'ulku_SPT5AA_PROseq_2.mat']);
% 
% disp('Loaded in UU SPT5 Anchor Away PROseq data')
% 
% %NOTES:
% %Correlation between these two sets is extremely good
% %The average levels in genes lie very close to the diagonal on visual
% %inspection
% xafter_start = 1000; 
% ybefore_start = 250;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes_start = zeros((xafter_start+ybefore_start),1) ;
%% SPT4/5 DMSO (Treating as repeats but strictly they are different)

% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ULKU_PRO_SEQ/BigWiG_spikein_normalised2019_PROandTEF-seq/AnchorAway/AA_PROCESSED_PROSEQ/SPT4AA_SPT5AA_DMSO_COMBO/';
% 
% PROseq_sample = importdata([dataloaddir 'ulku_SPT4AADMSO_PROseq_1.mat']);
% PROseq_sample2 = importdata([dataloaddir 'ulku_SPT5AADMSO_PROseq_1.mat']);
% 
% disp('Loaded in UU SPT5 Anchor Away PROseq data')
% 
% %NOTES:
% %Correlation between these two samples (which should be similar but are
% %technically distinct) is very good
% %The average levels in genes also lies reasonably close to the diagonal
% xafter_start = 1000; 
% ybefore_start = 250;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes_start = zeros((xafter_start+ybefore_start),1) ;
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
PRO_Seq_matrix_forward_start = zeros(gene_number_forward,length_pulled_out_start) ;
o = 0 ;
for i = 1:NoChromosomes

    %Isolating one chromosome%
    chromosomei_PROseq = PROseq_sample(i,1) ;
    chromosomei_genepos = gene_positions(i,1) ;
    
    
    %pulling out ybeforebp before the end of gene%
    chromosomei_startpoints = chromosomei_genepos{1,1}(:,1) - ybefore_start ;
    
    %start pos + 999 gives last 500 bases + 500 after%
    chromosomei_endpoints = chromosomei_startpoints + (length_pulled_out_start-1) ;
    
    %pulling out the individual genes%
    z = size(gene_positions{i,1}(:,1), 1) ;
   
    for k = 1:z 
        
        o = o + 1 ; 
        %add all genes together in single matrix for averaging later
        all_genes_start = all_genes_start + chromosomei_PROseq{1,1}(chromosomei_startpoints(k,1):chromosomei_endpoints(k,1),1) ;
        
        %construct a matrix with all gene profiles in
        PRO_Seq_matrix_forward_start(o, :) = chromosomei_PROseq{1,1}(chromosomei_startpoints(k,1):chromosomei_endpoints(k,1),1) ;
    
    end
end
    
%% reverse
PRO_Seq_matrix_reverse_start = zeros(gene_number_reverse,length_pulled_out_start) ;
w = 0 ;
fail_counter_reverse_start = 0;
for j = 1:NoChromosomes
   
    %Isolating one chromosome%
    chromosomej_PROseq = PROseq_sample(j,2) ;
    chromosomej_genepos = gene_positions(j,2) ;
    
    
    %pulling out the ybeforebp before end of genes%
    chromosomej_startpoints = chromosomej_genepos{1,1}(:,2)  ;
    %start points - x+y gives last x bases of reversed gene %
    chromosomej_endpoints = chromosomej_startpoints - (length_pulled_out_start - 1) ;
    
    %pulling out the individual genes%
    x = size(gene_positions{j,2}(:,1), 1) ;
    
    for k = 1:x
        
        
        w = w + 1  ;
       try
        %add to total signal matrix for all genes
        all_genes_start = all_genes_start + chromosomej_PROseq{1,1}(chromosomej_startpoints(k,1):-1:chromosomej_endpoints(k,1),1) ;
    
     
        %construct a matrix with all gene profiles in
        PRO_Seq_matrix_reverse_start(w, :) = chromosomej_PROseq{1,1}(chromosomej_startpoints(k,1):-1:chromosomej_endpoints(k,1),1);
       catch 
           fail_counter_reverse_start = fail_counter_reverse_start +1 
           
       end
       
    end
end


%% Look at metagene forward and reverse 

metagene_forward_start = zeros(1, length_pulled_out_start);
for k= 1:gene_number_forward
    metagene_forward_start = metagene_forward_start + PRO_Seq_matrix_forward_start(k, :) ;
end
metagene_forward_start = metagene_forward_start ./ gene_number_forward ;

metagene_reverse_start = zeros(1, length_pulled_out_start);
for k= 1:gene_number_reverse
    metagene_reverse_start = metagene_reverse_start + PRO_Seq_matrix_reverse_start(k, :) ;
end
metagene_reverse_start = metagene_reverse_start ./ gene_number_reverse ;

%% plot metagenes 
figure 
plot(metagene_forward_start')
title('metagene forward')
figure
plot(metagene_reverse_start')
title('metagene reverse')

%% combine reverse and forward matrices
PRO_Seq_Matrix_start = vertcat(PRO_Seq_matrix_forward_start,PRO_Seq_matrix_reverse_start) ; % edited removede 1 as only one proseq sample 

%% forward  repeat for biological repeat then average to get overall PROseq matrix 
% PRO_Seq_matrix_forward2 = zeros(gene_number_forward,length_pulled_out_start) ;
% o = 0 ;
% for i = 1:NoChromosomes
% 
%     %Isolating one chromosome%
%     chromosomei2_PROseq = PROseq_sample2(i,1) ;
%     chromosomei2_genepos = gene_positions(i,1) ;
%     
%     
%     %pulling out 250bp before the end of gene%
%     chromosomei2_startpoints = chromosomei2_genepos{1,1}(:,1) - ybefore_start ;
%     
%     %start pos + x+y gives last y bases + x after%
%     chromosomei2_endpoints = chromosomei2_startpoints + (length_pulled_out_start-1) ;
%     
%     %pulling out the individual genes%
%     z = size(gene_positions{i,1}(:,1), 1) ;
%    
%     for k = 1:z 
%         
%         o = o + 1 ; 
%         %add all genes together in single matrix for averaging later
%         all_genes_start = all_genes_start + chromosomei2_PROseq{1,1}(chromosomei2_startpoints(k,1):chromosomei2_endpoints(k,1),1) ;
%         
%         %construct a matrix with all gene profiles in
%         PRO_Seq_matrix_forward2(o, :) = chromosomei2_PROseq{1,1}(chromosomei2_startpoints(k,1):chromosomei2_endpoints(k,1),1) ;
%     
%     end
% end
%     
% %% reverse
% PRO_Seq_matrix_reverse_start2 = zeros(gene_number_reverse,length_pulled_out_start) ;
% w = 0 ;
% for j = 1:NoChromosomes
%    
%     %Isolating one chromosome%
%     chromosomej2_PROseq = PROseq_sample2(j,2) ;
%     chromosomej2_genepos = gene_positions(j,2) ;
%     
%     
%     %pulling out the 250bp before end of genes%
%     chromosomej2_startpoints = chromosomej2_genepos{1,1}(:,2) + ybefore_start  ;
%     %start points - 499 gives last 250 bases of reversed gene + 250 after%
%     chromosomej2_endpoints = chromosomej2_startpoints - (length_pulled_out_start - 1) ;
%                                                                                                                                                                                                                                                                                                                             
%     %pulling out the individual genes%
%     x = size(gene_positions{j,2}(:,1), 1) ;
%     
%     for k = 1:x
%         
%         
%         w = w + 1  ;                                            
% 
%         %add to total signal matrix for all genes
%         all_genes_start = all_genes_start + chromosomej2_PROseq{1,1}(chromosomej2_startpoints(k,1):-1:chromosomej2_endpoints(k,1),1) ;
%     
%      
%         %construct a matrix with all gene profiles in
%         PRO_Seq_matrix_reverse_start2(w, :) = chromosomej2_PROseq{1,1}(chromosomej2_startpoints(k,1):-1:chromosomej2_endpoints(k,1),1);
%     end
% end
% 
% %% Look at metagene forward and reverse 
% metagene_forward_start2 = zeros(1, length_pulled_out_start);
% for k= 1:gene_number_forward
%     metagene_forward_start2 = metagene_forward_start2 + PRO_Seq_matrix_forward_start(k, :) ;
% end
% metagene_forward_start2 = metagene_forward_start2 ./ gene_number_forward ;
% 
% metagene_reverse_start2 = zeros(1, length_pulled_out_start);
% for k= 1:gene_number_reverse
%     metagene_reverse_start2 = metagene_reverse_start2 + PRO_Seq_matrix_reverse_start(k, :) ;
% end
% metagene_reverse_start2 = metagene_reverse_start2 ./ gene_number_reverse ;
% 
%  %% plot metagenes 
% figure 
% plot(metagene_forward_start2)
% title('metagene foward repeat')
% figure
% plot(metagene_reverse_start2)
% title('metagene reverse repeat')
% 
% %% combine reverse and forward matrices
% PRO_Seq_Matrix_start2 = vertcat(PRO_Seq_matrix_forward2,PRO_Seq_matrix_reverse_start2) ;


%% Average the two samples 
% PRO_Seq_Matrix_start = (PRO_Seq_Matrix_start1 + PRO_Seq_Matrix_start2) ./ 2 ;

%% Remove negative values 
PROseq_matrix_no_negative_start = PRO_Seq_Matrix_start ; 
PROseq_matrix_no_negative_start(PRO_Seq_Matrix_start<0) = 0;


%g=size(PRO_Seq_Matrix(Logic_matrix_remove_negative));
% for h=1:g
%     Total_PRO
%     
% end
% PROseq_matrix_no_negative = PRO_Seq_Matrix(Logic_matrix_remove_negative)==0 
% %transpose to allow matrix multiplication
% no_negativity = Logic_matrix_remove_negative';
%generate 
% PROseq_matrix_no_negative = Total_PRO_Seq_Matrix*no_negativity ;

%% Remove low expression genes 

% %% Calculate mean across genome as 'background' 
% total_across_genome = 0;
% genome_size = 0 ;
% PROseq_composite = cell(size(PROseq_sample));
% for cctr = 1:1:NoChromosomes %chromosme counter 
%     for sctr = 1:1:2 %strand counter 
%         
%         PROseq_composite{cctr,sctr} = ...
%             nanmean( ...
%             [PROseq_sample{cctr,sctr},PROseq_sample2{cctr,sctr}],2 ...
%             );
%         
%     end
% end
% %calculate mean from composite cell array 
% for cctr = 1:NoChromosomes
%     for sctr = 1:2
%         total_across_genome = total_across_genome + sum(PROseq_composite{cctr, sctr}) ;
%         
%         genome_size = genome_size + sum(numel(PROseq_composite{cctr, sctr})) ;
%         
%         
%     end
% end
% 
% mean_across_genome = total_across_genome ./ genome_size ;
% %Combine samples
% 
% %% REMOVE VALUES BELOW MEAN
% 
% for g = size(PROseq_matrix_no_negative_start, 1):-1:1 
%     if nanmean(PROseq_matrix_no_negative_start(g, :), 2)<mean_across_genome %found above
%     
%        PROseq_matrix_no_negative_start(g, :) = [] ;
%     end      
% end
% 
%  

%% Normalise no negative PROseq Matrix
Normalised_PROseq_Matrix_start = Shape_normalisation_function(PROseq_matrix_no_negative_start); %only works if shape normalisation function is in current folder, need to find out how to specify path to function but for now will leave


% % remove NaNs
%  for k= size(Normalised_PROseq_Matrix_start, 1):-1:1 
%     
%     if  sum(isnan(Normalised_PROseq_Matrix_start(k, :)), 2)>0 %gives logical showing which genes contain any NaNs 
%         Normalised_PROseq_Matrix_start(k,:) = [] ;
% %     else 
% %         Normalised_PROseq_Matrix(k,:) = Normalised_PROseq_Matrix(k,:);
%     end
%     
%  end
 
 %% Find Normalised metagene 

 metagene_Normalised_start = zeros(1, length_pulled_out_start);
for k= 1:size(Normalised_PROseq_Matrix_start, 1)
    metagene_Normalised_start = metagene_Normalised_start + Normalised_PROseq_Matrix_start(k, :) ;
end
metagene_Normalised_start = metagene_Normalised_start ./ size(Normalised_PROseq_Matrix_start, 1) ;

%% plot normalised metagene 
%given no filtration often contains Nans here, not usually an issue as
%filtration occurs via copying 3end filtration later 
% figure 
% plot(metagene_Normalised_start)
% title('normalised metagene')


%% Bin the data prior to normalisation 
%% averaging every 10bps Bins
%make the matrix 
PROseq_Bins10bp_avg_start = zeros(size(PROseq_matrix_no_negative_start, 1), length_pulled_out_start/10);
for k = 1:size(PROseq_matrix_no_negative_start, 1)
    for x = 1:length_pulled_out_start/10
        Bin = mean(PROseq_matrix_no_negative_start(k,10*x-9:10*x));
        PROseq_Bins10bp_avg_start(k, x) = Bin  ;
    end
end


%% Normalise Binned PROseq
Normalised_PROseq_Bins_Matrix_start = Shape_normalisation_function(PROseq_Bins10bp_avg_start); %should probably integrate this into the code but i thought i was cool making a function
% %remove NaNs dont want to do this here as otherwise no lineup
% for k= size(Normalised_PROseq_Bins_Matrix_start, 1):-1:1 
%     
%     if  sum(isnan(Normalised_PROseq_Bins_Matrix_start(k, :)), 2)>0 %gives logical showing which genes contain any NaNs 
%         Normalised_PROseq_Bins_Matrix_start(k,:) = [] ;
% %     else not needed 
% %         Normalised_PROseq_Matrix(k,:) = Normalised_PROseq_Matrix(k,:);
%     end
%     
% end
%% NaN to zeros 
%Normalised_PROseq_Bins_Matrix(isnan(Normalised_PROseq_Bins_Matrix)) = 0; %converts
%Nans to zeros rather than removing
%% average gene calculation
% meta_gene_start = all_genes ./ (2*Total_gene_number) ;
% figure
% plot(meta_gene_start)
