%% load in annotations 
tic
%% Load in the processed annotation data

%NB when using alternative annotations, gene positions should always be 
%leftmost to rightmost irrespective of strand
%If this is not the case, need to reprocess the names here to conform

annotationloaddir = ... 
    /data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/
  
%'/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/TIF_SEQ_ANNOTATIONS/nature12121-s2/'; 
% annotationloaddir = ...
    % '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ENSEMBL/S_CEREVISIAE_ANNOTATIONS/';
    
extrasavepiece = '';


%loadfile1 = 'gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat';
%%original file with Andrew's annotations

%loadfile1 = 'Directly_downstream_removed_gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat'; 
%Original annotations with cramped end genes removed 

loadfile1 =

%loadfile1 = 'gene_positions_pelechano_TIFseq_nogene500.mat';
%TIF-seq annotations with cramped end genes removed 

% changed to specify new downstream removed genes
%above also has minimum gene size of 500 rather than 1000 as didn't seem
%relevant to exclude based on size given that I'm looking at the end

%loadfile2 = 'gene_ids_Saccharomyces_cerevisiae__R64_1_1__101.mat';
%original andrew's id annotations 

%should insert loadfile for Ids with cramped ends removed here 


loadfile2 = 

%loadfile2 = 'gene_names_pelechano_TIFseq_nogene500.mat';
%TIF-seq annotations

%loadfile3 = 'gene_positions_for_nogene_mean_Saccharomyces_cerevisiae__R64_1_1__101.mat';
%added loadfile3 to avoid having to integrate all components for background
%calculation into this code 
loadfile3 = 

%loadfile3 = 'gene_positions_pelechano_TIFseq_for_nogene_mean';

load([annotationloaddir loadfile1])
load([annotationloaddir loadfile2])
load([annotationloaddir loadfile3])
disp('Loaded annotation files.')

 %load gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat %commented
 %%out to load new annotations below 
%  load Directly_downstream_removed_gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat
 NoChromosomes = 24 ;
 % extra annotations can be loaded found in andrew's code 
%% Start with WT 

% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/notag_normalised/WT_nnt/';
% 
% NETseq_sample1 = importdata([dataloaddir 'ulku_WT_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'ulku_WT_NETseq_2.mat']);
% 
% disp('Loaded in UU WT NETseq data')
% 
% %NOTES:
% %The correlation between the two samples was very good
% %The average level at each gene fell onto the diagonal reasonably well
% 
% %leave things generic such that other data can be loaded into the same
% %variables and processed in the same way!
% 
% %Currently saving the output with a generic name to the load directory 
% % - this is fine when the original data is separated into different
% % directories (for data sets where they are all in one - consider splitting
% % up the data -- or at the very least, put the processed data into clearly
% % labellled subdirectories)
% xafter = 200;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;

%mutants to load can be found in andrew's code 
%% spt4 delete
% 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/notag_normalised/spt4del_nnt/';
% 
% NETseq_sample1 = importdata([dataloaddir 'ulku_spt4del_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'ulku_spt4del_NETseq_2.mat']);
% 
% disp('Loaded in UU spt4 delete NETseq data')
% 
% % NOTES:
% % The correlation between the two samples was very good
% % The average levels in the genes lie significantly off the diagonal - this
% % might indicate a levelling issue and a more careful average might be
% % needed in the longer run.
% xafter = 200;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;


%% dst1 delete

% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/notag_normalised/dst1del_nnt/';
% 
% NETseq_sample1 = importdata([dataloaddir 'ulku_dst1del_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'ulku_dst1del_NETseq_2.mat']);
% 
% disp('Loaded in UU dst1 delete NETseq data')
% 
% %NOTES:
% %The correlation between the two samples was very good
% %The average levels in each were reasonably close to the diagonal
% xafter = 200;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;

%% xrn1 delete

% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/notag_normalised/xrn1del_nnt/';
% 
% NETseq_sample1 = importdata([dataloaddir 'ulku_xrn1del_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'ulku_xrn1del_NETseq_2.mat']);
% 
% disp('Loaded in UU xrn1 delete NETseq data')
% % 
% % %NOTES:
% % %The correlation between the two samples was very good
% % %The average levels at genes were very close to the diagonal
% % 
% % 
% %  
% xafter = 200;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
%% SPT4 Anchor Away

% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/AnchorAway/AA_PROCESSED_NETSEQ/SPT4AA/';
% 
% NETseq_sample1 = importdata([dataloaddir 'ulku_SPT4AA_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'ulku_SPT4AA_NETseq_2.mat']);
% 
% disp('Loaded in UU SPT4 Anchor Away NETseq data')
% 
% %NOTES:
% %The correlation between the two samples is very good
% %The average levels ag genes were very close to the diagonal
% xafter = 200;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
 %% SPT5 Anchor Away
% 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/AnchorAway/AA_PROCESSED_NETSEQ/SPT5AA/';
% 
% NETseq_sample1 = importdata([dataloaddir 'ulku_SPT5AA_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'ulku_SPT5AA_NETseq_2.mat']);
% 
% disp('Loaded in UU SPT5 Anchor Away NETseq data')
% 
% %NOTES:
% %Correlation between these two sets is extremely good
% %The average levels in genes lie very close to the diagonal on visual
% %inspection
% xafter = 200;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
%% SPT4/5 DMSO (Treating as repeats but strictly they are different)
% 
dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/AnchorAway/AA_PROCESSED_NETSEQ/SPT4AA_SPT5AA_DMSO_COMBO/';

NETseq_sample1 = importdata([dataloaddir 'ulku_SPT4AADMSO_NETseq_1.mat']);
NETseq_sample2 = importdata([dataloaddir 'ulku_SPT5AADMSO_NETseq_1.mat']);

disp('Loaded in UU SPT5 Anchor Away NETseq data')

%NOTES:
%Correlation between these two samples (which should be similar but are
%technically distinct) is very good
%The average levels in genes also lies reasonably close to the diagonal
 xafter = 200;
ybefore = 100;
length_pulled_out = xafter+ybefore;
all_genes = zeros((xafter+ybefore),1) ;
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
NET_Seq_matrix_forward = zeros(gene_number_forward,(ybefore+xafter)) ;
o = 0 ;
for i = 1:NoChromosomes

    %Isolating one chromosome%
    chromosomei_NETseq = NETseq_sample1(i,1) ;
    chromosomei_genepos = gene_positions(i,1) ;
    
    
    %pulling out ybeforebp before the end of gene%
    chromosomei_startpoints = chromosomei_genepos{1,1}(:,2) - ybefore;
    
    %start pos + length pulled out minus 1 gives ybefore xafter end%
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
NET_Seq_matrix_reverse = zeros(gene_number_reverse,(ybefore+xafter)) ;
w = 0 ;
for j = 1:NoChromosomes
   
    %Isolating one chromosome%
    chromosomej_NETseq = NETseq_sample1(j,2) ;
    chromosomej_genepos = gene_positions(j,2) ;
    
    
    %pulling out the ybeforebp before end of genes%
    chromosomej_startpoints = chromosomej_genepos{1,1}(:,1) + ybefore ;
    %start points - x+y gives last x bases of reversed gene %
    chromosomej_endpoints = chromosomej_startpoints - (ybefore+xafter - 1) ;
    
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
title('metagene forward')
figure
plot(metagene_reverse')
title('metagene reverse')

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
    chromosomei2_startpoints = chromosomei2_genepos{1,1}(:,2) - ybefore ;
    
    %start pos + x+y gives last y bases + x after%
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
    chromosomej2_startpoints = chromosomej2_genepos{1,1}(:,1) + ybefore  ;
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
title('metagene forward repeat')
figure
plot(metagene_reverse2)
title('metagene reverse repeat')

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

%% Calculate mean across all non-gene areas of genome, better background.
%needs unfiltered gene positions before removal of any genes see
%ImportENSEMBL_S_Cerevisiae_find_areas_outside_genes 
NETseq_nogenes = NETseq_composite;
failed_removal_counter = 0 ;
for cctr = 1:1:NoChromosomes
    for sctr = 1:1:2
    for gctr = size(gene_positions_for_nogene_mean{cctr,sctr}, 1):-1:1 %need to input before removal of genes by overlaps etc 
        NETseq_nogenes{cctr, sctr}(gene_positions_for_nogene_mean{cctr, sctr}(gctr, 1):gene_positions_for_nogene_mean{cctr, sctr}(gctr, 2), :) = [];
        %try NETseq_nogenes(gene_positions_for_nogene_mean{cctr, sctr}(gctr, 1):gene_positions_for_nogene_mean{cctr, sctr}(gctr, 2), :) = [];
 %it worked without failing so didn't need try loop
        %catch
            %failed_removal_counter = failed_removal_counter +1 ; %if genes overlap will attempt to remove already removed area may fail
        %end
    end
   
    end
end

total_across_genome_nogene = 0;
genome_size_nogene = 0;
%Calculate mean from nogenes NETseq array
for cctr = 1:NoChromosomes
    for sctr = 1:2
        total_across_genome_nogene = total_across_genome_nogene + sum(NETseq_nogenes{cctr, sctr}) ;
        
        genome_size_nogene = genome_size_nogene + sum(numel(NETseq_nogenes{cctr, sctr})) ;
        
        
    end
end

%mean_across_genome_nogene = 0.2160 %answer from wild type NETseq data 
mean_across_genome_nogene = total_across_genome_nogene ./ genome_size ;
disp('done')

%% remove values below background mean 
list_remaining_genes = (1:1:size(NETseq_matrix_no_negative, 1))'; %Index to apply same selection at 5' end later 

for g = size(NETseq_matrix_no_negative, 1):-1:1
    %if nanmean(NETseq_matrix_no_negative(g, :), 2)<mean_across_genome
    %%found above = original filter
     if nanmean(NETseq_matrix_no_negative(g, :), 2)<mean_across_genome_nogene %found above new filter lose less signal   
        list_remaining_genes(g, :) =  [];
        
        NETseq_matrix_no_negative(g, :) = [] ;
    end
end

 

%% Normalise no negative NETseq Matrix
Normalised_NETseq_Matrix = Shape_normalisation_function(NETseq_matrix_no_negative); %only works if shape normalisation function is in current folder, need to find out how to specify path to function but for now will leave


% % remove NaNs
%  for k= size(Normalised_NETseq_Matrix, 1):-1:1 
%     
%     if  sum(isnan(Normalised_NETseq_Matrix(k, :)), 2)>0 %gives logical showing which genes contain any NaNs 
%         Normalised_NETseq_Matrix(k,:) = [] ;
% %     else 
% %         Normalised_NETseq_Matrix(k,:) = Normalised_NETseq_Matrix(k,:);
%     end
%     
%  end
 
 %% Find Normalised metagene 

 metagene_Normalised = zeros(1, length_pulled_out);
for k= 1:size(Normalised_NETseq_Matrix, 1)
    metagene_Normalised = metagene_Normalised + Normalised_NETseq_Matrix(k, :) ;
end
metagene_Normalised = metagene_Normalised ./ size(Normalised_NETseq_Matrix, 1) ;

%% plot normalised metagene 
figure 
plot(metagene_Normalised)
title('normalised metagene')
%% NaN to zeros 
% Normalised_NETseq_Matrix(isnan(Normalised_NETseq_Matrix)) = 0; converts
% Nans to zeros rather than removing

%% Bin the data prior to normalisation 
%% averaging every 10bps Bins
%make the matrix 
NETseq_Bins10bp_avg = zeros(size(NETseq_matrix_no_negative, 1), length_pulled_out/10);
for k = 1:size(NETseq_matrix_no_negative, 1)
    for x = 1:length_pulled_out/10
        Bin = mean(NETseq_matrix_no_negative(k,10*x-9:10*x));
        NETseq_Bins10bp_avg(k, x) = Bin  ;
    end
end


%% Normalise Binned NETseq
Normalised_NETseq_Bins_Matrix = Shape_normalisation_function(NETseq_Bins10bp_avg); %should probably integrate this into the code but i thought i was cool making a function
% %remove NaNs
% for k= size(Normalised_NETseq_Bins_Matrix, 1):-1:1 
%     
%     if  sum(isnan(Normalised_NETseq_Bins_Matrix(k, :)), 2)>0 %gives logical showing which genes contain any NaNs 
%         Normalised_NETseq_Bins_Matrix(k,:) = [] ;
% %     else not needed 
% %         Normalised_NETseq_Matrix(k,:) = Normalised_NETseq_Matrix(k,:);
%     end
%     
% end

%% average gene calculation
meta_gene = all_genes ./ (2*Total_gene_number) ;
figure
plot(meta_gene)

%%
toc
