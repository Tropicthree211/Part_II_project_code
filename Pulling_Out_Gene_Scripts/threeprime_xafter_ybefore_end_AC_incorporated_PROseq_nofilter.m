%% load in annotations 
tic
%% Load in the processed annotation data

%NB when using alternative annotations, gene positions should always be 
%leftmost to rightmost irrespective of strand
%If this is not the case, need to reprocess the names here to conform

annotationloaddir = ...
  '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/TIF_SEQ_ANNOTATIONS/nature12121-s2/'; 
% annotationloaddir = ...
    % '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ENSEMBL/S_CEREVISIAE_ANNOTATIONS/';
    
extrasavepiece = '';


%loadfile1 = 'gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat';
%%original file with Andrew's annotations

%loadfile1 = 'Directly_downstream_removed_gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat'; 
%Original annotations with cramped end genes removed 
loadfile1 = 'gene_positions_pelechano_TIFseq_no500paraconv.mat';


%loadfile1 = 'gene_positions_pelechano_TIFseq_nogene500.mat';
%TIF-seq annotations with cramed end genes removed 

% changed to specify new downstream removed genes
%above also has minimum gene size of 500 rather than 1000 as didn't seem
%relevant to exclude based on size given that I'm looking at the end

%loadfile2 = 'gene_ids_Saccharomyces_cerevisiae__R64_1_1__101.mat';
%original andrew's id annotations 

%should insert loadfile for Ids with cramped ends removed here 

loadfile2 = 'gene_names_pelechano_TIFseq_no500paraconv.mat';
%loadfile2 = 'gene_names_pelechano_TIFseq_nogene500.mat';
%TIF-seq annotations

%loadfile3 = 'gene_positions_for_nogene_mean_Saccharomyces_cerevisiae__R64_1_1__101.mat';
%added loadfile3 to avoid having to integrate all components for background
%calculation into this code 
loadfile3 = 'gene_positions_pelechano_TIFseq_for_nogene_mean';

load([annotationloaddir loadfile1])
load([annotationloaddir loadfile2])
load([annotationloaddir loadfile3])
disp('Loaded annotation files.')

 %load gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat %commented
 %%out to load new annotations below 
%  load Directly_downstream_removed_gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat
 NoChromosomes = 16 ;
 % extra annotations can be loaded found in andrew's code 
%% Start with WT 

% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/LIS_PRO_SEQ/LisLab_PROseq_data/';
% 
% PROseq_sample = importdata([dataloaddir 'Lis_PROseq_WT.mat']);
% 
% 
% disp('Loaded in WT PROseq data')
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
dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/LIS_PRO_SEQ/LisLab_PROseq_data/';
% 
 PROseq_sample = importdata([dataloaddir 'Lis_PROseq_spt4.mat']);

% 
disp('Loaded in LisLab spt4 delete PROseq data')

% NOTES:
% The correlation between the two samples was very good
% The average levels in the genes lie significantly off the diagonal - this
% might indicate a levelling issue and a more careful average might be
% needed in the longer run.
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
%% Remove genes that go outside the boundaries of the matrix % NEW for this script 
for cctr=1:1:NoChromosomes
    for sctr= 1:1:2
        for gctr = size(gene_positions{cctr, sctr}, 1):-1:1
            if gene_positions{cctr, sctr}(gctr, end) > size(PROseq_sample{cctr,sctr}, 1)
                
                gene_positions{cctr, sctr}(gctr, :) = [];
            end
        end
    end
end
%% Remove genes that go outside the boundaries of the matrix from the nogene mean positions matrix  % NEW for this script
for cctr=1:1:NoChromosomes
    for sctr= 1:1:2
        for gctr= size(gene_positions_for_nogene_mean{cctr, sctr}, 1):-1:1
            if gene_positions_for_nogene_mean{cctr, sctr}(gctr, end) > size(PROseq_sample{cctr,sctr}, 1)
                
                gene_positions_for_nogene_mean{cctr, sctr}(gctr, :) = [];
            end
        end
    end
end
%% forward
fail_counter_forward = 0;
PRO_Seq_matrix_forward = zeros(gene_number_forward,(ybefore+xafter)) ;
o = 0 ;
for i = 1:NoChromosomes

    %Isolating one chromosome%
    chromosomei_PROseq = PROseq_sample(i,1) ;
    chromosomei_genepos = gene_positions(i,1) ;
    
    
    %pulling out ybeforebp before the end of gene%
    chromosomei_startpoints = chromosomei_genepos{1,1}(:,2) - ybefore;
    
    %start pos + length pulled out minus 1 gives ybefore xafter end%
    chromosomei_endpoints = chromosomei_startpoints + (length_pulled_out-1) ;
    
    %pulling out the individual genes%
    z = size(gene_positions{i,1}(:,1), 1) ;
   
    for k = 1:z 
        
        o = o + 1 ; 
        try
        %add all genes together in single matrix for averaging later
        all_genes = all_genes + chromosomei_PROseq{1,1}(chromosomei_startpoints(k,1):chromosomei_endpoints(k,1),1) ;
        
        %construct a matrix with all gene profiles in
        PRO_Seq_matrix_forward(o, :) = chromosomei_PROseq{1,1}(chromosomei_startpoints(k,1):chromosomei_endpoints(k,1),1) ;
        catch 
            fail_counter_forward = fail_counter_forward +1
        end
       
    end
end
    
%% reverse
PRO_Seq_matrix_reverse = zeros(gene_number_reverse,(ybefore+xafter)) ;
w = 0 ;
fail_counter_reverse = 0 ;

for j = 1:NoChromosomes
   
    %Isolating one chromosome%
    chromosomej_PROseq = PROseq_sample(j,2) ;
    chromosomej_genepos = gene_positions(j,2) ;
    
    
    %pulling out the ybeforebp before end of genes%
    chromosomej_startpoints = chromosomej_genepos{1,1}(:,1) + ybefore ;
    %start points - x+y gives last x bases of reversed gene %
    chromosomej_endpoints = chromosomej_startpoints - (ybefore+xafter - 1) ;
    
    %pulling out the individual genes%
    x = size(gene_positions{j,2}(:,1), 1) ;
    
    for k = 1:x
        
        
        w = w + 1  ;
%        try
        %add to total signal matrix for all genes
        all_genes = all_genes + chromosomej_PROseq{1,1}(chromosomej_startpoints(k,1):-1:chromosomej_endpoints(k,1),1) ;
      
     
        %construct a matrix with all gene profiles in
        PRO_Seq_matrix_reverse(w, :) = chromosomej_PROseq{1,1}(chromosomej_startpoints(k,1):-1:chromosomej_endpoints(k,1),1);
    
%        catch
%            fail_counter_reverse = fail_counter_reverse + 1 
%        end
     end
    
end


%% Look at metagene forward and reverse 

metagene_forward = zeros(1, length_pulled_out);
for k= 1:gene_number_forward
    metagene_forward = metagene_forward + PRO_Seq_matrix_forward(k, :) ;
end
metagene_forward = metagene_forward ./ gene_number_forward ;

metagene_reverse = zeros(1, length_pulled_out);
for k= 1:gene_number_reverse
    metagene_reverse = metagene_reverse + PRO_Seq_matrix_reverse(k, :) ;
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
PRO_Seq_Matrix = vertcat(PRO_Seq_matrix_forward,PRO_Seq_matrix_reverse) ;

%% forward  repeat for biological repeat then average to get overall PROseq matrix 
% PRO_Seq_matrix_forward2 = zeros(gene_number_forward,length_pulled_out ) ;
% o = 0 ;
% for i = 1:NoChromosomes
% 
%     %Isolating one chromosome%
%     chromosomei2_PROseq = PROseq_sample2(i,1) ;
%     chromosomei2_genepos = gene_positions(i,1) ;
%     
%     
%     %pulling out 250bp before the end of gene%
%     chromosomei2_startpoints = chromosomei2_genepos{1,1}(:,2) - ybefore ;
%     
%     %start pos + x+y gives last y bases + x after%
%     chromosomei2_endpoints = chromosomei2_startpoints + (length_pulled_out-1) ;
%     
%     %pulling out the individual genes%
%     z = size(gene_positions{i,1}(:,1), 1) ;
%    
%     for k = 1:z 
%         
%         o = o + 1 ; 
%         %add all genes together in single matrix for averaging later
%         all_genes = all_genes + chromosomei2_PROseq{1,1}(chromosomei2_startpoints(k,1):chromosomei2_endpoints(k,1),1) ;
%         
%         %construct a matrix with all gene profiles in
%         PRO_Seq_matrix_forward2(o, :) = chromosomei2_PROseq{1,1}(chromosomei2_startpoints(k,1):chromosomei2_endpoints(k,1),1) ;
%     
%     end
% end
%     
% %% reverse
% PRO_Seq_matrix_reverse2 = zeros(gene_number_reverse,length_pulled_out ) ;
% w = 0 ;
% for j = 1:NoChromosomes
%    
%     %Isolating one chromosome%
%     chromosomej2_PROseq = PROseq_sample2(j,2) ;
%     chromosomej2_genepos = gene_positions(j,2) ;
%     
%     
%     %pulling out the 250bp before end of genes%
%     chromosomej2_startpoints = chromosomej2_genepos{1,1}(:,1) + ybefore  ;
%     %start points - 499 gives last 250 bases of reversed gene + 250 after%
%     chromosomej2_endpoints = chromosomej2_startpoints - (length_pulled_out - 1) ;
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
%         all_genes = all_genes + chromosomej2_PROseq{1,1}(chromosomej2_startpoints(k,1):-1:chromosomej2_endpoints(k,1),1) ;
%     
%      
%         %construct a matrix with all gene profiles in
%         PRO_Seq_matrix_reverse2(w, :) = chromosomej2_PROseq{1,1}(chromosomej2_startpoints(k,1):-1:chromosomej2_endpoints(k,1),1);
%     end
% end
% 
% %% Look at metagene forward and reverse 
% metagene_forward2 = zeros(1, length_pulled_out);
% for k= 1:gene_number_forward
%     metagene_forward2 = metagene_forward2 + PRO_Seq_matrix_forward(k, :) ;
% end
% metagene_forward2 = metagene_forward2 ./ gene_number_forward ;
% 
% metagene_reverse2 = zeros(1, length_pulled_out);
% for k= 1:gene_number_reverse
%     metagene_reverse2 = metagene_reverse2 + PRO_Seq_matrix_reverse(k, :) ;
% end
% metagene_reverse2 = metagene_reverse2 ./ gene_number_reverse ;
% 
%  %% plot metagenes 
% figure 
% plot(metagene_forward2)
% title('metagene forward repeat')
% figure
% plot(metagene_reverse2)
% title('metagene reverse repeat')
% 
% %% combine reverse and forward matrices
% PRO_Seq_Matrix2 = vertcat(PRO_Seq_matrix_forward2,PRO_Seq_matrix_reverse2) ;
% 
% 
% %% Average the two samples 
% PRO_Seq_Matrix = (PRO_Seq_Matrix1 + PRO_Seq_Matrix2) ./ 2 ;

%% Remove negative values 
PROseq_matrix_no_negative = PRO_Seq_Matrix ; 
PROseq_matrix_no_negative(PRO_Seq_Matrix<0) = 0;


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


%% Calculate mean across genome as 'background' 
% total_across_genome = 0;
% genome_size = 0 ;
% %below is for multiple sample averaging 
% % PROseq_composite = cell(size(PROseq_sample1));
% % for cctr = 1:1:NoChromosomes %chromosme counter 
% %     for sctr = 1:1:2 %strand counter 
% %         
% %         PROseq_composite{cctr,sctr} = ...
% %             nanmean( ...
% %             [PROseq_sample1{cctr,sctr},PROseq_sample2{cctr,sctr}],2 ...
% %             );
% %         
% %     end
% % end
% %calculate mean from composite cell array 
% for cctr = 1:NoChromosomes
%     for sctr = 1:2
%         total_across_genome = total_across_genome + sum(PROseq_sample{cctr, sctr}) ;% replaced PROseq_composite above as only one sample 
%         
%         genome_size = genome_size + sum(numel(PROseq_sample{cctr, sctr})) ;% ^same
%         
%         
%     end
% end
% 
% mean_across_genome = total_across_genome ./ genome_size ;
% 
% %% Calculate mean across all non-gene areas of genome, better background.
% %needs unfiltered gene positions before removal of any genes see
% %ImportENSEMBL_S_Cerevisiae_find_areas_outside_genes 
% PROseq_nogenes = PROseq_sample; %replaced PROseq_composite as no repeats 
% failed_removal_counter = 0 ;
% for cctr = 1:1:NoChromosomes
%     for sctr = 1:1:2
%     for gctr = size(gene_positions_for_nogene_mean{cctr,sctr}, 1):-1:1 %need to input before removal of genes by overlaps etc 
%         PROseq_nogenes{cctr, sctr}(gene_positions_for_nogene_mean{cctr, sctr}(gctr, 1):gene_positions_for_nogene_mean{cctr, sctr}(gctr, 2), :) = [];
% %         try PROseq_nogenes(gene_positions_for_nogene_mean{cctr, sctr}(gctr, 1):gene_positions_for_nogene_mean{cctr, sctr}(gctr, 2), :) = [];
%  %it worked without failing so didn't need try loop
% %         catch
% %             failed_removal_counter = failed_removal_counter +1 ; %if genes overlap 
%             %will attempt to remove already removed area may fail 
%             
%             %reinstated try loop to deal with gene out of bounds on
%             %chromosome 10 reverse strand 
% %         end
%     end
%    
%     end
% end
% 
% total_across_genome_nogene = 0;
% genome_size_nogene = 0;
% %Calculate mean from nogenes PROseq array
% for cctr = 1:NoChromosomes
%     for sctr = 1:2
%         total_across_genome_nogene = total_across_genome_nogene + sum(PROseq_nogenes{cctr, sctr}) ;
%         
%         genome_size_nogene = genome_size_nogene + sum(numel(PROseq_nogenes{cctr, sctr})) ;
%         
%         
%     end
% end
% 
% 
% mean_across_genome_nogene = total_across_genome_nogene ./ genome_size ;
% disp('done')
% 
% %% remove values below background mean 
% list_remaining_genes = (1:1:size(PROseq_matrix_no_negative, 1))'; %Index to apply same selection at 5' end later 
% 
% for g = size(PROseq_matrix_no_negative, 1):-1:1
%     %if nanmean(PROseq_matrix_no_negative(g, :), 2)<mean_across_genome
%     %%found above = original filter
%      if nanmean(PROseq_matrix_no_negative(g, :), 2)<mean_across_genome_nogene %found above new filter lose less signal   
%         list_remaining_genes(g, :) =  [];
%         
%         PROseq_matrix_no_negative(g, :) = [] ;
%     end
% end

 %commented out all the above to avoid filtering the data at all beyond the
 %annotation file so that filtering can occur later using WT-NETseq clusters. 

%% Normalise no negative PROseq Matrix
Normalised_PROseq_Matrix = Shape_normalisation_function(PROseq_matrix_no_negative); %only works if shape normalisation function is in current folder, need to find out how to specify path to function but for now will leave


% % remove NaNs
%  for k= size(Normalised_PROseq_Matrix, 1):-1:1 
%     
%     if  sum(isnan(Normalised_PROseq_Matrix(k, :)), 2)>0 %gives logical showing which genes contain any NaNs 
%         Normalised_PROseq_Matrix(k,:) = [] ;
% %     else 
% %         Normalised_PROseq_Matrix(k,:) = Normalised_PROseq_Matrix(k,:);
%     end
%     
%  end
 
 %% Find Normalised metagene 

 metagene_Normalised = zeros(1, length_pulled_out);
for k= 1:size(Normalised_PROseq_Matrix, 1)
    metagene_Normalised = metagene_Normalised + Normalised_PROseq_Matrix(k, :) ;
end
metagene_Normalised = metagene_Normalised ./ size(Normalised_PROseq_Matrix, 1) ;

%% plot normalised metagene 
figure 
plot(metagene_Normalised)
title('normalised metagene')
%% NaN to zeros 
% Normalised_PROseq_Matrix(isnan(Normalised_PROseq_Matrix)) = 0; converts
% Nans to zeros rather than removing

%% Bin the data prior to normalisation 
%% averaging every 10bps Bins
%make the matrix 
PROseq_Bins10bp_avg = zeros(size(PROseq_matrix_no_negative, 1), length_pulled_out/10);
for k = 1:size(PROseq_matrix_no_negative, 1)
    for x = 1:length_pulled_out/10
        Bin = mean(PROseq_matrix_no_negative(k,10*x-9:10*x));
        PROseq_Bins10bp_avg(k, x) = Bin  ;
    end
end


%% Normalise Binned PROseq
Normalised_ufdata_Bins_Matrix = Shape_normalisation_function(PROseq_Bins10bp_avg); %should probably integrate this into the code but i thought i was cool making a function
% %remove NaNs
% for k= size(Normalised_PROseq_Bins_Matrix, 1):-1:1 
%     
%     if  sum(isnan(Normalised_PROseq_Bins_Matrix(k, :)), 2)>0 %gives logical showing which genes contain any NaNs 
%         Normalised_PROseq_Bins_Matrix(k,:) = [] ;
% %     else not needed 
% %         Normalised_PROseq_Matrix(k,:) = Normalised_PROseq_Matrix(k,:);
%     end
%     
% end

%% average gene calculation
meta_gene = all_genes ./ (2*Total_gene_number) ;
figure
plot(meta_gene)
%%
save('Unfiltered_PROseqspt4del_TIF_nogene500pc100b200a.mat','Normalised_ufdata_Bins_Matrix')
unix('mv Unfiltered_PROseqspt4del_TIF_nogene500pc100b200a.mat Unfiltered_data')
%%
toc
