
%here to line ~750 is copy of loading and extraction into Normalised binned
%genes matrix 
%% load in annotations 
tic
%% Load in the processed annotation data

%NB when using alternative annotations, gene positions should always be 
%leftmost to rightmost irrespective of strand
%If this is not the case, need to reprocess the names here to conform
% 
% annotationloaddir = ...
%   '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/TIF_SEQ_ANNOTATIONS/nature12121-s2/'; 
% % annotationloaddir = ...
%     % '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ENSEMBL/S_CEREVISIAE_ANNOTATIONS/';
%     
% extrasavepiece = '';
% 
% 
% %loadfile1 = 'gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat';
% %%original file with Andrew's annotations
% 
% %loadfile1 = 'Directly_downstream_removed_gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat'; 
% %Original annotations with cramped end genes removed 
% loadfile1 = 'gene_positions_pelechano_TIFseq_no500paraconv.mat';
% annotations = 'nogene500pcTIFseq'
% %removed genes on opposing strands which converged 
% 
% 
% %loadfile1 = 'gene_positions_pelechano_TIFseq_nogene500.mat';
% %TIF-seq annotations with cramped end genes removed 
% 
% % changed to specify new downstream removed genes
% %above also has minimum gene size of 500 rather than 1000 as didn't seem
% %relevant to exclude based on size given that I'm looking at the end
% 
% %loadfile2 = 'gene_ids_Saccharomyces_cerevisiae__R64_1_1__101.mat';
% %original andrew's id annotations 
% 
% %should insert loadfile for Ids with cramped ends removed here 
% 
% loadfile2 ='gene_names_pelechano_TIFseq_no500paraconv.mat';
% %loadfile2 = 'gene_names_pelechano_TIFseq_nogene500.mat';
% %TIF-seq annotations
% 
% %loadfile3 = 'gene_positions_for_nogene_mean_Saccharomyces_cerevisiae__R64_1_1__101.mat';
% %added loadfile3 to avoid having to integrate all components for background
% %calculation into this code 
% loadfile3 = 'gene_positions_pelechano_TIFseq_for_nogene_mean';
% 
% load([annotationloaddir loadfile1])
% load([annotationloaddir loadfile2])
% load([annotationloaddir loadfile3])
% disp('Loaded annotation files.')
% 
%  %load gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat %commented
%  %%out to load new annotations below 
% %  load Directly_downstream_removed_gene_positions_Saccharomyces_cerevisiae__R64_1_1__101.mat
%  NoChromosomes = 16 ;
%  % extra annotations can be loaded found in andrew's code 

Current_modelling_job = 21
 %% Harry annotations
 NoChromosomes=16
 annotationloaddir = ...
    '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/OLD_ANNOTATIONS/Scer/';

extrasavepiece = 'ALTERNATIVE_ANNOTATIONS/HARRY/';


loadfile1 = 'harry_annotations.mat';

load([annotationloaddir loadfile1])


disp('Loaded annotation files.')

%This needs much extra processing - it's all in one big array for starters
% and the orientation needs checking amongst other things

%Split out the original file into two pieces to match the others

gene_ids = cell(NoChromosomes,2);
gene_positions = cell(NoChromosomes,2);

for cctr=1:1:NoChromosomes
   
    for sctr=1:1:2
       
        gene_ids{cctr,sctr} = harry_annotations{cctr,sctr+2};
        gene_positions{cctr,sctr} = harry_annotations{cctr,sctr};
        
    end
    
    %flip the orientation of the gnees on the second strand
    gene_positions{cctr,2} = fliplr(gene_positions{cctr,2});
    
end
annotations = 'HarrynoCUTs'
%% remove Cuts and Xuts
% 
for sctr = 1:1:2
    for cctr = 1:1:NoChromosomes
%         nongenelogical = zeros(size(gene_ids{cctr,sctr}(),2), 1);
        for gctr = size(gene_ids{cctr,sctr}, 2):-1:1
            splitname = split(gene_ids{cctr,sctr}{1,gctr}, '');
            if splitname{2} ~= 'Y'
                gene_ids{cctr,sctr}(:, gctr) = [] ;
                gene_positions{cctr,sctr}(gctr,:) = [];
                
            end
            
        end
        
        
    end
end
%% Glugal 15
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/Glugal_shift_raw/15min/'
% 
% NETseq_sample1 = importdata([dataloaddir 'Glugal_15_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'Glugal_15_NETseq_2.mat']);
% 
% disp('loaded in glugal shift 15 min data')
% xafter = 300;
% ybefore = 300;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant1 = 'Glugal15'
% datatype = 'NETseq'
%% Glugal 60
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/Glugal_shift_raw/60min/'
% 
% NETseq_sample1 = importdata([dataloaddir 'Glugal_60_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'Glugal_60_NETseq_2.mat']);
% 
% disp('loaded in glugal shift 60 min data')
% xafter = 300;
% ybefore = 300;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant1 = 'Glugal60'
% datatype = 'NETseq'
%% Start with WT PROseq

% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/LIS_PRO_SEQ/LisLab_PROseq_data/';
% 
% %NETseq_sample1 = importdata([dataloaddir 'Lis_PROseq_WT.mat']);
% %NETseq_sample2 = importdata([dataloaddir 'Lis_PROseq_WT.mat']); %duplicated the same dataset because shouldn't make a difference and can use same code that I have edited.
% NETseq_sample1 = importdata([dataloaddir 'LIS_PROseq_WT_amended.mat']); %increased number of zeros to reflect size of chromosomes and to allow use of Harry's annotations. 
% NETseq_sample2 = importdata([dataloaddir 'LIS_PROseq_WT_amended.mat']); %duplicated the same dataset because shouldn't make a difference and can use same code that I have edited.
% disp('Loaded in LIS WT PROseq data')
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
% mutant1 = 'WT'
% datatype = 'PROseq'
%mutants to load can be found in andrew's code 
%% spt4 delete
% 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/LIS_PRO_SEQ/LisLab_PROseq_data/'
% % 
%  NETseq_sample = importdata([dataloaddir 'Lis_PROseq_spt4.mat']);%  
% NETseq_sample1 = importdata([dataloaddir 'Lis_PROseq_spt4del_amended.mat']);   %increased number of zeros to reflect size of chromosomes and to allow use of Harry's annotations.
% NETseq_sample2 = importdata([dataloaddir 'Lis_PROseq_spt4del_amended.mat']);
%  NETseq_sample2 = importdata([dataloaddir 'ulku_spt4del_PROseq_2.mat']);

% isn't actually NETseq sample just calling it that so it plays nice with
% my code, should eventually move to generic naming of intermediates but
% functionally unimportant.
% % 
% disp('Loaded in LIS spt4 delete PROseq data')
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
%mutant1 = 'spt4del'
%datatype = 'PROseq'


%% Start with WT 

dataloaddir = '/data/FIXED_DATA/FIXED_TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/notag_normalised/WT_nnt/';

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
xafter = 300;
ybefore = 300;
length_pulled_out = xafter+ybefore;
all_genes = zeros((xafter+ybefore),1) ;
mutant1 = 'WT'
datatype = 'NETseq'
% % %mutants to load can be found in andrew's code 

%% set1 delete 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
% 
% NETseq_sample1 = importdata([dataloaddir 'weissman_set1del_NETseq.mat']);
% NETseq_sample2 = importdata([dataloaddir 'weissman_set1del_NETseq.mat']);
% 
% disp('Loaded in W set1del NETseq data')
% 
% xafter = 300;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant1 = 'set1del'
% datatype = 'NETseq'
%only one file for each no repeats but loading in two should mean the
%composite is the same as the original.
%% rco1 delete 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
% 
% NETseq_sample1 = importdata([dataloaddir 'weissman_rco1del_NETseq.mat']);
% NETseq_sample2 = importdata([dataloaddir 'weissman_rco1del_NETseq.mat']);
% 
% disp('Loaded in W rco1del NETseq data')
% 
% xafter = 300;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant1 = 'rco1del'
% datatype = 'NETseq'
% %only one file for each no repeats but loading in two should mean the
% %composite is the same as the original.

%% Weissman WT 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
% 
% NETseq_sample1 = importdata([dataloaddir 'weissman_WT_NETseq.mat']);
% NETseq_sample2 = importdata([dataloaddir 'weissman_WT_NETseq.mat']);
% 
% disp('Loaded in W WT NETseq data')
% 
% xafter = 300;
% ybefore = 300;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant1 = 'weissWT'
% datatype = 'NETseq'
%only one file for each no repeats but loading in two should mean the
%composite is the same as the original. 
%% set2 delete 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
% 
% NETseq_sample1 = importdata([dataloaddir 'weissman_set2del_NETseq.mat']);
% NETseq_sample2 = importdata([dataloaddir 'weissman_set2del_NETseq.mat']);
% 
% disp('Loaded in W set2del NETseq data')
% xafter = 300;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant1 = 'set2del'
% datatype = 'NETseq'
%% eaf3 delete 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
% 
% NETseq_sample1 = importdata([dataloaddir 'weissman_eaf3del_NETseq.mat']);
% NETseq_sample2 = importdata([dataloaddir 'weissman_eaf3del_NETseq.mat']);
% 
% disp('Loaded in W eaf3del NETseq data')
% 
% xafter = 300;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant1 = 'eaf3del'
% datatype = 'NETseq'
%% spt4 delete
% 
% dataloaddir = '/data/FIXED_DATA/FIXED_TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/notag_normalised/spt4del_nnt/';
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
% xafter = 300;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
%mutant1 = 'spt4del'
% datatype = 'NETseq'
%% dst1 delete

% dataloaddir = '/data/FIXED_DATA/FIXED_TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/notag_normalised/dst1del_nnt/';
% 
% NETseq_sample1 = importdata([dataloaddir 'ulku_dst1del_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'ulku_dst1del_NETseq_2.mat']);
% 
% disp('Loaded in UU dst1 delete NETseq data')
% 
% %NOTES:
% %The correlation between the two samples was very good
% %The average levels in each were reasonably close to the diagonal
% xafter = 300;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant1 = 'dst1del'
% datatype = 'NETseq'
%% xrn1 delete

% dataloaddir = '/data/FIXED_DATA/FIXED_TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/notag_normalised/xrn1del_nnt/';
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
% xafter = 300;
% ybefore = 300;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant1 = 'xrn1del'
% datatype = 'NETseq'
%% SPT4 Anchor Away

% dataloaddir = '/data/FIXED_DATA/FIXED_TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/AnchorAway/AA_PROCESSED_NETSEQ/SPT4AA/';
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
%mutant1 = 'spt4AA'
% datatype = 'NETseq'
 %% SPT5 Anchor Away
% 
% dataloaddir = '/data/FIXED_DATA/FIXED_TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/AnchorAway/AA_PROCESSED_NETSEQ/SPT5AA/';
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
% mutant1 = 'spt5AA'
% datatype = 'NETseq'
%% SPT4/5 DMSO (Treating as repeats but strictly they are different)
% 
% dataloaddir = '/data/FIXED_DATA/FIXED_TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/AnchorAway/AA_PROCESSED_NETSEQ/SPT4AA_SPT5AA_DMSO_COMBO/';
% 
% NETseq_sample1 = importdata([dataloaddir 'ulku_SPT4AADMSO_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'ulku_SPT5AADMSO_NETseq_1.mat']);
% 
% disp('Loaded in UU SPT5 Anchor Away NETseq data')
% 
% %NOTES:
% %Correlation between these two samples (which should be similar but are
% %technically distinct) is very good
% %The average levels in genes also lies reasonably close to the diagonal
%  xafter = 200;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant1 = 'spt4AA5AADMSO
% datatype = 'NETseq''

%% calculating total gene number

gene_number_forward = 0 ; %start with matrix of zeros
gene_number_reverse = 0 ;
for d = 1:NoChromosomes
    
    gene_number_forward = gene_number_forward + sum(length(gene_positions{d,1}(:,1)))  ;
    gene_number_reverse = gene_number_reverse + sum(length(gene_positions{d,2}(:,1))) ;
    %iteratively add up all the genes present
    
end
Total_gene_number = gene_number_forward + gene_number_reverse ;
%% forward strand
NET_Seq_matrix_forward = zeros(gene_number_forward,(ybefore+xafter)) ;
o = 0 ;
for cctr = 1:NoChromosomes

    %Isolating one chromosome%
    chromosomei_NETseq = NETseq_sample1(cctr,1) ;
    chromosomei_genepos = gene_positions(cctr,1) ;
    
    
    %pulling out ybeforebp before the end of gene%
    chromosomei_startpoints = chromosomei_genepos{1,1}(:,2) - ybefore;
    
    %start pos + length pulled out minus 1 gives ybefore xafter end%
    chromosomei_endpoints = chromosomei_startpoints + (length_pulled_out-1) ;
    
    %pulling out the individual genes%
    z = size(gene_positions{cctr,1}(:,1), 1) ;
   
    for gctr = 1:z 
        
        o = o + 1 ; 
        %add all genes together in single matrix for averaging later
        all_genes = all_genes + chromosomei_NETseq{1,1}(chromosomei_startpoints(gctr,1):chromosomei_endpoints(gctr,1),1) ;
        
        %construct a matrix with all gene profiles in
        NET_Seq_matrix_forward(o, :) = chromosomei_NETseq{1,1}(chromosomei_startpoints(gctr,1):chromosomei_endpoints(gctr,1),1) ;
    
    end
end
    
%% reverse strand
NET_Seq_matrix_reverse = zeros(gene_number_reverse,(ybefore+xafter)) ;
w = 0 ;
for cctr = 1:NoChromosomes
   
    %Isolating one chromosome%
    chromosomej_NETseq = NETseq_sample1(cctr,2) ;
    chromosomej_genepos = gene_positions(cctr,2) ;
    
    
    %pulling out the ybeforebp before end of genes%
    chromosomej_startpoints = chromosomej_genepos{1,1}(:,1) + ybefore ;
    %start points - x+y gives last x bases of reversed gene %
    chromosomej_endpoints = chromosomej_startpoints - (ybefore+xafter - 1) ;
    
    %pulling out the individual genes%
    x = size(gene_positions{cctr,2}(:,1), 1) ;
    
    for gctr = 1:x
        
        
        w = w + 1  ;
       
        %add to total signal matrix for all genes
        all_genes = all_genes + chromosomej_NETseq{1,1}(chromosomej_startpoints(gctr,1):-1:chromosomej_endpoints(gctr,1),1) ;
    
     
        %construct a matrix with all gene profiles in
        NET_Seq_matrix_reverse(w, :) = chromosomej_NETseq{1,1}(chromosomej_startpoints(gctr,1):-1:chromosomej_endpoints(gctr,1),1);
    end
end


%% Look at metagene forward and reverse 

metagene_forward = zeros(1, length_pulled_out);
for gctr= 1:gene_number_forward
    metagene_forward = metagene_forward + NET_Seq_matrix_forward(gctr, :) ;
end
metagene_forward = metagene_forward ./ gene_number_forward ;

metagene_reverse = zeros(1, length_pulled_out);
for gctr= 1:gene_number_reverse
    metagene_reverse = metagene_reverse + NET_Seq_matrix_reverse(gctr, :) ;
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

%% forward strand repeat for biological repeat then average to get overall NETseq matrix 
NET_Seq_matrix_forward2 = zeros(gene_number_forward,length_pulled_out ) ;
o = 0 ;
for cctr = 1:NoChromosomes

    %Isolating one chromosome%
    chromosomei2_NETseq = NETseq_sample2(cctr,1) ;
    chromosomei2_genepos = gene_positions(cctr,1) ;
    
    
    %pulling out 250bp before the end of gene%
    chromosomei2_startpoints = chromosomei2_genepos{1,1}(:,2) - ybefore ;
    
    %start pos + x+y gives last y bases + x after%
    chromosomei2_endpoints = chromosomei2_startpoints + (length_pulled_out-1) ;
    
    %pulling out the individual genes%
    z = size(gene_positions{cctr,1}(:,1), 1) ;
   
    for gctr = 1:z 
        
        o = o + 1 ; 
        %add all genes together in single matrix for averaging later
        all_genes = all_genes + chromosomei2_NETseq{1,1}(chromosomei2_startpoints(gctr,1):chromosomei2_endpoints(gctr,1),1) ;
        
        %construct a matrix with all gene profiles in
        NET_Seq_matrix_forward2(o, :) = chromosomei2_NETseq{1,1}(chromosomei2_startpoints(gctr,1):chromosomei2_endpoints(gctr,1),1) ;
    
    end
end
    
%% reverse strand
NET_Seq_matrix_reverse2 = zeros(gene_number_reverse,length_pulled_out ) ;
w = 0 ;
for cctr = 1:NoChromosomes
   
    %Isolating one chromosome%
    chromosomej2_NETseq = NETseq_sample2(cctr,2) ;
    chromosomej2_genepos = gene_positions(cctr,2) ;
    
    
    %pulling out the 250bp before end of genes%
    chromosomej2_startpoints = chromosomej2_genepos{1,1}(:,1) + ybefore  ;
    %start points - 499 gives last 250 bases of reversed gene + 250 after%
    chromosomej2_endpoints = chromosomej2_startpoints - (length_pulled_out - 1) ;
                                                                                                                                                                                                                                                                                                                            
    %pulling out the individual genes%
    x = size(gene_positions{cctr,2}(:,1), 1) ;
    
    for gctr = 1:x
        
        
        w = w + 1  ;                                            

        %add to total signal matrix for all genes
        all_genes = all_genes + chromosomej2_NETseq{1,1}(chromosomej2_startpoints(gctr,1):-1:chromosomej2_endpoints(gctr,1),1) ;
    
     
        %construct a matrix with all gene profiles in
        NET_Seq_matrix_reverse2(w, :) = chromosomej2_NETseq{1,1}(chromosomej2_startpoints(gctr,1):-1:chromosomej2_endpoints(gctr,1),1);
    end
end

%% Look at metagene forward and reverse 
metagene_forward2 = zeros(1, length_pulled_out);
for gctr= 1:gene_number_forward
    metagene_forward2 = metagene_forward2 + NET_Seq_matrix_forward(gctr, :) ;
end
metagene_forward2 = metagene_forward2 ./ gene_number_forward ;

metagene_reverse2 = zeros(1, length_pulled_out);
for gctr= 1:gene_number_reverse
    metagene_reverse2 = metagene_reverse2 + NET_Seq_matrix_reverse(gctr, :) ;
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
gene_positions_for_nogene_mean = gene_positions;
for cctr = 1:1:NoChromosomes
    for sctr = 1:1:2
    for gctr = size(gene_positions_for_nogene_mean{cctr,sctr}, 1):-1:1 %need to input before removal of genes by overlaps etc 
        if gctr == 1
            NETseq_nogenes{cctr, sctr}(gene_positions_for_nogene_mean{cctr, sctr}(gctr, 2):-1:gene_positions_for_nogene_mean{cctr, sctr}(gctr, 1), :) = [];
            continue;
            %should avoid trying to index with gctr = 0 which will fail my
            %code
        end
        if gene_positions_for_nogene_mean{cctr, sctr}(gctr, 1) < gene_positions_for_nogene_mean{cctr, sctr}((gctr -1) , 2)
            
        gene_positions_for_nogene_mean{cctr, sctr}((gctr -1), 2) = gene_positions_for_nogene_mean{cctr, sctr}(gctr , 1) -1;
        
        %added if statement should ensure not removing stuff that has already been removed in case of overlapping genes. 
        end
        NETseq_nogenes{cctr, sctr}(gene_positions_for_nogene_mean{cctr, sctr}(gctr, 2):-1:gene_positions_for_nogene_mean{cctr, sctr}(gctr, 1), :) = []; %
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
            %GM commented out the above 04/02/21 to quickly use Harry's
            %annotations rather than constructing a nogene version removing
            %all the genes. %pretty sure i can just define gene positions
            %for nogene mean as gene positions because no longer removing
            %any genes before this point
%mean_across_genome_nogene = 0.2160 %answer from wild type NETseq data 
 mean_across_genome_nogene = total_across_genome_nogene ./ genome_size ;
% disp('done')

%% remove values below background mean 
list_remaining_genes = (1:1:size(NETseq_matrix_no_negative, 1))'; %Index to apply same selection at 5' end later and on other datasets

for g = size(NETseq_matrix_no_negative, 1):-1:1
    %if nanmean(NETseq_matrix_no_negative(g, :), 2)<mean_across_genome
    %%found above = original filter %04/02/21 reverted to original see new file Harry 
     if nanmean(NETseq_matrix_no_negative(g, :), 2)<mean_across_genome_nogene %found above new filter lose less signal   
        list_remaining_genes(g, :) =  [];
        
        NETseq_matrix_no_negative(g, :) = [] ;
    end
end
save(append(mutant1,'_actively_expressing_genes_list', annotations, '.mat'),'list_remaining_genes') %must use same annotations file in both datasets if this is to work. 
%must remember to change this filename 
 

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
for gctr= 1:size(Normalised_NETseq_Matrix, 1)
    metagene_Normalised = metagene_Normalised + Normalised_NETseq_Matrix(gctr, :) ;
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
for gctr = 1:size(NETseq_matrix_no_negative, 1)
    for x = 1:length_pulled_out/10
        Bin = mean(NETseq_matrix_no_negative(gctr,10*x-9:10*x));
        NETseq_Bins10bp_avg(gctr, x) = Bin  ;
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

%% Constructing the mathematical model 
%presize matrix 
%modelparamsmatrix = zeros(4, 10000); %didn't need to presize no for loops
% random param generation 
rng('shuffle')
A = 0+(2-0)*rand(10000,1);
rng('shuffle')
B = -1 + (1+1)*rand(10000,1);
rng('shuffle')
a = 0 + 100*rand(10000,1);
rng('shuffle')
w = 0 +60*rand(10000,1) ;
%one param per column, each one alongside each other in 

%modelparamsmatrix = perms([A B a w]) ; tries to put every variable with
%every other too big 
%% Create the matrix containg many different permutations of ABaw
%

Arand = A(randperm(numel(A),10000));
Brand = B(randperm(numel(B),10000));
arand = a(randperm(numel(a),10000));
wrand = w(randperm(numel(w),10000));

ABapermw = [Arand Brand arand wrand]; 

% Arep = repmat(A, size(B/10,1), 1);
% Brep = repelem(B, size(A/10,1) , 1);
% ApermB = [Arep Brep];
% 
% ABrep = repmat(ApermB, size(a,1), 1);
% arep = repelem(a, size(ApermB ,1) , 1);
% ABperma = [ABrep arep];
% 
% ABarep = repmat(ABperma, size(w ,1), 1);
% wrep = repelem(w, size(ABperma ,1) , 1);
% ABapermw = [ABarep wrep];



%wrong code showing attempts
% y = perms(A)
% x = perms(B) wrong 
% z = perms(

% for x = 1:size(ApermB,1)
%  ApermB(x, 1) = Arep(x,1)
%  ApermB(x, 2) = Brep(x,1)
%  
% end %Don't need for loop I think
%ApermB = zeros(size(A,1)*size(B,1), 2)
 %ABperma = zeros(size(ApermB,1)*size(a,1), 1) 
%ABapermw = zeros(size(ApermB,1)*size(a,1), 1) 


%% Coding model using the parameters 
profile_matrix = zeros(size(ABapermw, 1), size(Normalised_NETseq_Bins_Matrix, 2));
for smplctr = 1:1: size(profile_matrix, 1)
    for x = 1:size(profile_matrix,2)
        if x<ABapermw(smplctr, 4)
           profile_matrix(smplctr, x) = ABapermw(smplctr, 1) +ABapermw(smplctr, 2);
        
        else
            
            profile_matrix(smplctr, x) = ABapermw(smplctr, 1).*exp((-ABapermw(smplctr, 3)).*(x-ABapermw(smplctr, 4))) + ABapermw(smplctr, 2) ;
        end
    end
end

%% compare each gene to each profile find minimum squared euclidean distance profile  
mindistprofiles = zeros(size(Normalised_NETseq_Bins_Matrix,1), size(Normalised_NETseq_Bins_Matrix,2));
mindistvalues = zeros(size(Normalised_NETseq_Bins_Matrix,1), 1);

dist_profile = zeros(size(profile_matrix, 1), 1);


for gctr = 1:1:size(Normalised_NETseq_Bins_Matrix, 1)
    
    for tmpctr = 1:1: size(profile_matrix, 1)
%         geneplusprofile = vertcat(Normalised_NETseq_Bins_Matrix(gctr, :), profile_matrix(tmpctr, :));
        
       
%         dist_profile(tmpctr ,1) = mean(pdist(geneplusprofile, 'squaredeuclidean'));
          dist_profile(tmpctr ,1) = mean(pdist2(Normalised_NETseq_Bins_Matrix(gctr, :),profile_matrix(tmpctr, :) , 'squaredeuclidean')); 
          %should run faster if not vertcating every run
    end
    [mindistvalues(gctr,1) , I] = min(dist_profile); 
    mindistprofiles(gctr, :) = profile_matrix(I, :);
   
end
%%
save(append('/home/orie3770/Modelling/', 'exp_modeled_UU', mutant1, 'NETseq_profiles', string(Current_modelling_job), '.mat'), ...
'mindistprofiles', 'mindistvalues') %change name before saving to avoid overwriting

toc

%% Pull out genes with modelled exponential drops in a certain range. 
% load from previous run can improve random model generation later just
% seeing if it works 
%load('/home/orie3770/Modelling/exp_modeled_UUGlugal15NETseq_profiles13.mat')
load(append('/home/orie3770/Modelling/', 'exp_modeled_UU', mutant1, 'NETseq_profiles', string(Current_modelling_job)))

%% find the drop point 
tic
expdropclasses = zeros(size(Normalised_NETseq_Bins_Matrix, 1),1 ); % matrix in which exponential drop location class will be stored.
for gctr =  1:1:size(Normalised_NETseq_Bins_Matrix, 1) % go through each gene in turn
    for pctr = 2:1:size(Normalised_NETseq_Bins_Matrix, 2) % go through each position in the gene
        if mindistprofiles(gctr,pctr) < mindistprofiles(gctr, (pctr-1)) && (2<=pctr) && (pctr<6)
            expdropclasses(gctr, 1) = 1 ;
            
        elseif mindistprofiles(gctr,pctr) < mindistprofiles(gctr, (pctr-1)) && (6<=pctr) && (pctr<10)
            expdropclasses(gctr, 1) = 2 ;            
            
        elseif mindistprofiles(gctr,pctr) < mindistprofiles(gctr, (pctr-1)) && (10<=pctr) && (pctr<14)
            expdropclasses(gctr, 1) = 3 ;
            
        elseif mindistprofiles(gctr,pctr) < mindistprofiles(gctr, (pctr-1)) &&(14<=pctr) && (pctr<18)
            expdropclasses(gctr, 1) = 4 ;
            
        elseif mindistprofiles(gctr,pctr) < mindistprofiles(gctr, (pctr-1)) &&(18<=pctr) && (pctr<22)
            expdropclasses(gctr, 1) = 5 ;
            
        elseif mindistprofiles(gctr,pctr) < mindistprofiles(gctr, (pctr-1)) && (22<=pctr) && (pctr<26)
            expdropclasses(gctr, 1) = 6 ;
            
        elseif mindistprofiles(gctr,pctr) < mindistprofiles(gctr, (pctr-1)) && (26<=pctr) && (pctr<30)
            expdropclasses(gctr, 1) = 7 ;
            
            %extend size of area with more elseifs
            
       elseif mindistprofiles(gctr,pctr) < mindistprofiles(gctr, (pctr-1)) && (30<=pctr) && (pctr<34)
            expdropclasses(gctr, 1) = 8 ;
       
        elseif mindistprofiles(gctr,pctr) < mindistprofiles(gctr, (pctr-1)) && (34<=pctr) && (pctr<38)
            expdropclasses(gctr, 1) = 9 ;  
       
        elseif mindistprofiles(gctr,pctr) < mindistprofiles(gctr, (pctr-1)) && (38<=pctr) && (pctr<42)
            expdropclasses(gctr, 1) = 10; 
            
        elseif mindistprofiles(gctr,pctr) < mindistprofiles(gctr, (pctr-1)) && (42<=pctr) && (pctr<46)
            expdropclasses(gctr, 1) = 11; 
            
        elseif mindistprofiles(gctr,pctr) < mindistprofiles(gctr, (pctr-1)) && (46<=pctr) && (pctr<50)
            expdropclasses(gctr, 1) = 12; 
            
        elseif mindistprofiles(gctr,pctr) < mindistprofiles(gctr, (pctr-1)) && (50<=pctr) && (pctr<54)
            expdropclasses(gctr, 1) = 13;     
            
        elseif mindistprofiles(gctr,pctr) < mindistprofiles(gctr, (pctr-1)) && (54<=pctr) && (pctr<58)
            expdropclasses(gctr, 1) = 14;   
          
        elseif mindistprofiles(gctr,pctr) < mindistprofiles(gctr, (pctr-1)) && (58<=pctr) && (pctr<60)
            expdropclasses(gctr, 1) = 15;   
        end
    end
end
toc                        
disp('drops found')
%% Use expdropclass to make 10 metagenes
%index to pull out the different genes in each class
%Don't need to do this 
%the information is all there can just work with indexes
% expdropclass1 = Normalised_NETseq_Bins_Matrix(expdropclasses == 1, :) ;                                       
% expdropclass2 = Normalised_NETseq_Bins_Matrix(expdropclasses == 2, :) ;
% expdropclass3 = Normalised_NETseq_Bins_Matrix(expdropclasses == 3, :) ;
% expdropclass4 = Normalised_NETseq_Bins_Matrix(expdropclasses == 4, :) ;
% expdropclass5 = Normalised_NETseq_Bins_Matrix(expdropclasses == 5, :) ;
% expdropclass6 = Normalised_NETseq_Bins_Matrix(expdropclasses == 6, :) ;
% expdropclass7 = Normalised_NETseq_Bins_Matrix(expdropclasses == 7, :) ;
% expdropclass8 = Normalised_NETseq_Bins_Matrix(expdropclasses == 8, :) ;
% expdropclass9 = Normalised_NETseq_Bins_Matrix(expdropclasses == 9, :) ;
% expdropclass10 = Normalised_NETseq_Bins_Matrix(expdropclasses == 10, :) ;
% expdropclass11 = Normalised_NETseq_Bins_Matrix(expdropclasses == 11, :) ;
% expdropclass12 = Normalised_NETseq_Bins_Matrix(expdropclasses == 12, :) ;
% expdropclass13 = Normalised_NETseq_Bins_Matrix(expdropclasses == 13, :) ;
% expdropclass14 = Normalised_NETseq_Bins_Matrix(expdropclasses == 14, :) ;
% expdropclass15 = Normalised_NETseq_Bins_Matrix(expdropclasses == 15, :) ;

%% make metagene for each class 
tic
unix(append('mkdir Current_Modelling_Job', string(Current_modelling_job)))
% metagene_exp1 = zeros(1, length_pulled_out/10);
% for gctr= 1:size(expdropclass1,1) 
%     metagene_exp1 = metagene_exp1 + expdropclass1(gctr, :) ;
% end
% metagene_exp1 = metagene_exp1 ./size(expdropclass1,1) ;
% figure 
% plot(metagene_exp1)
% %title('Normalised metagenes for exponential curve modelled clusters 1-4')
% title(append('exponential drop 80-40 pre polyA', datatype, mutant1, annotations))
% title('exponential drop 80-40 pre polyA')
% legend(string(size(expdropclass1,1)))
% saveas(gcf, 'modelling_metagene_80t040before', 'jpeg')
% %hold ON
% %
% metagene_exp2 = zeros(1, length_pulled_out/10);
% for gctr= 1:size(expdropclass2,1) 
%     metagene_exp2 = metagene_exp2 + expdropclass2(gctr, :) ;
% end
% metagene_exp2 = metagene_exp2 ./size(expdropclass2,1) ;
% figure 
% plot(metagene_exp2)
% title(append('exponential drop 40-0 pre polyA',datatype, mutant1, annotations) )
% %title('exponential drop 40-0 pre polyA')
% legend(string(size(expdropclass2,1)))
% saveas(gcf, 'modelling_metagene_40to0before', 'jpeg')
% %
% metagene_exp3 = zeros(1, length_pulled_out/10);
% for gctr= 1:size(expdropclass3,1) 
%     metagene_exp3 = metagene_exp3 + expdropclass3(gctr, :) ;
% end
% metagene_exp3 = metagene_exp3 ./size(expdropclass3,1) ;
% figure 
% plot(metagene_exp3)
% title(append('exponential drop 0-40 post polyA',datatype, mutant1, annotations ))
% %title('exponential drop 0-40 post polyA')
% legend(string(size(expdropclass3,1)))
% saveas(gcf, 'modelling_metagene_0to40after', 'jpeg')
% %
% metagene_exp4 = zeros(1, length_pulled_out/10);
% for gctr= 1:size(expdropclass4,1) 
%     metagene_exp4 = metagene_exp4 + expdropclass4(gctr, :) ;
% end
% metagene_exp4 = metagene_exp4 ./size(expdropclass4,1) ;
% figure 
% plot(metagene_exp4)
% title(append('exponential drop 40-80 post polyA',datatype, mutant1, annotations))
% %title('exponential drop 40-80 post polyA')
% legend(string(size(expdropclass4,1)))
% saveas(gcf, 'modelling_metagene_40to80after', 'jpeg')
% %
% metagene_exp5 = zeros(1, length_pulled_out/10);
% for gctr= 1:size(expdropclass5,1) 
%     metagene_exp5 = metagene_exp5 + expdropclass5(gctr, :) ;
% end
% metagene_exp5 = metagene_exp5 ./size(expdropclass5,1) ;
% figure 
% plot(metagene_exp5)
% title(append('exponential drop 80-120 post polyA', datatype, mutant1, annotations))
% %title('exponential drop 80-120 post polyA')
% legend(string(size(expdropclass5,1)))
% saveas(gcf, 'modelling_metagene_80to120after', 'jpeg')
% %hold ON
% %
% metagene_exp6 = zeros(1, length_pulled_out/10);
% for gctr= 1:size(expdropclass6,1) 
%     metagene_exp6 = metagene_exp6 + expdropclass6(gctr, :) ;
% end
% metagene_exp6 = metagene_exp6 ./size(expdropclass6,1) ;
% figure 
% plot(metagene_exp6)
% title(append('exponential drop 120-160 post polyA', datatype, mutant1, annotations))
% %title('exponential drop 120-160 post polyA')
% legend(string(size(expdropclass6,1)))
% saveas(gcf, 'modelling_metagene_120to160after', 'jpeg')
% %
% metagene_exp7 = zeros(1, length_pulled_out/10);
% for gctr= 1:size(expdropclass7,1) 
%     metagene_exp7 = metagene_exp7 + expdropclass7(gctr, :) ;
% end
% metagene_exp7 = metagene_exp7 ./size(expdropclass7,1) ;
% figure 
% plot(metagene_exp7)
% title(append('exponential drop 160-200 post polyA', datatype, mutant1, annotations ))
% %title('exponential drop 160-200 post polyA')
% legend(string(size(expdropclass7,1)))
% saveas(gcf, 'modelling_metagene_160to200after', 'jpeg')
% 
% metagene_exp8 = zeros(1, length_pulled_out/10);
% for gctr= 1:size(expdropclass8,1) 
%     metagene_exp8 = metagene_exp8 + expdropclass8(gctr, :) ;
% end
% metagene_exp8 = metagene_exp8 ./size(expdropclass8,1) ;
% figure 
% plot(metagene_exp8)
% title(append('exponential drop 200-240 post polyA', datatype, mutant1, annotations ))
% %title('exponential drop 200-240 post polyA')
% legend(string(size(expdropclass8,1)))
% saveas(gcf, 'modelling_metagene_200-240after', 'jpeg')
% 
% metagene_exp9 = zeros(1, length_pulled_out/10);
% for gctr= 1:size(expdropclass9,1) 
%     metagene_exp9 = metagene_exp9 +expdropclass9(gctr, :) ;
% end
% metagene_exp9 = metagene_exp9 ./size(expdropclass9,1) ;
% figure 
% plot(metagene_exp9)
% title(append('exponential drop 240-280 post polyA', datatype, mutant1, annotations ))
% % title('exponential drop 240-280 post polyA')
% legend(string(size(expdropclass9,1)))
% saveas(gcf, 'modelling_metagene_240-280after', 'jpeg')
% 
% metagene_exp10 = zeros(1, length_pulled_out/10);
% for gctr= 1:size(expdropclass10,1) 
%     metagene_exp10 = metagene_exp10 + expdropclass10(gctr, :) ;
% end
% metagene_exp10 = metagene_exp10 ./size(expdropclass10,1) ;
% figure 
% plot(metagene_exp10)
% title(append('exponential drop 280-300 post polyA', datatype, mutant1, annotations ))
% %title('exponential drop 280-300 post polyA')
% legend(string(size(expdropclass10,1)))
% saveas(gcf, 'modelling_metagene_280-300after', 'jpeg')

% below is code to pull out mutant data and group according to expdropclass of first mutant above.  
%% Glugal 15
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/Glugal_shift_raw/15min/'
% 
% NETseq_sample1 = importdata([dataloaddir 'Glugal_15_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'Glugal_15_NETseq_2.mat']);
% 
% disp('loaded in glugal shift 15 min data')
% xafter = 300;
% ybefore = 300;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant = 'Glugal15'
% datatype1 = 'NETseq'
%% Glugal 60
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/Glugal_shift_raw/60min/'
% 
% NETseq_sample1 = importdata([dataloaddir 'Glugal_60_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'Glugal_60_NETseq_2.mat']);
% 
% disp('loaded in glugal shift 60 min data')
% xafter = 300;
% ybefore = 300;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant = 'Glugal60'
% datatype1 = 'NETseq'
%% Start with WT PROseq

% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/LIS_PRO_SEQ/LisLab_PROseq_data/';
% 
% NETseq_sample1 = importdata([dataloaddir 'Lis_PROseq_WT.mat']);
% NETseq_sample2 = importdata([dataloaddir 'Lis_PROseq_WT.mat']); %duplicated the same dataset because shouldn't make a difference and can use same code that I have edited.
% NETseq_sample1 = importdata([dataloaddir 'LIS_PROseq_WT_amended.mat']); %increased number of zeros to reflect size of chromosomes and to allow use of Harry's annotations. 
% NETseq_sample2 = importdata([dataloaddir 'LIS_PROseq_WT_amended.mat']); %duplicated the same dataset because shouldn't make a difference and can use same code that I have edited.
% disp('Loaded in LIS WT PROseq data')
% 
% NOTES:
% The correlation between the two samples was very good
% The average level at each gene fell onto the diagonal reasonably well
% 
% leave things generic such that other data can be loaded into the same
% variables and processed in the same way!
% 
% Currently saving the output with a generic name to the load directory 
% - this is fine when the original data is separated into different
% directories (for data sets where they are all in one - consider splitting
% up the data -- or at the very least, put the processed data into clearly
% labellled subdirectories)
% xafter = 300;
% ybefore = 300;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant = 'WT'
% datatype1 = 'PROseq'
% mutants to load can be found in andrew's code 
%% spt4 delete
% 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/LIS_PRO_SEQ/LisLab_PROseq_data/'
% % 
%  NETseq_sample = importdata([dataloaddir 'Lis_PROseq_spt4.mat']);
% % NETseq_sample2 = importdata([dataloaddir 'ulku_spt4del_PROseq_2.mat']);
% isn't actually NETseq sample just calling it that so it plays nice with
% my code, should eventually move to generic naming of intermediates but
% functionally unimportant.
% % 
% disp('Loaded in LIS spt4 delete PROseq data')
% 
% % NOTES:
% % The correlation between the two samples was very good
% % The average levels in the genes lie significantly off the diagonal - this
% % might indicate a levelling issue and a more careful average might be
% % needed in the longer run.
% xafter = 300;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
%mutant = 'spt4del'
%datatype1 = 'PROseq'
%% Start with WT 

dataloaddir = '/data/FIXED_DATA/FIXED_TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/notag_normalised/WT_nnt/';

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
xafter = 300;
ybefore = 300;
length_pulled_out = xafter+ybefore;
all_genes = zeros((xafter+ybefore),1) ;
mutant = 'WT'
datatype1 = 'NETseq'
%mutants to load can be found in andrew's code 
%% set1 delete 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
% 
% NETseq_sample1 = importdata([dataloaddir 'weissman_set1del_NETseq.mat']);
% NETseq_sample2 = importdata([dataloaddir 'weissman_set1del_NETseq.mat']);
% 
% disp('Loaded in W set1del NETseq data')
% 
% xafter = 300;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant = 'set1del'
% datatype1 = 'NETseq'
%only one file for each no repeats but loading in two should mean the
%composite is the same as the original. 
%% rco1 delete 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
% 
% NETseq_sample1 = importdata([dataloaddir 'weissman_rco1del_NETseq.mat']);
% NETseq_sample2 = importdata([dataloaddir 'weissman_rco1del_NETseq.mat']);
% 
% disp('Loaded in W rco1del NETseq data')
% 
% xafter = 300;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant = 'rco1del'
% datatype1 = 'NETseq'
% %only one file for each no repeats but loading in two should mean the
% %composite is the same as the original.

%% Weissman WT 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
% 
% NETseq_sample1 = importdata([dataloaddir 'weissman_WT_NETseq.mat']);
% NETseq_sample2 = importdata([dataloaddir 'weissman_WT_NETseq.mat']);
% 
% disp('Loaded in W WT NETseq data')
% 
% xafter = 300;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant = 'weissWT'
% datatype1 = 'NETseq'
%only one file for each no repeats but loading in two should mean the
%composite is the same as the original. 
%% set2 delete 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
% 
% NETseq_sample1 = importdata([dataloaddir 'weissman_set2del_NETseq.mat']);
% NETseq_sample2 = importdata([dataloaddir 'weissman_set2del_NETseq.mat']);
% 
% disp('Loaded in W set2del NETseq data')
% xafter = 300;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant = 'set2del'
% datatype1 = 'NETseq'
%% eaf3 delete 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
% 
% NETseq_sample1 = importdata([dataloaddir 'weissman_eaf3del_NETseq.mat']);
% NETseq_sample2 = importdata([dataloaddir 'weissman_eaf3del_NETseq.mat']);
% 
% disp('Loaded in W eaf3del NETseq data')
% 
% xafter = 300;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant = 'eaf3del'
% datatype1 = 'NETseq'
%% spt4 delete

% dataloaddir = '/data/FIXED_DATA/FIXED_TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/notag_normalised/spt4del_nnt/';
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
% xafter = 300;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant = 'spt4del'
%datatype1 = 'NETseq'
%% dst1 delete

% dataloaddir = '/data/FIXED_DATA/FIXED_TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/notag_normalised/dst1del_nnt/';
% 
% NETseq_sample1 = importdata([dataloaddir 'ulku_dst1del_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'ulku_dst1del_NETseq_2.mat']);
% 
% disp('Loaded in UU dst1 delete NETseq data')
% 
% %NOTES:
% %The correlation between the two samples was very good
% %The average levels in each were reasonably close to the diagonal
% xafter = 300;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant = 'dst1del'
%datatype1 = 'NETseq'
%% xrn1 delete
% 
% dataloaddir = '/data/FIXED_DATA/FIXED_TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/notag_normalised/xrn1del_nnt/';
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
% xafter = 300;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant = 'xrn1del'
% datatype1 = 'NETseq'

%% SPT4 Anchor Away

% dataloaddir = '/data/FIXED_DATA/FIXED_TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/AnchorAway/AA_PROCESSED_NETSEQ/SPT4AA/';
% 
% NETseq_sample1 = importdata([dataloaddir 'ulku_SPT4AA_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'ulku_SPT4AA_NETseq_2.mat']);
% 
% disp('Loaded in UU SPT4 Anchor Away NETseq data')
% 
% %NOTES:
% %The correlation between the two samples is very good
% %The average levels ag genes were very close to the diagonal
% xafter = 300;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
%mutant = 'spt4AA'
%datatype1 = 'NETseq'

 %% SPT5 Anchor Away
% 
% dataloaddir = '/data/FIXED_DATA/FIXED_TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/AnchorAway/AA_PROCESSED_NETSEQ/SPT5AA/';
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
% xafter = 300;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
%mutant = 'spt5AA'
%datatype1 = 'NETseq'

%% SPT4/5 DMSO (Treating as repeats but strictly they are different)
% 
% dataloaddir = '/data/FIXED_DATA/FIXED_TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/AnchorAway/AA_PROCESSED_NETSEQ/SPT4AA_SPT5AA_DMSO_COMBO/';
% 
% NETseq_sample1 = importdata([dataloaddir 'ulku_SPT4AADMSO_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'ulku_SPT5AADMSO_NETseq_1.mat']);
% 
% disp('Loaded in UU SPT4/5 DMSO NETseq data')
% 
% %NOTES:
% %Correlation between these two samples (which should be similar but are
% %technically distinct) is very good
% %The average levels in genes also lies reasonably close to the diagonal
%  xafter = 300;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant = 'spt45DMSO'
%datatype1 = 'NETseq'

%% calculating total gene number

gene_number_forward = 0 ; %start with matrix of zeros
gene_number_reverse = 0 ;
for d = 1:NoChromosomes
    
    gene_number_forward = gene_number_forward + sum(length(gene_positions{d,1}(:,1)))  ;
    gene_number_reverse = gene_number_reverse + sum(length(gene_positions{d,2}(:,1))) ;
    %iteratively add up all the genes present
    
end
Total_gene_number = gene_number_forward + gene_number_reverse ;
%% forward strand
NET_Seq_matrix_forward = zeros(gene_number_forward,(ybefore+xafter)) ;
o = 0 ;
for cctr = 1:NoChromosomes

    %Isolating one chromosome%
    chromosomei_NETseq = NETseq_sample1(cctr,1) ;
    chromosomei_genepos = gene_positions(cctr,1) ;
    
    
    %pulling out ybeforebp before the end of gene%
    chromosomei_startpoints = chromosomei_genepos{1,1}(:,2) - ybefore;
    
    %start pos + length pulled out minus 1 gives ybefore xafter end%
    chromosomei_endpoints = chromosomei_startpoints + (length_pulled_out-1) ;
    
    %pulling out the individual genes%
    z = size(gene_positions{cctr,1}(:,1), 1) ;
   
    for gctr = 1:z 
        
        o = o + 1 ; 
        %add all genes together in single matrix for averaging later
        all_genes = all_genes + chromosomei_NETseq{1,1}(chromosomei_startpoints(gctr,1):chromosomei_endpoints(gctr,1),1) ;
        
        %construct a matrix with all gene profiles in
        NET_Seq_matrix_forward(o, :) = chromosomei_NETseq{1,1}(chromosomei_startpoints(gctr,1):chromosomei_endpoints(gctr,1),1) ;
    
    end
end
    
%% reverse strand
NET_Seq_matrix_reverse = zeros(gene_number_reverse,(ybefore+xafter)) ;
w = 0 ;
for cctr = 1:NoChromosomes
   
    %Isolating one chromosome%
    chromosomej_NETseq = NETseq_sample1(cctr,2) ;
    chromosomej_genepos = gene_positions(cctr,2) ;
    
    
    %pulling out the ybeforebp before end of genes%
    chromosomej_startpoints = chromosomej_genepos{1,1}(:,1) + ybefore ;
    %start points - x+y gives last x bases of reversed gene %
    chromosomej_endpoints = chromosomej_startpoints - (ybefore+xafter - 1) ;
    
    %pulling out the individual genes%
    x = size(gene_positions{cctr,2}(:,1), 1) ;
    
    for gctr = 1:x
        
        
        w = w + 1  ;
       
        %add to total signal matrix for all genes
        all_genes = all_genes + chromosomej_NETseq{1,1}(chromosomej_startpoints(gctr,1):-1:chromosomej_endpoints(gctr,1),1) ;
    
     
        %construct a matrix with all gene profiles in
        NET_Seq_matrix_reverse(w, :) = chromosomej_NETseq{1,1}(chromosomej_startpoints(gctr,1):-1:chromosomej_endpoints(gctr,1),1);
    end
end


%% Look at metagene forward and reverse 

metagene_forward = zeros(1, length_pulled_out);
for gctr= 1:gene_number_forward
    metagene_forward = metagene_forward + NET_Seq_matrix_forward(gctr, :) ;
end
metagene_forward = metagene_forward ./ gene_number_forward ;

metagene_reverse = zeros(1, length_pulled_out);
for gctr= 1:gene_number_reverse
    metagene_reverse = metagene_reverse + NET_Seq_matrix_reverse(gctr, :) ;
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

%% forward strand repeat for biological repeat then average to get overall NETseq matrix 
NET_Seq_matrix_forward2 = zeros(gene_number_forward,length_pulled_out ) ;
o = 0 ;
for cctr = 1:NoChromosomes

    %Isolating one chromosome%
    chromosomei2_NETseq = NETseq_sample2(cctr,1) ;
    chromosomei2_genepos = gene_positions(cctr,1) ;
    
    
    %pulling out 250bp before the end of gene%
    chromosomei2_startpoints = chromosomei2_genepos{1,1}(:,2) - ybefore ;
    
    %start pos + x+y gives last y bases + x after%
    chromosomei2_endpoints = chromosomei2_startpoints + (length_pulled_out-1) ;
    
    %pulling out the individual genes%
    z = size(gene_positions{cctr,1}(:,1), 1) ;
   
    for gctr = 1:z 
        
        o = o + 1 ; 
        %add all genes together in single matrix for averaging later
        all_genes = all_genes + chromosomei2_NETseq{1,1}(chromosomei2_startpoints(gctr,1):chromosomei2_endpoints(gctr,1),1) ;
        
        %construct a matrix with all gene profiles in
        NET_Seq_matrix_forward2(o, :) = chromosomei2_NETseq{1,1}(chromosomei2_startpoints(gctr,1):chromosomei2_endpoints(gctr,1),1) ;
    
    end
end
    
%% reverse strand
NET_Seq_matrix_reverse2 = zeros(gene_number_reverse,length_pulled_out ) ;
w = 0 ;
for cctr = 1:NoChromosomes
   
    %Isolating one chromosome%
    chromosomej2_NETseq = NETseq_sample2(cctr,2) ;
    chromosomej2_genepos = gene_positions(cctr,2) ;
    
    
    %pulling out the 250bp before end of genes%
    chromosomej2_startpoints = chromosomej2_genepos{1,1}(:,1) + ybefore  ;
    %start points - 499 gives last 250 bases of reversed gene + 250 after%
    chromosomej2_endpoints = chromosomej2_startpoints - (length_pulled_out - 1) ;
                                                                                                                                                                                                                                                                                                                            
    %pulling out the individual genes%
    x = size(gene_positions{cctr,2}(:,1), 1) ;
    
    for gctr = 1:x
        
        
        w = w + 1  ;                                            

        %add to total signal matrix for all genes
        all_genes = all_genes + chromosomej2_NETseq{1,1}(chromosomej2_startpoints(gctr,1):-1:chromosomej2_endpoints(gctr,1),1) ;
    
     
        %construct a matrix with all gene profiles in
        NET_Seq_matrix_reverse2(w, :) = chromosomej2_NETseq{1,1}(chromosomej2_startpoints(gctr,1):-1:chromosomej2_endpoints(gctr,1),1);
    end
end

%% Look at metagene forward and reverse 
metagene_forward2 = zeros(1, length_pulled_out);
for gctr= 1:gene_number_forward
    metagene_forward2 = metagene_forward2 + NET_Seq_matrix_forward(gctr, :) ;
end
metagene_forward2 = metagene_forward2 ./ gene_number_forward ;

metagene_reverse2 = zeros(1, length_pulled_out);
for gctr= 1:gene_number_reverse
    metagene_reverse2 = metagene_reverse2 + NET_Seq_matrix_reverse(gctr, :) ;
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


% %% Calculate mean across genome as 'background' 
% total_across_genome = 0;
% genome_size = 0 ;
% NETseq_composite = cell(size(NETseq_sample1));
% for cctr = 1:1:NoChromosomes %chromosme counter 
%     for sctr = 1:1:2 %strand counter 
%         
%         NETseq_composite{cctr,sctr} = ...
%             nanmean( ...
%             [NETseq_sample1{cctr,sctr},NETseq_sample2{cctr,sctr}],2 ...
%             );
%         
%     end
% end
% %calculate mean from composite cell array 
% for cctr = 1:NoChromosomes
%     for sctr = 1:2
%         total_across_genome = total_across_genome + sum(NETseq_composite{cctr, sctr}) ;
%         
%         genome_size = genome_size + sum(numel(NETseq_composite{cctr, sctr})) ;
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
% NETseq_nogenes = NETseq_composite;
% failed_removal_counter = 0 ;
% for cctr = 1:1:NoChromosomes
%     for sctr = 1:1:2
%     for gctr = size(gene_positions_for_nogene_mean{cctr,sctr}, 1):-1:1 %need to input before removal of genes by overlaps etc 
%         NETseq_nogenes{cctr, sctr}(gene_positions_for_nogene_mean{cctr, sctr}(gctr, 1):gene_positions_for_nogene_mean{cctr, sctr}(gctr, 2), :) = [];
%         %try NETseq_nogenes(gene_positions_for_nogene_mean{cctr, sctr}(gctr, 1):gene_positions_for_nogene_mean{cctr, sctr}(gctr, 2), :) = [];
%  %it worked without failing so didn't need try loop
%         %catch
%             %failed_removal_counter = failed_removal_counter +1 ; %if genes overlap will attempt to remove already removed area may fail
%         %end
%     end
%    
%     end
% end
% 
% total_across_genome_nogene = 0;
% genome_size_nogene = 0;
% %Calculate mean from nogenes NETseq array
% for cctr = 1:NoChromosomes
%     for sctr = 1:2
%         total_across_genome_nogene = total_across_genome_nogene + sum(NETseq_nogenes{cctr, sctr}) ;
%         
%         genome_size_nogene = genome_size_nogene + sum(numel(NETseq_nogenes{cctr, sctr})) ;
%         
%         
%     end
% end
% 
% %mean_across_genome_nogene = 0.2160 %answer from wild type NETseq data 
% mean_across_genome_nogene = total_across_genome_nogene ./ genome_size ;
% disp('done')
% 
% %% remove values below background mean 
% list_remaining_genes = (1:1:size(NETseq_matrix_no_negative, 1))'; %Index to apply same selection at 5' end later and on other datasets
% 
% for g = size(NETseq_matrix_no_negative, 1):-1:1
%     %if nanmean(NETseq_matrix_no_negative(g, :), 2)<mean_across_genome
%     %%found above = original filter
%      if nanmean(NETseq_matrix_no_negative(g, :), 2)<mean_across_genome_nogene %found above new filter lose less signal   
%         list_remaining_genes(g, :) =  [];
%         
%         NETseq_matrix_no_negative(g, :) = [] ;
%     end
% end
% save('WT_actively_expressing_genes_list.mat','list_remaining_genes') %must use same annotations file in both datasets if this is to work. 
 

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
for gctr= 1:size(Normalised_NETseq_Matrix, 1)
    metagene_Normalised = metagene_Normalised + Normalised_NETseq_Matrix(gctr, :) ;
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
for gctr = 1:size(NETseq_matrix_no_negative, 1)
    for x = 1:length_pulled_out/10
        Bin = mean(NETseq_matrix_no_negative(gctr,10*x-9:10*x));
        NETseq_Bins10bp_avg(gctr, x) = Bin  ;
    end
end


%% Normalise Binned NETseq
Normalised_ufdata_Bins_Matrix = Shape_normalisation_function(NETseq_Bins10bp_avg); %should probably integrate this into the code but i thought i was cool making a function
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
% meta_gene = all_genes ./ (2*Total_gene_number) ;
% figure
% plot(meta_gene)
filenameg = append('Unfiltered', datatype1, mutant, annotations,'300b300a', '.mat');
save(filenameg, 'Normalised_ufdata_Bins_Matrix')
unix(append('mv ', filenameg, ' Unfiltered_data'))

%% above is nofilter 3end pulling out script
toc


 %% Applying 3' end WT modelling to unfiltered 3end data 
        load(append(mutant1,'_actively_expressing_genes_list', annotations, '.mat'))
        ufdataloaddir = '/home/orie3770/Unfiltered_data/';
        %replace loadfilex with mutant to load 
        %loadfilex ='Unfiltered_PROseq_TIF_nogene500pc100b300a.mat';
        %loadfilex ='Unfiltered_NETseqspt4del_TIF_nogene500pc100b200a.mat';  
        %loadfilex = 'Unfiltered_NETseqdst1del_TIF_nogene500pc100b200a.mat';
        %loadfilex = 'Unfiltered_NETseqxrn1del_TIF_nogene500pc100b200a.mat';
        %loadfilex ='Unfiltered_NETseqspt4AA_TIF_nogene500pc100b200a.mat';
        %loadfilex = 'Unfiltered_PROseqspt4del_TIF_nogene500pc100b200a.mat';
        %loadfilex = 'Unfiltered_NETseqspt5AA_TIF_nogene500pc100b200a.mat';
        %loadfilex = 'Unfiltered_NETseqspt45DMSO_TIF_nogene500pc100b200a.mat';
        %loadfilex = 'Unfiltered_5_endset1delnogene500pcTIFseq.mat'; 
        %loadfilex = 'Unfiltered_5_endset2delnogene500pcTIFseq.mat';% need
        loadfilex = append('Unfiltered', datatype1, mutant, annotations,'300b300a', '.mat');
        %to redo clustering with this loaded in for set2 forgot to change 
        %loadfilex = 'Unfiltered_5_endeaf3delnogene500pcTIFseq.mat';
        %Can just load in different unfiltered files here, make sure
        %annotation stage is the same.
        load([ufdataloaddir loadfilex])
        %Apply same expression selection to uf data
        NETseq_3exprn_selected_ufdata = Normalised_ufdata_Bins_Matrix(list_remaining_genes, :);
        %Apply same sorting to PROseq 
        %NETseq_3class_sorted_ufdata = NETseq_3exprn_selected_ufdata(ix, :);

%% Attach gene names
%load(append(mutant1 ,'_actively_expressing_genes_list.mat', annotations)) %In current path or need to specify here

%load('WT_actively_expressing_genes_list_Harry_annotations_basic_filter.mat')
%take care to use correct gene annotations file.

%construct list of gene names in the same way as NET-seq matrix to apply
%filter to.
%%
gene_names_list = strings(Total_gene_number, 1) ;%empty presized list to add all the gene names into
%% forward names

w = 0;
for cctr = 1:1:NoChromosomes
    %x = size(gene_names{cctr,1} , 2); %labelled as names in TIF seq
    %annotations 
    x = size(gene_ids{cctr,1} , 2);
    for gctr = 1:1:x
        w = w+1;
        %gene_names_list(w,1) = gene_names{cctr,1}(1, gctr);
        %TIFseq labels in original annotations 
        gene_names_list(w,1) = gene_ids{cctr,1}(1, gctr);
    end
end

%% reverse names
% got error code in the middle and ended up duplicating 100 names. fixed
% now.
for cctr = 1:1:NoChromosomes
    %x = size(gene_names{cctr,2} , 2);
    x = size(gene_ids{cctr,2} , 2);
    for gctr = 1:1:x
        w = w+1;
        %gene_names_list(w,1) = gene_names{cctr,2}(1, gctr);
        gene_names_list(w,1) = gene_ids{cctr,2}(1, gctr);
        if w>Total_gene_number
            disp('too many gene names')
        end
    end
end
%%
filtered_gene_names = gene_names_list(list_remaining_genes);

 %% plot all expdropclasses overlayed with mutants
 
for n= 1:15

metagene=nanmean(Normalised_NETseq_Bins_Matrix(expdropclasses == n, :));
figure
plot(metagene)
%legend(append(mutant, string(size(Normalised_NETseq_Bins_Matrix(expdropclasses == n, :)))));
hold ON
mutmetagene=nanmean(NETseq_3exprn_selected_ufdata(expdropclasses == n, :));
plot(mutmetagene)

l1 = append(mutant1, datatype,' ', string(size(Normalised_NETseq_Bins_Matrix(expdropclasses == n), 1)))
l2 =append(mutant, datatype1,' ',  string(size(Normalised_NETseq_Bins_Matrix(expdropclasses == n), 1)))

legend(l1,l2) 

filename = append(mutant1, datatype, mutant, datatype1, 'expdropmodelled300b300a',string(n));

title(filename)
saveas(gcf, filename, 'jpeg')
saveas(gcf, filename, 'fig')

%save gene names alongside 
expdrop_gene_names = filtered_gene_names(expdropclasses == n);
directory = '/home/orie3770/';
fid = fopen(append(directory, filename , 'names','.txt' ) , 'w' );
fprintf(fid, '%s\n', expdrop_gene_names);
fclose(fid);
end



%% Compare mindist values between clusters 
mindistmeans = zeros(15,1) ;
for n = 1:15
mindistmeans(n,1) = mean(mindistvalues(expdropclasses == n));
end
figure
plot(mindistmeans)
title('mean minimum squared euclidean distance value')
saveas(gcf, append(filename, 'gof'), 'jpeg')
saveas(gcf, append(filename, 'gof'), 'fig')
%% Make metagene of just non-gene elements 
nonproteinmatrix = Normalised_NETseq_Bins_Matrix;
nongenelogical = zeros(size(filtered_gene_names));
for n = 1:size(filtered_gene_names, 1)
splitname = split(filtered_gene_names(n, 1), '');
if splitname(2) ~= 'Y' 
    nongenelogical(n,1) = 1 ; 
end

end

for n = size(filtered_gene_names, 1):-1:1 
if nongenelogical(n, 1) ==0
    nonproteinmatrix(n,:) = [] ;
    
end
 end   
    
    nongenemetagene = nanmean(nonproteinmatrix);

figure
plot(nongenemetagene)
title(append(mutant1, datatype, 'nonproteincodinggenes'))
legend(string(size(nonproteinmatrix, 1)))
saveas(gcf, append(mutant1, datatype, 'nonproteincodinggenes'), 'fig')

saveas(gcf, append(mutant1, datatype, 'nonproteincodinggenes'), 'jpeg')
unix(append('mv ',mutant1, datatype, 'nonproteincodinggenes','*', ' ./Current_Modelling_Job', string(Current_modelling_job)))
%%
unix(append('mv ',mutant1, datatype, mutant, datatype1, 'expdropmodelled300b300a' , '*', ' ./Current_Modelling_Job', string(Current_modelling_job)))
%
toc