% Have chagned every variable to end _start so this code can be run
% alongside the threeprime code without deleting each other's stuff

% Will have to alter the kmeans clustering code to accept _start variables
% done
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
 %% Harry Annotations
%  Nochromosomes =16;
%  annotationloaddir = ...
%     '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/OLD_ANNOTATIONS/Scer/';
% 
% extrasavepiece = 'ALTERNATIVE_ANNOTATIONS/HARRY/';
% 
% 
% loadfile1 = 'harry_annotations.mat';
% 
% load([annotationloaddir loadfile1])
% 
% 
% disp('Loaded annotation files.')
% 
% %This needs much extra processing - it's all in one big array for starters
% % and the orientation needs checking amongst other things
% 
% %Split out the original file into two pieces to match the others
% 
% gene_ids = cell(NoChromosomes,2);
% gene_positions = cell(NoChromosomes,2);
% 
% for cctr=1:1:NoChromosomes
%    
%     for sctr=1:1:2
%        
%         gene_ids{cctr,sctr} = harry_annotations{cctr,sctr+2};
%         gene_positions{cctr,sctr} = harry_annotations{cctr,sctr};
%         
%     end
%     
%     %flip the orientation of the gnees on the second strand
%     gene_positions{cctr,2} = fliplr(gene_positions{cctr,2});
%     
% end

%% Glugal 15
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/Glugal_shift_raw/15min/'
% 
% NETseq_sample1 = importdata([dataloaddir 'Glugal_15_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'Glugal_15_NETseq_2.mat']);
% 
% disp('loaded in glugal shift 15 min data')
% xafter_start = 1000;
% ybefore_start = 0;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes_start = zeros((xafter_start+ybefore_start),1) ;
% mutant = 'Glugal15'
% datatype1 = '5endNETseq'
%% Glugal 60
dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/Glugal_shift_raw/60min/'

NETseq_sample1 = importdata([dataloaddir 'Glugal_60_NETseq_1.mat']);
NETseq_sample2 = importdata([dataloaddir 'Glugal_60_NETseq_2.mat']);

disp('loaded in glugal shift 60 min data')
xafter_start = 1000;
ybefore_start = 0;
length_pulled_out_start = xafter_start+ybefore_start;
all_genes_start = zeros((xafter_start+ybefore_start),1) ;
mutant = 'Glugal60'
datatype1 = '5endNETseq'
%% Start with WT PROseq
% 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/LIS_PRO_SEQ/LisLab_PROseq_data/';
% 
% NETseq_sample1 = importdata([dataloaddir 'Lis_PROseq_WT.mat']);
% NETseq_sample2 = importdata([dataloaddir 'Lis_PROseq_WT.mat']); %duplicated the same dataset because shouldn't make a difference and can use same code that I have edited.
% 
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
% xafter_start = 1000;
% ybefore_start = 0;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes = zeros((xafter_start+ybefore_start),1) ;
% mutant1 = 'WT'
% datatype = '5endPROseq'
%mutants to load can be found in andrew's code 
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
% xafter_start = 1000;
% ybefore_start = 0;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes = zeros((xafter_start+ybefore_start),1) ;
% mutant1 = 'spt4del'
% datatype = '5endPROseq'

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
% 
% %mutants to load can be found in andrew's code 
% 
% 
%  
% xafter_start = 1000;
% ybefore_start = 0;
% length_pulled_out_start = (xafter_start+ybefore_start);
% all_genes_start = zeros(length_pulled_out_start,1) ;
%mutant = 'WT'
%%datatype1 = '5endNETseq'
%% rco1 delete 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
% 
% NETseq_sample1 = importdata([dataloaddir 'weissman_rco1del_NETseq.mat']);
% NETseq_sample2 = importdata([dataloaddir 'weissman_rco1del_NETseq.mat']);
% 
% disp('Loaded in W rco1del NETseq data')
% 
% xafter_start = 1000;
% ybefore_start = 0;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes_start = zeros((xafter_start+ybefore_start),1) ;
% mutant = 'rco1del'
% datatype1 = '5endNETseq'
% % %only one file for each no repeats but loading in two should mean the
% % %composite is the same as the original.

%% Weissman WT 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
% 
% NETseq_sample1 = importdata([dataloaddir 'weissman_WT_NETseq.mat']);
% NETseq_sample2 = importdata([dataloaddir 'weissman_WT_NETseq.mat']);
% 
% disp('Loaded in W WT NETseq data')
% 
% xafter_start = 1000;
% ybefore_start = 0;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes_start = zeros((xafter_start+ybefore_start),1) ;
% mutant = 'weissWT'
% datatype1 = '5endNETseq'
%only one file for each no repeats but loading in two should mean the
%composite is the same as the original. 
%% set1 delete 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
% 
% NETseq_sample1 = importdata([dataloaddir 'weissman_set1del_NETseq.mat']);
% NETseq_sample2 = importdata([dataloaddir 'weissman_set1del_NETseq.mat']);
% 
% disp('Loaded in W set1del NETseq data')
% 
% xafter_start = 1000;
% ybefore_start = 0;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes_start = zeros((xafter_start+ybefore_start),1) ;
% mutant ='set1del'
%datatype1 = '5endNETseq'
%only one file for each no repeats but loading in two should mean the
%composite is the same as the original. 
%% set2 delete 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
% 
% NETseq_sample1 = importdata([dataloaddir 'weissman_set2del_NETseq.mat']);
% NETseq_sample2 = importdata([dataloaddir 'weissman_set2del_NETseq.mat']);
% 
% disp('Loaded in W set2del NETseq data')
% xafter_start = 1000;
% ybefore_start = 0;
% length_pulled_out_start = xafter+ybefore;
% all_genes_start = zeros((xafter+ybefore),1) ;
% mutant = 'set2del'
% datatype1 = '5endNETseq'

%% eaf3 delete 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
% 
% NETseq_sample1 = importdata([dataloaddir 'weissman_eaf3del_NETseq.mat']);
% NETseq_sample2 = importdata([dataloaddir 'weissman_eaf3del_NETseq.mat']);
% 
% disp('Loaded in W eaf3del NETseq data')
% 
% xafter_start = 1000;
% ybefore_start = 0;
% length_pulled_out_start = xafter+ybefore;
% all_genes_start = zeros((xafter+ybefore),1) ;
% mutant = 'eaf3del'
% datatype1 = '5endNETseq'

%% spt4 delete

% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/notag_normalised/spt4del_nnt/';
% 
% NETseq_sample1 = importdata([dataloaddir 'ulku_spt4del_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'ulku_spt4del_NETseq_2.mat']);
% 
% disp('Loaded in UU spt4 delete NETseq data')

%NOTES:
%The correlation between the two samples was very good
%The average levels in the genes lie significantly off the diagonal - this
%might indicate a levelling issue and a more careful average might be
%needed in the longer run.
% xafter_start = 1000;
% ybefore_start = 250;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes_start = zeros((xafter_start+ybefore_start),1) ;
%mutant = 'spt4del'
%datatype1 = '5endNETseq'

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
% xafter_start = 1000;
% ybefore_start = 250;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes_start = zeros((xafter_start+ybefore_start),1) ;
% mutant = 'dst1del'
%% xrn1 delete

% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/notag_normalised/xrn1del_nnt/';
% 
% NETseq_sample1 = importdata([dataloaddir 'ulku_xrn1del_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'ulku_xrn1del_NETseq_2.mat']);
% 
% disp('Loaded in UU xrn1 delete NETseq data')
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
% mutant = 'xrn1del'
%datatype1 = '5endNETseq'
%% SPT4 Anchor Away
% 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/AnchorAway/AA_PROCESSED_NETSEQ/SPT4AA/';
% 
% NETseq_sample1 = importdata([dataloaddir 'ulku_SPT4AA_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'ulku_SPT4AA_NETseq_2.mat']);
% 
% disp('Loaded in UU SPT4 Anchor Away NETseq data')

%NOTES:
%The correlation between the two samples is very good
%The average levels ag genes were very close to the diagonal
% xafter_start = 1000;
% ybefore_start = 250;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes_start = zeros((xafter_start+ybefore_start),1) ;
% mutant = 'spt4AA'
%datatype1 = '5endNETseq'
%% SPT5 Anchor Away

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
% xafter_start = 1000; 
% ybefore_start = 250;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes_start = zeros((xafter_start+ybefore_start),1) ;
% mutant = 'spt5AA'
%datatype1 = '5endNETseq'
%% SPT4/5 DMSO (Treating as repeats but strictly they are different)

% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/AnchorAway/AA_PROCESSED_NETSEQ/SPT4AA_SPT5AA_DMSO_COMBO/';
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
% xafter_start = 1000; 
% ybefore_start = 250;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes_start = zeros((xafter_start+ybefore_start),1) ;
%mutant = 'spt45DMSO'
%datatype1 = '5endNETseq'
%% calculating total gene number

gene_number_forward = 0 ; %start with matrix of zeros
gene_number_reverse = 0 ;
for d = 1:NoChromosomes
    
    gene_number_forward = gene_number_forward + sum(length(gene_positions{d,1}(:,1)))  ;
    gene_number_reverse = gene_number_reverse + sum(length(gene_positions{d,2}(:,1))) ;
    %iteratively add up all the genes present
    
end
Total_gene_number_start = gene_number_forward + gene_number_reverse ;
%% forward
NET_Seq_matrix_forward_start = zeros(gene_number_forward,length_pulled_out_start) ;
o = 0 ;
for i = 1:NoChromosomes

    %Isolating one chromosome%
    chromosomei_NETseq = NETseq_sample1(i,1) ;
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
        all_genes_start = all_genes_start + chromosomei_NETseq{1,1}(chromosomei_startpoints(k,1):chromosomei_endpoints(k,1),1) ;
        
        %construct a matrix with all gene profiles in
        NET_Seq_matrix_forward_start(o, :) = chromosomei_NETseq{1,1}(chromosomei_startpoints(k,1):chromosomei_endpoints(k,1),1) ;
    
    end
end
    
%% reverse
NET_Seq_matrix_reverse_start = zeros(gene_number_reverse,length_pulled_out_start) ;
w = 0 ;
for j = 1:NoChromosomes
   
    %Isolating one chromosome%
    chromosomej_NETseq = NETseq_sample1(j,2) ;
    chromosomej_genepos = gene_positions(j,2) ;
    
    
    %pulling out the ybeforebp before end of genes%
    chromosomej_startpoints = chromosomej_genepos{1,1}(:,2)  ;
    %start points - x+y gives last x bases of reversed gene %
    chromosomej_endpoints = chromosomej_startpoints - (length_pulled_out_start - 1) ;
    
    %pulling out the individual genes%
    x = size(gene_positions{j,2}(:,1), 1) ;
    
    for k = 1:x
        
        
        w = w + 1  ;
       
        %add to total signal matrix for all genes
        all_genes_start = all_genes_start + chromosomej_NETseq{1,1}(chromosomej_startpoints(k,1):-1:chromosomej_endpoints(k,1),1) ;
    
     
        %construct a matrix with all gene profiles in
        NET_Seq_matrix_reverse_start(w, :) = chromosomej_NETseq{1,1}(chromosomej_startpoints(k,1):-1:chromosomej_endpoints(k,1),1);
    end
end


%% Look at metagene forward and reverse 

metagene_forward_start = zeros(1, length_pulled_out_start);
for k= 1:gene_number_forward
    metagene_forward_start = metagene_forward_start + NET_Seq_matrix_forward_start(k, :) ;
end
metagene_forward_start = metagene_forward_start ./ gene_number_forward ;

metagene_reverse_start = zeros(1, length_pulled_out_start);
for k= 1:gene_number_reverse
    metagene_reverse_start = metagene_reverse_start + NET_Seq_matrix_reverse_start(k, :) ;
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
NET_Seq_Matrix_start1 = vertcat(NET_Seq_matrix_forward_start,NET_Seq_matrix_reverse_start) ;

%% forward  repeat for biological repeat then average to get overall NETseq matrix 
NET_Seq_matrix_forward2 = zeros(gene_number_forward,length_pulled_out_start) ;
o = 0 ;
for i = 1:NoChromosomes

    %Isolating one chromosome%
    chromosomei2_NETseq = NETseq_sample2(i,1) ;
    chromosomei2_genepos = gene_positions(i,1) ;
    
    
    %pulling out 250bp before the end of gene%
    chromosomei2_startpoints = chromosomei2_genepos{1,1}(:,1) - ybefore_start ;
    
    %start pos + x+y gives last y bases + x after%
    chromosomei2_endpoints = chromosomei2_startpoints + (length_pulled_out_start-1) ;
    
    %pulling out the individual genes%
    z = size(gene_positions{i,1}(:,1), 1) ;
   
    for k = 1:z 
        
        o = o + 1 ; 
        %add all genes together in single matrix for averaging later
        all_genes_start = all_genes_start + chromosomei2_NETseq{1,1}(chromosomei2_startpoints(k,1):chromosomei2_endpoints(k,1),1) ;
        
        %construct a matrix with all gene profiles in
        NET_Seq_matrix_forward2(o, :) = chromosomei2_NETseq{1,1}(chromosomei2_startpoints(k,1):chromosomei2_endpoints(k,1),1) ;
    
    end
end
    
%% reverse
NET_Seq_matrix_reverse_start2 = zeros(gene_number_reverse,length_pulled_out_start) ;
w = 0 ;
for j = 1:NoChromosomes
   
    %Isolating one chromosome%
    chromosomej2_NETseq = NETseq_sample2(j,2) ;
    chromosomej2_genepos = gene_positions(j,2) ;
    
    
    %pulling out the 250bp before end of genes%
    chromosomej2_startpoints = chromosomej2_genepos{1,1}(:,2) + ybefore_start  ;
    %start points - 499 gives last 250 bases of reversed gene + 250 after%
    chromosomej2_endpoints = chromosomej2_startpoints - (length_pulled_out_start - 1) ;
                                                                                                                                                                                                                                                                                                                            
    %pulling out the individual genes%
    x = size(gene_positions{j,2}(:,1), 1) ;
    
    for k = 1:x
        
        
        w = w + 1  ;                                            

        %add to total signal matrix for all genes
        all_genes_start = all_genes_start + chromosomej2_NETseq{1,1}(chromosomej2_startpoints(k,1):-1:chromosomej2_endpoints(k,1),1) ;
    
     
        %construct a matrix with all gene profiles in
        NET_Seq_matrix_reverse_start2(w, :) = chromosomej2_NETseq{1,1}(chromosomej2_startpoints(k,1):-1:chromosomej2_endpoints(k,1),1);
    end
end

%% Look at metagene forward and reverse 
metagene_forward_start2 = zeros(1, length_pulled_out_start);
for k= 1:gene_number_forward
    metagene_forward_start2 = metagene_forward_start2 + NET_Seq_matrix_forward_start(k, :) ;
end
metagene_forward_start2 = metagene_forward_start2 ./ gene_number_forward ;

metagene_reverse_start2 = zeros(1, length_pulled_out_start);
for k= 1:gene_number_reverse
    metagene_reverse_start2 = metagene_reverse_start2 + NET_Seq_matrix_reverse_start(k, :) ;
end
metagene_reverse_start2 = metagene_reverse_start2 ./ gene_number_reverse ;

 %% plot metagenes 
figure 
plot(metagene_forward_start2)
title('metagene foward repeat')
figure
plot(metagene_reverse_start2)
title('metagene reverse repeat')

%% combine reverse and forward matrices
NET_Seq_Matrix_start2 = vertcat(NET_Seq_matrix_forward2,NET_Seq_matrix_reverse_start2) ;


%% Average the two samples 
NET_Seq_Matrix_start = (NET_Seq_Matrix_start1 + NET_Seq_Matrix_start2) ./ 2 ;

%% Remove negative values 
NETseq_matrix_no_negative_start = NET_Seq_Matrix_start ; 
NETseq_matrix_no_negative_start(NET_Seq_Matrix_start<0) = 0;


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
% %Combine samples
% 
% %% REMOVE VALUES BELOW MEAN
% 
% for g = size(NETseq_matrix_no_negative_start, 1):-1:1 
%     if nanmean(NETseq_matrix_no_negative_start(g, :), 2)<mean_across_genome %found above
%     
%        NETseq_matrix_no_negative_start(g, :) = [] ;
%     end      
% end
% 
%  

%% Normalise no negative NETseq Matrix
Normalised_NETseq_Matrix_start = Shape_normalisation_function(NETseq_matrix_no_negative_start); %only works if shape normalisation function is in current folder, need to find out how to specify path to function but for now will leave


% % remove NaNs
%  for k= size(Normalised_NETseq_Matrix_start, 1):-1:1 
%     
%     if  sum(isnan(Normalised_NETseq_Matrix_start(k, :)), 2)>0 %gives logical showing which genes contain any NaNs 
%         Normalised_NETseq_Matrix_start(k,:) = [] ;
% %     else 
% %         Normalised_NETseq_Matrix(k,:) = Normalised_NETseq_Matrix(k,:);
%     end
%     
%  end
 
 %% Find Normalised metagene 

 metagene_Normalised_start = zeros(1, length_pulled_out_start);
for k= 1:size(Normalised_NETseq_Matrix_start, 1)
    metagene_Normalised_start = metagene_Normalised_start + Normalised_NETseq_Matrix_start(k, :) ;
end
metagene_Normalised_start = metagene_Normalised_start ./ size(Normalised_NETseq_Matrix_start, 1) ;

%% plot normalised metagene 
%given no filtration often contains Nans here, not usually an issue as
%filtration occurs via copying 3end filtration later 
% figure 
% plot(metagene_Normalised_start)
% title('normalised metagene')


%% Bin the data prior to normalisation 
%% averaging every 10bps Bins
%make the matrix 
NETseq_Bins10bp_avg_start = zeros(size(NETseq_matrix_no_negative_start, 1), length_pulled_out_start/10);
for k = 1:size(NETseq_matrix_no_negative_start, 1)
    for x = 1:length_pulled_out_start/10
        Bin = mean(NETseq_matrix_no_negative_start(k,10*x-9:10*x));
        NETseq_Bins10bp_avg_start(k, x) = Bin  ;
    end
end


%% Normalise Binned NETseq
Normalised_ufdata_Bins_Matrix = Shape_normalisation_function(NETseq_Bins10bp_avg_start); %should probably integrate this into the code but i thought i was cool making a function
% %remove NaNs dont want to do this here as otherwise no lineup
% for k= size(Normalised_NETseq_Bins_Matrix_start, 1):-1:1 
%     
%     if  sum(isnan(Normalised_NETseq_Bins_Matrix_start(k, :)), 2)>0 %gives logical showing which genes contain any NaNs 
%         Normalised_NETseq_Bins_Matrix_start(k,:) = [] ;
% %     else not needed 
% %         Normalised_NETseq_Matrix(k,:) = Normalised_NETseq_Matrix(k,:);
%     end
%     
% end
%% NaN to zeros 
%Normalised_NETseq_Bins_Matrix(isnan(Normalised_NETseq_Bins_Matrix)) = 0; %converts
%Nans to zeros rather than removing
%% average gene calculation
% meta_gene_start = all_genes ./ (2*Total_gene_number) ;
% figure
% plot(meta_gene_start)
filenamey = append('Unfiltered', datatype1, mutant, annotations, '.mat');
save(filenamey, 'Normalised_ufdata_Bins_Matrix') 
unix(append('mv ', filenamey, ' Unfiltered_data'))