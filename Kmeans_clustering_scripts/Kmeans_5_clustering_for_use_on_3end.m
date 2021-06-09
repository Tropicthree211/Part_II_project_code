%% 
% load in annotations 
tic
Cluster_Job = 63;
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
annotations = 'HarrynoCuts'
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

%% calculating total gene number

gene_number_forward = 0 ; %start with matrix of zeros
gene_number_reverse = 0 ;
for d = 1:NoChromosomes
    
    gene_number_forward = gene_number_forward + sum(length(gene_positions{d,1}(:,1)))  ;
    gene_number_reverse = gene_number_reverse + sum(length(gene_positions{d,2}(:,1))) ;
    %iteratively add up all the genes present
    
end
Total_gene_number = gene_number_forward + gene_number_reverse ;
%% forward names
w = 0;
for cctr = 1:1:NoChromosomes
    x = size(gene_ids{cctr,1} , 2);
    for gctr = 1:1:x
        w = w+1;
        gene_names_list(w,1) = gene_ids{cctr,1}(1, gctr);
    end
end

%% reverse names
% got error code in the middle and ended up duplicating 100 names. fixed
% now.
for cctr = 1:1:NoChromosomes
    x = size(gene_ids{cctr,2} , 2);
    for gctr = 1:1:x
        w = w+1;
        gene_names_list(w,1) = gene_ids{cctr,2}(1, gctr);
        if w>Total_gene_number
            disp('too many gene names')
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
% xafter_start = 1000;
% ybefore_start = 0;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes_start = zeros((xafter_start+ybefore_start),1) ;
% mutant = 'Glugal15'
% datatype1 = '5endNETseq'
%% Glugal 60
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/Glugal_shift_raw/60min/'
% 
% NETseq_sample1 = importdata([dataloaddir 'Glugal_60_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'Glugal_60_NETseq_2.mat']);
% 
% disp('loaded in glugal shift 60 min data')
% xafter_start = 1000;
% ybefore_start = 0;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes_start = zeros((xafter_start+ybefore_start),1) ;
% mutant = 'Glugal60'
% datatype1 = '5endNETseq'
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

%mutants to load can be found in andrew's code 


 
xafter_start = 1000;
ybefore_start = 0;
length_pulled_out_start = (xafter_start+ybefore_start);
all_genes_start = zeros(length_pulled_out_start,1) ;
mutant = 'WT'
datatype1 = '5endNETseq'
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

% dataloaddir = '/data/FIXED_DATA/FIXED_TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/notag_normalised/spt4del_nnt/';
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
% xafter_start = 1000;
% ybefore_start = 250;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes_start = zeros((xafter_start+ybefore_start),1) ;
% mutant = 'dst1del'
%datatype1 = '5endNETseq'
%% xrn1 delete

% dataloaddir = '/data/FIXED_DATA/FIXED_TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/notag_normalised/xrn1del_nnt/';
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
% dataloaddir = '/data/FIXED_DATA/FIXED_TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/AnchorAway/AA_PROCESSED_NETSEQ/SPT4AA/';
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
% xafter_start = 1000; 
% ybefore_start = 250;
% length_pulled_out_start = xafter_start+ybefore_start;
% all_genes_start = zeros((xafter_start+ybefore_start),1) ;
% mutant = 'spt5AA'
%datatype1 = '5endNETseq'
%% SPT4/5 DMSO (Treating as repeats but strictly they are different)

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
    chromosomej_startpoints = chromosomej_genepos{1,1}(:,2) + ybefore_start ; %This ybefore start was missing, will have been why things looked out of line 
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
filenamei = append(datatype1, mutant, 'fwdmetagene');
saveas(gcf, filenamei, 'jpeg')
figure
plot(metagene_reverse_start')
title('metagene reverse')
filenamei = append(datatype1, mutant, 'rvsmetagene');
saveas(gcf, filenamei, 'jpeg')

%% combine reverse and forward matrices
NET_Seq_Matrix_start1 = vertcat(NET_Seq_matrix_forward_start,NET_Seq_matrix_reverse_start) ;

%% forward  repeat for biological repeat then average to get overall NETseq matrix 
NET_Seq_matrix_forward_start2 = zeros(gene_number_forward,length_pulled_out_start) ;
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
        NET_Seq_matrix_forward_start2(o, :) = chromosomei2_NETseq{1,1}(chromosomei2_startpoints(k,1):chromosomei2_endpoints(k,1),1) ;
    
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
    metagene_forward_start2 = metagene_forward_start2 + NET_Seq_matrix_forward_start2(k, :) ;
end
metagene_forward_start2 = metagene_forward_start2 ./ gene_number_forward ;

metagene_reverse_start2 = zeros(1, length_pulled_out_start);
for k= 1:gene_number_reverse
    metagene_reverse_start2 = metagene_reverse_start2 + NET_Seq_matrix_reverse_start2(k, :) ;
end
metagene_reverse_start2 = metagene_reverse_start2 ./ gene_number_reverse ;

 %% plot metagenes 
figure 
plot(metagene_forward_start2)
title('metagene foward repeat')
filenamei = append(datatype1, mutant, 'fwdrepmetagene');
saveas(gcf, filenamei, 'jpeg')

figure
plot(metagene_reverse_start2)
title('metagene reverse repeat')
filenamei = append(datatype1, mutant, 'rvsrepmetagene');
saveas(gcf, filenamei, 'jpeg')

%% combine reverse and forward matrices
NET_Seq_Matrix_start2 = vertcat(NET_Seq_matrix_forward_start2,NET_Seq_matrix_reverse_start2) ;


%% Average the two samples 
NET_Seq_Matrix_start = (NET_Seq_Matrix_start1 + NET_Seq_Matrix_start2) ./ 2 ;

%% Remove negative values 
NETseq_matrix_no_negative_start = NET_Seq_Matrix_start ; 
NETseq_matrix_no_negative_start(NET_Seq_Matrix_start<0) = 0;

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
list_remaining_genes = (1:1:size(NETseq_matrix_no_negative_start, 1))'; %Index to apply same selection at 5' end later and on other datasets

for g = size(NETseq_matrix_no_negative_start, 1):-1:1
    %if nanmean(NETseq_matrix_no_negative(g, :), 2)<mean_across_genome
    %%found above = original filter %04/02/21 reverted to original see new file Harry 
     if nanmean(NETseq_matrix_no_negative_start(g, :), 2)<mean_across_genome_nogene %found above new filter lose less signal   
        list_remaining_genes(g, :) =  [];
        
        NETseq_matrix_no_negative_start(g, :) = [] ;
    end
end
save(append(mutant,'_actively_expressing_genes_list', annotations, '.mat'),'list_remaining_genes') %must use same annotations file in both datasets if this is to work. 
%must remember to change this filename 
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
figure 
plot(metagene_Normalised_start)
title('5prime normalised metagene')

filenamei = append(datatype1, mutant, 'nmlsdmetagene');
saveas(gcf, filenamei, 'jpeg')

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
Normalised_NETseq_Bins_Matrix = Shape_normalisation_function(NETseq_Bins10bp_avg_start); %should probably integrate this into the code but i thought i was cool making a function

%16/04/21 have changed name no _start so this is clustered by the kmeans
%analysis.
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
%% Glugal 15
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/Glugal_shift_raw/15min/'
% 
% NETseq_sample1 = importdata([dataloaddir 'Glugal_15_NETseq_1.mat']);
% NETseq_sample2 = importdata([dataloaddir 'Glugal_15_NETseq_2.mat']);
% 
% disp('loaded in glugal shift 15 min data')
% xafter = 200;
% ybefore = 100;
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
% xafter = 200;
% ybefore = 100;
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
% xafter = 300;
% ybefore = 300;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant1 = 'WT'
% datatype = 'PROseq'
%mutants to load can be found in andrew's code 
%% spt4 delete
% 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/LIS_PRO_SEQ/LisLab_PROseq_data/'
% % 
%  %NETseq_sample = importdata([dataloaddir 'Lis_PROseq_spt4.mat']);%  
% NETseq_sample1 = importdata([dataloaddir 'LIS_PROseq_spt4del_amended.mat']);   %increased number of zeros to reflect size of chromosomes and to allow use of Harry's annotations.
% NETseq_sample2 = importdata([dataloaddir 'LIS_PROseq_spt4del_amended.mat']);
%  %NETseq_sample2 = importdata([dataloaddir 'ulku_spt4del_PROseq_2.mat']);
% 
% % isn't actually NETseq sample just calling it that so it plays nice with
% % my code, should eventually move to generic naming of intermediates but
% % functionally unimportant.
% % 
% disp('Loaded in LIS spt4 delete PROseq data')
% 
% % NOTES:
% % The correlation between the two samples was very good
% % The average levels in the genes lie significantly off the diagonal - this
% % might indicate a levelling issue and a more careful average might be
% % needed in the longer run.
% xafter = 300;
% ybefore = 300;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant1 = 'spt4del'
% datatype = 'PROseq'


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
% % 
%% set1 delete 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
% 
% NETseq_sample1 = importdata([dataloaddir 'weissman_set1del_NETseq.mat']);
% NETseq_sample2 = importdata([dataloaddir 'weissman_set1del_NETseq.mat']);
% 
% disp('Loaded in W set1del NETseq data')
% 
% xafter = 300;
% ybefore = 300;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant1 = 'set1del'
% datatype = 'NETseq'
% only one file for each no repeats but loading in two should mean the
% composite is the same as the original.
%% rco1 delete 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
% 
% NETseq_sample1 = importdata([dataloaddir 'weissman_rco1del_NETseq.mat']);
% NETseq_sample2 = importdata([dataloaddir 'weissman_rco1del_NETseq.mat']);
% 
% disp('Loaded in W rco1del NETseq data')
% 
% xafter = 200;
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
% %only one file for each no repeats but loading in two should mean the
%composite is the same as the original. 
%% set2 delete 
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
% 
% NETseq_sample1 = importdata([dataloaddir 'weissman_set2del_NETseq.mat']);
% NETseq_sample2 = importdata([dataloaddir 'weissman_set2del_NETseq.mat']);
% 
% disp('Loaded in W set2del NETseq data')
% xafter = 200;
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
% xafter = 200;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant1 = 'eaf3del'
% datatype = 'NETseq'
%% spt4 delete
% % 
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
% ybefore = 300;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant1 = 'spt4del'
% datatype = 'NETseq'
%% dst1 delete
% 
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
% ybefore = 300;
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
% xafter = 300;
% ybefore = 300;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant1 = 'spt4AA'
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
% xafter = 300;
% ybefore = 300;
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
%  xafter = 300;
% ybefore = 300;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant1 = 'spt4AA5AADMSO'
% datatype = 'NETseq'


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
filenamei = append(datatype, mutant1, 'forwardmetagene');
saveas(gcf, filenamei, 'jpeg')
figure
plot(metagene_reverse')
title('metagene reverse')
filenamei = append(datatype, mutant1, 'reversemetagene');
saveas(gcf, filenamei, 'jpeg')
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
    metagene_forward2 = metagene_forward2 + NET_Seq_matrix_forward2(gctr, :) ;
end
metagene_forward2 = metagene_forward2 ./ gene_number_forward ;

metagene_reverse2 = zeros(1, length_pulled_out);
for gctr= 1:gene_number_reverse
    metagene_reverse2 = metagene_reverse2 + NET_Seq_matrix_reverse2(gctr, :) ;
end
metagene_reverse2 = metagene_reverse2 ./ gene_number_reverse ;

 %% plot metagenes 
figure 
plot(metagene_forward2)
title('metagene forward repeat')
filenamei = append(datatype, mutant1, 'forwardrepeatmetagene');
saveas(gcf, filenamei, 'jpeg')
figure
plot(metagene_reverse2)
title('metagene reverse repeat')
filenamei = append(datatype, mutant1, 'reverserepeatmetagene');
saveas(gcf, filenamei, 'jpeg')
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
filenamei = append(datatype, mutant1, 'normalisedmetagene');
saveas(gcf, filenamei, 'jpeg')
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
meta_gene = all_genes ./ (2*Total_gene_number) ;
figure
plot(meta_gene)
filenamei = append(datatype, mutant1, 'metagene_allgenes');
saveas(gcf, filenamei, 'jpeg')

%% Bring information about gene size through
gene_sizes = cell(size(gene_positions));
for sctr = 1:1:2
    for cctr=1:1:NoChromosomes
        for gctr = 1:1:size(gene_positions{cctr,sctr},1)
            gene_sizes{cctr, sctr}(gctr,1) =  gene_positions{cctr,sctr}(gctr,2) - gene_positions{cctr,sctr}(gctr,1);
        end
    end
end
%convert cells with gene sizes into one matrix 
gene_sizes_matrix = zeros(Total_gene_number, 1) ;
w = 0;
for sctr = 1:1:2
    for cctr=1:1:NoChromosomes
        for gctr = 1:1:size(gene_positions{cctr,sctr},1)
            w = w+1;
            gene_sizes_matrix(w,1) =  gene_positions{cctr,sctr}(gctr,2) - gene_positions{cctr,sctr}(gctr,1);
        end
    end
end
%apply same filtering to gene size matrix 
gene_sizes_matrix_filtered = gene_sizes_matrix(list_remaining_genes,1) ; 
%% Bring information about gene expression through

gene_expression = cell(size(gene_positions));
for sctr = 1:1:2
    for cctr=1:1:NoChromosomes
        for gctr = 1:1:size(gene_positions{cctr,sctr},1)
            gene_expression{cctr, sctr}(gctr,1) =  sum(NETseq_composite{cctr, sctr}(gene_positions{cctr,sctr}(gctr,1):gene_positions{cctr,sctr}(gctr,2), 1))./gene_sizes{cctr, sctr}(gctr,1);
        end
    end
end
%convert cells with gene sizes into one matrix 
gene_expression_matrix = zeros(Total_gene_number, 1) ;
w = 0;
for sctr = 1:1:2
    for cctr=1:1:NoChromosomes
        for gctr = 1:1:size(gene_positions{cctr,sctr},1)
            w = w+1;
            gene_expression_matrix(w,1) =  sum(NETseq_composite{cctr, sctr}(gene_positions{cctr,sctr}(gctr,1):gene_positions{cctr,sctr}(gctr,2), 1))./gene_sizes{cctr, sctr}(gctr,1);
        end
    end
end
%apply same filtering to gene expression matrix 
gene_expression_matrix_filtered = gene_expression_matrix(list_remaining_genes,1) ;
%% Bring in info from Singleness measure 
load /home/orie3770/Modelling/Singlenessmeasure2individual.mat %make sure to use with NoCUTs data, has only 5579 values in. doesn't include Cuts
gene_singleness_matrix_filtered = singlenessbygenemeasure2(list_remaining_genes,1);
% might get issues with a bunch of NaNs in the data but we'll see how it
% comes out.

%% Bring in info from mTIF singleness measure 
load /home/orie3770/Modelling/mTIFderivedsinglenessmeasure.mat
gene_mTIFsingleness_matrix_filtered = singlenessbygenemeasuremTIFs(list_remaining_genes,1);
%%
 %save('/home/orie3770/Modelling/mTIFderivedsinglenessmeasure.mat', 'singlenessbygenemeasuremTIFs')

%%
toc


%% average gene calculation
% meta_gene_start = all_genes ./ (2*Total_gene_number) ;
% figure
% plot(meta_gene_start)
filenamey = append('Unfiltered', datatype, mutant1, annotations, '.mat');
save(filenamey, 'Normalised_ufdata_Bins_Matrix') 
unix(append('mv ', filenamey, ' Unfiltered_data'))

%%
toc

% Before running: change cluster job number 
%change mutant loaded in by pulling out genes code 
%change unfiltered data to be sorted by pulled out genes automated now

%Should try to automate to minimise the number of manual changes needed to save
%time 
klist=1:10;%the number of clusters to try
myfunc = @(X ,K)(kmeans(X, K));
clusters = append(datatype1, mutant, annotations)
%% evaluate clusters calinskiHarabasz
eva = evalclusters(Normalised_NETseq_Bins_Matrix ,myfunc,'CalinskiHarabasz','Klist',klist) 
figure 
plot(eva.CriterionValues)
title('Calinski Harabasz clustering evaluation')
xlabel('cluster number')
ylabel('result')
filenamen = append(clusters, string(ybefore), 'before', string(xafter), 'after', 'CHeval');
saveas(gcf, filenamen, 'jpeg')
%% evaluate clusters Silhouette
eva = evalclusters(Normalised_NETseq_Bins_Matrix ,myfunc,'silhouette','Klist',klist) 
figure
plot(eva.CriterionValues)
title('Silhouette clustering evaluation')
xlabel('cluster number')
ylabel('result')
filenamem = append(clusters, string(ybefore), 'before', string(xafter), 'after', 'Seval');
saveas(gcf, filenamem, 'jpeg')

%% elbow method implenented in code outside of MATLAB package 
dim=10;
% default number of test to get minimun under differnent random centroids
test_num=10;
distortion=zeros(dim,1);
for k_temp=1:dim(1)
    [~,~,sumd]=kmeans(Normalised_NETseq_Bins_Matrix,k_temp,'emptyaction','drop');
    distortion_temp=sum(sumd);
    % try differnet tests to find minimun disortion under k_temp clusters
    for test_count=2:test_num
        [~,~,sumd]=kmeans(Normalised_NETseq_Bins_Matrix,k_temp,'emptyaction','drop');
        distortion_temp=min(distortion_temp,sum(sumd));
    end
    distortion(k_temp,1)=distortion_temp;
end
variance=distortion(1:end-1)-distortion(2:end); %good way to pull out 
distortion_percent=cumsum(variance)/(distortion(1)-distortion(end));
figure
plot(distortion_percent,'b*--')
xlabel('Number of clusters')
ylabel('Percentage Distortion explained')
filenameo = append(clusters, string(ybefore), 'before', string(xafter), 'after', 'Elboweval');
saveas(gcf, filenameo, 'jpeg')

%% Eigenvalues
% 
% figure, sgtitle('set2')
% 
% subplot(1,2,1)
% 
% semilogy(diag(S),'k-o'), grid on
% 
% xlabel('r')
% 
% ylabel('\sigma_r')
% 
% set(gca)
% 
% subplot(1,2,2)
% 
% plot(cumsum(diag(S))./sum(diag(S)),'k-o') 
% 
% grid on
% 
% xlabel('r')
% 
% ylabel('Cumulative Energy')
% 
% set(gca)

%%
% variance_plot = zeros(1, 9);
% for K = 2:10 
%     clustering = kmeans(Normalised_NETseq_Bins_Matrix, K);
%     %Pull out each cluster 
%     clustern = Normalised_NETseq_Bins_Matrix(clustering == K, :) ;%might be wrong way round check
%     %calculate the variance of clustern, not sure if I need to do this for
%     %every single cluster and take an average or not.
%     
%     clusternvar = var(clustern, 1, [1 2]) ; % [1 2] here hopefully specifies the 1st and 2nd dimension of the cluster which will give the variance between 
%     %individual genes in the cluster rather than columns 
%     variance_plot(1,(K-1)) = clusternvar ;
% end   
% figure
% plot(variance_plot) %might need to transpose can't remember which dimension plot function likes. 
% title('Elbow method to determine optimal cluster number') 
% xlabel('Cluster number')
% ylabel('Variance')

%% perform the final kmeans cluster
filtered_gene_names = string(gene_names_list(list_remaining_genes,:));
JOB = append( 'Current_Cluster_Job', string(Cluster_Job));
unix(append('mkdir ', JOB))

clusters = append(datatype1, mutant, annotations)
% ufdata = 'NETseqspt4delxTIF';
%ufdata = 'NETseqdst1delxTIF';
%ufdata = 'NETseqxrn1delxTIF' ;
%ufdata = 'NETseqspt4AAxTIF';

%ufdata = 'NETseqspt5AA';
%ufdata = 'NETseqspt45DMSO';
%ufdata = 'NETseqset1del5endxTIF';
ufdata = append(datatype, mutant1, annotations) 
%ufdata = 'NETseqeaf3del5endxTIF';
%ufdata = 'PROseqWT'
%ufdata = 'PROseqspt4del';
 for K= 2:3
    for repeat = 1:3
        classes=kmeans(Normalised_NETseq_Bins_Matrix ,K);
        
        %% Pull out the classes from the clustering
        Class_sorted_NET_seq_matrix = zeros(size(Normalised_NETseq_Bins_Matrix, 1), size(Normalised_NETseq_Bins_Matrix, 2));
        Classmetagenes = zeros(size(Normalised_NETseq_Bins_Matrix, 2),K) ;
        classize_matrix = zeros(K, 1);
        classize_matrix0 = zeros(K+1, 1) ;
        classize_matrix_string = strings(K, 1);
        
        for n =1:K
            %sorted matrix for each class
            Classn = Normalised_NETseq_Bins_Matrix(classes==n, :) ;
            
            prevClasses = Normalised_NETseq_Bins_Matrix(classes<n, :);
            
            classize = size(Classn, 1);
            classize_matrix((n), 1) =  classize 
            classize_matrix0((n+1), 1) = classize; %this allows the 5' end clustering to work outside of the loop 
            
            prevclassizes = size(prevClasses,1);
            [sorted_classes, ix] = sort(classes);
            Class_sorted_NET_seq_matrix = Normalised_NETseq_Bins_Matrix(ix,:);
%             Class_sorted_NET_seq_matrix(prevclassizes+1:classize+prevclassizes,
%             :) = Classn ; old way of sorting use sort function allows
%             index to be pulled out that can sort 
            
            classize_matrix_string(n, 1) = string(classize_matrix(n, 1) ) ;
            
            
            
            %metagene for each class
            Classmetagenes(:, n) = sum(Class_sorted_NET_seq_matrix(prevclassizes+1:classize+prevclassizes, :))/classize ;
            
            %plot each metagene one at a time
            %     figure
            %     plot(Classmetagenes(:,n))
           
                classgenenames = filtered_gene_names(classes == (n));
                directory = '/home/orie3770/';
                fid = fopen(append(directory, mutant1 ,datatype, mutant, datatype1, string(n), 'cluster_of', string(K),'rep', string(repeat) , 'names','.txt' ) , 'w' );
                fprintf(fid, '%s\n', classgenenames);
                fclose(fid);
            
        end
        
        filename = append(clusters, string(ybefore_start), 'before', string(xafter_start), 'after', string(K), 'clusters', string(repeat));
        
        
       
        
        figure
        plot(Classmetagenes)
        title(filename) ;
        legend(classize_matrix_string)
        
        %% save as fig and jpeg
        saveas(gcf, filename, 'fig')
        
        saveas(gcf, filename, 'jpeg')
        disp('saved 3end clustering data') 
%Don't need to visualise WT
%         data again 
        
        %% Alternate sorting method Class_sorted_NET_seq_matrix = Class_sorted_NET_seq_matrix2 showing that ix contains original position info of genes.
%         [sorted_classes, ix] = sort(classes);
%         Class_sorted_NET_seq_matrix2 = Normalised_NETseq_Bins_Matrix(ix,:) ;
%         
%         Check_sorting_works = Normalised_NETseq_Bins_Matrix(ix, :) == Class_sorted_NET_seq_matrix ;
%         
%         clear original_order
%         original_order(ix, :) = Class_sorted_NET_seq_matrix ; % unsorts the above, so orginal_order == Normalised_NETseq_Bins_Matrix
        
        % %% Getting metagenes from sorted matrix
        % [p, index] = unique(sorted_classes)
        % Classmetagenes2 = zeros(size(Class_sorted_NET_seq_matrix2, 2),K);
        % for n=K:-1:1
        %
        %     Classmetagenes2(:, n) = sum(Class_sorted_NET_seq_matrix(size(sorted_classes,1):-1:index(n,1), :))/classize ;
        % end
        
        
%         
%         
        %% Applying 3' end WT clustering to 
        load(append(mutant,'_actively_expressing_genes_list', annotations, '.mat'))
        ufdataloaddir = '/home/orie3770/Unfiltered_data/';
        %replace loadfilex with mutant to load 
        %loadfilex ='Unfiltered_PROseq_TIF_nogene500pc100b200a.mat';
        %loadfilex ='Unfiltered_NETseqspt4del_TIF_nogene500pc100b200a.mat';  
        %loadfilex = 'Unfiltered_NETseqdst1del_TIF_nogene500pc100b200a.mat';
        %loadfilex = 'Unfiltered_NETseqxrn1del_TIF_nogene500pc100b200a.mat';
        %loadfilex ='Unfiltered_NETseqspt4AA_TIF_nogene500pc100b200a.mat';
        %loadfilex = 'Unfiltered_PROseqspt4del_TIF_nogene500pc100b200a.mat';
        %loadfilex = 'Unfiltered_NETseqspt5AA_TIF_nogene500pc100b200a.mat';
        %loadfilex = 'Unfiltered_NETseqspt45DMSO_TIF_nogene500pc100b200a.mat';
        %loadfilex = 'Unfiltered_5_endset1delnogene500pcTIFseq.mat'; 
        %loadfilex = 'Unfiltered_5_endset2delnogene500pcTIFseq.mat';% need
        loadfilex = append('Unfiltered', datatype, mutant1, annotations, '.mat'); %had to change to match switched around 
        %to redo clustering with this loaded in for set2 forgot to change 
        %loadfilex = 'Unfiltered_5_endeaf3delnogene500pcTIFseq.mat';
        %Can just load in different unfiltered files here, make sure
        %annotation stage is the same.
        load([ufdataloaddir loadfilex])
        %Apply same expression selection to PROseq data
        NETseq_3exprn_selected_ufdata = Normalised_ufdata_Bins_Matrix(list_remaining_genes, :);
        %Apply same sorting to PROseq 
        NETseq_3class_sorted_ufdata = NETseq_3exprn_selected_ufdata(ix, :);
        
        %% Find ufdata metagenes with NETseq clusters 
        
        Classmetagenesufdata = zeros(size(Normalised_ufdata_Bins_Matrix, 2), K);
        Cumulative_Classize = cumsum(classize_matrix0) ;
        for q = 2:K+1
            % Classmetagenes5prime_end(:, (q-1)) = sum(NETseq_matrix_3prime_Class_sorted_fiveprime_end( (classize_matrix((q-1),1)+1):(classize_matrix((q-1),1)+classize_matrix(q,1)), :)))/classize_matrix(q, 1) ;
            Classmetagenesufdata(:, (q-1)) = nansum(NETseq_3class_sorted_ufdata((Cumulative_Classize((q-1), 1) + 1): Cumulative_Classize(q , 1) , :))/classize_matrix0(q, 1) ;
            %had to use nansum due to NaNs, need to find source see if they're
            %impacting results
        end
        filename4 = append(clusters, ufdata, string(K), 'clusters', string(repeat)); % removed the following from filename4 string(ybefore), 'before', string(xafter), 'after',
        figure
        plot(Classmetagenesufdata)
        title(filename4);
        legend(classize_matrix_string)
        %% save as fig and jpeg
        saveas(gcf, filename4, 'fig')
        
        saveas(gcf, filename4, 'jpeg')
        
        disp('saved ufdata clusters')
        %% heatmap ufdata metagenes
        figure
        imagesc(NETseq_3class_sorted_ufdata)
        title(append(filename4, 'heatmap'))
        truesize
        filename5 = append(filename4, 'heatmap');
        saveas(gcf, filename5, 'jpeg')
        saveas(gcf, filename5, 'fig')
        disp('saved ufdata clustered heatmap')
        
                %% make boxplot of gene sizes within the clusters 
        filename8 = append(filename, 'sizes_boxplot');
        figure
        boxplot(gene_sizes_matrix_filtered, classes)
        title(filename8)
        ylabel('Gene Size')
        xlabel('cluster')
        %% Wilcoxon rank sum test to compare the medians of the two clusters in the boxplot
        size_ranksums = zeros(2, K)
        w = 0
        for q = 2:K
        [p, h] = ranksum(gene_sizes_matrix_filtered(classes==(q-1)), gene_sizes_matrix_filtered(classes == q));
        w=w+1
        size_ranksums(1, w) = p 
        size_ranksums(2, w) = h
        
        if q > 2 
            [q, i] = ranksum(gene_sizes_matrix_filtered(classes==(q-2)), gene_sizes_matrix_filtered(classes == q));
            w=w+1
            size_ranksums(1,w) = q 
            size_ranksums(2,w) = i
        end
        end
        directory = '/home/orie3770/';
                fid = fopen(append(directory, mutant, datatype1, string(K),'clusters' ,string(repeat),'repeat','size_ranksum.txt' ) , 'w' );
                fprintf(fid, '%f %f\n', size_ranksums );
                fclose(fid);
       

        %% save boxplot
        saveas(gcf, filename8, 'fig')
        saveas(gcf, filename8, 'jpeg')
        
        

        %% expression boxplot within clusters
        filenamet = append(filename, 'expression_boxplot');
        figure
        boxplot(gene_expression_matrix_filtered, classes)
        title(filenamet)
        ylabel('Expression level/gene size')
        xlabel('cluster')
        
         %% Wilcoxon rank sum test to compare the medians of the two clusters in the boxplot
        expression_ranksums = zeros(2, K)
        w = 0
        for q = 2:K
        [p, h] = ranksum(gene_expression_matrix_filtered(classes==(q-1)), gene_expression_matrix_filtered(classes == q));
        w=w+1
        expression_ranksums(1, w) = p 
        expression_ranksums(2, w) = h
        
        if q > 2 
            [q, i] = ranksum(gene_expression_matrix_filtered(classes==(q-2)), gene_expression_matrix_filtered(classes == q));
            w=w+1
            expression_ranksums(1,w) = q 
            expression_ranksums(2,w) = i
        end
        end
        directory = '/home/orie3770/';
                fid = fopen(append(directory, mutant, datatype1, string(K),'clusters',string(repeat),'repeat','expression_ranksum.txt' ) , 'w' );
                fprintf(fid, '%f %f\n', expression_ranksums );
                fclose(fid);
        %% save boxplot
        saveas(gcf, filenamet, 'fig')
        saveas(gcf, filenamet, 'jpeg')
        
        %% Boxplot of singleness within each of the clusters 
        filenamex = append(filename, 'singleness_boxplot');
        figure
        boxplot(gene_singleness_matrix_filtered, classes)
        title(filenamex)
        ylabel('singleness score')
        xlabel('cluster')
          %% Wilcoxon rank sum test to compare the medians of the two clusters in the boxplot
        singleness_ranksums = zeros(2, K)
        w = 0
        for q = 2:K
        [p, h] = ranksum(gene_singleness_matrix_filtered(classes==(q-1)), gene_singleness_matrix_filtered(classes == q));
        w=w+1
        singleness_ranksums(1, w) = p 
        singleness_ranksums(2, w) = h
        
        if q > 2 
            [q, i] = ranksum(gene_singleness_matrix_filtered(classes==(q-2)), gene_singleness_matrix_filtered(classes == q));
            w=w+1
            singleness_ranksums(1,w) = q 
            singleness_ranksums(2,w) = i
        end
        end
        directory = '/home/orie3770/';
                fid = fopen(append(directory, mutant, datatype1, string(K),'clusters',string(repeat),'repeat','singleness_ranksum.txt' ) , 'w' );
                fprintf(fid, '%f %f\n', singleness_ranksums );
                fclose(fid);
        %% save boxplot
        saveas(gcf, filenamex, 'fig')
        saveas(gcf, filenamex, 'jpeg')

         %% Boxplot of singleness within each of the clusters derived from S2 with similar ends and beginings collected.
        filenamez = append(filename, 'mTIFsingleness_boxplot');
        figure
        boxplot(gene_mTIFsingleness_matrix_filtered, classes)
        title(filenamez)
        ylabel('mTIFsingleness score')
        xlabel('cluster')
          %% Wilcoxon rank sum test to compare the medians of the two clusters in the boxplot
        mTIFsingleness_ranksums = zeros(2, K)
        w = 0
        for q = 2:K
        [p, h] = ranksum(gene_mTIFsingleness_matrix_filtered(classes==(q-1)), gene_mTIFsingleness_matrix_filtered(classes == q));
        w=w+1
        mTIFsingleness_ranksums(1, w) = p 
        mTIFsingleness_ranksums(2, w) = h
        
        if q > 2 
            [q, i] = ranksum(gene_mTIFsingleness_matrix_filtered(classes==(q-2)), gene_mTIFsingleness_matrix_filtered(classes == q));
            w=w+1
            mTIFsingleness_ranksums(1,w) = q 
            mTIFsingleness_ranksums(2,w) = i
        end
        end
        directory = '/home/orie3770/';
                fid = fopen(append(directory, mutant, datatype1, string(K),'clusters',string(repeat),'repeat','mTIFsingleness_ranksum.txt' ) , 'w' );
                fprintf(fid, '%f %f\n', mTIFsingleness_ranksums );
                fclose(fid);
        %% save boxplot
        saveas(gcf, filenamez, 'fig')
        saveas(gcf, filenamez, 'jpeg')

        
        %% Applying 3' end cluster sorting to 5' end data
%         %load in 5' end NETseqMatrix
%         load('Unfiltered_5_end_genes.mat') %need to either have in the current path or save it somewhere specific and specify the path here.Saves me having to run the 5' end every time.
%         %make sure to check the annotations files used match up. 
%         %Apply same selection to 5' end
%         NETseq_3primeexprnSelected_Normalised_Bins_Matrix_start = Normalised_NETseq_Bins_Matrix_start(list_remaining_genes, :);
%         %Apply same sorting to 5' end as above using ix
%         NETseq_matrix_3prime_Class_sorted_fiveprime_end = NETseq_3primeexprnSelected_Normalised_Bins_Matrix_start(ix,:);
        
%         %% Find metagenes at 5' end
%         Classmetagenes5prime_end = zeros(size(Normalised_NETseq_Bins_Matrix_start, 2), K);
%         Cumulative_Classize = cumsum(classize_matrix0) ;
%         for q = 2:K+1
%             % Classmetagenes5prime_end(:, (q-1)) = sum(NETseq_matrix_3prime_Class_sorted_fiveprime_end( (classize_matrix((q-1),1)+1):(classize_matrix((q-1),1)+classize_matrix(q,1)), :)))/classize_matrix(q, 1) ;
%             Classmetagenes5prime_end(:, (q-1)) = nansum(NETseq_matrix_3prime_Class_sorted_fiveprime_end((Cumulative_Classize((q-1), 1) + 1): Cumulative_Classize(q , 1) , :))/classize_matrix0(q, 1) ;
%             %had to use nansum due to NaNs, need to find source see if they're
%             %impacting results
%         end
%         filename2 = append(mutant, '5endsortedon3endclustering', string(ybefore), 'before', string(xafter), 'after', string(K), 'clusters', string(repeat));
%         figure
%         plot(Classmetagenes5prime_end)
%         title(filename2);
%         legend(classize_matrix_string)
%         %% save as fig and jpeg
%         saveas(gcf, filename2, 'fig')
%         
%         saveas(gcf, filename2, 'jpeg')
%         
%         disp('saved 5end clusters')
%% Attach gene names
%load(append(mutant,'_actively_expressing_genes_list', annotations, '.mat')) %In current path or need to specify here
%moved earlier to avoid issues.
%take care to use correct gene annotations file.

%construct list of gene names in the same way as NET-seq matrix to apply
%filter to.
%%

%save gene names alongside 



%% heatmap class metagenes
    figure
    imagesc(Class_sorted_NET_seq_matrix)
    title(append(filename, 'heatmap'))
    truesize
    filename3 = append(filename, 'heatmap');
    saveas(gcf, filename3, 'jpeg')
    saveas(gcf, filename3, 'fig')
    disp('saved NETseq clustered heatmap')

        end
end

% cluster_sorted_gene_names = filtered_gene_names(ix, :);
% 
% for clctr = 2:1:(K+1)
%     w=0;
%     for pctr = (Cumulative_Classize((clctr-1),1)+1):Cumulative_Classize(clctr,1)
%         w=w+1;
%     gene_names_clustered(w,(clctr-1)) = ...
%     cluster_sorted_gene_names(pctr,1);
%     end
% end
%% save gene names
% 
% save(append(datatype, mutant1, annotations, 'clustersortedgenenames',ufdata, string(K), 'clusters', string(repeat)), 'cluster_sorted_gene_names')
% save(append(datatype,mutant1, annotations, 'clusterednamesforGO',ufdata, string(K), 'clusters', string(repeat)), 'gene_names_clustered')
% disp('saved clustered gene names')
    

%% Gene size info into box plot of clusters
% cluster_sorted_gene_sizes = gene_sizes_filtered(ix,:);
% %prespecify size of gene_sizes_clustered
% gene_sizes_clustered = zeros(Total_gene_number,K);
% 
% for clctr = 2:1:(K+1)
%     w=0;
%     for pctr = (Cumulative_Classize((clctr-1),1)+1):Cumulative_Classize(clctr,1)
%         w=w+1;
%     gene_sizes_clustered(w,(clctr-1)) = ...
%     cluster_sorted_gene_sizes(pctr,1);
%     end
% end
% %Trim down gene size matrix to remove zeros.
% for clctr = 2:1:K
%     for gctr = 1:1:Total_gene_number
%         if gene_sizes_clustered(gctr, clctr) == 0 
%             gene_sizes_clustered(gctr,clctr) = [];
%         end
%     end
% end
%             
% % make boxplot of gene sizes within the clusters 
% filename4 = append(filename, 'boxplot');
% figure
% boxplot(gene_sizes_clustered)
% saveas(gcf, filename4, 'fig')
% saveas(gcf, filename4, 'jpeg')
% cluster_sorted_gene_sizes = gene_sizes_filtered(ix,:);
% %prespecify size of gene_sizes_clustered
% gene_sizes_clustered = zeros(Total_gene_number,K);
% 
% for clctr = 2:1:(K+1)
%     w=0;
%     for pctr = (Cumulative_Classize((clctr-1),1)+1):Cumulative_Classize(clctr,1)
%         w=w+1;
%     gene_sizes_clustered(w,(clctr-1)) = ...
%     cluster_sorted_gene_sizes(pctr,1);
%     end
% end
% %Trim down gene size matrix to remove zeros.
% for clctr = 2:1:K
%     for gctr = 1:1:Total_gene_number
%         if gene_sizes_clustered(gctr, clctr) == 0 
%             gene_sizes_clustered(gctr,clctr) = [];
%         end
%     end
% end
%             

    
 %% end
% unix(append('mv nogene500* ', JOB))
% unix(append('mv NETseq', mutant1, '* ', JOB))
% unix(append('mv PROseq', mutant1, '* ', JOB))
% unix(append('mv spt4delPROseq', mutant1, '* ', JOB))
unix(append('mv ', datatype1, mutant, '* ', JOB))
unix(append('mv ', mutant ,datatype1, '* ', JOB))
%unix(append('mv PROseq', mutant1, '* ', JOB))

%% Look for patterns in heatmap
% figure
% imagesc(NETseq_matrix_3prime_Class_sorted_fiveprime_end)
% truesize



