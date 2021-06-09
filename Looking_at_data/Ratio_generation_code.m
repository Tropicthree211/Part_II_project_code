%% Calculating ratios between different windows in the WT 
% Have pasted in the canonical processing pathway, will comment out
% normalisation part and calculate ratios. 
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
annotations = 'HarrynoCuts';
%% calculating total gene number

gene_number_forward = 0 ; %start with matrix of zeros
gene_number_reverse = 0 ;
for d = 1:NoChromosomes
    
    gene_number_forward = gene_number_forward + sum(length(gene_positions{d,1}(:,1)))  ;
    gene_number_reverse = gene_number_reverse + sum(length(gene_positions{d,2}(:,1))) ;
    %iteratively add up all the genes present
    
end
Total_gene_number = gene_number_forward + gene_number_reverse ;

%% remove Cuts and Xuts

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
%% Glugal 15/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA
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
% %mutants to load can be found in andrew's code 
% 
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
% %composite is the same as the original. 
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
% ybefore = 300;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant1 = 'spt4del'
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
% xafter = 200;
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
%list remaining genes is trimmed below in line with NETseq matrix hence if
%same filtration is applied before then same genes will be trimmed here

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
toc





%% Below is 5' pulling out data 
%going to use to pull out between 500 and 750 as a region to normalise to
%within the gene body. 


%% Harry annotations
%  NoChromosomes=16
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
% annotations = 'HarrynoCuts'
% %% remove Cuts and Xuts
% % 
% for sctr = 1:1:2
%     for cctr = 1:1:NoChromosomes
% %         nongenelogical = zeros(size(gene_ids{cctr,sctr}(),2), 1);
%         for gctr = size(gene_ids{cctr,sctr}, 2):-1:1
%             splitname = split(gene_ids{cctr,sctr}{1,gctr}, '');
%             if splitname{2} ~= 'Y'
%                 gene_ids{cctr,sctr}(:, gctr) = [] ;
%                 gene_positions{cctr,sctr}(gctr,:) = [];
%                 
%             end
%             
%         end
%         
%         
%     end
% end

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
% w = 0;
% for cctr = 1:1:NoChromosomes
%     x = size(gene_ids{cctr,1} , 2);
%     for gctr = 1:1:x
%         w = w+1;
%         gene_names_list(w,1) = gene_ids{cctr,1}(1, gctr);
%     end
% end
% 
% %% reverse names
% % got error code in the middle and ended up duplicating 100 names. fixed
% % now.
% for cctr = 1:1:NoChromosomes
%     x = size(gene_ids{cctr,2} , 2);
%     for gctr = 1:1:x
%         w = w+1;
%         gene_names_list(w,1) = gene_ids{cctr,2}(1, gctr);
%         if w>Total_gene_number
%             disp('too many gene names')
%         end
%     end
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


 
xafter_start = 750;
ybefore_start = -500;
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
% gene_positions_for_nogene_mean = gene_positions;
% for cctr = 1:1:NoChromosomes
%     for sctr = 1:1:2
%     for gctr = size(gene_positions_for_nogene_mean{cctr,sctr}, 1):-1:1 %need to input before removal of genes by overlaps etc 
%         if gctr == 1
%             NETseq_nogenes{cctr, sctr}(gene_positions_for_nogene_mean{cctr, sctr}(gctr, 2):-1:gene_positions_for_nogene_mean{cctr, sctr}(gctr, 1), :) = [];
%             continue;
%             %should avoid trying to index with gctr = 0 which will fail my
%             %code
%         end
%         if gene_positions_for_nogene_mean{cctr, sctr}(gctr, 1) < gene_positions_for_nogene_mean{cctr, sctr}((gctr -1) , 2)
%             
%         gene_positions_for_nogene_mean{cctr, sctr}((gctr -1), 2) = gene_positions_for_nogene_mean{cctr, sctr}(gctr , 1) -1;
%         
%         %added if statement should ensure not removing stuff that has already been removed in case of overlapping genes. 
%         end
%         NETseq_nogenes{cctr, sctr}(gene_positions_for_nogene_mean{cctr, sctr}(gctr, 2):-1:gene_positions_for_nogene_mean{cctr, sctr}(gctr, 1), :) = []; %
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
%             %GM commented out the above 04/02/21 to quickly use Harry's
%             %annotations rather than constructing a nogene version removing
%             %all the genes. %pretty sure i can just define gene positions
%             %for nogene mean as gene positions because no longer removing
%             %any genes before this point
% %mean_across_genome_nogene = 0.2160 %answer from wild type NETseq data 
%  mean_across_genome_nogene = total_across_genome_nogene ./ genome_size ;
% % disp('done')
% 
% %% remove values below background mean 
% list_remaining_genes = (1:1:size(NETseq_matrix_no_negative_start, 1))'; %Index to apply same selection at 5' end later and on other datasets
% 
% for g = size(NETseq_matrix_no_negative_start, 1):-1:1
%     %if nanmean(NETseq_matrix_no_negative(g, :), 2)<mean_across_genome
%     %%found above = original filter %04/02/21 reverted to original see new file Harry 
%      if nanmean(NETseq_matrix_no_negative_start(g, :), 2)<mean_across_genome_nogene %found above new filter lose less signal   
%         list_remaining_genes(g, :) =  [];
%         
%         NETseq_matrix_no_negative_start(g, :) = [] ;
%     end
% end
% save(append(mutant,'_actively_expressing_genes_list', annotations, '.mat'),'list_remaining_genes') %must use same annotations file in both datasets if this is to work. 
% %must remember to change this filename 
% %g=size(NET_Seq_Matrix(Logic_matrix_remove_negative));
% % for h=1:g
% %     Total_NET
% %     
% % end
% % NETseq_matrix_no_negative = NET_Seq_Matrix(Logic_matrix_remove_negative)==0 
% % %transpose to allow matrix multiplication
% % no_negativity = Logic_matrix_remove_negative';
% %generate 
% % NETseq_matrix_no_negative = Total_NET_Seq_Matrix*no_negativity ;

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

% %% Normalise no negative NETseq Matrix
% Normalised_NETseq_Matrix_start = Shape_normalisation_function(NETseq_matrix_no_negative_start); %only works if shape normalisation function is in current folder, need to find out how to specify path to function but for now will leave
% 
% 
% % % remove NaNs
% %  for k= size(Normalised_NETseq_Matrix_start, 1):-1:1 
% %     
% %     if  sum(isnan(Normalised_NETseq_Matrix_start(k, :)), 2)>0 %gives logical showing which genes contain any NaNs 
% %         Normalised_NETseq_Matrix_start(k,:) = [] ;
% % %     else 
% % %         Normalised_NETseq_Matrix(k,:) = Normalised_NETseq_Matrix(k,:);
% %     end
% %     
% %  end
%  
%  %% Find Normalised metagene 
% 
%  metagene_Normalised_start = zeros(1, length_pulled_out_start);
% for k= 1:size(Normalised_NETseq_Matrix_start, 1)
%     metagene_Normalised_start = metagene_Normalised_start + Normalised_NETseq_Matrix_start(k, :) ;
% end
% metagene_Normalised_start = metagene_Normalised_start ./ size(Normalised_NETseq_Matrix_start, 1) ;
% 
% 
% %% plot normalised metagene 
% %given no filtration often contains Nans here, not usually an issue as
% %filtration occurs via copying 3end filtration later 
% figure 
% plot(metagene_Normalised_start)
% title('5prime normalised metagene')
% 
% filenamei = append(datatype1, mutant, 'nmlsdmetagene');
% saveas(gcf, filenamei, 'jpeg')
% 
% %% Bin the data prior to normalisation 
% %% averaging every 10bps Bins
% %make the matrix 
% NETseq_Bins10bp_avg_start = zeros(size(NETseq_matrix_no_negative_start, 1), length_pulled_out_start/10);
% for k = 1:size(NETseq_matrix_no_negative_start, 1)
%     for x = 1:length_pulled_out_start/10
%         Bin = mean(NETseq_matrix_no_negative_start(k,10*x-9:10*x));
%         NETseq_Bins10bp_avg_start(k, x) = Bin  ;
%     end
% end
% 
% 
% %% Normalise Binned NETseq
% Normalised_NETseq_Bins_Matrix_start = Shape_normalisation_function(NETseq_Bins10bp_avg_start); %should probably integrate this into the code but i thought i was cool making a function
% 
% %16/04/21 have changed name no _start so this is clustered by the kmeans
% %analysis.
% % %remove NaNs dont want to do this here as otherwise no lineup
% % for k= size(Normalised_NETseq_Bins_Matrix_start, 1):-1:1 
% %     
% %     if  sum(isnan(Normalised_NETseq_Bins_Matrix_start(k, :)), 2)>0 %gives logical showing which genes contain any NaNs 
% %         Normalised_NETseq_Bins_Matrix_start(k,:) = [] ;
% % %     else not needed 
% % %         Normalised_NETseq_Matrix(k,:) = Normalised_NETseq_Matrix(k,:);
% %     end
% %     
% % end
% %% NaN to zeros 
% %Normalised_NETseq_Bins_Matrix(isnan(Normalised_NETseq_Bins_Matrix)) = 0; %converts
% %Nans to zeros rather than removing


%% Now pull out the average signal for each gene across the 3 windows at the 3' end and generate ratio to 500-750 signal


%Pull out the specific windows from the above code 
approach_window_matrix = zeros(size(NETseq_matrix_no_negative, 1), 100);
termination_window_matrix = zeros(size(NETseq_matrix_no_negative, 1), 100);
aftermath_window_matrix = zeros(size(NETseq_matrix_no_negative, 1), 100);
control_elongation_matrix = zeros(size(NETseq_matrix_no_negative, 1), 100);


for gctr = 1:1:size(NETseq_matrix_no_negative) 
    %need to ask jane about best exact window locations these are
    %placeholders to test the code
    approach_window_matrix(gctr, :) = NETseq_matrix_no_negative(gctr, 201:300);
    
    termination_window_matrix(gctr, :) =NETseq_matrix_no_negative(gctr, 301:400);
    
    aftermath_window_matrix(gctr, :) = NETseq_matrix_no_negative(gctr, 501:600);
    
    control_elongation_matrix(gctr,:) = NETseq_matrix_no_negative(gctr,1:100);
   
end


%Take means of the signal in each window 
approach_mean_matrix = zeros(size(NETseq_matrix_no_negative, 1),1);
termination_mean_matrix = zeros(size(NETseq_matrix_no_negative, 1),1);
aftermath_mean_matrix = zeros(size(NETseq_matrix_no_negative, 1),1);
control_mean_matrix = zeros(size(NETseq_matrix_no_negative, 1),1);



for gctr = 1:1:size(NETseq_matrix_no_negative) 

%Pretty confident allzero genes will have been filtered out by now so
%should be able to use normal mean. may have to change later. 
approach_mean_matrix(gctr,1) = mean(approach_window_matrix(gctr,:));

termination_mean_matrix(gctr,1) = mean(termination_window_matrix(gctr,:));

aftermath_mean_matrix(gctr,1) = mean(aftermath_window_matrix(gctr,:));

control_mean_matrix(gctr,1) = mean(control_elongation_matrix(gctr,:));

end 

%Divide through window mean by control mean to get ratio for eac window in
%each gene in 3 large matrices 
approach_ratio_matrix = zeros(size(NETseq_matrix_no_negative, 1),1);
termination_ratio_matrix = zeros(size(NETseq_matrix_no_negative, 1),1);
aftermath_ratio_matrix = zeros(size(NETseq_matrix_no_negative, 1),1);
control_ratio_matrix = zeros(size(NETseq_matrix_no_negative, 1),1);

for gctr = 1:1:size(NETseq_matrix_no_negative) 

    approach_ratio_matrix(gctr,1) = approach_mean_matrix(gctr,1) ./ control_mean_matrix(gctr,1);
    termination_ratio_matrix(gctr,1) = termination_mean_matrix(gctr,1)./ control_mean_matrix(gctr,1);
    aftermath_ratio_matrix(gctr,1) = aftermath_mean_matrix(gctr,1) ./control_mean_matrix(gctr,1);

    % can do unique on this to check tis all ones showing nothing has gone
    % awry.
    control_ratio_matrix(gctr,1) = control_mean_matrix(gctr,1) ./control_mean_matrix(gctr,1);
    
end

%% Getting termination ratio values to get more quantitative 







% %% Analysing the ratios
% figure
% boxplot(approach_ratio_matrix)
% 
% figure 
% boxplot(termination_ratio_matrix)
% 
% figure
% boxplot(aftermath_ratio_matrix)

% mean(


