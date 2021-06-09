%% Import Ulku's NETseq
%Including the newer AnchorAway data

%Adapted from ImportCramerNETseq_K562_AA_speedup.m

%Looking at the notag normalised version of the data 
%Not sure if the anchor away has been similarly treated - need to ask ulku

%90 min anchor away was done as a single repeat so will not be analysed
%in the initial pass but could consider looking at it later


%GM 18/02/21 adapted for glugal wig files to look at stressed termination
%observing whether similar to PRO-seq or if the second cluster is still
%present 

%% Start by testing the Glugal 15 NETseq

dataloaddir_start = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/Glugal_shift_raw/';
dataloaddir_end = '15min/';

dataloaddir_15 = [dataloaddir_start dataloaddir_end];

NoChromosomes = 16; %Read from the data (roman numerals used)
                  
Chromosome_Names = {'I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI'};

                    

%%

fileID1 = fopen([dataloaddir_15 'by_rpb3FLAG_gal15_repeat1_ntsub_unsmo_negstrand.wig']);
Glugal_15_NETseq_minus_1_readin = textscan(fileID1,'%s %f %f %f', ...
    'CommentStyle','#');
fclose(fileID1);
whos Glugal_15_NETseq_minus_1_readin



fileID2 = fopen([dataloaddir_15 'by_rpb3FLAG_gal15_repeat1_ntsub_unsmo_posstrand.wig']);
Glugal_15_NETseq_plus_1_readin = textscan(fileID2,'%s %f %f %f', ...
    'CommentStyle','#');
fclose(fileID2);
whos Glugal_15_NETseq_plus_1_readin

disp('Finished importing raw Ulku NETseq WT sample 1')



%%


fileID3 = fopen([dataloaddir_15 'by_rpb3FLAG_gal15_repeat2_ntsub_unsmo_negstrand.wig']);
Glugal_15_NETseq_minus_2_readin = textscan(fileID3,'%s %f %f %f', ...
    'CommentStyle','#');
fclose(fileID3);
whos Glugal_15_NETseq_minus_2_readin



fileID4 = fopen([dataloaddir_15 'by_rpb3FLAG_gal15_repeat2_ntsub_unsmo_posstrand.wig']);
Glugal_15_NETseq_plus_2_readin = textscan(fileID4,'%s %f %f %f', ...
    'CommentStyle','#');
fclose(fileID4);
whos Glugal_15_NETseq_plus_2_readin

disp('Finished importing raw Ulku NETseq WT sample 2')



%% Load in the Glugal60

dataloaddir_end = '60min/';

dataloaddir_60 = [dataloaddir_start dataloaddir_end];



fileID1 = fopen([dataloaddir_60 'by_rpb3FLAG_gal60_repeat1_ntsub_unsmo_negstrand.wig']);
Glugal_60_NETseq_minus_1_readin = textscan(fileID1,'%s %f %f %f', ...
    'CommentStyle','#');
fclose(fileID1);
whos Glugal_60_NETseq_minus_1_readin



fileID2 = fopen([dataloaddir_60 'by_rpb3FLAG_gal60_repeat1_ntsub_unsmo_posstrand.wig']);
Glugal_60_NETseq_plus_1_readin = textscan(fileID2,'%s %f %f %f', ...
    'CommentStyle','#');
fclose(fileID2);
whos Glugal_60_NETseq_plus_1_readin

disp('Finished importing raw Ulku NETseq dst1del sample 1')



fileID3 = fopen([dataloaddir_60 'by_rpb3FLAG_gal60_repeat2_ntsub_unsmo_negstrand.wig']);
Glugal_60_NETseq_minus_2_readin = textscan(fileID3,'%s %f %f %f', ...
    'CommentStyle','#');
fclose(fileID3);
whos Glugal_60_NETseq_minus_2_readin



fileID4 = fopen([dataloaddir_60 'by_rpb3FLAG_gal60_repeat2_ntsub_unsmo_posstrand.wig']);
Glugal_60_NETseq_plus_2_readin = textscan(fileID4,'%s %f %f %f', ...
    'CommentStyle','#');
fclose(fileID4);
whos Glugal_60_NETseq_plus_2_readin

disp('Finished importing raw Ulku NETseq dst1del sample 2')



% %% Load in the spt4 delete
% 
% dataloaddir_end = 'notag_normalised/spt4del_nnt/';
% 
% dataloaddir_spt4 = [dataloaddir_start dataloaddir_end];
% 
% 
% 
% fileID1 = fopen([dataloaddir_spt4 'spt4del_july18_nnt_negativestrand.wig']);
% ulku_spt4del_NETseq_minus_1_readin = textscan(fileID1,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID1);
% whos ulku_spt4del_NETseq_minus_1_readin
% 
% 
% 
% fileID2 = fopen([dataloaddir_spt4 'spt4del_july18_nnt_positivestrand.wig']);
% ulku_spt4del_NETseq_plus_1_readin = textscan(fileID2,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID2);
% whos ulku_spt4del_NETseq_plus_1_readin
% 
% disp('Finished importing raw Ulku NETseq spt4del sample 1')
% 
% 
% 
% fileID3 = fopen([dataloaddir_spt4 'spt4del_aug18_nnt_negativestrand.wig']);
% ulku_spt4del_NETseq_minus_2_readin = textscan(fileID3,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID3);
% whos ulku_spt4del_NETseq_minus_2_readin
% 
% 
% 
% fileID4 = fopen([dataloaddir_spt4 'spt4del_aug18_nnt_positivestrand.wig']);
% ulku_spt4del_NETseq_plus_2_readin = textscan(fileID4,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID4);
% whos ulku_spt4del_NETseq_plus_2_readin
% 
% disp('Finished importing raw Ulku NETseq spt4del sample 2')
% 
% 
% 
% %% Load in the xrn1 delete
% 
% dataloaddir_end = 'notag_normalised/xrn1del_nnt/';
% 
% dataloaddir_xrn1 = [dataloaddir_start dataloaddir_end];
% 
% 
% 
% fileID1 = fopen([dataloaddir_xrn1 'xrn1del_nnt_rpt1_negativestrand.wig']);
% ulku_xrn1del_NETseq_minus_1_readin = textscan(fileID1,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID1);
% whos ulku_xrn1del_NETseq_minus_1_readin
% 
% 
% 
% fileID2 = fopen([dataloaddir_xrn1 'xrn1del_nnt_rpt1_positivestrand.wig']);
% ulku_xrn1del_NETseq_plus_1_readin = textscan(fileID2,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID2);
% whos ulku_xrn1del_NETseq_plus_1_readin
% 
% disp('Finished importing raw Ulku NETseq xrn1del sample 1')
% 
% 
% 
% fileID3 = fopen([dataloaddir_xrn1 'xrn1del_nnt_rpt2_negativestrand.wig']);
% ulku_xrn1del_NETseq_minus_2_readin = textscan(fileID3,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID3);
% whos ulku_xrn1del_NETseq_minus_2_readin
% 
% 
% 
% fileID4 = fopen([dataloaddir_xrn1 'xrn1del_nnt_rpt2_positivestrand.wig']);
% ulku_xrn1del_NETseq_plus_2_readin = textscan(fileID4,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID4);
% whos ulku_xrn1del_NETseq_plus_2_readin
% 
% disp('Finished importing raw Ulku NETseq xrn1del sample 2')
% 
% 
% %% Load in the SPT4FRB Rapa 
% %Anchor Away SPT4 (I think)
% 
% dataloaddir_end = 'AnchorAway/';
% 
% dataloaddir_SPT4AA = [dataloaddir_start dataloaddir_end];
% 
% 
% 
% fileID1 = fopen([dataloaddir_SPT4AA 'Spt4FRB_Rapa_rpt1_negativestrand.wig']);
% ulku_SPT4AA_NETseq_minus_1_readin = textscan(fileID1,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID1);
% whos ulku_SPT4AA_NETseq_minus_1_readin
% 
% 
% 
% fileID2 = fopen([dataloaddir_SPT4AA 'Spt4FRB_Rapa_rpt1_positivestrand.wig']);
% ulku_SPT4AA_NETseq_plus_1_readin = textscan(fileID2,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID2);
% whos ulku_SPT4AA_NETseq_plus_1_readin
% 
% disp('Finished importing raw Ulku NETseq SPT4AA sample 1')
% 
% 
% 
% fileID3 = fopen([dataloaddir_SPT4AA 'Spt4FRB_Rapa_rpt2_negativestrand.wig']);
% ulku_SPT4AA_NETseq_minus_2_readin = textscan(fileID3,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID3);
% whos ulku_SPT4AA_NETseq_minus_2_readin
% 
% 
% 
% fileID4 = fopen([dataloaddir_SPT4AA 'Spt4FRB_Rapa_rpt2_positivestrand.wig']);
% ulku_SPT4AA_NETseq_plus_2_readin = textscan(fileID4,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID4);
% whos ulku_SPT4AA_NETseq_plus_2_readin
% 
% disp('Finished importing raw Ulku NETseq SPT4AA sample 2')
% 
% 
% %% Load in the SPT5FRB Rapa 
% %Anchor Away SPT5 (I think)
% 
% dataloaddir_end = 'AnchorAway/';
% 
% dataloaddir_SPT5AA = [dataloaddir_start dataloaddir_end];
% 
% 
% 
% fileID1 = fopen([dataloaddir_SPT5AA 'Spt5FRB_Rapa_rpt1_negativestrand.wig']);
% ulku_SPT5AA_NETseq_minus_1_readin = textscan(fileID1,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID1);
% whos ulku_SPT5AA_NETseq_minus_1_readin
% 
% 
% 
% fileID2 = fopen([dataloaddir_SPT5AA 'Spt5FRB_Rapa_rpt1_positivestrand.wig']);
% ulku_SPT5AA_NETseq_plus_1_readin = textscan(fileID2,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID2);
% whos ulku_SPT5AA_NETseq_plus_1_readin
% 
% disp('Finished importing raw Ulku NETseq SPT5AA sample 1')
% 
% 
% 
% fileID3 = fopen([dataloaddir_SPT5AA 'Spt5FRB_Rapa_rpt2_negativestrand.wig']);
% ulku_SPT5AA_NETseq_minus_2_readin = textscan(fileID3,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID3);
% whos ulku_SPT5AA_NETseq_minus_2_readin
% 
% 
% 
% fileID4 = fopen([dataloaddir_SPT5AA 'Spt5FRB_Rapa_rpt2_positivestrand.wig']);
% ulku_SPT5AA_NETseq_plus_2_readin = textscan(fileID4,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID4);
% whos ulku_SPT5AA_NETseq_plus_2_readin
% 
% disp('Finished importing raw Ulku NETseq SPT4AA sample 2')
% 
% 
% 
% %% Load in the SPT4FRB DMSO
% %Anchor Away SPT4 control (I think)
% 
% dataloaddir_end = 'AnchorAway/';
% 
% dataloaddir_SPT4AADMSO = [dataloaddir_start dataloaddir_end];
% 
% 
% 
% fileID1 = fopen([dataloaddir_SPT4AADMSO 'Spt4FRB_DMSO_negativestrand.wig']);
% ulku_SPT4AADMSO_NETseq_minus_1_readin = textscan(fileID1,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID1);
% whos ulku_SPT4AADMSO_NETseq_minus_1_readin
% 
% 
% 
% fileID2 = fopen([dataloaddir_SPT4AADMSO 'Spt4FRB_DMSO_positivestrand.wig']);
% ulku_SPT4AADMSO_NETseq_plus_1_readin = textscan(fileID2,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID2);
% whos ulku_SPT4AADMSO_NETseq_plus_1_readin
% 
% disp('Finished importing raw Ulku NETseq SPT4AA DMSO')
% 
% 
% %% Load in the SPT5FRB DMSO
% %Anchor Away SPT5 control (I think)
% 
% dataloaddir_end = 'AnchorAway/';
% 
% dataloaddir_SPT5AADMSO = [dataloaddir_start dataloaddir_end];
% 
% 
% 
% fileID1 = fopen([dataloaddir_SPT5AADMSO 'Spt5FRB_DMSO_negativestrand.wig']);
% ulku_SPT5AADMSO_NETseq_minus_1_readin = textscan(fileID1,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID1);
% whos ulku_SPT5AADMSO_NETseq_minus_1_readin
% 
% 
% 
% fileID2 = fopen([dataloaddir_SPT5AADMSO 'Spt5FRB_DMSO_positivestrand.wig']);
% ulku_SPT5AADMSO_NETseq_plus_1_readin = textscan(fileID2,'%s %f %f %f', ...
%     'CommentStyle','#');
% fclose(fileID2);
% whos ulku_SPT5AADMSO_NETseq_plus_1_readin
% 
% disp('Finished importing raw Ulku NETseq SPT5AA DMSO')
% 
% 
% 
% 
% %% NB the DMSO SPT4FRB and SPT5FRB may be combined
% %Only one repeat for each but they /it{should} be the same but also valid
% %to consider them each separately - they should identify any effect of the
% %FRB tag (I think).


%% Calculate the max point for each strand
%(might later need to extend it to the true points as could cause problems
%if annotated regions lie outside the range of measurement)
%Also, might need to extend it to all the different samples
%(in case additional signal is measured at the end of chromosomes in a
%specific sample)

%Do a separate loop for each sample so that adding additional samples is
%less hassle

%choosing to have the same extent for each stand for simplicity

max_chromosome_extent = zeros(NoChromosomes,1);
min_chromosome_loc = NaN*ones(NoChromosomes,1);

%WT sample 1

for chrctr=1:1:NoChromosomes

    %minus strand
    
    temp_chromo_locs = strfind(Glugal_15_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
    
    temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
    
    temp_max1 = max(max( ...
        [Glugal_15_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
        Glugal_15_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
    
    temp_min1 = min(min( ...
        [Glugal_15_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
        Glugal_15_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
    
    
    max_chromosome_extent(chrctr) = ...
        max(max_chromosome_extent(chrctr),temp_max1);
    min_chromosome_loc(chrctr) = ...
        min(min_chromosome_loc(chrctr),temp_min1);
    
    
    %plus strand
    
    temp_chromo_locs = strfind(Glugal_15_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
    
    temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
    
    temp_max1 = max(max( ...
        [Glugal_15_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
        Glugal_15_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
    
    temp_min1 = min(min( ...
        [Glugal_15_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
        Glugal_15_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
    
    
    max_chromosome_extent(chrctr) = ...
        max(max_chromosome_extent(chrctr),temp_max1);
    min_chromosome_loc(chrctr) = ...
        min(min_chromosome_loc(chrctr),temp_min1);
    
end


%WT sample 2

for chrctr=1:1:NoChromosomes

    %minus strand
    
    temp_chromo_locs = strfind(Glugal_15_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
    
    temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
    
    temp_max1 = max(max( ...
        [Glugal_15_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
        Glugal_15_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
    
    temp_min1 = min(min( ...
        [Glugal_15_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
        Glugal_15_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
    
    
    max_chromosome_extent(chrctr) = ...
        max(max_chromosome_extent(chrctr),temp_max1);
    min_chromosome_loc(chrctr) = ...
        min(min_chromosome_loc(chrctr),temp_min1);
    
    
    %plus strand
        
    temp_chromo_locs = strfind(Glugal_15_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
    
    temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
    
    temp_max1 = max(max( ...
        [Glugal_15_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
        Glugal_15_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
    
    temp_min1 = min(min( ...
        [Glugal_15_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
        Glugal_15_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
    
    
    max_chromosome_extent(chrctr) = ...
        max(max_chromosome_extent(chrctr),temp_max1);
    min_chromosome_loc(chrctr) = ...
        min(min_chromosome_loc(chrctr),temp_min1);
    
end




%dst1del sample 1

for chrctr=1:1:NoChromosomes

    %minus strand
    
    temp_chromo_locs = strfind(Glugal_60_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
    
    temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
    
    temp_max1 = max(max( ...
        [Glugal_60_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
        Glugal_60_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
    
    temp_min1 = min(min( ...
        [Glugal_60_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
        Glugal_60_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
    
    
    max_chromosome_extent(chrctr) = ...
        max(max_chromosome_extent(chrctr),temp_max1);
    min_chromosome_loc(chrctr) = ...
        min(min_chromosome_loc(chrctr),temp_min1);
    
    
    %plus strand
    
    temp_chromo_locs = strfind(Glugal_60_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
    
    temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
    
    temp_max1 = max(max( ...
        [Glugal_60_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
        Glugal_60_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
    
    temp_min1 = min(min( ...
        [Glugal_60_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
        Glugal_60_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
    
    
    max_chromosome_extent(chrctr) = ...
        max(max_chromosome_extent(chrctr),temp_max1);
    min_chromosome_loc(chrctr) = ...
        min(min_chromosome_loc(chrctr),temp_min1);
    
end


%dst1del sample 2

for chrctr=1:1:NoChromosomes

    %minus strand
    
    temp_chromo_locs = strfind(Glugal_60_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
    
    temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
    
    temp_max1 = max(max( ...
        [Glugal_60_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
        Glugal_60_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
    
    temp_min1 = min(min( ...
        [Glugal_60_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
        Glugal_60_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
    
    
    max_chromosome_extent(chrctr) = ...
        max(max_chromosome_extent(chrctr),temp_max1);
    min_chromosome_loc(chrctr) = ...
        min(min_chromosome_loc(chrctr),temp_min1);
    
    
    %plus strand
        
    temp_chromo_locs = strfind(Glugal_60_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
    
    temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
    
    temp_max1 = max(max( ...
        [Glugal_60_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
        Glugal_60_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
    
    temp_min1 = min(min( ...
        [Glugal_60_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
        Glugal_60_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
    
    
    max_chromosome_extent(chrctr) = ...
        max(max_chromosome_extent(chrctr),temp_max1);
    min_chromosome_loc(chrctr) = ...
        min(min_chromosome_loc(chrctr),temp_min1);
    
end




% %spt4del sample 1
% 
% for chrctr=1:1:NoChromosomes
% 
%     %minus strand
%     
%     temp_chromo_locs = strfind(ulku_spt4del_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_spt4del_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
%         ulku_spt4del_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_spt4del_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
%         ulku_spt4del_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
%     
%     %plus strand
%     
%     temp_chromo_locs = strfind(ulku_spt4del_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_spt4del_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
%         ulku_spt4del_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_spt4del_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
%         ulku_spt4del_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
% end
% 
% 
% %spt4del sample 2
% 
% for chrctr=1:1:NoChromosomes
% 
%     %minus strand
%     
%     temp_chromo_locs = strfind(ulku_spt4del_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_spt4del_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
%         ulku_spt4del_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_spt4del_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
%         ulku_spt4del_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
%     
%     %plus strand
%         
%     temp_chromo_locs = strfind(ulku_spt4del_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_spt4del_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
%         ulku_spt4del_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_spt4del_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
%         ulku_spt4del_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
% end
% 
% 
% 
% 
% %xrn1del sample 1
% 
% for chrctr=1:1:NoChromosomes
% 
%     %minus strand
%     
%     temp_chromo_locs = strfind(ulku_xrn1del_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_xrn1del_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
%         ulku_xrn1del_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_xrn1del_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
%         ulku_xrn1del_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
%     
%     %plus strand
%     
%     temp_chromo_locs = strfind(ulku_xrn1del_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_xrn1del_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
%         ulku_xrn1del_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_xrn1del_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
%         ulku_xrn1del_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
% end
% 
% 
% %xrn1del sample 2
% 
% for chrctr=1:1:NoChromosomes
% 
%     %minus strand
%     
%     temp_chromo_locs = strfind(ulku_xrn1del_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_xrn1del_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
%         ulku_xrn1del_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_xrn1del_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
%         ulku_xrn1del_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
%     
%     %plus strand
%         
%     temp_chromo_locs = strfind(ulku_xrn1del_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_xrn1del_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
%         ulku_xrn1del_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_xrn1del_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
%         ulku_xrn1del_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
% end
% 
% 
% 
% 
% 
% %SPT4AA sample 1
% 
% for chrctr=1:1:NoChromosomes
% 
%     %minus strand
%     
%     temp_chromo_locs = strfind(ulku_SPT4AA_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_SPT4AA_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
%         ulku_SPT4AA_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_SPT4AA_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
%         ulku_SPT4AA_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
%     
%     %plus strand
%     
%     temp_chromo_locs = strfind(ulku_SPT4AA_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_SPT4AA_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
%         ulku_SPT4AA_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_SPT4AA_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
%         ulku_SPT4AA_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
% end
% 
% 
% %SPT4AA sample 2
% 
% for chrctr=1:1:NoChromosomes
% 
%     %minus strand
%     
%     temp_chromo_locs = strfind(ulku_SPT4AA_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_SPT4AA_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
%         ulku_SPT4AA_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_SPT4AA_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
%         ulku_SPT4AA_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
%     
%     %plus strand
%         
%     temp_chromo_locs = strfind(ulku_SPT4AA_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_SPT4AA_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
%         ulku_SPT4AA_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_SPT4AA_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
%         ulku_SPT4AA_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
% end
% 
% 
% 
% 
% 
% %SPT5AA sample 1
% 
% for chrctr=1:1:NoChromosomes
% 
%     %minus strand
%     
%     temp_chromo_locs = strfind(ulku_SPT5AA_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_SPT5AA_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
%         ulku_SPT5AA_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_SPT5AA_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
%         ulku_SPT5AA_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
%     
%     %plus strand
%     
%     temp_chromo_locs = strfind(ulku_SPT5AA_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_SPT5AA_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
%         ulku_SPT5AA_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_SPT5AA_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
%         ulku_SPT5AA_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
% end
% 
% 
% %SPT5AA sample 2
% 
% for chrctr=1:1:NoChromosomes
% 
%     %minus strand
%     
%     temp_chromo_locs = strfind(ulku_SPT5AA_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_SPT5AA_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
%         ulku_SPT5AA_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_SPT5AA_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
%         ulku_SPT5AA_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
%     
%     %plus strand
%         
%     temp_chromo_locs = strfind(ulku_SPT5AA_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_SPT5AA_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
%         ulku_SPT5AA_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_SPT5AA_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
%         ulku_SPT5AA_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
% end
% 
% 
% 
% %SPT4AA DMSO
% 
% for chrctr=1:1:NoChromosomes
% 
%     %minus strand
%     
%     temp_chromo_locs = strfind(ulku_SPT4AADMSO_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_SPT4AADMSO_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
%         ulku_SPT4AADMSO_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_SPT4AADMSO_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
%         ulku_SPT4AADMSO_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
%     
%     %plus strand
%     
%     temp_chromo_locs = strfind(ulku_SPT4AADMSO_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_SPT4AADMSO_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
%         ulku_SPT4AADMSO_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_SPT4AADMSO_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
%         ulku_SPT4AADMSO_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
% end
% 
% 
% 
% %SPT5AA DMSO
% 
% for chrctr=1:1:NoChromosomes
% 
%     %minus strand
%     
%     temp_chromo_locs = strfind(ulku_SPT5AADMSO_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_SPT5AADMSO_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
%         ulku_SPT5AADMSO_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_SPT5AADMSO_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
%         ulku_SPT5AADMSO_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
%     
%     %plus strand
%     
%     temp_chromo_locs = strfind(ulku_SPT5AADMSO_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [ulku_SPT5AADMSO_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
%         ulku_SPT5AADMSO_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [ulku_SPT5AADMSO_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
%         ulku_SPT5AADMSO_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
%     
%     
%     max_chromosome_extent(chrctr) = ...
%         max(max_chromosome_extent(chrctr),temp_max1);
%     min_chromosome_loc(chrctr) = ...
%         min(min_chromosome_loc(chrctr),temp_min1);
%     
% end



%%

Glugal_15_NETseq_1 = cell(NoChromosomes,2);
Glugal_15_NETseq_2 = cell(NoChromosomes,2);

for chrctr=1:NoChromosomes
    
    for strandctr = 1:2
   
        Glugal_15_NETseq_1{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
        Glugal_15_NETseq_2{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
    
    end
    
end


Glugal_60_NETseq_1 = cell(NoChromosomes,2);
Glugal_60_NETseq_2 = cell(NoChromosomes,2);

for chrctr=1:NoChromosomes
    
    for strandctr = 1:2
   
        Glugal_60_NETseq_1{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
        Glugal_60_NETseq_2{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
    
    end
    
end


Ulku_spt4del_NETseq_1 = cell(NoChromosomes,2);
Ulku_spt4del_NETseq_2 = cell(NoChromosomes,2);

for chrctr=1:NoChromosomes
    
    for strandctr = 1:2
   
        Ulku_spt4del_NETseq_1{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
        Ulku_spt4del_NETseq_2{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
    
    end
    
end


Ulku_xrn1del_NETseq_1 = cell(NoChromosomes,2);
Ulku_xrn1del_NETseq_2 = cell(NoChromosomes,2);

for chrctr=1:NoChromosomes
    
    for strandctr = 1:2
   
        Ulku_xrn1del_NETseq_1{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
        Ulku_xrn1del_NETseq_2{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
    
    end
    
end


Ulku_SPT4AA_NETseq_1 = cell(NoChromosomes,2);
Ulku_SPT4AA_NETseq_2 = cell(NoChromosomes,2);

for chrctr=1:NoChromosomes
    
    for strandctr = 1:2
   
        Ulku_SPT4AA_NETseq_1{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
        Ulku_SPT4AA_NETseq_2{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
    
    end
    
end


Ulku_SPT5AA_NETseq_1 = cell(NoChromosomes,2);
Ulku_SPT5AA_NETseq_2 = cell(NoChromosomes,2);

for chrctr=1:NoChromosomes
    
    for strandctr = 1:2
   
        Ulku_SPT5AA_NETseq_1{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
        Ulku_SPT5AA_NETseq_2{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
    
    end
    
end


Ulku_SPT4AADMSO_NETseq_1 = cell(NoChromosomes,2);
Ulku_SPT5AADMSO_NETseq_1 = cell(NoChromosomes,2);

for chrctr=1:NoChromosomes
    
    for strandctr = 1:2
   
        Ulku_SPT4AADMSO_NETseq_1{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
        Ulku_SPT5AADMSO_NETseq_1{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
    
    end
    
end



%% Process the WT data

%NB I am going from the minimum chromosome location +1 to the end location
%might mean things are out by a basepair depending on how things were
%supposed to be interpreted - referencing in the original file starts from
%0 so I am assuming this is the correct way

for chrctr=1:1:NoChromosomes

    %sample1
    
    %minus strand
    tempmatches = strfind(Glugal_15_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
    
    tempmatches_processed = ~cellfun(@isempty,tempmatches);
    
    temp_starts = Glugal_15_NETseq_minus_1_readin{2}(tempmatches_processed);
    temp_ends = Glugal_15_NETseq_minus_1_readin{3}(tempmatches_processed);
    temp_vals = Glugal_15_NETseq_minus_1_readin{4}(tempmatches_processed);
    
    
    for ictr=1:1:sum(tempmatches_processed)
       
        Glugal_15_NETseq_1{chrctr,2}(...
            (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
        
        if mod(ictr,10000) == 0
            disp(['Completed ' num2str(ictr) ' of ' ...
                num2str(length(temp_starts)) ...
                ' 1 minus iterations for chromosome' num2str(chrctr) '.'])
        end
        
    end
    
    %plus strand
    tempmatches = strfind(Glugal_15_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
    
    tempmatches_processed = ~cellfun(@isempty,tempmatches);
    
    temp_starts = Glugal_15_NETseq_plus_1_readin{2}(tempmatches_processed);
    temp_ends = Glugal_15_NETseq_plus_1_readin{3}(tempmatches_processed);
    temp_vals = Glugal_15_NETseq_plus_1_readin{4}(tempmatches_processed);
    
    
    for ictr=1:1:sum(tempmatches_processed)
       
        Glugal_15_NETseq_1{chrctr,1}(...
            (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
        
        if mod(ictr,10000) == 0
            disp(['Completed ' num2str(ictr) ' of ' ...
                num2str(length(temp_starts)) ...
                ' 1 plus iterations for chromosome' num2str(chrctr) '.'])
        end
        
    end
    
    
    
    %sample2
    
    %minus strand
    tempmatches = strfind(Glugal_15_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
    
    tempmatches_processed = ~cellfun(@isempty,tempmatches);
    
    temp_starts = Glugal_15_NETseq_minus_2_readin{2}(tempmatches_processed);
    temp_ends = Glugal_15_NETseq_minus_2_readin{3}(tempmatches_processed);
    temp_vals = Glugal_15_NETseq_minus_2_readin{4}(tempmatches_processed);
    
    
    for ictr=1:1:sum(tempmatches_processed)
       
        Glugal_15_NETseq_2{chrctr,2}(...
            (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
        
        if mod(ictr,10000) == 0
            disp(['Completed ' num2str(ictr) ' of ' ...
                num2str(length(temp_starts)) ...
                ' 2 minus iterations for chromosome' num2str(chrctr) '.'])
        end
        
    end
    
    
    %plus strand
    tempmatches = strfind(Glugal_15_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
    
    tempmatches_processed = ~cellfun(@isempty,tempmatches);
    
    temp_starts = Glugal_15_NETseq_plus_2_readin{2}(tempmatches_processed);
    temp_ends = Glugal_15_NETseq_plus_2_readin{3}(tempmatches_processed);
    temp_vals = Glugal_15_NETseq_plus_2_readin{4}(tempmatches_processed);
    
    
    for ictr=1:1:sum(tempmatches_processed)
       
        Glugal_15_NETseq_2{chrctr,1}(...
            (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
        
        if mod(ictr,10000) == 0
            disp(['Completed ' num2str(ictr) ' of ' ...
                num2str(length(temp_starts)) ...
                ' 2 plus iterations for chromosome' num2str(chrctr) '.'])
        end
        
    end
    
    
end


%%

save([dataloaddir_15 'Glugal_15_NETseq_1.mat'],'Glugal_15_NETseq_1','-v7.3')


%%

save([dataloaddir_15 'Glugal_15_NETseq_2.mat'],'Glugal_15_NETseq_2','-v7.3')






%% Process the dst1del data

%NB I am going from the minimum chromosome location +1 to the end location
%might mean things are out by a basepair depending on how things were
%supposed to be interpreted - referencing in the original file starts from
%0 so I am assuming this is the correct way

for chrctr=1:1:NoChromosomes

    %sample1
    
    %minus strand
    tempmatches = strfind(Glugal_60_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
    
    tempmatches_processed = ~cellfun(@isempty,tempmatches);
    
    temp_starts = Glugal_60_NETseq_minus_1_readin{2}(tempmatches_processed);
    temp_ends = Glugal_60_NETseq_minus_1_readin{3}(tempmatches_processed);
    temp_vals = Glugal_60_NETseq_minus_1_readin{4}(tempmatches_processed);
    
    
    for ictr=1:1:sum(tempmatches_processed)
       
        Glugal_60_NETseq_1{chrctr,2}(...
            (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
        
        if mod(ictr,10000) == 0
            disp(['Completed ' num2str(ictr) ' of ' ...
                num2str(length(temp_starts)) ...
                ' 1 minus iterations for chromosome' num2str(chrctr) '.'])
        end
        
    end
    
    %plus strand
    tempmatches = strfind(Glugal_60_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
    
    tempmatches_processed = ~cellfun(@isempty,tempmatches);
    
    temp_starts = Glugal_60_NETseq_plus_1_readin{2}(tempmatches_processed);
    temp_ends = Glugal_60_NETseq_plus_1_readin{3}(tempmatches_processed);
    temp_vals = Glugal_60_NETseq_plus_1_readin{4}(tempmatches_processed);
    
    
    for ictr=1:1:sum(tempmatches_processed)
       
        Glugal_60_NETseq_1{chrctr,1}(...
            (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
        
        if mod(ictr,10000) == 0
            disp(['Completed ' num2str(ictr) ' of ' ...
                num2str(length(temp_starts)) ...
                ' 1 plus iterations for chromosome' num2str(chrctr) '.'])
        end
        
    end
    
    
    
    %sample2
    
    %minus strand
    tempmatches = strfind(Glugal_60_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
    
    tempmatches_processed = ~cellfun(@isempty,tempmatches);
    
    temp_starts = Glugal_60_NETseq_minus_2_readin{2}(tempmatches_processed);
    temp_ends = Glugal_60_NETseq_minus_2_readin{3}(tempmatches_processed);
    temp_vals = Glugal_60_NETseq_minus_2_readin{4}(tempmatches_processed);
    
    
    for ictr=1:1:sum(tempmatches_processed)
       
        Glugal_60_NETseq_2{chrctr,2}(...
            (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
        
        if mod(ictr,10000) == 0
            disp(['Completed ' num2str(ictr) ' of ' ...
                num2str(length(temp_starts)) ...
                ' 2 minus iterations for chromosome' num2str(chrctr) '.'])
        end
        
    end
    
    
    %plus strand
    tempmatches = strfind(Glugal_60_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
    
    tempmatches_processed = ~cellfun(@isempty,tempmatches);
    
    temp_starts = Glugal_60_NETseq_plus_2_readin{2}(tempmatches_processed);
    temp_ends = Glugal_60_NETseq_plus_2_readin{3}(tempmatches_processed);
    temp_vals = Glugal_60_NETseq_plus_2_readin{4}(tempmatches_processed);
    
    
    for ictr=1:1:sum(tempmatches_processed)
       
        Glugal_60_NETseq_2{chrctr,1}(...
            (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
        
        if mod(ictr,10000) == 0
            disp(['Completed ' num2str(ictr) ' of ' ...
                num2str(length(temp_starts)) ...
                ' 2 plus iterations for chromosome' num2str(chrctr) '.'])
        end
        
    end
    
    
end


%%

save([dataloaddir_60 'Glugal_60_NETseq_1.mat'],'Glugal_60_NETseq_1','-v7.3')


%%

save([dataloaddir_60 'Glugal_60_NETseq_2.mat'],'Glugal_60_NETseq_2','-v7.3')












% %% Process the spt4del data
% 
% %STILL TO BE UPDATED
% 
% %NB I am going from the minimum chromosome location +1 to the end location
% %might mean things are out by a basepair depending on how things were
% %supposed to be interpreted - referencing in the original file starts from
% %0 so I am assuming this is the correct way
% 
% for chrctr=1:1:NoChromosomes
% 
%     %sample1
%     
%     %minus strand
%     tempmatches = strfind(ulku_spt4del_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_spt4del_NETseq_minus_1_readin{2}(tempmatches_processed);
%     temp_ends = ulku_spt4del_NETseq_minus_1_readin{3}(tempmatches_processed);
%     temp_vals = ulku_spt4del_NETseq_minus_1_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_spt4del_NETseq_1{chrctr,2}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 1 minus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     %plus strand
%     tempmatches = strfind(ulku_spt4del_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_spt4del_NETseq_plus_1_readin{2}(tempmatches_processed);
%     temp_ends = ulku_spt4del_NETseq_plus_1_readin{3}(tempmatches_processed);
%     temp_vals = ulku_spt4del_NETseq_plus_1_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_spt4del_NETseq_1{chrctr,1}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 1 plus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     
%     
%     %sample2
%     
%     %minus strand
%     tempmatches = strfind(ulku_spt4del_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_spt4del_NETseq_minus_2_readin{2}(tempmatches_processed);
%     temp_ends = ulku_spt4del_NETseq_minus_2_readin{3}(tempmatches_processed);
%     temp_vals = ulku_spt4del_NETseq_minus_2_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_spt4del_NETseq_2{chrctr,2}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 2 minus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     
%     %plus strand
%     tempmatches = strfind(ulku_spt4del_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_spt4del_NETseq_plus_2_readin{2}(tempmatches_processed);
%     temp_ends = ulku_spt4del_NETseq_plus_2_readin{3}(tempmatches_processed);
%     temp_vals = ulku_spt4del_NETseq_plus_2_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_spt4del_NETseq_2{chrctr,1}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 2 plus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     
% end
% 
% 
% %%
% 
% save([dataloaddir_spt4 'ulku_spt4del_NETseq_1.mat'],'Ulku_spt4del_NETseq_1','-v7.3')
% 
% 
% %%
% 
% save([dataloaddir_spt4 'ulku_spt4del_NETseq_2.mat'],'Ulku_spt4del_NETseq_2','-v7.3')
% 
% 
% 
% 
% 
% %% Process the xrn1del data
% 
% %NEED TO FINISH UPDATING THIS!
% 
% %NB I am going from the minimum chromosome location +1 to the end location
% %might mean things are out by a basepair depending on how things were
% %supposed to be interpreted - referencing in the original file starts from
% %0 so I am assuming this is the correct way
% 
% for chrctr=1:1:NoChromosomes
% 
%     %sample1
%     
%     %minus strand
%     tempmatches = strfind(ulku_xrn1del_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_xrn1del_NETseq_minus_1_readin{2}(tempmatches_processed);
%     temp_ends = ulku_xrn1del_NETseq_minus_1_readin{3}(tempmatches_processed);
%     temp_vals = ulku_xrn1del_NETseq_minus_1_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_xrn1del_NETseq_1{chrctr,2}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 1 minus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     %plus strand
%     tempmatches = strfind(ulku_xrn1del_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_xrn1del_NETseq_plus_1_readin{2}(tempmatches_processed);
%     temp_ends = ulku_xrn1del_NETseq_plus_1_readin{3}(tempmatches_processed);
%     temp_vals = ulku_xrn1del_NETseq_plus_1_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_xrn1del_NETseq_1{chrctr,1}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 1 plus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     
%     
%     %sample2
%     
%     %minus strand
%     tempmatches = strfind(ulku_xrn1del_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_xrn1del_NETseq_minus_2_readin{2}(tempmatches_processed);
%     temp_ends = ulku_xrn1del_NETseq_minus_2_readin{3}(tempmatches_processed);
%     temp_vals = ulku_xrn1del_NETseq_minus_2_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_xrn1del_NETseq_2{chrctr,2}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 2 minus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     
%     %plus strand
%     tempmatches = strfind(ulku_xrn1del_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_xrn1del_NETseq_plus_2_readin{2}(tempmatches_processed);
%     temp_ends = ulku_xrn1del_NETseq_plus_2_readin{3}(tempmatches_processed);
%     temp_vals = ulku_xrn1del_NETseq_plus_2_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_xrn1del_NETseq_2{chrctr,1}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 2 plus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     
% end
% 
% 
% %%
% 
% save([dataloaddir_xrn1 'ulku_xrn1del_NETseq_1.mat'],'Ulku_xrn1del_NETseq_1','-v7.3')
% 
% 
% %%
% 
% save([dataloaddir_xrn1 'ulku_xrn1del_NETseq_2.mat'],'Ulku_xrn1del_NETseq_2','-v7.3')
% 
% 
% 
% %% Process the SPT4AA data
% 
% %NB I am going from the minimum chromosome location +1 to the end location
% %might mean things are out by a basepair depending on how things were
% %supposed to be interpreted - referencing in the original file starts from
% %0 so I am assuming this is the correct way
% 
% for chrctr=1:1:NoChromosomes
% 
%     %sample1
%     
%     %minus strand
%     tempmatches = strfind(ulku_SPT4AA_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_SPT4AA_NETseq_minus_1_readin{2}(tempmatches_processed);
%     temp_ends = ulku_SPT4AA_NETseq_minus_1_readin{3}(tempmatches_processed);
%     temp_vals = ulku_SPT4AA_NETseq_minus_1_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_SPT4AA_NETseq_1{chrctr,2}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 1 minus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     %plus strand
%     tempmatches = strfind(ulku_SPT4AA_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_SPT4AA_NETseq_plus_1_readin{2}(tempmatches_processed);
%     temp_ends = ulku_SPT4AA_NETseq_plus_1_readin{3}(tempmatches_processed);
%     temp_vals = ulku_SPT4AA_NETseq_plus_1_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_SPT4AA_NETseq_1{chrctr,1}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 1 plus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     
%     
%     %sample2
%     
%     %minus strand
%     tempmatches = strfind(ulku_SPT4AA_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_SPT4AA_NETseq_minus_2_readin{2}(tempmatches_processed);
%     temp_ends = ulku_SPT4AA_NETseq_minus_2_readin{3}(tempmatches_processed);
%     temp_vals = ulku_SPT4AA_NETseq_minus_2_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_SPT4AA_NETseq_2{chrctr,2}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 2 minus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     
%     %plus strand
%     tempmatches = strfind(ulku_SPT4AA_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_SPT4AA_NETseq_plus_2_readin{2}(tempmatches_processed);
%     temp_ends = ulku_SPT4AA_NETseq_plus_2_readin{3}(tempmatches_processed);
%     temp_vals = ulku_SPT4AA_NETseq_plus_2_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_SPT4AA_NETseq_2{chrctr,1}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 2 plus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     
% end
% 
% 
% %%
% 
% save([dataloaddir_SPT4AA 'ulku_SPT4AA_NETseq_1.mat'],'Ulku_SPT4AA_NETseq_1','-v7.3')
% 
% 
% %%
% 
% save([dataloaddir_SPT4AA 'ulku_SPT4AA_NETseq_2.mat'],'Ulku_SPT4AA_NETseq_2','-v7.3')
% 
% 
% 
% 
% %% Process the SPT5AA data
% 
% %NB I am going from the minimum chromosome location +1 to the end location
% %might mean things are out by a basepair depending on how things were
% %supposed to be interpreted - referencing in the original file starts from
% %0 so I am assuming this is the correct way
% 
% for chrctr=1:1:NoChromosomes
% 
%     %sample1
%     
%     %minus strand
%     tempmatches = strfind(ulku_SPT5AA_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_SPT5AA_NETseq_minus_1_readin{2}(tempmatches_processed);
%     temp_ends = ulku_SPT5AA_NETseq_minus_1_readin{3}(tempmatches_processed);
%     temp_vals = ulku_SPT5AA_NETseq_minus_1_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_SPT5AA_NETseq_1{chrctr,2}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 1 minus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     %plus strand
%     tempmatches = strfind(ulku_SPT5AA_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_SPT5AA_NETseq_plus_1_readin{2}(tempmatches_processed);
%     temp_ends = ulku_SPT5AA_NETseq_plus_1_readin{3}(tempmatches_processed);
%     temp_vals = ulku_SPT5AA_NETseq_plus_1_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_SPT5AA_NETseq_1{chrctr,1}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 1 plus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     
%     
%     %sample2
%     
%     %minus strand
%     tempmatches = strfind(ulku_SPT5AA_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_SPT5AA_NETseq_minus_2_readin{2}(tempmatches_processed);
%     temp_ends = ulku_SPT5AA_NETseq_minus_2_readin{3}(tempmatches_processed);
%     temp_vals = ulku_SPT5AA_NETseq_minus_2_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_SPT5AA_NETseq_2{chrctr,2}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 2 minus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     
%     %plus strand
%     tempmatches = strfind(ulku_SPT5AA_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_SPT5AA_NETseq_plus_2_readin{2}(tempmatches_processed);
%     temp_ends = ulku_SPT5AA_NETseq_plus_2_readin{3}(tempmatches_processed);
%     temp_vals = ulku_SPT5AA_NETseq_plus_2_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_SPT5AA_NETseq_2{chrctr,1}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 2 plus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     
% end
% 
% 
% %%
% 
% save([dataloaddir_SPT5AA 'ulku_SPT5AA_NETseq_1.mat'],'Ulku_SPT5AA_NETseq_1','-v7.3')
% 
% 
% %%
% 
% save([dataloaddir_SPT5AA 'ulku_SPT5AA_NETseq_2.mat'],'Ulku_SPT5AA_NETseq_2','-v7.3')
% 
% 
% 
% 
% %% Process the SPT4AADMSO data
% 
% %NB I am going from the minimum chromosome location +1 to the end location
% %might mean things are out by a basepair depending on how things were
% %supposed to be interpreted - referencing in the original file starts from
% %0 so I am assuming this is the correct way
% 
% for chrctr=1:1:NoChromosomes
% 
%     %sample1
%     
%     %minus strand
%     tempmatches = strfind(ulku_SPT4AADMSO_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_SPT4AADMSO_NETseq_minus_1_readin{2}(tempmatches_processed);
%     temp_ends = ulku_SPT4AADMSO_NETseq_minus_1_readin{3}(tempmatches_processed);
%     temp_vals = ulku_SPT4AADMSO_NETseq_minus_1_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_SPT4AADMSO_NETseq_1{chrctr,2}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 1 minus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     %plus strand
%     tempmatches = strfind(ulku_SPT4AADMSO_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_SPT4AADMSO_NETseq_plus_1_readin{2}(tempmatches_processed);
%     temp_ends = ulku_SPT4AADMSO_NETseq_plus_1_readin{3}(tempmatches_processed);
%     temp_vals = ulku_SPT4AADMSO_NETseq_plus_1_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_SPT4AADMSO_NETseq_1{chrctr,1}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 1 plus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     
%     
% end
% 
% 
% %%
% 
% save([dataloaddir_SPT4AADMSO 'ulku_SPT4AADMSO_NETseq_1.mat'],'Ulku_SPT4AADMSO_NETseq_1','-v7.3')
% 
% 
% 
% 
% %% Process the SPT5AADMSO data
% 
% %NB I am going from the minimum chromosome location +1 to the end location
% %might mean things are out by a basepair depending on how things were
% %supposed to be interpreted - referencing in the original file starts from
% %0 so I am assuming this is the correct way
% 
% for chrctr=1:1:NoChromosomes
% 
%     %sample1
%     
%     %minus strand
%     tempmatches = strfind(ulku_SPT5AADMSO_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_SPT5AADMSO_NETseq_minus_1_readin{2}(tempmatches_processed);
%     temp_ends = ulku_SPT5AADMSO_NETseq_minus_1_readin{3}(tempmatches_processed);
%     temp_vals = ulku_SPT5AADMSO_NETseq_minus_1_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_SPT5AADMSO_NETseq_1{chrctr,2}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 1 minus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     %plus strand
%     tempmatches = strfind(ulku_SPT5AADMSO_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = ulku_SPT5AADMSO_NETseq_plus_1_readin{2}(tempmatches_processed);
%     temp_ends = ulku_SPT5AADMSO_NETseq_plus_1_readin{3}(tempmatches_processed);
%     temp_vals = ulku_SPT5AADMSO_NETseq_plus_1_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         Ulku_SPT5AADMSO_NETseq_1{chrctr,1}(...
%             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
%         
%         if mod(ictr,10000) == 0
%             disp(['Completed ' num2str(ictr) ' of ' ...
%                 num2str(length(temp_starts)) ...
%                 ' 1 plus iterations for chromosome' num2str(chrctr) '.'])
%         end
%         
%     end
%     
%     
%     
% end
% 
% 
% %%
% 
% save([dataloaddir_SPT5AADMSO 'ulku_SPT5AADMSO_NETseq_1.mat'],'Ulku_SPT5AADMSO_NETseq_1','-v7.3')
% 
% 
% 
% %% NOTES
% 
% %Still needs to be run to pick out obvious bugs
% 
% %Also needs to be check for any persisting silly errors like writing the
% %wrong data to the save file or processing the wrong data in a given
% %section