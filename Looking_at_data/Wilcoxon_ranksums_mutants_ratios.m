%% Compare mutants to WT with wilcoxon rank sum test



%load in the WT ratio matrices
dataloaddirr = '/home/orie3770/Current_Cluster_Job90/';
loadfile = 'NETseqWTratios.mat';

load([dataloaddirr loadfile])


%rename the WT variables so not overwrittewn by mutant variables when
%loaded in later.
ratios_matrix_WT = ratios_matrix;
approach_ratio_matrix_WT = approach_ratio_matrix;
termination_ratio_matrix_WT = termination_ratio_matrix;
aftermath_ratio_matrix_WT = aftermath_ratio_matrix;
Removetopbottom10pclogical_WT = Removetopbottom10pclogical;
Removetopbottom10pclogical2_WT = Removetopbottom10pclogical2;
Removetopbottom10pclogical3_WT = Removetopbottom10pclogical3;


%load in the WT clusters
classes2clust = zeros(size(approach_ratio_matrix));
classes3clust = zeros(size(approach_ratio_matrix));
for K = 2:3
    
    loadfile2 = char(append('NETseqWT',string(K), 'clusters1', 'clustering.mat'));
    
    load([dataloaddirr loadfile2])
    if K == 2 % added if clauses so that loading in the next dataset will not remove the previous.
        classes2clust = classes ;
    end
    
    if K == 3
        classes3clust = classes ;
        
    end
end


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
% mutant = 'Glugal15'
% datatype1 = 'NETseq'
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
% mutant = 'Glugal60'
% datatype1 = 'NETseq'
%% Start with WT PROseq

% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/LIS_PRO_SEQ/LisLab_PROseq_data/';
%
% % NETseq_sample1 = importdata([dataloaddir 'Lis_PROseq_WT.mat']);
% % NETseq_sample2 = importdata([dataloaddir 'Lis_PROseq_WT.mat']); %duplicated the same dataset because shouldn't make a difference and can use same code that I have edited.
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
% mutant = 'WT'
% datatype1 = 'PROseq'
% % mutants to load can be found in andrew's code
%% spt4 delete
% %
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/LIS_PRO_SEQ/LisLab_PROseq_data/'
% %
%  NETseq_sample1 = importdata([dataloaddir 'LIS_PROseq_spt4del_amended.mat']);
%  NETseq_sample2 = importdata([dataloaddir 'LIS_PROseq_spt4del_amended.mat']);
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
% mutant = 'spt4del'
% datatype1 = 'PROseq'
%% Start with WT
%
% dataloaddir = '/data/FIXED_DATA/FIXED_TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/notag_normalised/WT_nnt/';
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
% xafter = 300;
% ybefore = 300;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant = 'WT'
% datatype1 = 'NETseq'
% % mutants to load can be found in andrew's code
%% set1 delete
% dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
%
% NETseq_sample1 = importdata([dataloaddir 'weissman_set1del_NETseq.mat']);
% NETseq_sample2 = importdata([dataloaddir 'weissman_set1del_NETseq.mat']);
%
% disp('Loaded in W set1del NETseq data')
%
% xafter = 200;
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
% xafter = 200;
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
% ybefore = 300;
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
% xafter = 200;
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
% xafter = 200;
% ybefore = 100;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant = 'eaf3del'
% datatype1 = 'NETseq'
%% spt4 delete

dataloaddir = '/data/FIXED_DATA/FIXED_TRANSCRIPTION_DATA/ULKU_NET_SEQ/BigWiG_spikein_normalised2019_NETandTEF-seq/notag_normalised/spt4del_nnt/';

NETseq_sample1 = importdata([dataloaddir 'ulku_spt4del_NETseq_1.mat']);
NETseq_sample2 = importdata([dataloaddir 'ulku_spt4del_NETseq_2.mat']);

disp('Loaded in UU spt4 delete NETseq data')

% NOTES:
% The correlation between the two samples was very good
% The average levels in the genes lie significantly off the diagonal - this
% might indicate a levelling issue and a more careful average might be
% needed in the longer run.
xafter = 300;
ybefore = 300;
length_pulled_out = xafter+ybefore;
all_genes = zeros((xafter+ybefore),1) ;
mutant = 'spt4del'
datatype1 = 'NETseq'
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
% ybefore = 300;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant = 'dst1del'
% datatype1 = 'NETseq'
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
% ybefore = 300;
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
% ybefore = 300;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant = 'spt4AA'
% datatype1 = 'NETseq'

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
% mutant = 'spt5AA'
% datatype1 = 'NETseq'

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
% ybefore = 300;
% length_pulled_out = xafter+ybefore;
% all_genes = zeros((xafter+ybefore),1) ;
% mutant = 'spt45DMSO'
% datatype1 = 'NETseq'


%% load in mutant ratios
dataloadirr2 = '/home/orie3770/Current_Cluster_Job91/'%change X to mutant of interest ratios.


loadfile3 = append(datatype1,mutant,'ratios.mat');
load([dataloadirr2, loadfile3])
%%
for clctr = 2:3
    %%
    repeat = 1;
    %will need to replace ratio matrices with avg_ratio_matrix - its alright saved everything
    % PrePolyA window
    PrePolyA_ratio_ranksums = zeros(2, clctr);
    colnames = strings(clctr,1);
    w = 0;
    %%
    
    if clctr == 2
        [p, h] = ranksum(approach_ratio_matrix_WT(classes2clust==(clctr-1)& Removetopbottom10pclogical_WT == 1 ), approach_ratio_matrix(classes2clust == (clctr-1) & Removetopbottom10pclogical == 1 ));
        w=w+1;
        PrePolyA_ratio_ranksums(1, w) = p; %might need to switch around orientation of these to get nice table.
        PrePolyA_ratio_ranksums(2, w) = h;
        colnames(clctr-1,1) = append('WT',string(clctr-1),'vs ', mutant, string(clctr-1));
        
        
        [p, h] = ranksum(approach_ratio_matrix_WT(classes2clust==(clctr)& Removetopbottom10pclogical_WT == 1 ), approach_ratio_matrix(classes2clust == (clctr) & Removetopbottom10pclogical == 1 ));
        w=w+1;
        PrePolyA_ratio_ranksums(1, w) = p; %might need to switch around orientation of these to get nice table.
        PrePolyA_ratio_ranksums(2, w) = h;
        colnames(clctr,1) = append('WT',string(clctr),'vs ', mutant, string(clctr));
    end
    % need separate if clause for clctr == 3 because need to change
    % clustering used.
    %%
    if clctr == 3
        [p, h] = ranksum(approach_ratio_matrix_WT(classes3clust==(clctr-1)& Removetopbottom10pclogical_WT == 1 ), approach_ratio_matrix(classes3clust == (clctr-1) & Removetopbottom10pclogical == 1 ));
        w=w+1;
        PrePolyA_ratio_ranksums(1, w) = p; %might need to switch around orientation of these to get nice table.
        PrePolyA_ratio_ranksums(2, w) = h;
        colnames(clctr-1,1) = append('WT',string(clctr-1),'vs ', mutant, string(clctr-1));
        
        
        [p, h] = ranksum(approach_ratio_matrix_WT(classes3clust==(clctr)& Removetopbottom10pclogical_WT == 1 ), approach_ratio_matrix(classes3clust == (clctr) & Removetopbottom10pclogical == 1 ));
        w=w+1;
        PrePolyA_ratio_ranksums(1, w) = p; %might need to switch around orientation of these to get nice table.
        PrePolyA_ratio_ranksums(2, w) = h;
        colnames(clctr,1) = append('WT',string(clctr-1),'vs ', mutant, string(clctr-1));
        
        if clctr ==3
            [q, i] = ranksum(approach_ratio_matrix_WT(classes3clust==(clctr-2) & Removetopbottom10pclogical_WT == 1), approach_ratio_matrix(classes3clust == (clctr-2) & Removetopbottom10pclogical == 1));
            w=w+1;
            PrePolyA_ratio_ranksums(1,w) = q ;
            PrePolyA_ratio_ranksums(2,w) = i;
            colnames(clctr-2,1) = append('WT',string(clctr-2),' vs ', mutant, string(clctr-2));
        end
    end
    
    %% make table instead of text file
    rownames = ["p value", "result"];
    
    f = figure ;
    z = uitable(f, 'Data', PrePolyA_ratio_ranksums, 'ColumnName', colnames, 'RowName', rownames)
    z.Position(3:4) = z.Extent(3:4);
    
    
    
    filename = append(mutant, datatype1,string(clctr),'WTclusters',string(repeat));
    %% save table and textfile
    filenamew = append(filename, 'PrePolyA_ranksums');
    saveas(gcf, filenamew, 'fig')
    saveas(gcf,filenamew, 'jpeg')
    
    %save as text file as before as well no headers though.
    directory = '/home/orie3770/';
    fid = fopen(append(directory, mutant, datatype1, string(clctr),'clusters',string(repeat),'repeat','PrePolyA_ranksum.txt' ) , 'w' );
    fprintf(fid,'%s %s\n' ,  rownames)
    fprintf(fid, '%f %f\n', PrePolyA_ratio_ranksums );
    fclose(fid);
    
    
    %Post PolyA window
    
    PostPolyA_ratio_ranksums = zeros(2, clctr);
    colnames = strings(clctr,1);
    w = 0;
    
    
    if clctr ==2
        [p, h] = ranksum(termination_ratio_matrix(classes2clust==(clctr-1)& Removetopbottom10pclogical2_WT == 1 ), termination_ratio_matrix(classes2clust == (clctr-1) & Removetopbottom10pclogical2 == 1 ));
        w=w+1;
        PostPolyA_ratio_ranksums(1, w) = p; %might need to switch around orientation of these to get nice table.
        PostPolyA_ratio_ranksums(2, w) = h;
        colnames(clctr-1,1) = append('WT',string(clctr-1),' vs ', mutant, string(clctr));
        
        [p, h] = ranksum(termination_ratio_matrix(classes2clust==(clctr)& Removetopbottom10pclogical2_WT == 1 ), termination_ratio_matrix(classes2clust == clctr & Removetopbottom10pclogical2 == 1 ));
        w=w+1;
        PostPolyA_ratio_ranksums(1, w) = p; %might need to switch around orientation of these to get nice table.
        PostPolyA_ratio_ranksums(2, w) = h;
        colnames(clctr,1) = append('WT',string(clctr),' vs ', mutant, string(clctr));
    end
    
    if clctr ==3
        [p, h] = ranksum(termination_ratio_matrix(classes3clust==(clctr-1)& Removetopbottom10pclogical2_WT == 1 ), termination_ratio_matrix(classes3clust == (clctr-1) & Removetopbottom10pclogical2 == 1 ));
        w=w+1;
        PostPolyA_ratio_ranksums(1, w) = p; %might need to switch around orientation of these to get nice table.
        PostPolyA_ratio_ranksums(2, w) = h;
        colnames(clctr-1,w) = append('WT',string(clctr-1),' vs ', mutant, string(clctr));
        
        [p, h] = ranksum(termination_ratio_matrix(classes3clust==(clctr)& Removetopbottom10pclogical2_WT == 1 ), termination_ratio_matrix(classes3clust == clctr & Removetopbottom10pclogical2 == 1 ));
        w=w+1;
        PostPolyA_ratio_ranksums(1, w) = p; %might need to switch around orientation of these to get nice table.
        PostPolyA_ratio_ranksums(2, w) = h;
        colnames(clctr,w) = append('WT',string(clctr),' vs ', mutant, string(clctr));
        if clctr ==3
            [q, i] = ranksum(termination_ratio_matrix(classes3clust==(clctr-2) & Removetopbottom10pclogical2 == 1), termination_ratio_matrix(classes3clust == (clctr-2) & Removetopbottom10pclogical2 == 1));
            w=w+1;
            PostPolyA_ratio_ranksums(1,w) = q ;
            PostPolyA_ratio_ranksums(2,w) = i;
            colnames(clctr-2,w) = append('WT',string(clctr-2),' vs ', mutant, string(clctr-2));
        end
    end
    
    %% make table instead of text file
    rownames = ["p value", "result"];
    
    f = figure ;
    z = uitable(f, 'Data', PostPolyA_ratio_ranksums, 'ColumnName', colnames, 'RowName', rownames)
    z.Position(3:4) = z.Extent(3:4);
    %% save table and textfile
    filenamew = append(filename, 'PostPolyA_ranksums');
    saveas(gcf, filenamew, 'fig')
    saveas(gcf,filenamew, 'jpeg')
    
    %save as text file as before as well no headers though.
    directory = '/home/orie3770/';
    fid = fopen(append(directory, mutant, datatype1, string(clctr),'clusters',string(repeat),'repeat','PostPolyA_ranksum.txt' ) , 'w' );
    fprintf(fid,'%s %s\n' ,  rownames)
    fprintf(fid, '%f %f\n', PostPolyA_ratio_ranksums );
    fclose(fid);
    
    % Post Termination window
    
    PostTerm_ratio_ranksums = zeros(2, clctr);
    colnames = strings(clctr,1);
    w = 0;
    
    if clctr ==2
        [p, h] = ranksum(aftermath_ratio_matrix(classes2clust==(clctr-1)& Removetopbottom10pclogical3_WT == 1 ), aftermath_ratio_matrix(classes2clust == (clctr-1) & Removetopbottom10pclogical3 == 1 ));
        w=w+1;
        PostTerm_ratio_ranksums(1, w) = p; %might need to switch around orientation of these to get nice table.
        PostTerm_ratio_ranksums(2, w) = h;
        colnames(clctr-1,1) = append('WT',string(clctr-1),' vs ', mutant, string(clctr-1));
        
        
        [p, h] = ranksum(aftermath_ratio_matrix(classes2clust==(clctr)& Removetopbottom10pclogical3_WT == 1 ), aftermath_ratio_matrix(classes2clust == (clctr) & Removetopbottom10pclogical3 == 1 ));
        w=w+1;
        PostTerm_ratio_ranksums(1, w) = p; %might need to switch around orientation of these to get nice table.
        PostTerm_ratio_ranksums(2, w) = h;
        colnames(clctr,1) = append('WT',string(clctr),' vs ', mutant, string(clctr));
    end
    
    if clctr ==3
        [p, h] = ranksum(aftermath_ratio_matrix(classes2clust==(clctr-1)& Removetopbottom10pclogical3_WT == 1 ), aftermath_ratio_matrix(classes2clust == (clctr-1) & Removetopbottom10pclogical3 == 1 ));
        w=w+1;
        PostTerm_ratio_ranksums(1, w) = p; %might need to switch around orientation of these to get nice table.
        PostTerm_ratio_ranksums(2, w) = h;
        colnames(clctr-1,1) = append('WT',string(clctr-1),' vs ', mutant, string(clctr-1));
        
        
        [p, h] = ranksum(aftermath_ratio_matrix(classes2clust==(clctr)& Removetopbottom10pclogical3_WT == 1 ), aftermath_ratio_matrix(classes2clust == (clctr) & Removetopbottom10pclogical3 == 1 ));
        w=w+1;
        PostTerm_ratio_ranksums(1, w) = p; %might need to switch around orientation of these to get nice table.
        PostTerm_ratio_ranksums(2, w) = h;
        colnames(clctr,1) = append('WT',string(clctr),' vs ', mutant, string(clctr));
        if clctr ==3
            [q, i] = ranksum(aftermath_ratio_matrix(classes3clust==(clctr-2) & Removetopbottom10pclogical3_WT == 1), aftermath_ratio_matrix(classes3clust == (clctr-2) & Removetopbottom10pclogical3 == 1));
            w=w+1;
            PostTerm_ratio_ranksums(1,w) = q ;
            PostTerm_ratio_ranksums(2,w) = i;
            colnames(clctr-2,1) = append('WT',string(clctr-1),' vs ', mutant, string(clctr));
        end
    end
    
    %% make table instead of text file
    rownames = ["p value", "result"];
    
    f = figure ;
    z = uitable(f, 'Data', PostTerm_ratio_ranksums, 'ColumnName', colnames, 'RowName', rownames)
    z.Position(3:4) = z.Extent(3:4);
    %% save table and textfile
    filenamew = append(filename, 'PostTerm_ranksums');
    saveas(gcf, filenamew, 'fig')
    saveas(gcf,filenamew, 'jpeg')
    
    %save as text file as before as well no headers though.
    directory = '/home/orie3770/';
    fid = fopen(append(directory, mutant, datatype1, string(clctr),'clusters',string(repeat),'repeat','PostTerm_ranksum.txt' ) , 'w' );
    fprintf(fid, '%s %s\n' , rownames)
    fprintf(fid, '%f %f\n', PostTerm_ratio_ranksums );
    fclose(fid);
    
    
    
    
    
end
%% Make place to put all the ranksums
unix('mkdir Mutant_Ratio_Wilcoxon_Ranksums')
unix('cd Mutant_Ratio_Wilcoxon_Ranksums')
unix(append('mkdir ', mutant, datatype1))
unix('cd ..')

unix(append('mv *ranksum* ./Mutant_Ratio_Wilcoxon_Ranksums/', mutant, datatype1))


