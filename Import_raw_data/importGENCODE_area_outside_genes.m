%% Import gene lists from the GENCODE data

%Getting this right will allow me to update everything when new 
%annotations are released

%%

dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/GENCODE/';

%Start with the basic one - I really only need gene starts and ends right
%now

GENCODE_raw = importdata([dataloaddir 'gencode.v33.basic.annotation.gtf']);

%%

GENCODE_start = 6;
%as the first few lines contain header information
%there doesn't seem to be a footer - could change for future sets so check
GENCODE_raw = {GENCODE_raw{GENCODE_start:length(GENCODE_raw)}};
GENCODE_raw = GENCODE_raw';

%% Select just the genes

GENCODE_genes_logical = zeros(length(GENCODE_raw),1);
GENCODE_chrs = zeros(length(GENCODE_raw),1);
GENCODE_firsts = zeros(length(GENCODE_raw),1);
GENCODE_lasts = zeros(length(GENCODE_raw),1);
GENCODE_strands = zeros(length(GENCODE_raw),1);
GENCODE_source = zeros(length(GENCODE_raw),1);
GENCODE_six = zeros(length(GENCODE_raw),1);
GENCODE_eight = zeros(length(GENCODE_raw),1);

GENCODE_names = cell(length(GENCODE_raw),1);

%% Chromosome Loop

formatspec = '%s %s %s %s %s %s %s %s %s';



for ictr=1:1:length(GENCODE_raw) %short start for testing!

    temptextscan = textscan(GENCODE_raw{ictr},formatspec,'Delimiter','\t');
    %1 = chromosome
    %2 = source of annotation (HAVANA is the manual team)
    %3 = gene (or something else)
    %4 = one gene end
    %5 = other gene end
    %6 = ? (it's a . in the early entries)
    %7 = plus or minus strand?
    %8 = ? (it's a . in the early entries)
    

    switch temptextscan{1}{1}
        
        case 'chr1'
            GENCODE_chrs(ictr) = 1;
        case 'chr2'
            GENCODE_chrs(ictr) = 2;
        case 'chr3'
            GENCODE_chrs(ictr) = 3;
        case 'chr4'
            GENCODE_chrs(ictr) = 4;
        case 'chr5'
            GENCODE_chrs(ictr) = 5;
        case 'chr6'
            GENCODE_chrs(ictr) = 6;
        case 'chr7'
            GENCODE_chrs(ictr) = 7;
        case 'chr8'
            GENCODE_chrs(ictr) = 8;
        case 'chr9'
            GENCODE_chrs(ictr) = 9;
        case 'chr10'
            GENCODE_chrs(ictr) = 10;
        case 'chr11'
            GENCODE_chrs(ictr) = 11;
        case 'chr12'
            GENCODE_chrs(ictr) = 12;
        case 'chr13'
            GENCODE_chrs(ictr) = 13;
        case 'chr14'
            GENCODE_chrs(ictr) = 14;
        case 'chr15'
            GENCODE_chrs(ictr) = 15;
        case 'chr16'
            GENCODE_chrs(ictr) = 16;
        case 'chr17'
            GENCODE_chrs(ictr) = 17;
        case 'chr18'
            GENCODE_chrs(ictr) = 18;
        case 'chr19'
            GENCODE_chrs(ictr) = 19;
        case 'chr20'
            GENCODE_chrs(ictr) = 20;
        case 'chr21'
            GENCODE_chrs(ictr) = 21;
        case 'chr22'
            GENCODE_chrs(ictr) = 22;
        case 'chrX'
            GENCODE_chrs(ictr) = 23;
        case 'chrY'
            GENCODE_chrs(ictr) = 24;
        case 'chrM'
            GENCODE_chrs(ictr) = -1;
        otherwise
            GENCODE_chrs(ictr) = -2;
    end
    
    if strcmp(temptextscan{3},'gene')
        GENCODE_genes_logical(ictr) = 1;
    end
    
    GENCODE_firsts(ictr) = str2double(temptextscan{4});
    GENCODE_lasts(ictr) = str2double(temptextscan{5});
    
    if strcmp(temptextscan{7},'+')
        GENCODE_strands(ictr) = 1;
    elseif strcmp(temptextscan{7},'-')
        GENCODE_strands(ictr) = 0;
    else
        GENCODE_strands(ictr) = -1;
    end
    %put the -1 case in for subsequent testing
    
    switch temptextscan{2}{1}
        
        case 'HAVANA'
            GENCODE_source(ictr) = 1;
        case 'ENSEMBL'
            GENCODE_source(ictr) = 0;
        otherwise
            GENCODE_source(ictr) = -1;
    end
    
    
    
end


%% Notes

%It seems that lasts is always greater than firsts in location
%Which makes it difficult to tell exactly how the strands are oriented
%may have to do some flipping when pulling out the relevant NETseq stuff
%here!

%It seems that some genes come from a manual source (HAVANA) and the others
%come from an automated source so should I be using only the manual ones?


%It might be possible to extract some further information out of the last
%string

%% Quick Test of textscan

%might be better to use textscan than import data

%formatspec = '%s %s %s %s %s %s %s %s %s';
%temptextscan = textscan(GENCODE_raw{1},formatspec,'Delimiter','\t');

%It is way better than the alternative but can still be used on the
%imported cell array!

%% How many genes come from Havana




%% Reduce this down to a list of genes

gene_positions_for_nogene_mean = cell(24,2);
%24 chromosomes 2 strands 

for cctr = 1:1:24
    
    chr_logical = GENCODE_chrs == cctr;
    
    %+strand

    temp_fwd_locs = ...
        logical(chr_logical.*GENCODE_genes_logical.*(GENCODE_strands==1));
    
    gene_positions_for_nogene_mean{cctr,1} = ...
        [GENCODE_firsts(temp_fwd_locs) , ...
        GENCODE_lasts(temp_fwd_locs)];
    
    %-strand
    
    temp_rev_locs = ...
        logical(chr_logical.*GENCODE_genes_logical.*(GENCODE_strands==0));
    
    gene_positions_for_nogene_mean{cctr,2} = ...
        [GENCODE_firsts(temp_rev_locs) , ...
        GENCODE_lasts(temp_rev_locs)];
    
    
end


%Could remove overlapping genes at this point!

%% Remove overlapping genes from this set

%Are the genes in the order in which they appear on the chromosome?
% for cctr=1:1:24
%     figure
%     
%     templocs = logical((GENCODE_chrs == cctr).*GENCODE_genes_logical);
%     
%     histogram(diff(GENCODE_firsts(templocs)))
%     
% end

%From the above it seems that the genes are reasonably in order but the
%non-gene annotations might not be

% double check with some sums

% for cctr=1:1:24
%    
%     templocs = logical((GENCODE_chrs == cctr).*GENCODE_genes_logical);
%     
%     sum(diff(GENCODE_firsts(templocs))<0)
%     
% end

%The above should also be indpendent of strand so it seems that the genes
%are ordered by where they appear and not strand in the original set - 
%seems logical enough

% Should I later limit to protein-conding TUs?



    
%Assume that the genes are in order of their 'first' positions
%so only need to check each gene for if it overlaps with the next one
%actually this means it is ok only to check for the end of one gene
%overlapping with the next ones
%but what about genes within genes? - I don't think Tom removed these but I
%think I might!
%Copying the logic from Tom's code remove_overlapping.m
%Although Tom seems to check for all genes above - perhaps this is
%necessary for cases in which more than two genes overlap on the same 
%strand?
%Tom also seems to be checking for genes that overlap on the other strand
% - should I also remove these?


% remove_genes = cell(24,2);
% for cctr = 1:1:24
%    for sctr = 1:2
%       remove_genes{cctr,sctr} = zeros(size(gene_positions{cctr,sctr},1),1);
%    end
% end
% 
% 
% for cctr = 1:1:24
%     
%     strand = 1; %forward
%     
%     for gctr=1:size(gene_positions{cctr,strand})
%         for ictr = gctr+1:size(gene_positions{cctr,strand},1)
%             if( (gene_positions{cctr,strand}(gctr,2) >= gene_positions{cctr,strand}(ictr,1)) && ...
%                     (gene_positions{cctr,strand}(gctr,1) <= gene_positions{cctr,strand}(ictr,2)) )
%                 remove_genes{cctr,strand}(gctr) = 1;
%                 remove_genes{cctr,strand}(ictr) = 1;
%             end
%         end
%         %Also check for overlaps with the opposing strand here
%         for ictr = 1:size(gene_positions{cctr,3-strand},1)
%             if( (gene_positions{cctr,strand}(gctr,2) >= gene_positions{cctr,3-strand}(ictr,1)) && ...
%                     (gene_positions{cctr,strand}(gctr,1) <= gene_positions{cctr,3-strand}(ictr,2)) )
%                 remove_genes{cctr,strand}(gctr) = 1;
%                 remove_genes{cctr,3-strand}(ictr) = 1;
%             end
%         end
%         %Here I differ from Tom in that the gencode that I downloaded has
%         %the 'leftmost' and 'rightmost' points of genes not the start and
%         %end - so there is no 'flipping' to be done in the comparison
%         
%     end
%     
%     
%     strand = 2; %reverse
%     for gctr=1:size(gene_positions{cctr,strand})
%         for ictr = gctr+1:size(gene_positions{cctr,strand},1)
%             if( (gene_positions{cctr,strand}(gctr,2) >= gene_positions{cctr,strand}(ictr,1)) && ...
%                     (gene_positions{cctr,strand}(gctr,1) <= gene_positions{cctr,strand}(ictr,2)) )
%                 remove_genes{cctr,strand}(gctr) = 1;
%                 remove_genes{cctr,strand}(ictr) = 1;
%             end
%         end
%         %Also check for overlaps with the opposing strand here
%         for ictr = 1:size(gene_positions{cctr,3-strand},1)
%             if( (gene_positions{cctr,strand}(gctr,2) >= gene_positions{cctr,3-strand}(ictr,1)) && ...
%                     (gene_positions{cctr,strand}(gctr,1) <= gene_positions{cctr,3-strand}(ictr,2)) )
%                 remove_genes{cctr,strand}(gctr) = 1;
%                 remove_genes{cctr,3-strand}(ictr) = 1;
%             end
%         end
%         %Here I differ from Tom in that the gencode that I downloaded has
%         %the 'leftmost' and 'rightmost' points of genes not the start and
%         %end - so there is no 'flipping' to be done in the comparison
%         
%     end
%     
%     %Should I really be checking for all overlapping things even if they
%     %are not genes?
%     
%    
% end
% 
% % I am also going to go ahead and remove genes that are shorter than 1000
% %bp at this point
% 
% for cctr=1:1:24
%    
%     for sctr=1:1:2
%         
%         for gctr=1:1:size(gene_positions{cctr,sctr},1)
%             
%             if (abs(gene_positions{cctr,sctr}(gctr,2) - gene_positions{cctr,sctr}(gctr,1)) < 1000)
%                 remove_genes{cctr,sctr}(gctr) = 1;
%             end
%             
%         end
%         
%     end
%     
% end
% 
% 
% 
% %Now remove the actual genes
% %NB doing this in a different way than Tom
% %going from right to left such that don't need to check if the structure
% %has gotten too small
% %Also not 100% sure how Tom's algorithm for this works
% for cctr = 1:1:24
%     for sctr=1:2
%        for ictr=size(remove_genes{cctr,sctr},1):-1:1
%            if remove_genes{cctr,sctr}(ictr)==1
%                gene_positions{cctr,sctr}(ictr,:) = [];
%            end
%        end
%     end
% end
% 
% %Currently not tracking gene names - will have to go back to here and also
% %remove gene names or give each gene position some kind of identifier
% 
% %I had missed a loop over genes when removing those that were too small
% %Now fixed
% 
% 
% %I should check at the end that no overlaps remain:
% %% Separate section as this could take quite a while
% 
% %determine the largest location for each chromosome
% chr_rough_sizes = zeros(24,1);
% rough_chromosomes = cell(24,1);
% 
% for cctr=1:1:24
% 
%     chr_rough_sizes(cctr) = ...
%         max( ...
%         [max(max(gene_positions{cctr,1})) ...
%         max(max(gene_positions{cctr,2}))] );
%     
%     rough_chromosomes{cctr} = zeros(chr_rough_sizes(cctr),1);
%     
% end
% 
% %I think the matrix I was envisioning for this will be simply too big to 
% %store in memory - actually it might just work I made one chromosome and it
% %worked
% 
% for cctr = 1:1:24
%     
%     for sctr=1:1:2
%        
%         for gctr = 1:1:size(gene_positions{cctr,sctr},1)
%         
%             rough_chromosomes{cctr}( ...
%                 gene_positions{cctr,sctr}(gctr,1):gene_positions{cctr,sctr}(gctr,2) ) = ...
%                 rough_chromosomes{cctr}( ...
%                 gene_positions{cctr,sctr}(gctr,1):gene_positions{cctr,sctr}(gctr,2) ) + 1;
%             
%         end
%         
%     end
%     
% end
% 
% overlaps = zeros(24,1);
% 
% for ictr=1:1:24
%    
%     overlaps(ictr) = sum(rough_chromosomes{ictr}>1);
%     
% end
% 
% %Seems to have worked - there are no detected overlaps
% 
% %I could also remove genes that are considered to be too small here!
% 
% %% Quick comparison with previously saved GENCODE stuff
% 
% previous_annotations_dir = ...
%     '/home/bioc1242/R_ANGEL/TOM_BROWN/Elongation_clustering/HeLa_analysis/gene_annotations/';
% 
% Previous_GENCODE = ...
%     load([previous_annotations_dir 'HeLa_gencodeAnnotations.mat']);
% 
% %Structure of this is a little different and I am not sure how to compare
% %the two
% 
%% Save the imported and processed gene positions

save([dataloaddir 'GENCODE_basic_processed_positions_for_nogene_mean.mat'],'gene_positions_for_nogene_mean')

% %Consider saving other things as they might become useful later on
% %but just use this for now
% 
% 
% %% Getting gene names also
% 
% %It turns out that I will need to extract the gene names to identify
% %things in some of the analysis - develop that here!
% 
% 
% %%
% 
% nameformatspec = '';
% 
% for ictr=1:1:max(semicolonno)
%    
%     nameformatspec = [nameformatspec '%s '];
%     
% end
% 
% wrapper = @(x) strfind(x,'gene_name');
% 
% tempnamematches = cellfun(wrapper,tempnamescan{1},'UniformOutput',false);
% 
% %this works but now to extract the one with the success (if there is one)
% 
% wrapper2 = @(x) ~isempty(x);
% 
% tempnamelocs = cellfun(wrapper2,tempnamematches);
% 
% tempgenename = tempnamescan{1}{tempnamelocs}(12:end-1);
% 
% 
% %% Extract all gene names
% 
% %It might also be useful to get the ids (in fact this might well be more
% %important than the names for gene ontology
% 
% %holder for the gene names
% GENCODE_names = cell(length(GENCODE_raw),1);
% %holder for the gene ids
% GENCODE_ids = cell(length(GENCODE_raw),1);
% 
% %set up some compound functions for use with cellfun
% wrapper = @(x) strfind(x,'gene_name');
% wrapper2 = @(x) ~isempty(x);
% wrapper3 = @(x) strfind(x,'gene_id');
% 
% for ictr=1:1:length(GENCODE_raw)
%    
%     temptextscan = textscan(GENCODE_raw{ictr},formatspec,'Delimiter','\t');
%     
%     tempnamescan = textscan(temptextscan{9}{1},nameformatspec, ...
%         'CollectOutput',true,'Delimiter',';');
%     
%     tempnamematches = cellfun(wrapper,tempnamescan{1},'UniformOutput',false);
%     
%     tempnamelocs = cellfun(wrapper2,tempnamematches);
%     
%     tempgenename = tempnamescan{1}{tempnamelocs}(12:end-1);
%     
%     GENCODE_names{ictr} = tempgenename;
%     
%     tempidmatches = cellfun(wrapper3,tempnamescan{1},'UniformOutput',false);
%     
%     tempidlocs = cellfun(wrapper2,tempidmatches);
%     
%     tempgeneid = tempnamescan{1}{tempidlocs}(10:end-1);
%     
%     GENCODE_ids{ictr} = tempgeneid;
%     
% end
% 
% 
% %% reduce the gene names and ids down as for the gene positions
% 
% gene_names = cell(24,2);
% gene_ids = cell(24,2);
% 
% for cctr = 1:1:24
%     
%     chr_logical = GENCODE_chrs == cctr;
%     
%     %+strand
%     
%     temp_fwd_locs = ...
%         logical(chr_logical.*GENCODE_genes_logical.*(GENCODE_strands==1));
%     
%     gene_names{cctr,1} = ...
%         GENCODE_names(temp_fwd_locs);
% 
%     gene_ids{cctr,1} = ...
%         GENCODE_ids(temp_fwd_locs);
%     
%     
%     %-strand
%     
%     temp_rev_locs = ...
%         logical(chr_logical.*GENCODE_genes_logical.*(GENCODE_strands==0));
%     
%     gene_names{cctr,2} = ...
%         GENCODE_names(temp_rev_locs);
%     
%     gene_ids{cctr,2} = ...
%         GENCODE_ids(temp_rev_locs);
%     
% end
% 
% %% Further reduce ids and names by removal of non-conformers
% 
% for cctr = 1:1:24
%     
%     for sctr=1:2
%         
%         for ictr=size(remove_genes{cctr,sctr},1):-1:1
%             if remove_genes{cctr,sctr}(ictr)==1
%                 gene_names{cctr,sctr}(ictr,:) = [];
%                 gene_ids{cctr,sctr}(ictr,:) = [];
%             end
%         end
%         
%     end
%         
% end

%This appears to have worked but it was suspiciously painless



%% Save the imported and processed gene names and ids

% save([dataloaddir 'GENCODE_basic_processed_names.mat'],'gene_names')
% save([dataloaddir 'GENCODE_basic_processed_ids.mat'],'gene_ids')


%Consider saving other things as they might become useful later on
%but just use this for now