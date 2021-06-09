%% Import gene lists from the Ensembl Cerevisiae annotations

%data seems to be largely in the same format as the GENCODE - should
%double check that there is the same number of columns
%Also might be able to do things in a slightly better/faster way than
%before - with textscan and explicit acceptance of comment lines (

%NB this script does some things more efficiently than previous imports

tic


%% Load in the data as raw


dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/ENSEMBL/S_CEREVISIAE_ANNOTATIONS/';

dataloadfile = 'Saccharomyces_cerevisiae.R64-1-1.101.gtf';

ENSEMBL_raw = importdata([dataloaddir dataloadfile]);


%% Test for the character(s) used to denote comments


commlocs = strfind(ENSEMBL_raw,'#');

commno = sum(~cellfun(@isempty,commlocs));

%this is the right character (or at least it does not appear other than in
%the first five lines of the file - means it is probably OK to use textscan
%to read in the file in a way which will require less processing down the
%line

%% Also test for the number of columns

tablocs = strfind(ENSEMBL_raw,sprintf('\t'));

tabno = unique(cellfun(@numel,tablocs));

%it seeems that there are only 0 and 8 tabs in each of the lines read in - 
%this means it is safe to assume the same number of columns for all 
%non-comment lines



%% Open the raw data

formatspec = '%s %s %s %f %f %s %s %s %s';


fileID1 = fopen([dataloaddir dataloadfile]);

ensembl_sc_readin = textscan(fileID1,formatspec, ...
    'Delimiter','\t', ...
    'CommentStyle','#');

fclose(fileID1);
whos ensembl_sc_readin

%% Test for the contents of each column

col1cont = unique(ensembl_sc_readin{1}(:));
%16 chromosomes (+mitochondrial) in roman numerals

col2cont = unique(ensembl_sc_readin{2}(:));
%sgd - source of the annotation (appears to be the only source)

col3cont = unique(ensembl_sc_readin{3}(:));
%CDS, exon, five_prime_utr, gene, start_codon, stop_codon, transcript
% - I think I should just concentrate on 'gene' to start with

%4 and 5 contain the start and end locations - I should check the 
%orientation

col_4_5_orientation = sum((ensembl_sc_readin{4}-ensembl_sc_readin{5}) > 0);
%always 5>4 (except 4 lines when they are the same)


col6cont = unique(ensembl_sc_readin{6}(:));
%always '.'

col7cont = unique(ensembl_sc_readin{7}(:));
% '+' or '-' (I guess this is the forward or reverse strand)

col8cont = unique(ensembl_sc_readin{8}(:));
% '.', '0', '1' or '2'

col9cont = unique(ensembl_sc_readin{9}(:));
%~35000 distinct entries here - contains a lot of peripheral information
% most relevant is 'gene_id' at the start which directly precedes the 
%gene id (which I do need to extract)
%should later think about extracting which (if any) genes have a 'gene
%biotype' that is 'noncoding'



%%

chromosome_labels = ...
    {'I','II','III','IV','V','VI','VII','VIII','IX','X', ...
    'XI','XII','XIII','XIV','XV','XVI'};

NoChromosomes = 16;

%minGeneSize = 500; %GM changed min gene size as don't need large genes to visualise termination event and can get more signal


%% Reduce the data down to a list of genes



gene_positions_for_nogene_mean = cell(NoChromosomes,2);
%16 chromosomes, 2strands

%Currently ignoring the mitochondrial stuff - consider looking at it later

gene_logical = strcmp(ensembl_sc_readin{3},'gene');

fwd_logical = strcmp(ensembl_sc_readin{7},'+');

rev_logical = strcmp(ensembl_sc_readin{7},'-');

for cctr = 1:1:NoChromosomes
    
    %chr_logical = strfind(ensembl_sc_readin{1},chromosome_labels{cctr});
    %above didn't work
    chr_logical = strcmp(ensembl_sc_readin{1},chromosome_labels{cctr});
    
    %+strand
    
    temp_fwd_locs = logical(chr_logical.*gene_logical.*fwd_logical);
    
    gene_positions_for_nogene_mean{cctr,1} = ...
        [ensembl_sc_readin{4}(temp_fwd_locs), ...
        ensembl_sc_readin{5}(temp_fwd_locs)];
    
    
    %-strand
    
    temp_rev_locs = logical(chr_logical.*gene_logical.*rev_logical);
    
    gene_positions_for_nogene_mean{cctr,2} = ...
        [ensembl_sc_readin{4}(temp_rev_locs), ...
        ensembl_sc_readin{5}(temp_rev_locs)];
    
    
    
end


%GM - Above is all I needed to pull out gene positions for the nogene mean, have
%commented out rest of script rather than delete for the unlikely scenario
%where everything is lost bar this file.


%could remove overlapping genes at this point (and ones that are too short)



% %% Remove overlapping genes and ones < minGeneSize bp - (normally 1000bp)
% %or a different length if set earlier
% 
% 
% %Are the genes in the order in which they appear on the chromosome
% 
% % for cctr=1:1:NoChromosomes
% %     figure
% %     
% %     templocs = ...
% %         logical( ...
% %         strcmp(ensembl_sc_readin{1},chromosome_labels{cctr}).* ...
% %         strcmp(ensembl_sc_readin{3},'gene') ...
% %         );
% %     
% %     histogram(diff(ensembl_sc_readin{4}(templocs)))
% %    
% %     
% % end
% 
% % double check with some sums
% 
% % for cctr=1:1:NoChromosomes
% % 
% %     templocs = ...
% %         logical( ...
% %         strcmp(ensembl_sc_readin{1},chromosome_labels{cctr}).* ...
% %         strcmp(ensembl_sc_readin{3},'gene') ...
% %         );
% % 
% %     sum(diff(ensembl_sc_readin{4}(templocs))<0)
% % 
% % end
% 
% %The above checks worked out
% 
% %There are some additional checks that could be done here
% 
% 
% 
% %Assume that the genes are in order of their 'first' positions
% %so only need to check each gene for if it overlaps with the next one
% %actually this means it is ok only to check for the end of one gene
% %overlapping with the next ones
% %but what about genes within genes? - I don't think Tom removed these but I
% %think I might!
% %Copying the logic from Tom's code remove_overlapping.m
% %Although Tom seems to check for all genes above - perhaps this is
% %necessary for cases in which more than two genes overlap on the same 
% %strand?
% %Tom also seems to be checking for genes that overlap on the other strand
% % - should I also remove these?
% 
% 
% remove_genes = cell(NoChromosomes,2);
% for cctr=1:1:NoChromosomes
%     for sctr=1:2
%         remove_genes{cctr,sctr} = zeros(size(gene_positions_for_nogene_mean{cctr,sctr},1),1);
%     end
% end
% 
% 
% for cctr=1:1:NoChromosomes
%    
%     strand = 1; %forward
%     
%     for gctr=1:1:size(gene_positions_for_nogene_mean{cctr,strand},1)
%        
%         for ictr = (gctr+1):size(gene_positions_for_nogene_mean{cctr,strand},1)
%             if( (gene_positions_for_nogene_mean{cctr,strand}(gctr,2) >= gene_positions_for_nogene_mean{cctr,strand}(ictr,1)-500) && ...%GM added -500 to remove genes with another gene within 500bp from end 
%                     (gene_positions_for_nogene_mean{cctr,strand}(gctr,1) <= gene_positions_for_nogene_mean{cctr,strand}(ictr,2)) )
%                 remove_genes{cctr,strand}(gctr) = 1;
%                 remove_genes{cctr,strand}(ictr) = 1;
%             end            
%         end
%         %Also check for overlaps with the opposing strand here
%         for ictr = 1:size(gene_positions_for_nogene_mean{cctr,3-strand},1)
%            if( (gene_positions_for_nogene_mean{cctr,strand}(gctr,2) >= gene_positions_for_nogene_mean{cctr,3-strand}(ictr,1)) && ...
%                    (gene_positions_for_nogene_mean{cctr,strand}(gctr,1) <= gene_positions_for_nogene_mean{cctr,3-strand}(ictr,2)) )
%                remove_genes{cctr,strand}(gctr) = 1;
%                remove_genes{cctr,3-strand}(ictr) = 1;
%            end
%         end
%         
%     end
%     
%     
%     strand = 2; %reverse
%     
%     for gctr=1:1:size(gene_positions_for_nogene_mean{cctr,strand},1)
%        
%         for ictr = (gctr+1):size(gene_positions_for_nogene_mean{cctr,strand},1)
%             if( (gene_positions_for_nogene_mean{cctr,strand}(gctr,2) >= gene_positions_for_nogene_mean{cctr,strand}(ictr,1)-500) && ...%GM -500 removes genes with another gene close downstream
%                     (gene_positions_for_nogene_mean{cctr,strand}(gctr,1) <= gene_positions_for_nogene_mean{cctr,strand}(ictr,2)) )
%                 remove_genes{cctr,strand}(gctr) = 1;
%                 remove_genes{cctr,strand}(ictr) = 1;
%             end            
%         end
%         %Also check for overlaps with the opposing strand here
%         for ictr = 1:size(gene_positions_for_nogene_mean{cctr,3-strand},1)
%            if( (gene_positions_for_nogene_mean{cctr,strand}(gctr,2) >= gene_positions_for_nogene_mean{cctr,3-strand}(ictr,1)) && ...
%                    (gene_positions_for_nogene_mean{cctr,strand}(gctr,1) <= gene_positions_for_nogene_mean{cctr,3-strand}(ictr,2)) )
%                remove_genes{cctr,strand}(gctr) = 1;
%                remove_genes{cctr,3-strand}(ictr) = 1;
%            end
%         end
%         
%     end
%     
%     
%     %I should consider also checking for overlap with non-coding
%     %transcripts etc.
%     %For the ensembl set it does seem to be detecing some ncRNA and others
%     % so this might be OK
%     
% end
% 
% 
% %Also remove genes that are shorter than the specified minimum size here
% 
% for cctr=1:1:NoChromosomes
%     
%    for sctr=1:2
%        
%        for gctr=1:1:size(gene_positions_for_nogene_mean{cctr,sctr},1)
%        
%            if (abs(gene_positions_for_nogene_mean{cctr,sctr}(gctr,2) - gene_positions_for_nogene_mean{cctr,sctr}(gctr,1)) < minGeneSize)
%                remove_genes{cctr,sctr}(gctr) = 1;
%            end
%            
%        end
%    end
%     
%     
% end
% 
% %Now remove the actual genes
% %NB doing this in a different way than Tom
% %going from right to left such that don't need to check if the structure
% %has gotten too small
% %Also not 100% sure how Tom's algorithm for this works
% for cctr=1:1:NoChromosomes
%     for sctr=1:2
%        for ictr=size(remove_genes{cctr,sctr},1):-1:1
%            if remove_genes{cctr,sctr}(ictr)==1
%                gene_positions_for_nogene_mean{cctr,sctr}(ictr,:) = [];
%            end
%        end
%     end
% end
% 
% %Currently not tracking gene names - will have to go back to here and also
% %remove gene names or give each gene position some kind of identifier
% 
% 
% %% I should check at the end that no overlaps remain:
% 
% % Separate section as this could take quite a while
% 
% %determine the largest location for each chromosome
% chr_rough_sizes = zeros(NoChromosomes,1);
% rough_chromosomes = cell(NoChromosomes,1);
% 
% for cctr=1:1:NoChromosomes
% 
%     chr_rough_sizes(cctr) = ...
%         max( ...
%         [max(max(gene_positions_for_nogene_mean{cctr,1})) ...
%         max(max(gene_positions_for_nogene_mean{cctr,2}))] );
%     
%     rough_chromosomes{cctr} = zeros(chr_rough_sizes(cctr),1);
%     
% end
% 
% %I think the matrix I was envisioning for this will be simply too big to 
% %store in memory - actually it might just work I made one chromosome and it
% %worked
% 
% for cctr = 1:1:NoChromosomes
%     
%     for sctr=1:1:2
%        
%         for gctr = 1:1:size(gene_positions_for_nogene_mean{cctr,sctr},1)
%         
%             rough_chromosomes{cctr}( ...
%                 gene_positions_for_nogene_mean{cctr,sctr}(gctr,1):gene_positions_for_nogene_mean{cctr,sctr}(gctr,2) ) = ...
%                 rough_chromosomes{cctr}( ...
%                 gene_positions_for_nogene_mean{cctr,sctr}(gctr,1):gene_positions_for_nogene_mean{cctr,sctr}(gctr,2) ) + 1;
%             
%         end
%         
%     end
%     
% end
% 
% overlaps = zeros(NoChromosomes,1);
% 
% for ictr=1:1:NoChromosomes
%    
%     overlaps(ictr) = sum(rough_chromosomes{ictr}>1);
%     
% end
% 
% %Seems to have worked - there are no detected overlaps
% 
% %I could also remove genes that are considered to be too small here!
% %Think I already did that above! - limiting to 1000bp
% 
% 
%% Save the imported and processed gene positions

  save([dataloaddir 'gene_positions_for_nogene_mean_Saccharomyces_cerevisiae__R64_1_1__101.mat'],'gene_positions_for_nogene_mean')
%changed filename so won't overide canonical andrew file 
% 
% %% Also extract the gene names to identify things in the analysis
% 
% %% See if all the genes have the same number of description pieces
% 
% semicolonno = zeros(length(ensembl_sc_readin{9}),1);
% 
% for ictr=1:1:length(ensembl_sc_readin{9})
%     
%     semicolonno(ictr) = length(strfind(ensembl_sc_readin{9}{ictr},';'));
%     
% end
% 
% disp(['Unique numbers of delimiters (;) is ' num2str(length(unique(semicolonno)))])
% 
% %look at just the elements identified as genes to see if these have
% %different values.
% 
% disp(['Unique numbers of delimiters (;) for elements identified as genes is ' ...
%     num2str(length(unique(semicolonno(gene_logical))))])
% 
% %This was a worthwhile exercise as some have 3 and some 4 ';' delimiters
% % Need to be careful of this
% 
% %Seems that some of them have a gene id and a separate gene name - 
% %it's probably best to just extract the gene ids for now
% 
% 
% nameformatspec = '';
% 
% for ictr=1:1:max(semicolonno)
%     
%     nameformatspec = [nameformatspec '%s '];
%     
% end
% 
% %will be an extra space at the end of nameformatspec - hopefully this will
% %not affect things
% 
% %% Select out the gene IDs
% 
% %% Initial test
% 
% for ictr=1:1:length(ensembl_sc_readin{9})
%     
%     %temptextscan = textscan(GENCODE_raw{ictr},formatspec,'Delimiter','\t');
%     
%     tempnamescan = textscan(ensembl_sc_readin{9}{ictr},nameformatspec, ...
%         'CollectOutput',true,'Delimiter',';');
%     
% end
% 
% 
% %% test finding the cell with gene_id in it
% 
% wrapper = @(x) strfind(x,'gene_id');
% 
% tempnamematches = cellfun(wrapper,tempnamescan{1},'UniformOutput',false);
% 
% %this works but now to extract the one with the success (if there is one)
% 
% wrapper2 = @(x) ~isempty(x);
% 
% tempnamelocs = cellfun(wrapper2,tempnamematches);
% 
% %tempgenename = tempnamescan{1}{tempnamelocs}(12:end-1);
% 
% tempgeneid = tempnamescan{1}{tempnamelocs}((length('gene_id')+3):end-1);
% 
% 
% %Might not need to do it this way but it works so probably best to just
% %stick with it for now
% 
% 
% %% Extract all gene ids
% 
% 
% ENSEMBL_ids = cell(length(ensembl_sc_readin{9}),1);
% %ENSEMBL_names = cell(length(ensembl_sc_readin{9}),1);
% %perhaps collect names in the future
% 
% wrapper = @(x) strfind(x,'gene_id');
% wrapper2 = @(x) ~isempty(x);
% 
% for ictr=1:1:length(ensembl_sc_readin{9})
%    
%     tempnamescan = textscan(ensembl_sc_readin{9}{ictr},nameformatspec, ...
%         'CollectOutput',true,'Delimiter',';');
%     
%     tempnamematches = cellfun(wrapper,tempnamescan{1},'UniformOutput',false);
%     
%     tempnamelocs = cellfun(wrapper2,tempnamematches);
%     
%     tempgeneid = tempnamescan{1}{tempnamelocs}((length('gene_id')+3):end-1);
%     
%     ENSEMBL_ids{ictr} = tempgeneid;
%     
% end
% 
% 
% %% Reduce the gene ids down as for the gene positions
% 
% gene_ids = cell(NoChromosomes,2);
% 
% for cctr=1:1:NoChromosomes
%    
%     chr_logical = strcmp(ensembl_sc_readin{1},chromosome_labels{cctr});
%     
%     %+strand
%     
%     temp_fwd_locs = logical(chr_logical.*gene_logical.*fwd_logical);
%     
%     gene_ids{cctr,1} = ...
%         ENSEMBL_ids(temp_fwd_locs);
%     
%     
%     %-strand
%     
%     temp_rev_locs = logical(chr_logical.*gene_logical.*rev_logical);
%     
%     gene_ids{cctr,2} = ...
%         ENSEMBL_ids(temp_rev_locs);
%     
% end
% 
% 
% %% Further reduce IDs by removal of non-conformers
% 
% for cctr=1:1:NoChromosomes
%     
%    for sctr=1:2
%       
%        for ictr=size(remove_genes{cctr,sctr},1):-1:1
%            if remove_genes{cctr,sctr}(ictr) == 1
%                gene_ids{cctr,sctr}(ictr,:) = [];
%                %gene_names{cctr,sctr}(ictr,:) = [];
%            end
%        end
%        
%        
%    end
%     
%     
% end
% 
% 
% %% Save the imported and processed gene ids
% 
% % save([dataloaddir 'gene_ids_Saccharomyces_cerevisiae__R64_1_1__101.mat'],'gene_ids')
% 
% 
% %% Also extract the 'gene biotype' information
% 
% 
% ENSEMBL_biotypes = cell(length(ensembl_sc_readin{9}),1);
% %ENSEMBL_names = cell(length(ensembl_sc_readin{9}),1);
% %perhaps collect names in the future
% 
% wrapper3 = @(x) strfind(x,'gene_biotype');
% wrapper4 = @(x) ~isempty(x);
% 
% for ictr=1:1:length(ensembl_sc_readin{9})
%    
%     tempnamescan = textscan(ensembl_sc_readin{9}{ictr},nameformatspec, ...
%         'CollectOutput',true,'Delimiter',';');
%     
%     tempnamematches = cellfun(wrapper3,tempnamescan{1},'UniformOutput',false);
%     
%     tempnamelocs = cellfun(wrapper4,tempnamematches);
%     
%     tempgenebiotype = tempnamescan{1}{tempnamelocs}((length('gene_biotype')+3):end-1);
%     
%     ENSEMBL_biotypes{ictr} = tempgenebiotype;
%     
% end
% 
% %% Have a look at the different biotypes
% 
% all_biotypes_found = unique(ENSEMBL_biotypes);
% 
% %How many of each type are there?
% 
% biotype_nos = zeros(length(all_biotypes_found),1);
% gene_biotype_nos = zeros(length(all_biotypes_found),1);
% 
% for ictr=1:1:length(all_biotypes_found)
%     
%     biotype_nos(ictr) = sum(~cellfun(@isempty,strfind(ENSEMBL_biotypes,all_biotypes_found{ictr})));
%     %this gives the total number (but probably counts things for repeats of
%     %'gene' due to exons, trancripts etc. - I should also pick out how many
%     %are there from 
%     
%     gene_biotype_nos(ictr) = ...
%         sum(~cellfun(@isempty,strfind(ENSEMBL_biotypes,all_biotypes_found{ictr})).* ...
%         gene_logical);
%     
% end
% 
% 
% 
% %% Reduce the gene biotypes down as for the gene positions
% 
% gene_biotypes = cell(NoChromosomes,2);
% 
% for cctr=1:1:NoChromosomes
%    
%     chr_logical = strcmp(ensembl_sc_readin{1},chromosome_labels{cctr});
%     
%     %+strand
%     
%     temp_fwd_locs = logical(chr_logical.*gene_logical.*fwd_logical);
%     
%     gene_biotypes{cctr,1} = ...
%         ENSEMBL_biotypes(temp_fwd_locs);
%     
%     
%     %-strand
%     
%     temp_rev_locs = logical(chr_logical.*gene_logical.*rev_logical);
%     
%     gene_biotypes{cctr,2} = ...
%         ENSEMBL_biotypes(temp_rev_locs);
%     
% end
% 
% 
% %% Further reduce biotypes by removal of non-conformers
% 
% for cctr=1:1:NoChromosomes
%     
%    for sctr=1:2
%       
%        for ictr=size(remove_genes{cctr,sctr},1):-1:1
%            if remove_genes{cctr,sctr}(ictr) == 1
%                gene_biotypes{cctr,sctr}(ictr,:) = [];
%                %gene_names{cctr,sctr}(ictr,:) = [];
%            end
%        end
%        
%        
%    end
%     
%     
% end
% 
% 
% 
% %% Save the imported and processed gene biotypes
% 
% % save([dataloaddir 'gene_biotypes_Saccharomyces_cerevisiae__R64_1_1__101.mat'],'gene_biotypes')
% 
% 
% 
% %% Notes
% 
% %THE 'GENE' TYPE IN THE EARLY ENSEMBL LIST IS ACTUALLY PICKING UP
% %NON-PROTEIN-CODING TUs INCLUDING rRNA, tRNA AND ncRNA 
% %I SHOULD CONSIDER REMOVING THESE LATER
% 
% %%
% 
% toc
