%% Import TIFseq annotations

%Follow what harry described and select the most prominent transript 
%(for each unique gene id) as the starting point

%% Load in the raw data
%(might later do a faster/better textscan import)

dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/TIF_SEQ_ANNOTATIONS/nature12121-s2/';

dataloadfile = 'S2_tcd_mTIFAnno.txt';

TIFseq_raw = importdata([dataloaddir dataloadfile]);


%% See how many lines start with chr (like the first line)

descriplocs = strfind(TIFseq_raw,'chr');
%nb won't give good results if any other text contains the string 'chr'

descripno = sum(~cellfun(@isempty,descriplocs));

%Came out as one - so there is only a single header (which can be safely
%and easily ignored when loading in the data)


%% See if the fields can be regarded as being separated by tabs
%or if it is just spaces

tablocs = strfind(TIFseq_raw,sprintf('\t'));

tabno = unique(cellfun(@numel,tablocs));

%0 detected 

%% How about how many spaces are present

spacelocs = strfind(TIFseq_raw,sprintf(' '));

spaceno = unique(cellfun(@numel,spacelocs));

%between 7 and 11 so will possibly have to get creative with how to 
%nicely read things in

no_readin_els = max(spaceno)+1;

%% Open the raw data with a textscan

%The data in the fields is:
%chr strand t5 t3 ypd gal type name

%So the bit I thought was the final sentence is in fact two pieces -
%hopefully the final space always separates the final column! Otherwise
%things could get a bit tricky!
%This also means the name will appear in different cells depending on the
%line - what a pain.

formatspec = '%f %s %f %f %f %f %s %s %s %s %s %s';
%see if this works with space as a delimiter (the spaces in the final piece
%of text are not delimiters, unfortunately)
%Just had to pad with as many strings as corresponded to the maximum
%detected number of spaces in a line - try recombinging them at the end

fileID1 = fopen([dataloaddir dataloadfile]);

TIFseq_sc_readin = textscan(fileID1,formatspec, ...
    'HeaderLines',1, ...
    'CommentStyle','#');
%'Delimiter','\t', ...
fclose(fileID1);
whos TIFseq_sc_readin



%% Test for the contents of each column

%double check really as this data file was fairly clearly labelled

col1cont = unique(TIFseq_sc_readin{1}(:));
%16 chromosomes in standard number format

col2cont = unique(TIFseq_sc_readin{2}(:));
%'+' or '-' - presumably the strand

col3cont = unique(TIFseq_sc_readin{3}(:));
%positions - presumably the start positoins

col4cont = unique(TIFseq_sc_readin{4}(:));
%positions - presumably the end positions
%oddly far fewer unique values than the starts - the 3' of trancscripts is
%less variable perhaps?

col5cont = unique(TIFseq_sc_readin{5}(:));
%how many transcripts detected in ypd (strangely low number of unique
%entries) - 249399 entries > 0, 1775 unique entries - could be due to many with
%low reads (e.g. 1:100) and only a few with much higher reads - does seem
%to be the case

col6cont = unique(TIFseq_sc_readin{6}(:));
%how many transcripts in gal - similar story about distributions as for ypd
%(although not fully tested)

col7cont = unique(TIFseq_sc_readin{7}(:));
%the beginning of the final text piece with multiple elements of
%information

col8cont = unique(TIFseq_sc_readin{8}(:));

col9cont = unique(TIFseq_sc_readin{9}(:));

col10cont = unique(TIFseq_sc_readin{10}(:));

col11cont = unique(TIFseq_sc_readin{11}(:));

col12cont = unique(TIFseq_sc_readin{12}(:));


%% Check which orientation the start and ends are in

TIFseq_orients_1 = sum((TIFseq_sc_readin{4}(:)-TIFseq_sc_readin{3}(:)) > 0);
TIFseq_orients_2 = sum((TIFseq_sc_readin{3}(:)-TIFseq_sc_readin{4}(:)) > 0);

%% Attempt to reassemble the TIF names

popd_cols = zeros(no_readin_els,1);

TIF_names_recon = cell(size(TIFseq_sc_readin{1},1),1);

TIF_type_recon = cell(size(TIFseq_sc_readin{1},1),1);

for ictr=1:1:size(TIFseq_sc_readin{1},1)
    
    popd_cols(1:6,1) = 1;
    
    for jctr=7:1:no_readin_els
        
        popd_cols(jctr,1) = ~isempty(TIFseq_sc_readin{jctr}{ictr});
        
    end
    
    TIF_names_recon{ictr} = TIFseq_sc_readin{find(popd_cols,1,'last')}{ictr};
    
    
    TIF_type_recon{ictr,1} = TIFseq_sc_readin{7}{ictr};
    for jctr = 8:1:(find(popd_cols,1,'last')-1)
    
        TIF_type_recon{ictr,1} = ...
            horzcat(TIF_type_recon{ictr,1},' ',TIFseq_sc_readin{jctr}{ictr});
    
    end
    
    %NEED TO ALSO EXTRACT AND RECOMBINE THE 'TYPE' STRINGS!!!!!
    
end

%Am I pulling out the CUT/XUT and SUT identifiers correctly?
%Do they even all have identifiers?

%% Test if all CUT/XUT and SUTs have identifiers of some kind


cutxut_logical = strcmp(TIF_type_recon,'CUT/XUT');

cutxut_names = TIF_names_recon(cutxut_logical);

%all seem to have a name of some kind - there are multiple repeats


sut_logical = strcmp(TIF_type_recon,'SUT');

sut_names = TIF_names_recon(sut_logical);

%again, all seem to have a name of some kind


%% How many of the transcripts are not named

NAlocs = strfind(TIF_names_recon,'NA');

NAno = sum(~cellfun(@isempty,NAlocs));

%~34K in the first sample tested

%Will require care when finding the transcripts such that the maximum of
%'N/A' is not selected as a single gene TU location!

%% How many of the transcripts are 'Covering one intact ORF'

COIOlocs = strfind(TIF_type_recon,'Covering one intact ORF');

COIOno = sum(~cellfun(@isempty,COIOlocs));

%~185K in the first sample tested



%% How many unnamed transcripts also cover one intact ORF

NA_COIO_overlap = ...
    sum((~cellfun(@isempty,NAlocs)).*(~cellfun(@isempty,COIOlocs)));

%~9000 for the first sample tested


%% How many TIFseq transcripts are observed in YPD

TIFseq_inypd = TIFseq_sc_readin{5} > 0;

TIFseq_inypd_no = sum(TIFseq_inypd);



%% How many of the transcripts 'Overlap 5' of one ORF' ?

TIFseq_types = unique(TIF_type_recon);
%The one I am interested in is usually the 7th
%I am going to copy it straight out of this, as the name contains inverted
%commas for the 5' (prime) part which could mess up the interpretation of
%strings

OL5PRIMElocs = strcmp(TIF_type_recon,TIFseq_types{7});

OL5PRIMEno = sum(OL5PRIMElocs);

%very much fewer of these


%%

%no need to specify chromosome labels in this case as they are stored as
%numerical data

NoChromosomes = 16;

minGeneSize = 1000; %GM changed to 500 for more signal


%% Reduce the data down to a list of genes

%In this case will also have to select the most populous transcript for
%each unique ID
%Also will require a very different loop structure to get things working
%properly

unique_IDs = unique(TIF_names_recon);

gene_positions = cell(NoChromosomes,2);
%16 chromosomes, 2 strands
gene_names = cell(NoChromosomes,2);
gene_types = cell(NoChromosomes,2);
%Get these at the same time as the positions in this script - makes even
%more sense to extract at the same time for this annotation set

%No mitochondrial stuff present in this data set


chromosome_counter = zeros(NoChromosomes,2);

for ictr=1:1:length(unique_IDs)
    
    %ignore the entries if the ID is 'NA'
    if ~strcmp(unique_IDs{ictr},'NA')
        
        temp_locs = strcmp(TIF_names_recon,unique_IDs{ictr});
        
        temp_chr = unique(TIFseq_sc_readin{1}(temp_locs));
        
        if (numel(temp_chr) > 1)
            disp(['Error detected. ' unique_IDs{ictr} ...
                ' seems to be present on multiple chromosomes.' ])
        end
        
        temp_strand_char = unique(TIFseq_sc_readin{2}(temp_locs));
        
        if (numel(temp_strand_char) > 1)
            disp(['Error detected. ' unique_IDs{ictr} ...
                ' seems to be present on multiple strands.' ])
        end
        
        if strcmp(temp_strand_char,'+')
            temp_strand = 1;
        else
            temp_strand = 2;
        end
        
        
        %Needs to have an entry in the ypd column or don't use this transcript
        
        temp_ypd_vals = ...
            TIFseq_sc_readin{5}(temp_locs);
        
        
        if (sum(temp_ypd_vals) > 0)
            
            temp_starts = ...
                TIFseq_sc_readin{3}(temp_locs);
            
            temp_ends = ...
                TIFseq_sc_readin{4}(temp_locs);
            
            temp_names = TIF_names_recon(temp_locs);
            temp_types = TIF_type_recon(temp_locs);
            
            
            temp_most_occ_ypd = ...
                find(temp_ypd_vals == max(temp_ypd_vals));
            
            %NB if there are multiple most represented transcripts - simply
            %pick the first one (a mean would be strange)
            
            chromosome_counter(temp_chr,temp_strand) = ...
                chromosome_counter(temp_chr,temp_strand) + 1;
            
            if (temp_strand == 1)
                gene_positions{temp_chr,temp_strand}(chromosome_counter(temp_chr,temp_strand),:) = ...
                    [temp_starts(temp_most_occ_ypd(1)), ...
                    temp_ends(temp_most_occ_ypd(1))];
            else
                gene_positions{temp_chr,temp_strand}(chromosome_counter(temp_chr,temp_strand),:) = ...
                    [temp_ends(temp_most_occ_ypd(1)), ...
                    temp_starts(temp_most_occ_ypd(1))];
            end
            %NB to keep conformity with previous data sets, the start and end
            %of each gene are presented in order of appearance (from left to
            %right) so the start appears second for the reverse strand!
            
            
            
            
            gene_names{temp_chr,temp_strand}{chromosome_counter(temp_chr,temp_strand)} = ...
                temp_names(temp_most_occ_ypd(1));
            gene_types{temp_chr,temp_strand}{chromosome_counter(temp_chr,temp_strand)} = ...
                temp_types(temp_most_occ_ypd(1));
            
            
            
            
            
            %error checking:
            if (temp_starts(temp_most_occ_ypd(1)) == 0)
                disp(['There was a transcript with start 0 ' num2str(ictr)])
            end
            
            
        end
        
        
    end
    
    
end


%% To remove overlaps - TUs must be sorted by chromosomal position

%sorting by the leftmost position should be sufficient
%should probably do some checking that nothing is going wrong at some point

gene_positions_old = gene_positions;
gene_names_old = gene_names;
gene_types_old = gene_types;

for cctr=1:1:NoChromosomes
    
   for sctr=1:2
      
       [new_starts,new_positions_idx] = ...
           sort(gene_positions{cctr,sctr}(:,1));
       
       new_ends = gene_positions{cctr,sctr}(new_positions_idx,2);
       
       new_names = gene_names{cctr,sctr}(new_positions_idx);
       
       new_types = gene_types{cctr,sctr}(new_positions_idx);
      
       
       gene_positions{cctr,sctr} = [new_starts, new_ends];
       
       gene_names{cctr,sctr} = new_names;
       
       gene_types{cctr,sctr} = new_types;
       
       %for checking
       
       if ( sum( ...
               (gene_positions{cctr,sctr}(2:end,1) - ...
               gene_positions{cctr,sctr}(1:end-1,2)) < 0 ) ~= 0)
           disp([num2str(sum( ...
               (gene_positions{cctr,sctr}(2:end,1) - ...
               gene_positions{cctr,sctr}(1:end-1,2)) < 0 )) ...
               ' overlap(s) on chr ' num2str(cctr) ...
               ' strand ' num2str(sctr)])
           
       end
       
       
       if ( sum( ...
               (gene_positions{cctr,sctr}(:,2) - ...
               gene_positions{cctr,sctr}(:,1)) < 0 ) ~= 0)
           disp(['Ordering oddness on chr ' num2str(cctr) ...
               'strand' num2str(sctr)])
           
       end
       
       
       
   end
    
end


%% Remove overlapping genes and ones < minGeneSize bp - (normally 1000bp) GM changed to 500 
%GM duplicated and edited below to remove genes with another gene within
%500bp
%or a different length if set earlier


%Are the genes in the order in which they appear on the chromosome?
%This was not the case for the TIFseq annotations but they were then sorted
%and checked for appropriate ordering!

remove_genes = cell(NoChromosomes,2);
for cctr=1:1:NoChromosomes
    for sctr=1:2
        remove_genes{cctr,sctr} = zeros(size(gene_positions{cctr,sctr},1),1);
    end
end


for cctr=1:1:NoChromosomes
   
    strand = 1; %forward
    
    for gctr=1:1:size(gene_positions{cctr,strand},1)
       
        for ictr = (gctr+1):size(gene_positions{cctr,strand},1)
            if( (gene_positions{cctr,strand}(gctr,2) >= gene_positions{cctr,strand}(ictr,1)) && ...   
                    (gene_positions{cctr,strand}(gctr,1) <= gene_positions{cctr,strand}(ictr,2)) )
                remove_genes{cctr,strand}(gctr) = 1;
                remove_genes{cctr,strand}(ictr) = 1;
            end            
        end
        %Also check for overlaps with the opposing strand here
        for ictr = 1:size(gene_positions{cctr,3-strand},1)
           if( (gene_positions{cctr,strand}(gctr,2) >= gene_positions{cctr,3-strand}(ictr,1)) && ...
                   (gene_positions{cctr,strand}(gctr,1) <= gene_positions{cctr,3-strand}(ictr,2)) )
               remove_genes{cctr,strand}(gctr) = 1;
               remove_genes{cctr,3-strand}(ictr) = 1;
           end
        end
        
    end
    
    
    strand = 2; %reverse
    
    for gctr=1:1:size(gene_positions{cctr,strand},1)
       
        for ictr = (gctr+1):size(gene_positions{cctr,strand},1)
            if( (gene_positions{cctr,strand}(gctr,2) >= gene_positions{cctr,strand}(ictr,1)) && ...
                    (gene_positions{cctr,strand}(gctr,1) <= gene_positions{cctr,strand}(ictr,2)) )
                remove_genes{cctr,strand}(gctr) = 1;
                remove_genes{cctr,strand}(ictr) = 1;
            end            
        end
        %Also check for overlaps with the opposing strand here
        for ictr = 1:size(gene_positions{cctr,3-strand},1)
           if( (gene_positions{cctr,strand}(gctr,2) >= gene_positions{cctr,3-strand}(ictr,1)) && ...
                   (gene_positions{cctr,strand}(gctr,1) <= gene_positions{cctr,3-strand}(ictr,2)) )
               remove_genes{cctr,strand}(gctr) = 1;
               remove_genes{cctr,3-strand}(ictr) = 1;
           end
        end
        
    end
    
    
    %I should consider also checking for overlap with non-coding
    %transcripts etc.
    %For the ensembl set it does seem to be detecing some ncRNA and others
    % so this might be OK (for the ENSEMBL case, might stilll be necessary
    % with other annotation sets)
    
    
end


%Also remove genes that are shorter than the specified minimum size here

for cctr=1:1:NoChromosomes
    
   for sctr=1:2
       
       for gctr=1:1:size(gene_positions{cctr,sctr},1)
       
           if (abs(gene_positions{cctr,sctr}(gctr,2) - gene_positions{cctr,sctr}(gctr,1)) < minGeneSize)
               remove_genes{cctr,sctr}(gctr) = 1;
           end
           
       end
   end
    
    
end


%Now remove the actual genes
%NB doing this in a different way than Tom
%going from right to left such that don't need to check if the structure
%has gotten too small
%Also not 100% sure how Tom's algorithm for this works
for cctr=1:1:NoChromosomes
    for sctr=1:2
       for ictr=size(remove_genes{cctr,sctr},1):-1:1
           if remove_genes{cctr,sctr}(ictr)==1
               gene_positions{cctr,sctr}(ictr,:) = [];
           end
       end
    end
end


no_genes_aftertrimming = 0;

for cctr=1:1:NoChromosomes

    for sctr=1:2
       
        no_genes_aftertrimming = no_genes_aftertrimming + ...
            size(gene_positions{cctr,sctr},1);
        
    end

end
    
%In this version tracked gene IDs and types from the beginning
%These should also have the overlaps removed

for cctr=1:1:NoChromosomes
    for sctr=1:2
       for ictr=size(remove_genes{cctr,sctr},1):-1:1
           if remove_genes{cctr,sctr}(ictr)==1
               gene_names{cctr,sctr}(ictr) = [];
               gene_types{cctr,sctr}(ictr) = [];
           end
       end
    end
end

%% Remove overlapping genes and ones < minGeneSize bp - (normally 1000bp) GM changed to 500
%or a different length if set earlier


%Are the genes in the order in which they appear on the chromosome?
%This was not the case for the TIFseq annotations but they were then sorted
%and checked for appropriate ordering!

remove_genes = cell(NoChromosomes,2);
for cctr=1:1:NoChromosomes
    for sctr=1:2
        remove_genes{cctr,sctr} = zeros(size(gene_positions{cctr,sctr},1),1);
    end
end


for cctr=1:1:NoChromosomes
   
    strand = 1; %forward
    
    for gctr=1:1:size(gene_positions{cctr,strand},1)
       
        for ictr = (gctr+1):size(gene_positions{cctr,strand},1)
            if( (gene_positions{cctr,strand}(gctr,2) >= (gene_positions{cctr,strand}(ictr,1) - 500)) && ... %GM -500 allows genes with  
                    (gene_positions{cctr,strand}(gctr,1) <= gene_positions{cctr,strand}(ictr,2)) )
                remove_genes{cctr,strand}(gctr) = 1;
                remove_genes{cctr,strand}(ictr) = 0; %GM changed this value to zero so only the gene with cramped end is removed.
            end            
        end
        %Also check for overlaps with the opposing strand here
        for ictr = 1:size(gene_positions{cctr,3-strand},1)
           if( (gene_positions{cctr,strand}(gctr,2) >= gene_positions{cctr,3-strand}(ictr,1)) && ...
                   (gene_positions{cctr,strand}(gctr,1) <= gene_positions{cctr,3-strand}(ictr,2)) )
               remove_genes{cctr,strand}(gctr) = 1;
               remove_genes{cctr,3-strand}(ictr) = 1;
           end
        end
        
    end
    
    
    strand = 2; %reverse
    
    for gctr=1:1:size(gene_positions{cctr,strand},1)
       
        for ictr = (gctr+1):size(gene_positions{cctr,strand},1)
            if( (gene_positions{cctr,strand}(gctr,2) >= (gene_positions{cctr,strand}(ictr,1) - 500)) && ... %GM -500 to remove cramped end genes
                    (gene_positions{cctr,strand}(gctr,1) <= gene_positions{cctr,strand}(ictr,2)) )
                remove_genes{cctr,strand}(gctr) = 1;
                remove_genes{cctr,strand}(ictr) = 0;
            end            
        end
        %Also check for overlaps with the opposing strand here
        for ictr = 1:size(gene_positions{cctr,3-strand},1)
           if( (gene_positions{cctr,strand}(gctr,2) >= gene_positions{cctr,3-strand}(ictr,1)) && ...
                   (gene_positions{cctr,strand}(gctr,1) <= gene_positions{cctr,3-strand}(ictr,2)) )
               remove_genes{cctr,strand}(gctr) = 1;
               remove_genes{cctr,3-strand}(ictr) = 1;
           end
        end
        
    end
    
    
    %I should consider also checking for overlap with non-coding
    %transcripts etc.
    %For the ensembl set it does seem to be detecing some ncRNA and others
    % so this might be OK (for the ENSEMBL case, might stilll be necessary
    % with other annotation sets)
    
    
end


%Also remove genes that are shorter than the specified minimum size here

for cctr=1:1:NoChromosomes
    
   for sctr=1:2
       
       for gctr=1:1:size(gene_positions{cctr,sctr},1)
       
           if (abs(gene_positions{cctr,sctr}(gctr,2) - gene_positions{cctr,sctr}(gctr,1)) < minGeneSize)
               remove_genes{cctr,sctr}(gctr) = 1;
           end
           
       end
   end
    
    
end


%Now remove the actual genes
%NB doing this in a different way than Tom
%going from right to left such that don't need to check if the structure
%has gotten too small
%Also not 100% sure how Tom's algorithm for this works
for cctr=1:1:NoChromosomes
    for sctr=1:2
       for ictr=size(remove_genes{cctr,sctr},1):-1:1
           if remove_genes{cctr,sctr}(ictr)==1
               gene_positions{cctr,sctr}(ictr,:) = [];
           end
       end
    end
end


no_genes_aftertrimming = 0;

for cctr=1:1:NoChromosomes

    for sctr=1:2
       
        no_genes_aftertrimming = no_genes_aftertrimming + ...
            size(gene_positions{cctr,sctr},1);
        
    end

end
    
%In this version tracked gene IDs and types from the beginning
%These should also have the overlaps removed

for cctr=1:1:NoChromosomes
    for sctr=1:2
       for ictr=size(remove_genes{cctr,sctr},1):-1:1
           if remove_genes{cctr,sctr}(ictr)==1
               gene_names{cctr,sctr}(ictr) = [];
               gene_types{cctr,sctr}(ictr) = [];
           end
       end
    end
end



%% I should check at the end that no overlaps remain:

% Separate section as this could take quite a while

%determine the largest location for each chromosome
chr_rough_sizes = zeros(NoChromosomes,1);
rough_chromosomes = cell(NoChromosomes,1);

for cctr=1:1:NoChromosomes

    chr_rough_sizes(cctr) = ...
        max( ...
        [max(max(gene_positions{cctr,1})) ...
        max(max(gene_positions{cctr,2}))] );
    
    rough_chromosomes{cctr} = zeros(chr_rough_sizes(cctr),1);
    
end

%I think the matrix I was envisioning for this will be simply too big to 
%store in memory - actually it might just work I made one chromosome and it
%worked

for cctr = 1:1:NoChromosomes
    
    for sctr=1:1:2
       
        for gctr = 1:1:size(gene_positions{cctr,sctr},1)
        
            rough_chromosomes{cctr}( ...
                gene_positions{cctr,sctr}(gctr,1):gene_positions{cctr,sctr}(gctr,2) ) = ...
                rough_chromosomes{cctr}( ...
                gene_positions{cctr,sctr}(gctr,1):gene_positions{cctr,sctr}(gctr,2) ) + 1;
            
        end
        
    end
    
end

overlaps = zeros(NoChromosomes,1);

for ictr=1:1:NoChromosomes
   
    overlaps(ictr) = sum(rough_chromosomes{ictr}>1);
    
end

%Seems to have worked - there are no detected overlaps

%I could also remove genes that are considered to be too small here!
%Think I already did that above! - limiting to 1000bp


%% Save the imported and processed gene positions, names and types

save([dataloaddir 'gene_positions_pelechano_TIFseq_nogene500_min1000.mat'],'gene_positions')
save([dataloaddir 'gene_names_pelechano_TIFseq_nogene500_min1000.mat'],'gene_names')
save([dataloaddir 'gene_types_pelechano_TIFseq_nogene500_min1000.mat'],'gene_types')


%% Notes:

%Initially TUs do not appear in the order they appear on the chromosome
%Maybe CUT/XUTs go before genes for example.
%In order to properly remove overlaps, I think they will have to be sorted!
%They have now been sorted and it appears to be working correctly

%Might need to move gene_names to gene_ids for conformity with previous
%annotation sets




