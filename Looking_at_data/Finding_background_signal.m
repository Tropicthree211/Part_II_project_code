
NETseq_nogenes = NETseq_composite;
failed_removal_counter = 0 ;
for cctr = 1:1:NoChromosomes
    for gctr = size(gene_positions{cctr,: }, 1):-1:1 %need to input before removal of genes by overlaps etc 
        try NETseq_nogenes(gene_positions{cctr, :}(gctr, 1):gene_positions{cctr, :}(gctr, 2), :) = [];
            
        catch
            failed_removal_counter = failed_removal_counter +1 ; %if genes overlap will attempt to remove already removed area may fail
        end
    end
end
Background_signal = nanmean(NETseq_nogenes);      


%prototype code for finding background in genome with no genes, finished
%code in file threeprime_xafter_ybefore _end_AC_incorporated.m