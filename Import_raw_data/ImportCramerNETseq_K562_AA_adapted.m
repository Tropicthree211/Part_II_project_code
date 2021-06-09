%% Import the Cramer K562 mNETseq and RNAseq data
%Adapted from ImportCramerNETseq_AA.m

%Perhaps import the RNAseq in a separate script for crude parallelisation

dataloaddir = '/data/angeldata/TRANSCRIPTION_DATA/CRAMER_NET_SEQ/K562/';



%%

cramer_K562_mNETseq_minus_1_readin = ...
    importdata([dataloaddir 'GSM3518117_N_tCTD_1_CGGAAT_1_K562_wildtype.minus.wig']);

cramer_K562_mNETseq_plus_1_readin = ...
    importdata([dataloaddir 'GSM3518117_N_tCTD_1_CGGAAT_1_K562_wildtype.plus.wig']);


disp('Finished importing raw mNETseq sample 1')



%% 


cramer_K562_mNETseq_1 = cell(24,2);

for strandctr=1:2
   
    cramer_K562_mNETseq_1{1,strandctr} = sparse(zeros(248956422,1));
    cramer_K562_mNETseq_1{2,strandctr} = sparse(zeros(242193529,1));
    cramer_K562_mNETseq_1{3,strandctr} = sparse(zeros(198295559,1));
    cramer_K562_mNETseq_1{4,strandctr} = sparse(zeros(190214555,1));
    cramer_K562_mNETseq_1{5,strandctr} = sparse(zeros(181538259,1));
    cramer_K562_mNETseq_1{6,strandctr} = sparse(zeros(170805979,1));
    cramer_K562_mNETseq_1{7,strandctr} = sparse(zeros(159345973,1));
    cramer_K562_mNETseq_1{8,strandctr} = sparse(zeros(145138636,1));
    cramer_K562_mNETseq_1{9,strandctr} = sparse(zeros(138394717,1));
    cramer_K562_mNETseq_1{10,strandctr} = sparse(zeros(133797422,1));
    cramer_K562_mNETseq_1{11,strandctr} = sparse(zeros(135086622,1));
    cramer_K562_mNETseq_1{12,strandctr} = sparse(zeros(133275309,1));
    cramer_K562_mNETseq_1{13,strandctr} = sparse(zeros(114364328,1));
    cramer_K562_mNETseq_1{14,strandctr} = sparse(zeros(107043718,1));
    cramer_K562_mNETseq_1{15,strandctr} = sparse(zeros(101991189,1));
    cramer_K562_mNETseq_1{16,strandctr} = sparse(zeros(90338345,1));
    cramer_K562_mNETseq_1{17,strandctr} = sparse(zeros(83257441,1));
    cramer_K562_mNETseq_1{18,strandctr} = sparse(zeros(80373285,1));
    cramer_K562_mNETseq_1{19,strandctr} = sparse(zeros(58617616,1));
    cramer_K562_mNETseq_1{20,strandctr} = sparse(zeros(64444167,1));
    cramer_K562_mNETseq_1{21,strandctr} = sparse(zeros(46709983,1));
    cramer_K562_mNETseq_1{22,strandctr} = sparse(zeros(50818468,1));
    cramer_K562_mNETseq_1{23,strandctr} = sparse(zeros(156040895,1));
    cramer_K562_mNETseq_1{24,strandctr} = sparse(zeros(57227415,1));
    
end


%should I be calculating the max point for each cell from the data itself?
%Something for the future perhaps

%%

cramer_K562_mNETseq_minus_1_readin.textdata = ...
    cramer_K562_mNETseq_minus_1_readin.textdata(2:end);

%No longer sure this step is entirely necessary - it seems similar text
%pieces are repeated throughout breaking up the results into sections on 
%chromosomes
%This step simply removes the first of these
%The text data is only one larger than the full data though
%why aren't the data lines for the other 'separators' empty?

%Checked this and it is fine and the removal is a necessary step
%The separators after the first one have NaNs in the data part
%I have checked around one of the separators to show that the correct
%data was aligned with the correct text.


ow_counter = 0;

for ictr=1:1:length(cramer_K562_mNETseq_minus_1_readin.textdata)
   
    switch cramer_K562_mNETseq_minus_1_readin.textdata{ictr}
        
        case 'chr1'
            chr=1;
        case 'chr2'
            chr=2;
        case 'chr3'
            chr=3;
        case 'chr4'
            chr=4;
        case 'chr5'
            chr=5;
        case 'chr6'
            chr=6;
        case 'chr7'
            chr=7;
        case 'chr8'
            chr=8;
        case 'chr9'
            chr=9;
        case 'chr10'
            chr=10;
        case 'chr11'
            chr=11;
        case 'chr12'
            chr=12;
        case 'chr13'
            chr=13;
        case 'chr14'
            chr=14;
        case 'chr15'
            chr=15;
        case 'chr16'
            chr=16;
        case 'chr17'
            chr=17;
        case 'chr18'
            chr=18;
        case 'chr19'
            chr=19;
        case 'chr20'
            chr=20;
        case 'chr21'
            chr=21;
        case 'chr22'
            chr=22;
        case 'chrX'
            chr=23;
        case 'chrY'
            chr=24;
        otherwise
            ow_counter = ow_counter + 1;
            continue;
    end
    
    %cramer_Raji_NETseq_1a{chr,2}( ...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,1):...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,2)-1) = ...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,3);
    
    %cramer_Raji_NETseq_1a{chr,2}( ...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,1)+1:...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,2) ) = ...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,3);
    %cramer data starts from zero! - Need to make sure what the actual
    %positions are here!
    
    cramer_K562_mNETseq_1{chr,2}( ...
        cramer_K562_mNETseq_minus_1_readin.data(ictr,1)+1:...
        cramer_K562_mNETseq_minus_1_readin.data(ictr,2) ) = ...
        cramer_K562_mNETseq_minus_1_readin.data(ictr,3);
    %cramer data starts from zero1! - Still need to make sure what the 
    %actual positions are here!
    
    if mod(ictr,10000) == 0
        disp(['Completed ' num2str(ictr) ' of ' ...
            num2str(length(cramer_K562_mNETseq_minus_1_readin.textdata)) ...
            ' 1 minus iterations.'])
    end
    
end

%%

cramer_K562_mNETseq_plus_1_readin.textdata = ...
    cramer_K562_mNETseq_plus_1_readin.textdata(2:end);

%No longer sure this step is entirely necessary - it seems similar text
%pieces are repeated throughout breaking up the results into sections on 
%chromosomes
%This step simply removes the first of these
%The text data is only one larger than the full data though
%why aren't the data lines for the other 'separators' empty?

%Checked this and it is fine and the removal is a necessary step
%The separators after the first one have NaNs in the data part
%I have checked around one of the separators to show that the correct
%data was aligned with the correct text.


ow2_counter = 0;

for ictr=1:1:length(cramer_K562_mNETseq_plus_1_readin.textdata)
   
    switch cramer_K562_mNETseq_plus_1_readin.textdata{ictr}
        
        case 'chr1'
            chr=1;
        case 'chr2'
            chr=2;
        case 'chr3'
            chr=3;
        case 'chr4'
            chr=4;
        case 'chr5'
            chr=5;
        case 'chr6'
            chr=6;
        case 'chr7'
            chr=7;
        case 'chr8'
            chr=8;
        case 'chr9'
            chr=9;
        case 'chr10'
            chr=10;
        case 'chr11'
            chr=11;
        case 'chr12'
            chr=12;
        case 'chr13'
            chr=13;
        case 'chr14'
            chr=14;
        case 'chr15'
            chr=15;
        case 'chr16'
            chr=16;
        case 'chr17'
            chr=17;
        case 'chr18'
            chr=18;
        case 'chr19'
            chr=19;
        case 'chr20'
            chr=20;
        case 'chr21'
            chr=21;
        case 'chr22'
            chr=22;
        case 'chrX'
            chr=23;
        case 'chrY'
            chr=24;
        otherwise
            ow2_counter = ow2_counter + 1;
            continue;
    end
    
    %cramer_Raji_NETseq_1a{chr,2}( ...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,1):...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,2)-1) = ...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,3);
    
    %cramer_Raji_NETseq_1a{chr,2}( ...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,1)+1:...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,2) ) = ...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,3);
    %cramer data starts from zero! - Need to make sure what the actual
    %positions are here!
    
    cramer_K562_mNETseq_1{chr,1}( ...
        cramer_K562_mNETseq_plus_1_readin.data(ictr,1)+1:...
        cramer_K562_mNETseq_plus_1_readin.data(ictr,2) ) = ...
        cramer_K562_mNETseq_plus_1_readin.data(ictr,3);
    %cramer data starts from zero1! - Still need to make sure what the 
    %actual positions are here!
    
    if mod(ictr,10000) == 0
        disp(['Completed ' num2str(ictr) ' of ' ...
            num2str(length(cramer_K562_mNETseq_minus_1_readin.textdata)) ...
            ' 1 plus iterations.'])
    end
    
end

%%

save([dataloaddir 'cramer_K562_mNETseq_1.mat'],'cramer_K562_mNETseq_1','-v7.3')


%%

cramer_K562_mNETseq_minus_2_readin = ...
    importdata([dataloaddir 'GSM3518118_N_tCTD_2_CTAGCT_2_K562_wildtype.minus.wig']);

cramer_K562_mNETseq_plus_2_readin = ...
    importdata([dataloaddir 'GSM3518118_N_tCTD_2_CTAGCT_2_K562_wildtype.plus.wig']);


disp('Finished importing raw mNETseq sample 2')


%% 


cramer_K562_mNETseq_2 = cell(24,2);

for strandctr=1:2
   
    cramer_K562_mNETseq_2{1,strandctr} = sparse(zeros(248956422,1));
    cramer_K562_mNETseq_2{2,strandctr} = sparse(zeros(242193529,1));
    cramer_K562_mNETseq_2{3,strandctr} = sparse(zeros(198295559,1));
    cramer_K562_mNETseq_2{4,strandctr} = sparse(zeros(190214555,1));
    cramer_K562_mNETseq_2{5,strandctr} = sparse(zeros(181538259,1));
    cramer_K562_mNETseq_2{6,strandctr} = sparse(zeros(170805979,1));
    cramer_K562_mNETseq_2{7,strandctr} = sparse(zeros(159345973,1));
    cramer_K562_mNETseq_2{8,strandctr} = sparse(zeros(145138636,1));
    cramer_K562_mNETseq_2{9,strandctr} = sparse(zeros(138394717,1));
    cramer_K562_mNETseq_2{10,strandctr} = sparse(zeros(133797422,1));
    cramer_K562_mNETseq_2{11,strandctr} = sparse(zeros(135086622,1));
    cramer_K562_mNETseq_2{12,strandctr} = sparse(zeros(133275309,1));
    cramer_K562_mNETseq_2{13,strandctr} = sparse(zeros(114364328,1));
    cramer_K562_mNETseq_2{14,strandctr} = sparse(zeros(107043718,1));
    cramer_K562_mNETseq_2{15,strandctr} = sparse(zeros(101991189,1));
    cramer_K562_mNETseq_2{16,strandctr} = sparse(zeros(90338345,1));
    cramer_K562_mNETseq_2{17,strandctr} = sparse(zeros(83257441,1));
    cramer_K562_mNETseq_2{18,strandctr} = sparse(zeros(80373285,1));
    cramer_K562_mNETseq_2{19,strandctr} = sparse(zeros(58617616,1));
    cramer_K562_mNETseq_2{20,strandctr} = sparse(zeros(64444167,1));
    cramer_K562_mNETseq_2{21,strandctr} = sparse(zeros(46709983,1));
    cramer_K562_mNETseq_2{22,strandctr} = sparse(zeros(50818468,1));
    cramer_K562_mNETseq_2{23,strandctr} = sparse(zeros(156040895,1));
    cramer_K562_mNETseq_2{24,strandctr} = sparse(zeros(57227415,1));
    
end


%should I be calculating the max point for each cell from the data itself?
%Something for the future perhaps


%%

cramer_K562_mNETseq_minus_2_readin.textdata = ...
    cramer_K562_mNETseq_minus_2_readin.textdata(2:end);



ow3_counter = 0;

for ictr=1:1:length(cramer_K562_mNETseq_minus_2_readin.textdata)
   
    switch cramer_K562_mNETseq_minus_2_readin.textdata{ictr}
        
        case 'chr1'
            chr=1;
        case 'chr2'
            chr=2;
        case 'chr3'
            chr=3;
        case 'chr4'
            chr=4;
        case 'chr5'
            chr=5;
        case 'chr6'
            chr=6;
        case 'chr7'
            chr=7;
        case 'chr8'
            chr=8;
        case 'chr9'
            chr=9;
        case 'chr10'
            chr=10;
        case 'chr11'
            chr=11;
        case 'chr12'
            chr=12;
        case 'chr13'
            chr=13;
        case 'chr14'
            chr=14;
        case 'chr15'
            chr=15;
        case 'chr16'
            chr=16;
        case 'chr17'
            chr=17;
        case 'chr18'
            chr=18;
        case 'chr19'
            chr=19;
        case 'chr20'
            chr=20;
        case 'chr21'
            chr=21;
        case 'chr22'
            chr=22;
        case 'chrX'
            chr=23;
        case 'chrY'
            chr=24;
        otherwise
            ow3_counter = ow3_counter + 1;
            continue;
    end
    
    %cramer_Raji_NETseq_1a{chr,2}( ...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,1):...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,2)-1) = ...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,3);
    
    %cramer_Raji_NETseq_1a{chr,2}( ...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,1)+1:...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,2) ) = ...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,3);
    %cramer data starts from zero! - Need to make sure what the actual
    %positions are here!
    
    cramer_K562_mNETseq_2{chr,2}( ...
        cramer_K562_mNETseq_minus_2_readin.data(ictr,1)+1:...
        cramer_K562_mNETseq_minus_2_readin.data(ictr,2) ) = ...
        cramer_K562_mNETseq_minus_2_readin.data(ictr,3);
    %cramer data starts from zero1! - Still need to make sure what the 
    %actual positions are here!
    
    if mod(ictr,10000) == 0
        disp(['Completed ' num2str(ictr) ' of ' ...
            num2str(length(cramer_K562_mNETseq_minus_2_readin.textdata)) ...
            ' 2 minus iterations.'])
    end
    
end



%%

cramer_K562_mNETseq_plus_2_readin.textdata = ...
    cramer_K562_mNETseq_plus_2_readin.textdata(2:end);

%No longer sure this step is entirely necessary - it seems similar text
%pieces are repeated throughout breaking up the results into sections on 
%chromosomes
%This step simply removes the first of these
%The text data is only one larger than the full data though
%why aren't the data lines for the other 'separators' empty?

%Checked this and it is fine and the removal is a necessary step
%The separators after the first one have NaNs in the data part
%I have checked around one of the separators to show that the correct
%data was aligned with the correct text.


ow4_counter = 0;

for ictr=1:1:length(cramer_K562_mNETseq_plus_2_readin.textdata)
   
    switch cramer_K562_mNETseq_plus_2_readin.textdata{ictr}
        
        case 'chr1'
            chr=1;
        case 'chr2'
            chr=2;
        case 'chr3'
            chr=3;
        case 'chr4'
            chr=4;
        case 'chr5'
            chr=5;
        case 'chr6'
            chr=6;
        case 'chr7'
            chr=7;
        case 'chr8'
            chr=8;
        case 'chr9'
            chr=9;
        case 'chr10'
            chr=10;
        case 'chr11'
            chr=11;
        case 'chr12'
            chr=12;
        case 'chr13'
            chr=13;
        case 'chr14'
            chr=14;
        case 'chr15'
            chr=15;
        case 'chr16'
            chr=16;
        case 'chr17'
            chr=17;
        case 'chr18'
            chr=18;
        case 'chr19'
            chr=19;
        case 'chr20'
            chr=20;
        case 'chr21'
            chr=21;
        case 'chr22'
            chr=22;
        case 'chrX'
            chr=23;
        case 'chrY'
            chr=24;
        otherwise
            ow4_counter = ow4_counter + 1;
            continue;
    end
    
    %cramer_Raji_NETseq_1a{chr,2}( ...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,1):...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,2)-1) = ...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,3);
    
    %cramer_Raji_NETseq_1a{chr,2}( ...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,1)+1:...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,2) ) = ...
    %    cramer_Raji_NETseq_minus_1a_readin.data(ictr,3);
    %cramer data starts from zero! - Need to make sure what the actual
    %positions are here!
    
    cramer_K562_mNETseq_2{chr,1}( ...
        cramer_K562_mNETseq_plus_2_readin.data(ictr,1)+1:...
        cramer_K562_mNETseq_plus_2_readin.data(ictr,2) ) = ...
        cramer_K562_mNETseq_plus_2_readin.data(ictr,3);
    %cramer data starts from zero1! - Still need to make sure what the 
    %actual positions are here!
    
    if mod(ictr,10000) == 0
        disp(['Completed ' num2str(ictr) ' of ' ...
            num2str(length(cramer_K562_mNETseq_minus_1_readin.textdata)) ...
            ' 2 plus iterations.'])
    end
    
end



%%
save([dataloaddir 'cramer_K562_mNETseq_2.mat'],'cramer_K562_mNETseq_2','-v7.3')
