%% Import the Cramer NETseq data
%This version written by Andrew (to hopefully deal with bigwig files)
%Adapted from Tom Brown's convert_churchmanHelaFiles.m in
%Elongation_clustering/HeLa_analysis/raw_data/

dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/CRAMER_NET_SEQ/';


%%

cramer_Raji_NETseq_minus_1a_readin = ...
    importdata([dataloaddir 'GSM2728735_N_tCTDeDMSO_8_DMSON1_1_Raji_wildtype.minus.wig']);
cramer_Raji_NETseq_plus_1a_readin = ...
    importdata([dataloaddir 'GSM2728735_N_tCTDeDMSO_8_DMSON1_1_Raji_wildtype.plus.wig']);

disp('Finished importing raw data 1a.')
%Trying the DMSO one first 
%The churchman data was bedgraph so this was never going to work!
%should work for Phil's data though

%It's possible I just have to convert bigwig to wig and then import the
%data!

%Data has now been converted to wig which is hopefully more
%straightforward to process


%%


cramer_Raji_NETseq_1a = cell(24,2);

for i=1:2
    cramer_Raji_NETseq_1a{1,i} = sparse(zeros(248956422,1));
    cramer_Raji_NETseq_1a{2,i} = sparse(zeros(242193529,1));
    cramer_Raji_NETseq_1a{3,i} = sparse(zeros(198295559,1));
    cramer_Raji_NETseq_1a{4,i} = sparse(zeros(190214555,1));
    cramer_Raji_NETseq_1a{5,i} = sparse(zeros(181538259,1));
    cramer_Raji_NETseq_1a{6,i} = sparse(zeros(170805979,1));
    cramer_Raji_NETseq_1a{7,i} = sparse(zeros(159345973,1));
    cramer_Raji_NETseq_1a{8,i} = sparse(zeros(145138636,1));
    cramer_Raji_NETseq_1a{9,i} = sparse(zeros(138394717,1));
    cramer_Raji_NETseq_1a{10,i} = sparse(zeros(133797422,1));
    cramer_Raji_NETseq_1a{11,i} = sparse(zeros(135086622,1));
    cramer_Raji_NETseq_1a{12,i} = sparse(zeros(133275309,1));
    cramer_Raji_NETseq_1a{13,i} = sparse(zeros(114364328,1));
    cramer_Raji_NETseq_1a{14,i} = sparse(zeros(107043718,1));
    cramer_Raji_NETseq_1a{15,i} = sparse(zeros(101991189,1));
    cramer_Raji_NETseq_1a{16,i} = sparse(zeros(90338345,1));
    cramer_Raji_NETseq_1a{17,i} = sparse(zeros(83257441,1));
    cramer_Raji_NETseq_1a{18,i} = sparse(zeros(80373285,1));
    cramer_Raji_NETseq_1a{19,i} = sparse(zeros(58617616,1));
    cramer_Raji_NETseq_1a{20,i} = sparse(zeros(64444167,1));
    cramer_Raji_NETseq_1a{21,i} = sparse(zeros(46709983,1));
    cramer_Raji_NETseq_1a{22,i} = sparse(zeros(50818468,1));
    cramer_Raji_NETseq_1a{23,i} = sparse(zeros(156040895,1));
    cramer_Raji_NETseq_1a{24,i} = sparse(zeros(57227415,1));
    
end

%%

cramer_Raji_NETseq_minus_1a_readin.textdata = ...
    cramer_Raji_NETseq_minus_1a_readin.textdata(2:end);

ow_counter = 0;

for ictr=1:1:length(cramer_Raji_NETseq_minus_1a_readin.textdata)
   
    switch cramer_Raji_NETseq_minus_1a_readin.textdata{ictr}
        
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
    
    cramer_Raji_NETseq_1a{chr,2}( ...
        cramer_Raji_NETseq_minus_1a_readin.data(ictr,1)+1:...
        cramer_Raji_NETseq_minus_1a_readin.data(ictr,2) ) = ...
        cramer_Raji_NETseq_minus_1a_readin.data(ictr,3);
    %cramer data starts from zero! - Need to make sure what the actual
    %positions are here!
    
    if mod(ictr,10000) == 0
        disp(['Completed ' num2str(ictr/10000) ' 1a minus iterations.'])
    end
    
end

%% Some quick(ish) checking

%In the initial test set there were 32919 'otherwise' trackings in the 
%switch statement (i.e. not chr?) - I should have a quick look to see what
%they are (although a quick look is largely impossible)

%possibly some are chrM? I saw that in the GENCODE stuff, I think


%% Now do the plus strand


cramer_Raji_NETseq_plus_1a_readin.textdata = ...
    cramer_Raji_NETseq_plus_1a_readin.textdata(2:end);

ow2_counter = 0;

for ictr=1:1:length(cramer_Raji_NETseq_plus_1a_readin.textdata)
   
    switch cramer_Raji_NETseq_plus_1a_readin.textdata{ictr}
        
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
    
    cramer_Raji_NETseq_1a{chr,1}( ...
        cramer_Raji_NETseq_plus_1a_readin.data(ictr,1)+1:...
        cramer_Raji_NETseq_plus_1a_readin.data(ictr,2) ) = ...
        cramer_Raji_NETseq_plus_1a_readin.data(ictr,3);
    %cramer data starts from zero! - Need to make sure what the actual
    %positions are here!
    
    if mod(ictr,10000) == 0
        disp(['Completed ' num2str(ictr/10000) ' 1a plus iterations.'])
    end
    
end


%%

save([dataloaddir 'cramer_Raji_NETseq_1a.mat'],'cramer_Raji_NETseq_1a','-v7.3')



%%
cramer_Raji_NETseq_minus_1b_readin = ...
    importdata([dataloaddir 'GSM2728736_N_tCTDeDMSO_8_DMSON2_1_Raji_wildtype.minus.wig']);
cramer_Raji_NETseq_plus_1b_readin = ...
    importdata([dataloaddir 'GSM2728736_N_tCTDeDMSO_8_DMSON2_1_Raji_wildtype.plus.wig']);

disp('Finished importing raw data 1b.')
%Trying the DMSO one first 
%The churchman data was bedgraph so this was never going to work!
%should work for Phil's data though

%It's possible I just have to convert bigwig to wig and then import the
%data!

%Data has now been converted to wig which is hopefully more
%straightforward to process


%%


cramer_Raji_NETseq_1b = cell(24,2);

for i=1:2
    cramer_Raji_NETseq_1b{1,i} = sparse(zeros(248956422,1));
    cramer_Raji_NETseq_1b{2,i} = sparse(zeros(242193529,1));
    cramer_Raji_NETseq_1b{3,i} = sparse(zeros(198295559,1));
    cramer_Raji_NETseq_1b{4,i} = sparse(zeros(190214555,1));
    cramer_Raji_NETseq_1b{5,i} = sparse(zeros(181538259,1));
    cramer_Raji_NETseq_1b{6,i} = sparse(zeros(170805979,1));
    cramer_Raji_NETseq_1b{7,i} = sparse(zeros(159345973,1));
    cramer_Raji_NETseq_1b{8,i} = sparse(zeros(145138636,1));
    cramer_Raji_NETseq_1b{9,i} = sparse(zeros(138394717,1));
    cramer_Raji_NETseq_1b{10,i} = sparse(zeros(133797422,1));
    cramer_Raji_NETseq_1b{11,i} = sparse(zeros(135086622,1));
    cramer_Raji_NETseq_1b{12,i} = sparse(zeros(133275309,1));
    cramer_Raji_NETseq_1b{13,i} = sparse(zeros(114364328,1));
    cramer_Raji_NETseq_1b{14,i} = sparse(zeros(107043718,1));
    cramer_Raji_NETseq_1b{15,i} = sparse(zeros(101991189,1));
    cramer_Raji_NETseq_1b{16,i} = sparse(zeros(90338345,1));
    cramer_Raji_NETseq_1b{17,i} = sparse(zeros(83257441,1));
    cramer_Raji_NETseq_1b{18,i} = sparse(zeros(80373285,1));
    cramer_Raji_NETseq_1b{19,i} = sparse(zeros(58617616,1));
    cramer_Raji_NETseq_1b{20,i} = sparse(zeros(64444167,1));
    cramer_Raji_NETseq_1b{21,i} = sparse(zeros(46709983,1));
    cramer_Raji_NETseq_1b{22,i} = sparse(zeros(50818468,1));
    cramer_Raji_NETseq_1b{23,i} = sparse(zeros(156040895,1));
    cramer_Raji_NETseq_1b{24,i} = sparse(zeros(57227415,1));
    
end

%%

cramer_Raji_NETseq_minus_1b_readin.textdata = ...
    cramer_Raji_NETseq_minus_1b_readin.textdata(2:end);

ow_counter = 0;

for ictr=1:1:length(cramer_Raji_NETseq_minus_1b_readin.textdata)
   
    switch cramer_Raji_NETseq_minus_1b_readin.textdata{ictr}
        
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
    
    %cramer_Raji_NETseq_1b{chr,2}( ...
    %    cramer_Raji_NETseq_minus_1b_readin.data(ictr,1):...
    %    cramer_Raji_NETseq_minus_1b_readin.data(ictr,2)-1) = ...
    %    cramer_Raji_NETseq_minus_1b_readin.data(ictr,3);
    
    cramer_Raji_NETseq_1b{chr,2}( ...
        cramer_Raji_NETseq_minus_1b_readin.data(ictr,1)+1:...
        cramer_Raji_NETseq_minus_1b_readin.data(ictr,2) ) = ...
        cramer_Raji_NETseq_minus_1b_readin.data(ictr,3);
    %cramer data starts from zero! - Need to make sure what the actual
    %positions are here!
    
    if mod(ictr,10000) == 0
        disp(['Completed ' num2str(ictr/10000) ' 1b minus iterations.'])
    end
    
end

%% Some quick(ish) checking

%In the initial test set there were 32919 'otherwise' trackings in the 
%switch statement (i.e. not chr?) - I should have a quick look to see what
%they are (although a quick look is largely impossible)


%% Now do the plus strand


cramer_Raji_NETseq_plus_1b_readin.textdata = ...
    cramer_Raji_NETseq_plus_1b_readin.textdata(2:end);

ow2_counter = 0;

for ictr=1:1:length(cramer_Raji_NETseq_plus_1b_readin.textdata)
   
    switch cramer_Raji_NETseq_plus_1b_readin.textdata{ictr}
        
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
    
    %cramer_Raji_NETseq_1b{chr,2}( ...
    %    cramer_Raji_NETseq_minus_1b_readin.data(ictr,1):...
    %    cramer_Raji_NETseq_minus_1b_readin.data(ictr,2)-1) = ...
    %    cramer_Raji_NETseq_minus_1b_readin.data(ictr,3);
    
    cramer_Raji_NETseq_1b{chr,1}( ...
        cramer_Raji_NETseq_plus_1b_readin.data(ictr,1)+1:...
        cramer_Raji_NETseq_plus_1b_readin.data(ictr,2) ) = ...
        cramer_Raji_NETseq_plus_1b_readin.data(ictr,3);
    %cramer data starts from zero! - Need to make sure what the actual
    %positions are here!
    
    if mod(ictr,10000) == 0
        disp(['Completed ' num2str(ictr/10000) ' 1b plus iterations.'])
    end
    
end


%%

save([dataloaddir 'cramer_Raji_NETseq_1b.mat'],'cramer_Raji_NETseq_1b','-v7.3')

%%

cramer_Raji_NETseq_minus_2_readin = ...
    importdata([dataloaddir 'GSM2728737_N_tCTDeDMSO_10_DMSON3_1_Raji_wildtype.minus.wig']);
cramer_Raji_NETseq_plus_2_readin = ...
    importdata([dataloaddir 'GSM2728737_N_tCTDeDMSO_10_DMSON3_1_Raji_wildtype.plus.wig']);

disp('Finished importing raw data.')
%Trying the DMSO one first 
%The churchman data was bedgraph so this was never going to work!
%should work for Phil's data though

%It's possible I just have to convert bigwig to wig and then import the
%data!

%Data has now been converted to wig which is hopefully more
%straightforward to process


%%


cramer_Raji_NETseq_2 = cell(24,2);

for i=1:2
    cramer_Raji_NETseq_2{1,i} = sparse(zeros(248956422,1));
    cramer_Raji_NETseq_2{2,i} = sparse(zeros(242193529,1));
    cramer_Raji_NETseq_2{3,i} = sparse(zeros(198295559,1));
    cramer_Raji_NETseq_2{4,i} = sparse(zeros(190214555,1));
    cramer_Raji_NETseq_2{5,i} = sparse(zeros(181538259,1));
    cramer_Raji_NETseq_2{6,i} = sparse(zeros(170805979,1));
    cramer_Raji_NETseq_2{7,i} = sparse(zeros(159345973,1));
    cramer_Raji_NETseq_2{8,i} = sparse(zeros(145138636,1));
    cramer_Raji_NETseq_2{9,i} = sparse(zeros(138394717,1));
    cramer_Raji_NETseq_2{10,i} = sparse(zeros(133797422,1));
    cramer_Raji_NETseq_2{11,i} = sparse(zeros(135086622,1));
    cramer_Raji_NETseq_2{12,i} = sparse(zeros(133275309,1));
    cramer_Raji_NETseq_2{13,i} = sparse(zeros(114364328,1));
    cramer_Raji_NETseq_2{14,i} = sparse(zeros(107043718,1));
    cramer_Raji_NETseq_2{15,i} = sparse(zeros(101991189,1));
    cramer_Raji_NETseq_2{16,i} = sparse(zeros(90338345,1));
    cramer_Raji_NETseq_2{17,i} = sparse(zeros(83257441,1));
    cramer_Raji_NETseq_2{18,i} = sparse(zeros(80373285,1));
    cramer_Raji_NETseq_2{19,i} = sparse(zeros(58617616,1));
    cramer_Raji_NETseq_2{20,i} = sparse(zeros(64444167,1));
    cramer_Raji_NETseq_2{21,i} = sparse(zeros(46709983,1));
    cramer_Raji_NETseq_2{22,i} = sparse(zeros(50818468,1));
    cramer_Raji_NETseq_2{23,i} = sparse(zeros(156040895,1));
    cramer_Raji_NETseq_2{24,i} = sparse(zeros(57227415,1));
    
end

%%

cramer_Raji_NETseq_minus_2_readin.textdata = ...
    cramer_Raji_NETseq_minus_2_readin.textdata(2:end);

ow_counter = 0;

for ictr=1:1:length(cramer_Raji_NETseq_minus_2_readin.textdata)
   
    switch cramer_Raji_NETseq_minus_2_readin.textdata{ictr}
        
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
    
    %cramer_Raji_NETseq_2{chr,2}( ...
    %    cramer_Raji_NETseq_minus_2_readin.data(ictr,1):...
    %    cramer_Raji_NETseq_minus_2_readin.data(ictr,2)-1) = ...
    %    cramer_Raji_NETseq_minus_2_readin.data(ictr,3);
    
    cramer_Raji_NETseq_2{chr,2}( ...
        cramer_Raji_NETseq_minus_2_readin.data(ictr,1)+1:...
        cramer_Raji_NETseq_minus_2_readin.data(ictr,2) ) = ...
        cramer_Raji_NETseq_minus_2_readin.data(ictr,3);
    %cramer data starts from zero! - Need to make sure what the actual
    %positions are here!
    
    if mod(ictr,10000) == 0
        disp(['Completed ' num2str(ictr/10000) ' 2 minus iterations.'])
    end
    
end

%% Some quick(ish) checking

%In the initial test set there were 32919 'otherwise' trackings in the 
%switch statement (i.e. not chr?) - I should have a quick look to see what
%they are (although a quick look is largely impossible)


%% Now do the plus strand


cramer_Raji_NETseq_plus_2_readin.textdata = ...
    cramer_Raji_NETseq_plus_2_readin.textdata(2:end);

ow2_counter = 0;

for ictr=1:1:length(cramer_Raji_NETseq_plus_2_readin.textdata)
   
    switch cramer_Raji_NETseq_plus_2_readin.textdata{ictr}
        
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
    
    %cramer_Raji_NETseq_2{chr,2}( ...
    %    cramer_Raji_NETseq_minus_2_readin.data(ictr,1):...
    %    cramer_Raji_NETseq_minus_2_readin.data(ictr,2)-1) = ...
    %    cramer_Raji_NETseq_minus_2_readin.data(ictr,3);
    
    cramer_Raji_NETseq_2{chr,1}( ...
        cramer_Raji_NETseq_plus_2_readin.data(ictr,1)+1:...
        cramer_Raji_NETseq_plus_2_readin.data(ictr,2) ) = ...
        cramer_Raji_NETseq_plus_2_readin.data(ictr,3);
    %cramer data starts from zero! - Need to make sure what the actual
    %positions are here!
    
    if mod(ictr,10000) == 0
        disp(['Completed ' num2str(ictr/10000) ' 2 plus iterations.'])
    end
    
end


%%

save([dataloaddir 'cramer_Raji_NETseq_2.mat'],'cramer_Raji_NETseq_2','-v7.3')




%% Notes

%The way Tom treated data may differ from the source file type
%Just make sure that everything is compatible at the end point


%The chromosome lengths appear to differ (from script to script)
% Perhaps this is why the try is necessary? Are some reads going over the
% edge?
%Perhaps I can check with Phil what the most appropriate annotation is

%Perhaps check what the largest thing is in the original data file

%Also will need to sort out the starting from zero problem - otherwise all
%my reads will be 1bp off (and could just go over the end of the
%chromosome)
