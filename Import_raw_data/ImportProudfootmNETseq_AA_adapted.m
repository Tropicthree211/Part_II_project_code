%% Import Proudfoot mNETseq data
% This version written by Andrew (to hopefully deal with bigwig files)
% Adapted from Andrew's import cramer NETseq script 
% (which was in turn adapted from one of Tom's scripts)
% Note it deals only with wig files - the conversion from bigwig
% is most easily done externally

dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/PROUDFOOT_NET_SEQ/';

%%

proudfoot_HeLa_mNETseq_minus_readin = ...
    importdata([dataloaddir 'GSM1474229_ANET_CMA301_rep1_R.wig']);

proudfoot_HeLa_mNETseq_plus_readin = ...
    importdata([dataloaddir 'GSM1474229_ANET_CMA301_rep1_F.wig']);    



%%

proudfoot_HeLa_mNETseq = cell(24,2);

for ictr=1:2
    proudfoot_HeLa_mNETseq{1,ictr} = sparse(zeros(248956422,1));
    proudfoot_HeLa_mNETseq{2,ictr} = sparse(zeros(242193529,1));
    proudfoot_HeLa_mNETseq{3,ictr} = sparse(zeros(198295559,1));
    proudfoot_HeLa_mNETseq{4,ictr} = sparse(zeros(190214555,1));
    proudfoot_HeLa_mNETseq{5,ictr} = sparse(zeros(181538259,1));
    proudfoot_HeLa_mNETseq{6,ictr} = sparse(zeros(170805979,1));
    proudfoot_HeLa_mNETseq{7,ictr} = sparse(zeros(159345973,1));
    proudfoot_HeLa_mNETseq{8,ictr} = sparse(zeros(145138636,1));
    proudfoot_HeLa_mNETseq{9,ictr} = sparse(zeros(138394717,1));
    proudfoot_HeLa_mNETseq{10,ictr} = sparse(zeros(133797422,1));
    proudfoot_HeLa_mNETseq{11,ictr} = sparse(zeros(135086622,1));
    proudfoot_HeLa_mNETseq{12,ictr} = sparse(zeros(133275309,1));
    proudfoot_HeLa_mNETseq{13,ictr} = sparse(zeros(114364328,1));
    proudfoot_HeLa_mNETseq{14,ictr} = sparse(zeros(107043718,1));
    proudfoot_HeLa_mNETseq{15,ictr} = sparse(zeros(101991189,1));
    proudfoot_HeLa_mNETseq{16,ictr} = sparse(zeros(90338345,1));
    proudfoot_HeLa_mNETseq{17,ictr} = sparse(zeros(83257441,1));
    proudfoot_HeLa_mNETseq{18,ictr} = sparse(zeros(80373285,1));
    proudfoot_HeLa_mNETseq{19,ictr} = sparse(zeros(58617616,1));
    proudfoot_HeLa_mNETseq{20,ictr} = sparse(zeros(64444167,1));
    proudfoot_HeLa_mNETseq{21,ictr} = sparse(zeros(46709983,1));
    proudfoot_HeLa_mNETseq{22,ictr} = sparse(zeros(50818468,1));
    proudfoot_HeLa_mNETseq{23,ictr} = sparse(zeros(156040895,1));
    proudfoot_HeLa_mNETseq{24,ictr} = sparse(zeros(57227415,1));
end

%%

proudfoot_HeLa_mNETseq_minus_readin.textdata = ...
    proudfoot_HeLa_mNETseq_minus_readin.textdata(2:end);

ow_counter = 0;

for ictr=1:1:length(proudfoot_HeLa_mNETseq_minus_readin.textdata)
    
    switch proudfoot_HeLa_mNETseq_minus_readin.textdata{ictr}
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
    
    proudfoot_HeLa_mNETseq{chr,2}( ...
        proudfoot_HeLa_mNETseq_minus_readin.data(ictr,1): ...
        proudfoot_HeLa_mNETseq_minus_readin.data(ictr,2) ) = ...
        proudfoot_HeLa_mNETseq_minus_readin.data(ictr,3);
    %does this data start from zero? - I guess a warning/error will
    %be thrown if so - it did come up with the Cramer data set
    
    if mod(ictr,10000) == 0
        disp(['Completed ' num2str(ictr/10000) ' R strand iterations.'])
    end
       
    
end



%%

%%

%save([dataloaddir 'proudfoot_temphalf_reverse.mat'],'proudfoot_HeLa_mNETseq','-v7.3')




%% Now do the forward strand



proudfoot_HeLa_mNETseq_plus_readin.textdata = ...
    proudfoot_HeLa_mNETseq_plus_readin.textdata(2:end);

ow2_counter = 0;

for ictr=1:1:length(proudfoot_HeLa_mNETseq_plus_readin.textdata)
    
    switch proudfoot_HeLa_mNETseq_plus_readin.textdata{ictr}
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
    
    proudfoot_HeLa_mNETseq{chr,1}( ...
        proudfoot_HeLa_mNETseq_plus_readin.data(ictr,1): ...
        proudfoot_HeLa_mNETseq_plus_readin.data(ictr,2) ) = ...
        proudfoot_HeLa_mNETseq_plus_readin.data(ictr,3);
    %does this data start from zero? - I guess a warning/error will
    %be thrown if so - it did come up with the Cramer data set
    
    if mod(ictr,10000) == 0
        disp(['Completed ' num2str(ictr/10000) ' F strand iterations.'])
    end
       
    
end


%%

save([dataloaddir 'proudfoot_HeLa_mNETseq.mat'],'proudfoot_HeLa_mNETseq','-v7.3')

