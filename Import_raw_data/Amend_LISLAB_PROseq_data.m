LIS_Lab_PROseq_WT_amended{1,1} = zeros(230218,1)
LIS_Lab_PROseq_WT_amended{1,2} = zeros(230218,1)
LIS_Lab_PROseq_WT_amended{2,1} = zeros(813184,1)
LIS_Lab_PROseq_WT_amended{2,2} = zeros(813184,1)
LIS_Lab_PROseq_WT_amended{3,1} = zeros(316620,1)
LIS_Lab_PROseq_WT_amended{3,2} = zeros(316620,1)
LIS_Lab_PROseq_WT_amended{4,1} = zeros(1531933,1)
LIS_Lab_PROseq_WT_amended{4,2} = zeros(1531933,1)
LIS_Lab_PROseq_WT_amended{5,1} = zeros(576874,1)
LIS_Lab_PROseq_WT_amended{5,2} = zeros(576874,1)
LIS_Lab_PROseq_WT_amended{6,1} = zeros(270161,1)
LIS_Lab_PROseq_WT_amended{6,2} = zeros(270161,1)
LIS_Lab_PROseq_WT_amended{7,1} = zeros(1090940,1)
LIS_Lab_PROseq_WT_amended{7,2} = zeros(1090940,1)
LIS_Lab_PROseq_WT_amended{8,1} = zeros(562643,1)
LIS_Lab_PROseq_WT_amended{8,2} = zeros(562643,1)
LIS_Lab_PROseq_WT_amended{9,1} = zeros(439888,1)
LIS_Lab_PROseq_WT_amended{9,2} = zeros(439888,1)
LIS_Lab_PROseq_WT_amended{10,1} = zeros(745751,1)
LIS_Lab_PROseq_WT_amended{10,2} = zeros(745751,1)
LIS_Lab_PROseq_WT_amended{11,1} = zeros(666816,1)
LIS_Lab_PROseq_WT_amended{11,2} = zeros(666816,1)
LIS_Lab_PROseq_WT_amended{12,1} = zeros(1078177,1)
LIS_Lab_PROseq_WT_amended{12,2} = zeros(1078177,1)
LIS_Lab_PROseq_WT_amended{13,1} = zeros(924431,1)
LIS_Lab_PROseq_WT_amended{13,2} = zeros(924431,1)
LIS_Lab_PROseq_WT_amended{14,1} = zeros(784333,1)
LIS_Lab_PROseq_WT_amended{14,2} = zeros(784333,1)
LIS_Lab_PROseq_WT_amended{15,1} = zeros(1091291,1)
LIS_Lab_PROseq_WT_amended{15,2} = zeros(1091291,1)
LIS_Lab_PROseq_WT_amended{16,1} = zeros(948066,1)
LIS_Lab_PROseq_WT_amended{16,2} = zeros(948066,1)

dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/LIS_PRO_SEQ/LisLab_PROseq_data/';

PROseq_sample = importdata([dataloaddir 'Lis_PROseq_WT.mat']);
%%
for cctr = 1:1:16
    for sctr = 1:1:2
        for pctr = 1:1:size(PROseq_sample{cctr,sctr})
        if PROseq_sample{cctr,sctr}(pctr,1)  ~= 0 
            LIS_Lab_PROseq_WT_amended{cctr,sctr}(pctr,1) = PROseq_sample{cctr,sctr}(pctr,1);
        end
        end
    end
end

%% TEST 
LIS_Lab_PROseq_WT_amended{14,1}(781576, 1)
for cctr = 1:1:16
    for sctr = 1:1:2
        for pctr = 1:1:size(PROseq_sample{cctr,sctr})
        if PROseq_sample{cctr,sctr}(pctr,1)  ~= LIS_Lab_PROseq_WT_amended{cctr,sctr}(pctr,1)
            disp('something is wrong, I can feel it')
        end
        end
    end
end

%% save output
save(append(dataloaddir, 'LIS_PROseq_WT_amended.mat'), 'LIS_Lab_PROseq_spt4del_amended')

%% Make frame based on chromosome size on the SGD
LIS_Lab_PROseq_spt4del_amended{1,1} = zeros(230218,1)
LIS_Lab_PROseq_spt4del_amended{1,2} = zeros(230218,1)
LIS_Lab_PROseq_spt4del_amended{2,1} = zeros(813184,1)
LIS_Lab_PROseq_spt4del_amended{2,2} = zeros(813184,1)
LIS_Lab_PROseq_spt4del_amended{3,1} = zeros(316620,1)
LIS_Lab_PROseq_spt4del_amended{3,2} = zeros(316620,1)
LIS_Lab_PROseq_spt4del_amended{4,1} = zeros(1531933,1)
LIS_Lab_PROseq_spt4del_amended{4,2} = zeros(1531933,1)
LIS_Lab_PROseq_spt4del_amended{5,1} = zeros(576874,1)
LIS_Lab_PROseq_spt4del_amended{5,2} = zeros(576874,1)
LIS_Lab_PROseq_spt4del_amended{6,1} = zeros(270161,1)
LIS_Lab_PROseq_spt4del_amended{6,2} = zeros(270161,1)
LIS_Lab_PROseq_spt4del_amended{7,1} = zeros(1090940,1)
LIS_Lab_PROseq_spt4del_amended{7,2} = zeros(1090940,1)
LIS_Lab_PROseq_spt4del_amended{8,1} = zeros(562643,1)
LIS_Lab_PROseq_spt4del_amended{8,2} = zeros(562643,1)
LIS_Lab_PROseq_spt4del_amended{9,1} = zeros(439888,1)
LIS_Lab_PROseq_spt4del_amended{9,2} = zeros(439888,1)
LIS_Lab_PROseq_spt4del_amended{10,1} = zeros(745751,1)
LIS_Lab_PROseq_spt4del_amended{10,2} = zeros(745751,1)
LIS_Lab_PROseq_spt4del_amended{11,1} = zeros(666816,1)
LIS_Lab_PROseq_spt4del_amended{11,2} = zeros(666816,1)
LIS_Lab_PROseq_spt4del_amended{12,1} = zeros(1078177,1)
LIS_Lab_PROseq_spt4del_amended{12,2} = zeros(1078177,1)
LIS_Lab_PROseq_spt4del_amended{13,1} = zeros(924431,1)
LIS_Lab_PROseq_spt4del_amended{13,2} = zeros(924431,1)
LIS_Lab_PROseq_spt4del_amended{14,1} = zeros(784333,1)
LIS_Lab_PROseq_spt4del_amended{14,2} = zeros(784333,1)
LIS_Lab_PROseq_spt4del_amended{15,1} = zeros(1091291,1)
LIS_Lab_PROseq_spt4del_amended{15,2} = zeros(1091291,1)
LIS_Lab_PROseq_spt4del_amended{16,1} = zeros(948066,1)
LIS_Lab_PROseq_spt4del_amended{16,2} = zeros(948066,1)

dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/LIS_PRO_SEQ/LisLab_PROseq_data/';

PROseq_sample1 = importdata([dataloaddir 'Lis_PROseq_spt4.mat']);
%%
for cctr = 1:1:16
    for sctr = 1:1:2
        for pctr = 1:1:size(PROseq_sample1{cctr,sctr})
        if PROseq_sample1{cctr,sctr}(pctr,1)  ~= 0 
            LIS_Lab_PROseq_spt4del_amended{cctr,sctr}(pctr,1) = PROseq_sample1{cctr,sctr}(pctr,1);
        end
        end
    end
end

%% TEST 
LIS_Lab_PROseq_spt4del_amended{14,1}(781576, 1)
for cctr = 1:1:16
    for sctr = 1:1:2
        for pctr = 1:1:size(PROseq_sample1{cctr,sctr})
        if PROseq_sample1{cctr,sctr}(pctr,1)  ~= LIS_Lab_PROseq_spt4del_amended{cctr,sctr}(pctr,1)
            disp('something is wrong, I can feel it')
        end
        end
    end
end

%%
save(append(dataloaddir, 'LIS_PROseq_spt4del_amended.mat'), 'LIS_Lab_PROseq_spt4del_amended')