%% Import Weissman's NETseq

%GM 02/02/21 Have commented out sections for mutants and replaced the ulku_WT with the
%set1 weissman data. For weissman other data can replace some of the
%commented mutants. Hopefully the file structure is the same otherwise I'll
%have to get really into the nitty gritty of the code. 


%Including the newer AnchorAway data

%Adapted from ImportCramerNETseq_K562_AA_speedup.m

%Looking at the notag normalised version of the data 
%Not sure if the anchor away has been similarly treated - need to ask ulku

%90 min anchor away was done as a single repeat so will not be analysed
%in the initial pass but could consider looking at it later

%% Start by testing the WT NETseq

dataloaddir = '/data/STUDENT_PLAY_AREA/TRANSCRIPTION_DATA/WEISSMAN_NET_SEQ/';
%dataloaddir_end = '';




NoChromosomes = 16; %Read from the data (roman numerals used)GM think no roman numerals used in this file  
                  
Chromosome_Names = {'chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI'};

                    

 %% Initial textscan to show only 16 chromosome labels 
% 
% fileID1 = fopen([dataloaddir 'GSM617033_SET1D_minus.wig']);
% weissman_set1del_NETseq_minus_1_readin = textscan(fileID1, '%s', ... %Initial textscan with only one string specified
%     'CommentStyle','track')
% fclose(fileID1);
% whos weissman_set1del_NETseq_minus_1_readin
%%
%strfind(weissman_set1del_NETseq_minus_1_readin{2}{:}, Chromosome_Names)
% %% Find the boundaries between chromosomes
% Chromosomepos = zeros((NoChromosomes+1), 1) ; 
% cctr = 0 ;
% for k = 1:1:size(weissman_set1del_NETseq_minus_1_readin{1,1})
%     
%     if weissman_set1del_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'v' 
%     %if weissman_set1del_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'c'
%         cctr = cctr +1;
%         Chromosomepos(cctr, 1) = k
%     end
%     
% end
%%
% Chromosomepos((NoChromosomes+1),1) = size(weissman_set1del_NETseq_minus_1_readin{1,1}, 1) + 1 
% %pull out data from whole textscan by using a for loop at every other
% %position?
% for pctr = (Chromosomepos(cctr,1)+2):2:(Chromosomepos((cctr+1),1)-1)
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
% end
    
    
%% second textscan to find chromosome positions 
fileID1 = fopen([dataloaddir 'GSM617033_SET1D_minus.wig']);
weissman_set1del_NETseq_minus_1_readin = textscan(fileID1, '%s %s', ... %two strings 
    'CommentStyle','track')
fclose(fileID1);
whos weissman_set1del_NETseq_minus_1_readin
%%
%strfind(weissman_set1del_NETseq_minus_1_readin{2}{:}, Chromosome_Names)
%% Find the boundaries between chromosomes
Chromosomepos = zeros((NoChromosomes+1), 1) ; 
cctr = 0 ;
for k = 1:1:size(weissman_set1del_NETseq_minus_1_readin{1,1})
    
    if weissman_set1del_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'v' 
    %if weissman_set1del_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'c'
        cctr = cctr +1;
        Chromosomepos(cctr, 1) = k
    end
   
    
end
Chromosomepos((NoChromosomes+1),1) = size(weissman_set1del_NETseq_minus_1_readin{1,1}(:), 1) + 1
%% Looped textscan over each chromosome between the strings to pull out floating values 
% weissman_set1del_NETseq_sample_nozeros = cell(16,2);
% % for cctr = 1:1:NoChromosomes
% %     for sctr = 1:1:2 
% %preallocate the size of the cell array before adding in data, done via for
% %loop hopefully won't be too slow 
%        for cctr = 1:1:NoChromosomes
%            weissman_set1del_NETseq_sample_nozeros{cctr,2} = zeros(size((Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)))
%        end
%  %this gives correctly sized numeric cell array to receive data. Need to make data numeric to input.    
%    %w=0
   
  %% following should target and convert all the data between the chromosome labels from string to numeric
   %converts positional data
   numeric_weissman_set1del_NETseq_minus_1_readin = zeros(size(weissman_set1del_NETseq_minus_1_readin{1,1}, 1), 2);
for cctr = 1:1:NoChromosomes 
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_set1del_NETseq_minus_1_readin(pctr,1) = str2double(weissman_set1del_NETseq_minus_1_readin{1,1}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end
%converts signal data 
for cctr = 1:1:NoChromosomes
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_set1del_NETseq_minus_1_readin(pctr, 2) = str2double(weissman_set1del_NETseq_minus_1_readin{1,2}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end   
%% check there are only 16 zeros in the dataset where the labels were 
sum(numeric_weissman_set1del_NETseq_minus_1_readin == 0) %It has worked!!! 


%If this works will now need to input the zeros where chromosome labels
%were before now have zeros, should be okay because I can still selectively
%target between the zeros by pctr. 

%% Make the sample structure with the right number of zeros 
weissman_set1del_NETseq_sample = cell(16,2) %need to give the correct number of zeros according to the SGD? 
%Then can input the data into the zeros to create the og NETseq signal.
weissman_set1del_NETseq_sample{1,1} = zeros(230218,1)
weissman_set1del_NETseq_sample{1,2} = zeros(230218,1)
weissman_set1del_NETseq_sample{2,1} = zeros(813184,1)
weissman_set1del_NETseq_sample{2,2} = zeros(813184,1)
weissman_set1del_NETseq_sample{3,1} = zeros(316620,1)
weissman_set1del_NETseq_sample{3,2} = zeros(316620,1)
weissman_set1del_NETseq_sample{4,1} = zeros(1531933,1)
weissman_set1del_NETseq_sample{4,2} = zeros(1531933,1)
weissman_set1del_NETseq_sample{5,1} = zeros(576874,1)
weissman_set1del_NETseq_sample{5,2} = zeros(576874,1)
weissman_set1del_NETseq_sample{6,1} = zeros(270161,1)
weissman_set1del_NETseq_sample{6,2} = zeros(270161,1)
weissman_set1del_NETseq_sample{7,1} = zeros(1090940,1)
weissman_set1del_NETseq_sample{7,2} = zeros(1090940,1)
weissman_set1del_NETseq_sample{8,1} = zeros(562643,1)
weissman_set1del_NETseq_sample{8,2} = zeros(562643,1)
weissman_set1del_NETseq_sample{9,1} = zeros(439888,1)
weissman_set1del_NETseq_sample{9,2} = zeros(439888,1)
weissman_set1del_NETseq_sample{10,1} = zeros(745751,1)
weissman_set1del_NETseq_sample{10,2} = zeros(745751,1)
weissman_set1del_NETseq_sample{11,1} = zeros(666816,1)
weissman_set1del_NETseq_sample{11,2} = zeros(666816,1)
weissman_set1del_NETseq_sample{12,1} = zeros(1078177,1)
weissman_set1del_NETseq_sample{12,2} = zeros(1078177,1)
weissman_set1del_NETseq_sample{13,1} = zeros(924431,1)
weissman_set1del_NETseq_sample{13,2} = zeros(924431,1)
weissman_set1del_NETseq_sample{14,1} = zeros(784333,1)
weissman_set1del_NETseq_sample{14,2} = zeros(784333,1)
weissman_set1del_NETseq_sample{15,1} = zeros(1091291,1)
weissman_set1del_NETseq_sample{15,2} = zeros(1091291,1)
weissman_set1del_NETseq_sample{16,1} = zeros(948066,1)
weissman_set1del_NETseq_sample{16,2} = zeros(948066,1)

%should create correctly sized matrix in which to insert the NETseq data. 

%% Inserting the NETseq signal into the correctly sized array of zeros.
for cctr = 1:1:NoChromosomes
    
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
        weissman_set1del_NETseq_sample{cctr, 2}(numeric_weissman_set1del_NETseq_minus_1_readin(pctr,1), 1) =...
            numeric_weissman_set1del_NETseq_minus_1_readin(pctr,2);
    end
end

%% Repeat above for plus strand 

% Initial textscan to show only 16 chromosome labels 

% fileID1 = fopen([dataloaddir 'GSM617033_SET1D_plus.wig']);
% weissman_set1del_NETseq_plus_1_readin = textscan(fileID1, '%s', ... %Initial textscan with only one string specified
%     'CommentStyle','track')
% fclose(fileID1);
% whos weissman_set1del_NETseq_plus_1_readin
% 
% %% Find the boundaries between chromosomes
% Chromosomepos = zeros((NoChromosomes+1), 1) ; 
% cctr = 0 ;
% for k = 1:1:size(weissman_set1del_NETseq_plus_1_readin{1,1})
%     
%     if weissman_set1del_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'v' 
%     %if weissman_set1del_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'c'
%         cctr = cctr +1;
%         Chromosomepos(cctr, 1) = k
%     end
%     
% end
%%
% Chromosomepos((NoChromosomes+1),1) = size(weissman_set1del_NETseq_plus_1_readin{1,1}, 1) + 1 
% %pull out data from whole textscan by using a for loop at every other
% %position?
% for pctr = (Chromosomepos(cctr,1)+2):2:(Chromosomepos((cctr+1),1)-1)
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
% end
    
    
%% second textscan to find chromosome positions 
fileID2 = fopen([dataloaddir 'GSM617033_SET1D_plus.wig']);
weissman_set1del_NETseq_plus_1_readin = textscan(fileID2, '%s %s', ... %two strings 
    'CommentStyle','track')
fclose(fileID2);
whos weissman_set1del_NETseq_plus_1_readin
%%
%strfind(weissman_set1del_NETseq_plus_1_readin{2}{:}, Chromosome_Names)
%% Find the boundaries between chromosomes
Chromosomepos = zeros((NoChromosomes+1), 1) ; 
cctr = 0 ;
for k = 1:1:size(weissman_set1del_NETseq_plus_1_readin{1,1})
    
    if weissman_set1del_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'v' 
    %if weissman_set1del_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'c'
        cctr = cctr +1;
        Chromosomepos(cctr, 1) = k
    end
   
    
end
Chromosomepos((NoChromosomes+1),1) = size(weissman_set1del_NETseq_plus_1_readin{1,1}(:), 1) + 1
%% Looped textscan over each chromosome between the strings to pull out floating values 
% weissman_set1del_NETseq_sample_nozeros = cell(16,2);
% % for cctr = 1:1:NoChromosomes
% %     for sctr = 1:1:2 
% %preallocate the size of the cell array before adding in data, done via for
% %loop hopefully won't be too slow 
%        for cctr = 1:1:NoChromosomes
%            weissman_set1del_NETseq_sample_nozeros{cctr,2} = zeros(size((Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)))
%        end
%  %this gives correctly sized numeric cell array to receive data. Need to make data numeric to input.    
%    %w=0
   
  %% following should target and convert all the data between the chromosome labels from string to numeric
   %converts positional data
   numeric_weissman_set1del_NETseq_plus_1_readin = zeros(size(weissman_set1del_NETseq_plus_1_readin{1,1}, 1), 2);
for cctr = 1:1:NoChromosomes 
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_set1del_NETseq_plus_1_readin(pctr,1) = str2double(weissman_set1del_NETseq_plus_1_readin{1,1}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end
%converts signal data 
for cctr = 1:1:NoChromosomes
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_set1del_NETseq_plus_1_readin(pctr, 2) = str2double(weissman_set1del_NETseq_plus_1_readin{1,2}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end   
%% check there are only 16 zeros in the dataset where the labels were 
sum(numeric_weissman_set1del_NETseq_plus_1_readin == 0) %It has worked!!! 


%If this works will now need to input the zeros where chromosome labels
%were before now have zeros, should be okay because I can still selectively
%target between the zeros by pctr. 


%% Inserting the NETseq signal into the correctly sized array of zeros.
for cctr = 1:1:NoChromosomes
    
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
        weissman_set1del_NETseq_sample{cctr, 1}(numeric_weissman_set1del_NETseq_plus_1_readin(pctr,1), 1) =...
            numeric_weissman_set1del_NETseq_plus_1_readin(pctr,2);
    end
end

%% Save the output 

save(append(dataloaddir, 'weissman_set1del_NETseq.mat'), 'weissman_set1del_NETseq_sample')
%%
        




















%% Initial textscan to show only 16 chromosome labels SET2DEL

fileID1 = fopen([dataloaddir 'GSM617032_SET2D_minus.wig']);
weissman_set2del_NETseq_minus_1_readin = textscan(fileID1, '%s', ... %Initial textscan with only one string specified
    'CommentStyle','track')
fclose(fileID1);
whos weissman_set2del_NETseq_minus_1_readin
%%
%strfind(weissman_set2del_NETseq_minus_1_readin{2}{:}, Chromosome_Names)
%% Find the boundaries between chromosomes
Chromosomepos = zeros((NoChromosomes+1), 1) ; 
cctr = 0 ;
for k = 1:1:size(weissman_set2del_NETseq_minus_1_readin{1,1})
    
    if weissman_set2del_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'v' 
    %if weissman_set2del_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'c'
        cctr = cctr +1;
        Chromosomepos(cctr, 1) = k
    end
    
end
%%
% Chromosomepos((NoChromosomes+1),1) = size(weissman_set2del_NETseq_minus_1_readin{1,1}, 1) + 1 
% %pull out data from whole textscan by using a for loop at every other
% %position?
% for pctr = (Chromosomepos(cctr,1)+2):2:(Chromosomepos((cctr+1),1)-1)
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
% end
    
    
%% second textscan to find chromosome positions 
fileID1 = fopen([dataloaddir 'GSM617032_SET2D_minus.wig']);
weissman_set2del_NETseq_minus_1_readin = textscan(fileID1, '%s %s', ... %two strings 
    'CommentStyle','track')
fclose(fileID1);
whos weissman_set2del_NETseq_minus_1_readin
%%
%strfind(weissman_set2del_NETseq_minus_1_readin{2}{:}, Chromosome_Names)
%% Find the boundaries between chromosomes
Chromosomepos = zeros((NoChromosomes+1), 1) ; 
cctr = 0 ;
for k = 1:1:size(weissman_set2del_NETseq_minus_1_readin{1,1})
    
    if weissman_set2del_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'v' 
    %if weissman_set2del_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'c'
        cctr = cctr +1;
        Chromosomepos(cctr, 1) = k
    end
   
    
end
Chromosomepos((NoChromosomes+1),1) = size(weissman_set2del_NETseq_minus_1_readin{1,1}(:), 1) + 1
%% Looped textscan over each chromosome between the strings to pull out floating values 
% weissman_set2del_NETseq_sample_nozeros = cell(16,2);
% % for cctr = 1:1:NoChromosomes
% %     for sctr = 1:1:2 
% %preallocate the size of the cell array before adding in data, done via for
% %loop hopefully won't be too slow 
%        for cctr = 1:1:NoChromosomes
%            weissman_set2del_NETseq_sample_nozeros{cctr,2} = zeros(size((Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)))
%        end
%  %this gives correctly sized numeric cell array to receive data. Need to make data numeric to input.    
%    %w=0
   
  %% following should target and convert all the data between the chromosome labels from string to numeric
   %converts positional data
   numeric_weissman_set2del_NETseq_minus_1_readin = zeros(size(weissman_set2del_NETseq_minus_1_readin{1,1}, 1), 2);
for cctr = 1:1:NoChromosomes 
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_set2del_NETseq_minus_1_readin(pctr,1) = str2double(weissman_set2del_NETseq_minus_1_readin{1,1}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end
%converts signal data 
for cctr = 1:1:NoChromosomes
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_set2del_NETseq_minus_1_readin(pctr, 2) = str2double(weissman_set2del_NETseq_minus_1_readin{1,2}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end   
%% check there are only 16 zeros in the dataset where the labels were 
sum(numeric_weissman_set2del_NETseq_minus_1_readin == 0) %It has worked!!! 


%If this works will now need to input the zeros where chromosome labels
%were before now have zeros, should be okay because I can still selectively
%target between the zeros by pctr. 

%% Make the sample structure with the right number of zeros 
weissman_set2del_NETseq_sample = cell(16,2) %need to give the correct number of zeros according to the SGD? 
%Then can input the data into the zeros to create the og NETseq signal.
weissman_set2del_NETseq_sample{1,1} = zeros(230218,1)
weissman_set2del_NETseq_sample{1,2} = zeros(230218,1)
weissman_set2del_NETseq_sample{2,1} = zeros(813184,1)
weissman_set2del_NETseq_sample{2,2} = zeros(813184,1)
weissman_set2del_NETseq_sample{3,1} = zeros(316620,1)
weissman_set2del_NETseq_sample{3,2} = zeros(316620,1)
weissman_set2del_NETseq_sample{4,1} = zeros(1531933,1)
weissman_set2del_NETseq_sample{4,2} = zeros(1531933,1)
weissman_set2del_NETseq_sample{5,1} = zeros(576874,1)
weissman_set2del_NETseq_sample{5,2} = zeros(576874,1)
weissman_set2del_NETseq_sample{6,1} = zeros(270161,1)
weissman_set2del_NETseq_sample{6,2} = zeros(270161,1)
weissman_set2del_NETseq_sample{7,1} = zeros(1090940,1)
weissman_set2del_NETseq_sample{7,2} = zeros(1090940,1)
weissman_set2del_NETseq_sample{8,1} = zeros(562643,1)
weissman_set2del_NETseq_sample{8,2} = zeros(562643,1)
weissman_set2del_NETseq_sample{9,1} = zeros(439888,1)
weissman_set2del_NETseq_sample{9,2} = zeros(439888,1)
weissman_set2del_NETseq_sample{10,1} = zeros(745751,1)
weissman_set2del_NETseq_sample{10,2} = zeros(745751,1)
weissman_set2del_NETseq_sample{11,1} = zeros(666816,1)
weissman_set2del_NETseq_sample{11,2} = zeros(666816,1)
weissman_set2del_NETseq_sample{12,1} = zeros(1078177,1)
weissman_set2del_NETseq_sample{12,2} = zeros(1078177,1)
weissman_set2del_NETseq_sample{13,1} = zeros(924431,1)
weissman_set2del_NETseq_sample{13,2} = zeros(924431,1)
weissman_set2del_NETseq_sample{14,1} = zeros(784333,1)
weissman_set2del_NETseq_sample{14,2} = zeros(784333,1)
weissman_set2del_NETseq_sample{15,1} = zeros(1091291,1)
weissman_set2del_NETseq_sample{15,2} = zeros(1091291,1)
weissman_set2del_NETseq_sample{16,1} = zeros(948066,1)
weissman_set2del_NETseq_sample{16,2} = zeros(948066,1)

%should create correctly sized matrix in which to insert the NETseq data. 

%% Inserting the NETseq signal into the correctly sized array of zeros.
for cctr = 1:1:NoChromosomes
    
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
        weissman_set2del_NETseq_sample{cctr, 2}(numeric_weissman_set2del_NETseq_minus_1_readin(pctr,1), 1) =...
            numeric_weissman_set2del_NETseq_minus_1_readin(pctr,2);
    end
end

%% Repeat above for plus strand 

% Initial textscan to show only 16 chromosome labels 

fileID1 = fopen([dataloaddir 'GSM617032_SET2D_plus.wig']);
weissman_set2del_NETseq_plus_1_readin = textscan(fileID1, '%s', ... %Initial textscan with only one string specified
    'CommentStyle','track')
fclose(fileID1);
whos weissman_set2del_NETseq_plus_1_readin

%% Find the boundaries between chromosomes
Chromosomepos = zeros((NoChromosomes+1), 1) ; 
cctr = 0 ;
for k = 1:1:size(weissman_set2del_NETseq_plus_1_readin{1,1})
    
    if weissman_set2del_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'v' 
    %if weissman_set2del_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'c'
        cctr = cctr +1;
        Chromosomepos(cctr, 1) = k
    end
    
end
%%
% Chromosomepos((NoChromosomes+1),1) = size(weissman_set2del_NETseq_plus_1_readin{1,1}, 1) + 1 
% %pull out data from whole textscan by using a for loop at every other
% %position?
% for pctr = (Chromosomepos(cctr,1)+2):2:(Chromosomepos((cctr+1),1)-1)
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
% end
    
    
%% second textscan to find chromosome positions 
fileID2 = fopen([dataloaddir 'GSM617032_SET2D_plus.wig']);
weissman_set2del_NETseq_plus_1_readin = textscan(fileID2, '%s %s', ... %two strings 
    'CommentStyle','track')
fclose(fileID2);
whos weissman_set2del_NETseq_plus_1_readin
%%
%strfind(weissman_set2del_NETseq_plus_1_readin{2}{:}, Chromosome_Names)
%% Find the boundaries between chromosomes
Chromosomepos = zeros((NoChromosomes+1), 1) ; 
cctr = 0 ;
for k = 1:1:size(weissman_set2del_NETseq_plus_1_readin{1,1})
    
    if weissman_set2del_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'v' 
    %if weissman_set2del_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'c'
        cctr = cctr +1;
        Chromosomepos(cctr, 1) = k
    end
   
    
end
Chromosomepos((NoChromosomes+1),1) = size(weissman_set2del_NETseq_plus_1_readin{1,1}(:), 1) + 1
%% Looped textscan over each chromosome between the strings to pull out floating values 
% weissman_set2del_NETseq_sample_nozeros = cell(16,2);
% % for cctr = 1:1:NoChromosomes
% %     for sctr = 1:1:2 
% %preallocate the size of the cell array before adding in data, done via for
% %loop hopefully won't be too slow 
%        for cctr = 1:1:NoChromosomes
%            weissman_set2del_NETseq_sample_nozeros{cctr,2} = zeros(size((Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)))
%        end
%  %this gives correctly sized numeric cell array to receive data. Need to make data numeric to input.    
%    %w=0
   
  %% following should target and convert all the data between the chromosome labels from string to numeric
   %converts positional data
   numeric_weissman_set2del_NETseq_plus_1_readin = zeros(size(weissman_set2del_NETseq_plus_1_readin{1,1}, 1), 2);
for cctr = 1:1:NoChromosomes 
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_set2del_NETseq_plus_1_readin(pctr,1) = str2double(weissman_set2del_NETseq_plus_1_readin{1,1}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end
%converts signal data 
for cctr = 1:1:NoChromosomes
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_set2del_NETseq_plus_1_readin(pctr, 2) = str2double(weissman_set2del_NETseq_plus_1_readin{1,2}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end   
%% check there are only 16 zeros in the dataset where the labels were 
sum(numeric_weissman_set2del_NETseq_plus_1_readin == 0) %It has worked!!! 


%If this works will now need to input the zeros where chromosome labels
%were before now have zeros, should be okay because I can still selectively
%target between the zeros by pctr. 


%% Inserting the NETseq signal into the correctly sized array of zeros.
for cctr = 1:1:NoChromosomes
    
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
        weissman_set2del_NETseq_sample{cctr, 1}(numeric_weissman_set2del_NETseq_plus_1_readin(pctr,1), 1) =...
            numeric_weissman_set2del_NETseq_plus_1_readin(pctr,2);
    end
end

%% Save the output 

save(append(dataloaddir, 'weissman_set2del_NETseq.mat'), 'weissman_set2del_NETseq_sample')
%%








  %% Initial textscan to show only 16 chromosome labels EAF3DEL

fileID1 = fopen([dataloaddir 'GSM617031_EAF3D_minus.wig']);
weissman_eaf3del_NETseq_minus_1_readin = textscan(fileID1, '%s', ... %Initial textscan with only one string specified
    'CommentStyle','track')
fclose(fileID1);
whos weissman_eaf3del_NETseq_minus_1_readin
%%
%strfind(weissman_eaf3del_NETseq_minus_1_readin{2}{:}, Chromosome_Names)
%% Find the boundaries between chromosomes
Chromosomepos = zeros((NoChromosomes+1), 1) ; 
cctr = 0 ;
for k = 1:1:size(weissman_eaf3del_NETseq_minus_1_readin{1,1})
    
    if weissman_eaf3del_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'v' 
    %if weissman_eaf3del_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'c'
        cctr = cctr +1;
        Chromosomepos(cctr, 1) = k
    end
    
end
%%
% Chromosomepos((NoChromosomes+1),1) = size(weissman_eaf3del_NETseq_minus_1_readin{1,1}, 1) + 1 
% %pull out data from whole textscan by using a for loop at every other
% %position?
% for pctr = (Chromosomepos(cctr,1)+2):2:(Chromosomepos((cctr+1),1)-1)
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
% end
    
    
%% second textscan to find chromosome positions 
fileID1 = fopen([dataloaddir 'GSM617031_EAF3D_minus.wig']);
weissman_eaf3del_NETseq_minus_1_readin = textscan(fileID1, '%s %s', ... %two strings 
    'CommentStyle','track')
fclose(fileID1);
whos weissman_eaf3del_NETseq_minus_1_readin
%%
%strfind(weissman_eaf3del_NETseq_minus_1_readin{2}{:}, Chromosome_Names)
%% Find the boundaries between chromosomes
Chromosomepos = zeros((NoChromosomes+1), 1) ; 
cctr = 0 ;
for k = 1:1:size(weissman_eaf3del_NETseq_minus_1_readin{1,1})
    
    if weissman_eaf3del_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'v' 
    %if weissman_eaf3del_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'c'
        cctr = cctr +1;
        Chromosomepos(cctr, 1) = k
    end
   
    
end
Chromosomepos((NoChromosomes+1),1) = size(weissman_eaf3del_NETseq_minus_1_readin{1,1}(:), 1) + 1
%% Looped textscan over each chromosome between the strings to pull out floating values 
% weissman_eaf3del_NETseq_sample_nozeros = cell(16,2);
% % for cctr = 1:1:NoChromosomes
% %     for sctr = 1:1:2 
% %preallocate the size of the cell array before adding in data, done via for
% %loop hopefully won't be too slow 
%        for cctr = 1:1:NoChromosomes
%            weissman_eaf3del_NETseq_sample_nozeros{cctr,2} = zeros(size((Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)))
%        end
%  %this gives correctly sized numeric cell array to receive data. Need to make data numeric to input.    
%    %w=0
   
  %% following should target and convert all the data between the chromosome labels from string to numeric
   %converts positional data
   numeric_weissman_eaf3del_NETseq_minus_1_readin = zeros(size(weissman_eaf3del_NETseq_minus_1_readin{1,1}, 1), 2);
for cctr = 1:1:NoChromosomes 
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_eaf3del_NETseq_minus_1_readin(pctr,1) = str2double(weissman_eaf3del_NETseq_minus_1_readin{1,1}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end
%converts signal data 
for cctr = 1:1:NoChromosomes
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_eaf3del_NETseq_minus_1_readin(pctr, 2) = str2double(weissman_eaf3del_NETseq_minus_1_readin{1,2}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end   
%% check there are only 16 zeros in the dataset where the labels were 
sum(numeric_weissman_eaf3del_NETseq_minus_1_readin == 0) %It has worked!!! 


%If this works will now need to input the zeros where chromosome labels
%were before now have zeros, should be okay because I can still selectively
%target between the zeros by pctr. 

%% Make the sample structure with the right number of zeros 
weissman_eaf3del_NETseq_sample = cell(16,2) %need to give the correct number of zeros according to the SGD? 
%Then can input the data into the zeros to create the og NETseq signal.
weissman_eaf3del_NETseq_sample{1,1} = zeros(230218,1)
weissman_eaf3del_NETseq_sample{1,2} = zeros(230218,1)
weissman_eaf3del_NETseq_sample{2,1} = zeros(813184,1)
weissman_eaf3del_NETseq_sample{2,2} = zeros(813184,1)
weissman_eaf3del_NETseq_sample{3,1} = zeros(316620,1)
weissman_eaf3del_NETseq_sample{3,2} = zeros(316620,1)
weissman_eaf3del_NETseq_sample{4,1} = zeros(1531933,1)
weissman_eaf3del_NETseq_sample{4,2} = zeros(1531933,1)
weissman_eaf3del_NETseq_sample{5,1} = zeros(576874,1)
weissman_eaf3del_NETseq_sample{5,2} = zeros(576874,1)
weissman_eaf3del_NETseq_sample{6,1} = zeros(270161,1)
weissman_eaf3del_NETseq_sample{6,2} = zeros(270161,1)
weissman_eaf3del_NETseq_sample{7,1} = zeros(1090940,1)
weissman_eaf3del_NETseq_sample{7,2} = zeros(1090940,1)
weissman_eaf3del_NETseq_sample{8,1} = zeros(562643,1)
weissman_eaf3del_NETseq_sample{8,2} = zeros(562643,1)
weissman_eaf3del_NETseq_sample{9,1} = zeros(439888,1)
weissman_eaf3del_NETseq_sample{9,2} = zeros(439888,1)
weissman_eaf3del_NETseq_sample{10,1} = zeros(745751,1)
weissman_eaf3del_NETseq_sample{10,2} = zeros(745751,1)
weissman_eaf3del_NETseq_sample{11,1} = zeros(666816,1)
weissman_eaf3del_NETseq_sample{11,2} = zeros(666816,1)
weissman_eaf3del_NETseq_sample{12,1} = zeros(1078177,1)
weissman_eaf3del_NETseq_sample{12,2} = zeros(1078177,1)
weissman_eaf3del_NETseq_sample{13,1} = zeros(924431,1)
weissman_eaf3del_NETseq_sample{13,2} = zeros(924431,1)
weissman_eaf3del_NETseq_sample{14,1} = zeros(784333,1)
weissman_eaf3del_NETseq_sample{14,2} = zeros(784333,1)
weissman_eaf3del_NETseq_sample{15,1} = zeros(1091291,1)
weissman_eaf3del_NETseq_sample{15,2} = zeros(1091291,1)
weissman_eaf3del_NETseq_sample{16,1} = zeros(948066,1)
weissman_eaf3del_NETseq_sample{16,2} = zeros(948066,1)

%should create correctly sized matrix in which to insert the NETseq data. 

%% Inserting the NETseq signal into the correctly sized array of zeros.
for cctr = 1:1:NoChromosomes
    
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
        weissman_eaf3del_NETseq_sample{cctr, 2}(numeric_weissman_eaf3del_NETseq_minus_1_readin(pctr,1), 1) =...
            numeric_weissman_eaf3del_NETseq_minus_1_readin(pctr,2);
    end
end

%% Repeat above for plus strand 

% Initial textscan to show only 16 chromosome labels 

fileID1 = fopen([dataloaddir 'GSM617031_EAF3D_plus.wig']);
weissman_eaf3del_NETseq_plus_1_readin = textscan(fileID1, '%s', ... %Initial textscan with only one string specified
    'CommentStyle','track')
fclose(fileID1);
whos weissman_eaf3del_NETseq_plus_1_readin

%% Find the boundaries between chromosomes
Chromosomepos = zeros((NoChromosomes+1), 1) ; 
cctr = 0 ;
for k = 1:1:size(weissman_eaf3del_NETseq_plus_1_readin{1,1})
    
    if weissman_eaf3del_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'v' 
    %if weissman_eaf3del_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'c'
        cctr = cctr +1;
        Chromosomepos(cctr, 1) = k
    end
    
end
%%
% Chromosomepos((NoChromosomes+1),1) = size(weissman_eaf3del_NETseq_plus_1_readin{1,1}, 1) + 1 
% %pull out data from whole textscan by using a for loop at every other
% %position?
% for pctr = (Chromosomepos(cctr,1)+2):2:(Chromosomepos((cctr+1),1)-1)
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
% end
    
    
%% second textscan to find chromosome positions 
fileID2 = fopen([dataloaddir 'GSM617031_EAF3D_plus.wig']);
weissman_eaf3del_NETseq_plus_1_readin = textscan(fileID2, '%s %s', ... %two strings 
    'CommentStyle','track')
fclose(fileID2);
whos weissman_eaf3del_NETseq_plus_1_readin
%%
%strfind(weissman_eaf3del_NETseq_plus_1_readin{2}{:}, Chromosome_Names)
%% Find the boundaries between chromosomes
Chromosomepos = zeros((NoChromosomes+1), 1) ; 
cctr = 0 ;
for k = 1:1:size(weissman_eaf3del_NETseq_plus_1_readin{1,1})
    
    if weissman_eaf3del_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'v' 
    %if weissman_eaf3del_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'c'
        cctr = cctr +1;
        Chromosomepos(cctr, 1) = k
    end
   
    
end
Chromosomepos((NoChromosomes+1),1) = size(weissman_eaf3del_NETseq_plus_1_readin{1,1}(:), 1) + 1
%% Looped textscan over each chromosome between the strings to pull out floating values 
% weissman_eaf3del_NETseq_sample_nozeros = cell(16,2);
% % for cctr = 1:1:NoChromosomes
% %     for sctr = 1:1:2 
% %preallocate the size of the cell array before adding in data, done via for
% %loop hopefully won't be too slow 
%        for cctr = 1:1:NoChromosomes
%            weissman_eaf3del_NETseq_sample_nozeros{cctr,2} = zeros(size((Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)))
%        end
%  %this gives correctly sized numeric cell array to receive data. Need to make data numeric to input.    
%    %w=0
   
  %% following should target and convert all the data between the chromosome labels from string to numeric
   %converts positional data
   numeric_weissman_eaf3del_NETseq_plus_1_readin = zeros(size(weissman_eaf3del_NETseq_plus_1_readin{1,1}, 1), 2);
for cctr = 1:1:NoChromosomes 
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_eaf3del_NETseq_plus_1_readin(pctr,1) = str2double(weissman_eaf3del_NETseq_plus_1_readin{1,1}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end
%converts signal data 
for cctr = 1:1:NoChromosomes
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_eaf3del_NETseq_plus_1_readin(pctr, 2) = str2double(weissman_eaf3del_NETseq_plus_1_readin{1,2}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end   
%% check there are only 16 zeros in the dataset where the labels were 
sum(numeric_weissman_eaf3del_NETseq_plus_1_readin == 0) %It has worked!!! 


%If this works will now need to input the zeros where chromosome labels
%were before now have zeros, should be okay because I can still selectively
%target between the zeros by pctr. 


%% Inserting the NETseq signal into the correctly sized array of zeros.
for cctr = 1:1:NoChromosomes
    
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
        weissman_eaf3del_NETseq_sample{cctr, 1}(numeric_weissman_eaf3del_NETseq_plus_1_readin(pctr,1), 1) =...
            numeric_weissman_eaf3del_NETseq_plus_1_readin(pctr,2);
    end
end

%% Save the output 

save(append(dataloaddir, 'weissman_eaf3del_NETseq.mat'), 'weissman_eaf3del_NETseq_sample')      
        



 %% Initial textscan to show only 16 chromosome labels 

fileID1 = fopen([dataloaddir 'GSM617027_WT_NC_minus.wig']);
weissman_WT_NETseq_minus_1_readin = textscan(fileID1, '%s', ... %Initial textscan with only one string specified
    'CommentStyle','track')
fclose(fileID1);
whos weissman_WT_NETseq_minus_1_readin
%%
%strfind(weissman_WT_NETseq_minus_1_readin{2}{:}, Chromosome_Names)
% %% Find the boundaries between chromosomes
% Chromosomepos = zeros((NoChromosomes+1), 1) ; 
% cctr = 0 ;
% for k = 1:1:size(weissman_WT_NETseq_minus_1_readin{1,1})
%     
%     if weissman_WT_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'v' 
%     %if weissman_WT_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'c'
%         cctr = cctr +1;
%         Chromosomepos(cctr, 1) = k
%     end
%     
% end
%%
% Chromosomepos((NoChromosomes+1),1) = size(weissman_WT_NETseq_minus_1_readin{1,1}, 1) + 1 
% %pull out data from whole textscan by using a for loop at every other
% %position?
% for pctr = (Chromosomepos(cctr,1)+2):2:(Chromosomepos((cctr+1),1)-1)
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
% end
    
    
%% second textscan to find chromosome positions 
fileID1 = fopen([dataloaddir 'GSM617027_WT_NC_minus.wig']);
weissman_WT_NETseq_minus_1_readin = textscan(fileID1, '%s %s', ... %two strings 
    'CommentStyle','track')
fclose(fileID1);
whos weissman_WT_NETseq_minus_1_readin
%%
%strfind(weissman_WT_NETseq_minus_1_readin{2}{:}, Chromosome_Names)
%% Find the boundaries between chromosomes
Chromosomepos = zeros((NoChromosomes+1), 1) ; 
cctr = 0 ;
for k = 1:1:size(weissman_WT_NETseq_minus_1_readin{1,1})
    
    if weissman_WT_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'v' 
    %if weissman_WT_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'c'
        cctr = cctr +1;
        Chromosomepos(cctr, 1) = k
    end
   
    
end
Chromosomepos((NoChromosomes+1),1) = size(weissman_WT_NETseq_minus_1_readin{1,1}(:), 1) + 1
%% Looped textscan over each chromosome between the strings to pull out floating values 
% weissman_WT_NETseq_sample_nozeros = cell(16,2);
% % for cctr = 1:1:NoChromosomes
% %     for sctr = 1:1:2 
% %preallocate the size of the cell array before adding in data, done via for
% %loop hopefully won't be too slow 
%        for cctr = 1:1:NoChromosomes
%            weissman_WT_NETseq_sample_nozeros{cctr,2} = zeros(size((Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)))
%        end
%  %this gives correctly sized numeric cell array to receive data. Need to make data numeric to input.    
%    %w=0
   
  %% following should target and convert all the data between the chromosome labels from string to numeric
   %converts positional data
   numeric_weissman_WT_NETseq_minus_1_readin = zeros(size(weissman_WT_NETseq_minus_1_readin{1,1}, 1), 2);
for cctr = 1:1:NoChromosomes 
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_WT_NETseq_minus_1_readin(pctr,1) = str2double(weissman_WT_NETseq_minus_1_readin{1,1}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end
%converts signal data 
for cctr = 1:1:NoChromosomes
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_WT_NETseq_minus_1_readin(pctr, 2) = str2double(weissman_WT_NETseq_minus_1_readin{1,2}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end   
%% check there are only 16 zeros in the dataset where the labels were 
sum(numeric_weissman_WT_NETseq_minus_1_readin == 0) %It has worked!!! 


%If this works will now need to input the zeros where chromosome labels
%were before now have zeros, should be okay because I can still selectively
%target between the zeros by pctr. 

%% Make the sample structure with the right number of zeros 
weissman_WT_NETseq_sample = cell(16,2) %need to give the correct number of zeros according to the SGD? 
%Then can input the data into the zeros to create the og NETseq signal.
weissman_WT_NETseq_sample{1,1} = zeros(230218,1)
weissman_WT_NETseq_sample{1,2} = zeros(230218,1)
weissman_WT_NETseq_sample{2,1} = zeros(813184,1)
weissman_WT_NETseq_sample{2,2} = zeros(813184,1)
weissman_WT_NETseq_sample{3,1} = zeros(316620,1)
weissman_WT_NETseq_sample{3,2} = zeros(316620,1)
weissman_WT_NETseq_sample{4,1} = zeros(1531933,1)
weissman_WT_NETseq_sample{4,2} = zeros(1531933,1)
weissman_WT_NETseq_sample{5,1} = zeros(576874,1)
weissman_WT_NETseq_sample{5,2} = zeros(576874,1)
weissman_WT_NETseq_sample{6,1} = zeros(270161,1)
weissman_WT_NETseq_sample{6,2} = zeros(270161,1)
weissman_WT_NETseq_sample{7,1} = zeros(1090940,1)
weissman_WT_NETseq_sample{7,2} = zeros(1090940,1)
weissman_WT_NETseq_sample{8,1} = zeros(562643,1)
weissman_WT_NETseq_sample{8,2} = zeros(562643,1)
weissman_WT_NETseq_sample{9,1} = zeros(439888,1)
weissman_WT_NETseq_sample{9,2} = zeros(439888,1)
weissman_WT_NETseq_sample{10,1} = zeros(745751,1)
weissman_WT_NETseq_sample{10,2} = zeros(745751,1)
weissman_WT_NETseq_sample{11,1} = zeros(666816,1)
weissman_WT_NETseq_sample{11,2} = zeros(666816,1)
weissman_WT_NETseq_sample{12,1} = zeros(1078177,1)
weissman_WT_NETseq_sample{12,2} = zeros(1078177,1)
weissman_WT_NETseq_sample{13,1} = zeros(924431,1)
weissman_WT_NETseq_sample{13,2} = zeros(924431,1)
weissman_WT_NETseq_sample{14,1} = zeros(784333,1)
weissman_WT_NETseq_sample{14,2} = zeros(784333,1)
weissman_WT_NETseq_sample{15,1} = zeros(1091291,1)
weissman_WT_NETseq_sample{15,2} = zeros(1091291,1)
weissman_WT_NETseq_sample{16,1} = zeros(948066,1)
weissman_WT_NETseq_sample{16,2} = zeros(948066,1)

%should create correctly sized matrix in which to insert the NETseq data. 

%% Inserting the NETseq signal into the correctly sized array of zeros.
for cctr = 1:1:NoChromosomes
    
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
        weissman_WT_NETseq_sample{cctr, 2}(numeric_weissman_WT_NETseq_minus_1_readin(pctr,1), 1) =...
            numeric_weissman_WT_NETseq_minus_1_readin(pctr,2);
    end
end

%% Repeat above for plus strand 

% Initial textscan to show only 16 chromosome labels 

% fileID1 = fopen([dataloaddir 'GSM617027_WT_NC_plus.wig']);
% weissman_WT_NETseq_plus_1_readin = textscan(fileID1, '%s', ... %Initial textscan with only one string specified
%     'CommentStyle','track')
% fclose(fileID1);
% whos weissman_WT_NETseq_plus_1_readin
% 
% %% Find the boundaries between chromosomes
% Chromosomepos = zeros((NoChromosomes+1), 1) ; 
% cctr = 0 ;
% for k = 1:1:size(weissman_WT_NETseq_plus_1_readin{1,1})
%     
%     if weissman_WT_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'v' 
%     %if weissman_WT_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'c'
%         cctr = cctr +1;
%         Chromosomepos(cctr, 1) = k
%     end
%     
% end
%%
% Chromosomepos((NoChromosomes+1),1) = size(weissman_WT_NETseq_plus_1_readin{1,1}, 1) + 1 
% %pull out data from whole textscan by using a for loop at every other
% %position?
% for pctr = (Chromosomepos(cctr,1)+2):2:(Chromosomepos((cctr+1),1)-1)
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
% end
    
    
%% second textscan to find chromosome positions 
fileID2 = fopen([dataloaddir 'GSM617027_WT_NC_plus.wig']);
weissman_WT_NETseq_plus_1_readin = textscan(fileID2, '%s %s', ... %two strings 
    'CommentStyle','track')
fclose(fileID2);
whos weissman_WT_NETseq_plus_1_readin
%%
%strfind(weissman_WT_NETseq_plus_1_readin{2}{:}, Chromosome_Names)
%% Find the boundaries between chromosomes
Chromosomepos = zeros((NoChromosomes+1), 1) ; 
cctr = 0 ;
for k = 1:1:size(weissman_WT_NETseq_plus_1_readin{1,1})
    
    if weissman_WT_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'v' 
    %if weissman_WT_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'c'
        cctr = cctr +1;
        Chromosomepos(cctr, 1) = k
    end
   
    
end
Chromosomepos((NoChromosomes+1),1) = size(weissman_WT_NETseq_plus_1_readin{1,1}(:), 1) + 1
%% Looped textscan over each chromosome between the strings to pull out floating values 
% weissman_WT_NETseq_sample_nozeros = cell(16,2);
% % for cctr = 1:1:NoChromosomes
% %     for sctr = 1:1:2 
% %preallocate the size of the cell array before adding in data, done via for
% %loop hopefully won't be too slow 
%        for cctr = 1:1:NoChromosomes
%            weissman_WT_NETseq_sample_nozeros{cctr,2} = zeros(size((Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)))
%        end
%  %this gives correctly sized numeric cell array to receive data. Need to make data numeric to input.    
%    %w=0
   
  %% following should target and convert all the data between the chromosome labels from string to numeric
   %converts positional data
   numeric_weissman_WT_NETseq_plus_1_readin = zeros(size(weissman_WT_NETseq_plus_1_readin{1,1}, 1), 2);
for cctr = 1:1:NoChromosomes 
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_WT_NETseq_plus_1_readin(pctr,1) = str2double(weissman_WT_NETseq_plus_1_readin{1,1}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end
%converts signal data 
for cctr = 1:1:NoChromosomes
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_WT_NETseq_plus_1_readin(pctr, 2) = str2double(weissman_WT_NETseq_plus_1_readin{1,2}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end   
%% check there are only 16 zeros in the dataset where the labels were 
sum(numeric_weissman_WT_NETseq_plus_1_readin == 0) %It has worked!!! 


%If this works will now need to input the zeros where chromosome labels
%were before now have zeros, should be okay because I can still selectively
%target between the zeros by pctr. 


%% Inserting the NETseq signal into the correctly sized array of zeros.
for cctr = 1:1:NoChromosomes
    
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
        weissman_WT_NETseq_sample{cctr, 1}(numeric_weissman_WT_NETseq_plus_1_readin(pctr,1), 1) =...
            numeric_weissman_WT_NETseq_plus_1_readin(pctr,2);
    end
end

%% Save the output 

save(append(dataloaddir, 'weissman_WT_NETseq.mat'), 'weissman_WT_NETseq_sample')
        








 %% Initial textscan to show only 16 chromosome labels 
% 
fileID1 = fopen([dataloaddir 'GSM617029_RCO1D_minus.wig']);
weissman_rco1del_NETseq_minus_1_readin = textscan(fileID1, '%s', ... %Initial textscan with only one string specified
    'CommentStyle','track')
fclose(fileID1);
whos weissman_rco1del_NETseq_minus_1_readin
%%
%strfind(weissman_rco1del_NETseq_minus_1_readin{2}{:}, Chromosome_Names)
% %% Find the boundaries between chromosomes
% Chromosomepos = zeros((NoChromosomes+1), 1) ; 
% cctr = 0 ;
% for k = 1:1:size(weissman_rco1del_NETseq_minus_1_readin{1,1})
%     
%     if weissman_rco1del_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'v' 
%     %if weissman_rco1del_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'c'
%         cctr = cctr +1;
%         Chromosomepos(cctr, 1) = k
%     end
%     
% end
%%
% Chromosomepos((NoChromosomes+1),1) = size(weissman_rco1del_NETseq_minus_1_readin{1,1}, 1) + 1 
% %pull out data from whole textscan by using a for loop at every other
% %position?
% for pctr = (Chromosomepos(cctr,1)+2):2:(Chromosomepos((cctr+1),1)-1)
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
% end
    
    
%% second textscan to find chromosome positions 
fileID1 = fopen([dataloaddir 'GSM617029_RCO1D_minus.wig']);
weissman_rco1del_NETseq_minus_1_readin = textscan(fileID1, '%s %s', ... %two strings 
    'CommentStyle','track')
fclose(fileID1);
whos weissman_rco1del_NETseq_minus_1_readin
%%
%strfind(weissman_rco1del_NETseq_minus_1_readin{2}{:}, Chromosome_Names)
%% Find the boundaries between chromosomes
Chromosomepos = zeros((NoChromosomes+1), 1) ; 
cctr = 0 ;
for k = 1:1:size(weissman_rco1del_NETseq_minus_1_readin{1,1})
    
    if weissman_rco1del_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'v' 
    %if weissman_rco1del_NETseq_minus_1_readin{1,1}{k, 1}(1,1) == 'c'
        cctr = cctr +1;
        Chromosomepos(cctr, 1) = k
    end
   
    
end
Chromosomepos((NoChromosomes+1),1) = size(weissman_rco1del_NETseq_minus_1_readin{1,1}(:), 1) + 1
%% Looped textscan over each chromosome between the strings to pull out floating values 
% weissman_rco1del_NETseq_sample_nozeros = cell(16,2);
% % for cctr = 1:1:NoChromosomes
% %     for sctr = 1:1:2 
% %preallocate the size of the cell array before adding in data, done via for
% %loop hopefully won't be too slow 
%        for cctr = 1:1:NoChromosomes
%            weissman_rco1del_NETseq_sample_nozeros{cctr,2} = zeros(size((Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)))
%        end
%  %this gives correctly sized numeric cell array to receive data. Need to make data numeric to input.    
%    %w=0
   
  %% following should target and convert all the data between the chromosome labels from string to numeric
   %converts positional data
   numeric_weissman_rco1del_NETseq_minus_1_readin = zeros(size(weissman_rco1del_NETseq_minus_1_readin{1,1}, 1), 2);
for cctr = 1:1:NoChromosomes 
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_rco1del_NETseq_minus_1_readin(pctr,1) = str2double(weissman_rco1del_NETseq_minus_1_readin{1,1}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end
%converts signal data 
for cctr = 1:1:NoChromosomes
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_rco1del_NETseq_minus_1_readin(pctr, 2) = str2double(weissman_rco1del_NETseq_minus_1_readin{1,2}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end   
%% check there are only 16 zeros in the dataset where the labels were 
sum(numeric_weissman_rco1del_NETseq_minus_1_readin == 0) %It has worked!!! 


%If this works will now need to input the zeros where chromosome labels
%were before now have zeros, should be okay because I can still selectively
%target between the zeros by pctr. 

%% Make the sample structure with the right number of zeros 
weissman_rco1del_NETseq_sample = cell(16,2) %need to give the correct number of zeros according to the SGD? 
%Then can input the data into the zeros to create the og NETseq signal.
weissman_rco1del_NETseq_sample{1,1} = zeros(230218,1)
weissman_rco1del_NETseq_sample{1,2} = zeros(230218,1)
weissman_rco1del_NETseq_sample{2,1} = zeros(813184,1)
weissman_rco1del_NETseq_sample{2,2} = zeros(813184,1)
weissman_rco1del_NETseq_sample{3,1} = zeros(316620,1)
weissman_rco1del_NETseq_sample{3,2} = zeros(316620,1)
weissman_rco1del_NETseq_sample{4,1} = zeros(1531933,1)
weissman_rco1del_NETseq_sample{4,2} = zeros(1531933,1)
weissman_rco1del_NETseq_sample{5,1} = zeros(576874,1)
weissman_rco1del_NETseq_sample{5,2} = zeros(576874,1)
weissman_rco1del_NETseq_sample{6,1} = zeros(270161,1)
weissman_rco1del_NETseq_sample{6,2} = zeros(270161,1)
weissman_rco1del_NETseq_sample{7,1} = zeros(1090940,1)
weissman_rco1del_NETseq_sample{7,2} = zeros(1090940,1)
weissman_rco1del_NETseq_sample{8,1} = zeros(562643,1)
weissman_rco1del_NETseq_sample{8,2} = zeros(562643,1)
weissman_rco1del_NETseq_sample{9,1} = zeros(439888,1)
weissman_rco1del_NETseq_sample{9,2} = zeros(439888,1)
weissman_rco1del_NETseq_sample{10,1} = zeros(745751,1)
weissman_rco1del_NETseq_sample{10,2} = zeros(745751,1)
weissman_rco1del_NETseq_sample{11,1} = zeros(666816,1)
weissman_rco1del_NETseq_sample{11,2} = zeros(666816,1)
weissman_rco1del_NETseq_sample{12,1} = zeros(1078177,1)
weissman_rco1del_NETseq_sample{12,2} = zeros(1078177,1)
weissman_rco1del_NETseq_sample{13,1} = zeros(924431,1)
weissman_rco1del_NETseq_sample{13,2} = zeros(924431,1)
weissman_rco1del_NETseq_sample{14,1} = zeros(784333,1)
weissman_rco1del_NETseq_sample{14,2} = zeros(784333,1)
weissman_rco1del_NETseq_sample{15,1} = zeros(1091291,1)
weissman_rco1del_NETseq_sample{15,2} = zeros(1091291,1)
weissman_rco1del_NETseq_sample{16,1} = zeros(948066,1)
weissman_rco1del_NETseq_sample{16,2} = zeros(948066,1)

%should create correctly sized matrix in which to insert the NETseq data. 

%% Inserting the NETseq signal into the correctly sized array of zeros.
for cctr = 1:1:NoChromosomes
    
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
        weissman_rco1del_NETseq_sample{cctr, 2}(numeric_weissman_rco1del_NETseq_minus_1_readin(pctr,1), 1) =...
            numeric_weissman_rco1del_NETseq_minus_1_readin(pctr,2);
    end
end

%% Repeat above for plus strand 

% Initial textscan to show only 16 chromosome labels 

% fileID1 = fopen([dataloaddir 'GSM617029_RCO1D_plus.wig']);
% weissman_rco1del_NETseq_plus_1_readin = textscan(fileID1, '%s', ... %Initial textscan with only one string specified
%     'CommentStyle','track')
% fclose(fileID1);
% whos weissman_rco1del_NETseq_plus_1_readin
% 
% %% Find the boundaries between chromosomes
% Chromosomepos = zeros((NoChromosomes+1), 1) ; 
% cctr = 0 ;
% for k = 1:1:size(weissman_rco1del_NETseq_plus_1_readin{1,1})
%     
%     if weissman_rco1del_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'v' 
%     %if weissman_rco1del_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'c'
%         cctr = cctr +1;
%         Chromosomepos(cctr, 1) = k
%     end
%     
% end
%%
% Chromosomepos((NoChromosomes+1),1) = size(weissman_rco1del_NETseq_plus_1_readin{1,1}, 1) + 1 
% %pull out data from whole textscan by using a for loop at every other
% %position?
% for pctr = (Chromosomepos(cctr,1)+2):2:(Chromosomepos((cctr+1),1)-1)
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
% end
    
    
%% second textscan to find chromosome positions 
fileID2 = fopen([dataloaddir 'GSM617029_RCO1D_plus.wig']);
weissman_rco1del_NETseq_plus_1_readin = textscan(fileID2, '%s %s', ... %two strings 
    'CommentStyle','track')
fclose(fileID2);
whos weissman_rco1del_NETseq_plus_1_readin
%%
%strfind(weissman_rco1del_NETseq_plus_1_readin{2}{:}, Chromosome_Names)
%% Find the boundaries between chromosomes
Chromosomepos = zeros((NoChromosomes+1), 1) ; 
cctr = 0 ;
for k = 1:1:size(weissman_rco1del_NETseq_plus_1_readin{1,1})
    
    if weissman_rco1del_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'v' 
    %if weissman_rco1del_NETseq_plus_1_readin{1,1}{k, 1}(1,1) == 'c'
        cctr = cctr +1;
        Chromosomepos(cctr, 1) = k
    end
   
    
end
Chromosomepos((NoChromosomes+1),1) = size(weissman_rco1del_NETseq_plus_1_readin{1,1}(:), 1) + 1
%% Looped textscan over each chromosome between the strings to pull out floating values 
% weissman_rco1del_NETseq_sample_nozeros = cell(16,2);
% % for cctr = 1:1:NoChromosomes
% %     for sctr = 1:1:2 
% %preallocate the size of the cell array before adding in data, done via for
% %loop hopefully won't be too slow 
%        for cctr = 1:1:NoChromosomes
%            weissman_rco1del_NETseq_sample_nozeros{cctr,2} = zeros(size((Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)))
%        end
%  %this gives correctly sized numeric cell array to receive data. Need to make data numeric to input.    
%    %w=0
   
  %% following should target and convert all the data between the chromosome labels from string to numeric
   %converts positional data
   numeric_weissman_rco1del_NETseq_plus_1_readin = zeros(size(weissman_rco1del_NETseq_plus_1_readin{1,1}, 1), 2);
for cctr = 1:1:NoChromosomes 
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_rco1del_NETseq_plus_1_readin(pctr,1) = str2double(weissman_rco1del_NETseq_plus_1_readin{1,1}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end
%converts signal data 
for cctr = 1:1:NoChromosomes
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
       % w=w+1;
       format longg
       numeric_weissman_rco1del_NETseq_plus_1_readin(pctr, 2) = str2double(weissman_rco1del_NETseq_plus_1_readin{1,2}{pctr,1}) ;
        %textscan(fileID1((Chromosomepos(cctr,1)+1):(Chromosomepos((cctr+1),1)-1)), '%f %f', 'CommentStyle', 'track')
    end
end   
%% check there are only 16 zeros in the dataset where the labels were 
sum(numeric_weissman_rco1del_NETseq_plus_1_readin == 0) %It has worked!!! 


%If this works will now need to input the zeros where chromosome labels
%were before now have zeros, should be okay because I can still selectively
%target between the zeros by pctr. 


%% Inserting the NETseq signal into the correctly sized array of zeros.
for cctr = 1:1:NoChromosomes
    
    for pctr = (Chromosomepos(cctr,1)+1):1:(Chromosomepos((cctr+1),1)-1)
        weissman_rco1del_NETseq_sample{cctr, 1}(numeric_weissman_rco1del_NETseq_plus_1_readin(pctr,1), 1) =...
            numeric_weissman_rco1del_NETseq_plus_1_readin(pctr,2);
    end
end

%% Save the output 

save(append(dataloaddir, 'weissman_rco1del_NETseq.mat'), 'weissman_rco1del_NETseq_sample')

    %weissman_WT_NETseq_plus_1_readin =   str2double(weissman_WT_NETseq_plus_1_readin);

% weissman_WT_NETseq_plus_1_rkeadin2 = textscan(fileID1,'%f %f', ...
%     'CommentStyle','track');
% fclose(fileID1);
% whos weissman_WT_NETseq_plus_1_readin

% fileID2 = fopen([dataloaddir 'GSM617033_SET1D_plus.wig']);
% weissman_WT_NETseq_plus_1_readin = textscan(fileID2,'%s %s', ...
%     'CommentStyle','track');
% fclose(fileID2);
% whos weissman_WT_NETseq_plus_1_readin
% 
% disp('Finished importing raw weissman NETseq WT sample 1')



% %%
% 
% 
% % fileID3 = fopen([dataloaddir_WT 'Rpb3_aug18_nnt_negativestrand.wig']);
% % weissman_WT_NETseq_plus_2_readin = textscan(fileID3,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID3);
% % whos weissman_WT_NETseq_plus_2_readin
% % 
% % 
% % 
% % fileID4 = fopen([dataloaddir_WT 'Rpb3_aug18_nnt_positivestrand.wig']);
% % weissman_WT_NETseq_plus_2_readin = textscan(fileID4,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID4);
% % whos weissman_WT_NETseq_plus_2_readin
% % 
% % disp('Finished importing raw Ulku NETseq WT sample 2')
% 
% 
% 
% %% Load in the dst1 delete
% 
% % dataloaddir_end = 'notag_normalised/dst1del_nnt/';
% % 
% % dataloaddir_dst1 = [dataloaddir_start dataloaddir_end];
% % 
% % 
% % 
% % fileID1 = fopen([dataloaddir_dst1 'dst1de_nnt_rp1_negativestrand.wig']);
% % ulku_dst1del_NETseq_minus_1_readin = textscan(fileID1,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID1);
% % whos ulku_dst1del_NETseq_minus_1_readin
% % 
% % 
% % 
% % fileID2 = fopen([dataloaddir_dst1 'dst1de_nnt_rp1_positivestrand.wig']);
% % ulku_dst1del_NETseq_plus_1_readin = textscan(fileID2,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID2);
% % whos ulku_dst1del_NETseq_plus_1_readin
% % 
% % disp('Finished importing raw Ulku NETseq dst1del sample 1')
% % 
% % 
% % 
% % fileID3 = fopen([dataloaddir_dst1 'dst1de_nnt_rp2_negativestrand.wig']);
% % ulku_dst1del_NETseq_minus_2_readin = textscan(fileID3,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID3);
% % whos ulku_dst1del_NETseq_minus_2_readin
% % 
% % 
% % 
% % fileID4 = fopen([dataloaddir_dst1 'dst1de_nnt_rp2_positivestrand.wig']);
% % ulku_dst1del_NETseq_plus_2_readin = textscan(fileID4,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID4);
% % whos ulku_dst1del_NETseq_plus_2_readin
% % 
% % disp('Finished importing raw Ulku NETseq dst1del sample 2')
% % 
% % 
% % 
% % %% Load in the spt4 delete
% % 
% % dataloaddir_end = 'notag_normalised/spt4del_nnt/';
% % 
% % dataloaddir_spt4 = [dataloaddir_start dataloaddir_end];
% % 
% % 
% % 
% % fileID1 = fopen([dataloaddir_spt4 'spt4del_july18_nnt_negativestrand.wig']);
% % ulku_spt4del_NETseq_minus_1_readin = textscan(fileID1,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID1);
% % whos ulku_spt4del_NETseq_minus_1_readin
% % 
% % 
% % 
% % fileID2 = fopen([dataloaddir_spt4 'spt4del_july18_nnt_positivestrand.wig']);
% % ulku_spt4del_NETseq_plus_1_readin = textscan(fileID2,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID2);
% % whos ulku_spt4del_NETseq_plus_1_readin
% % 
% % disp('Finished importing raw Ulku NETseq spt4del sample 1')
% % 
% % 
% % 
% % fileID3 = fopen([dataloaddir_spt4 'spt4del_aug18_nnt_negativestrand.wig']);
% % ulku_spt4del_NETseq_minus_2_readin = textscan(fileID3,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID3);
% % whos ulku_spt4del_NETseq_minus_2_readin
% % 
% % 
% % 
% % fileID4 = fopen([dataloaddir_spt4 'spt4del_aug18_nnt_positivestrand.wig']);
% % ulku_spt4del_NETseq_plus_2_readin = textscan(fileID4,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID4);
% % whos ulku_spt4del_NETseq_plus_2_readin
% % 
% % disp('Finished importing raw Ulku NETseq spt4del sample 2')
% % 
% % 
% % 
% % %% Load in the xrn1 delete
% % 
% % dataloaddir_end = 'notag_normalised/xrn1del_nnt/';
% % 
% % dataloaddir_xrn1 = [dataloaddir_start dataloaddir_end];
% % 
% % 
% % 
% % fileID1 = fopen([dataloaddir_xrn1 'xrn1del_nnt_rpt1_negativestrand.wig']);
% % ulku_xrn1del_NETseq_minus_1_readin = textscan(fileID1,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID1);
% % whos ulku_xrn1del_NETseq_minus_1_readin
% % 
% % 
% % 
% % fileID2 = fopen([dataloaddir_xrn1 'xrn1del_nnt_rpt1_positivestrand.wig']);
% % ulku_xrn1del_NETseq_plus_1_readin = textscan(fileID2,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID2);
% % whos ulku_xrn1del_NETseq_plus_1_readin
% % 
% % disp('Finished importing raw Ulku NETseq xrn1del sample 1')
% % 
% % 
% % 
% % fileID3 = fopen([dataloaddir_xrn1 'xrn1del_nnt_rpt2_negativestrand.wig']);
% % ulku_xrn1del_NETseq_minus_2_readin = textscan(fileID3,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID3);
% % whos ulku_xrn1del_NETseq_minus_2_readin
% % 
% % 
% % 
% % fileID4 = fopen([dataloaddir_xrn1 'xrn1del_nnt_rpt2_positivestrand.wig']);
% % ulku_xrn1del_NETseq_plus_2_readin = textscan(fileID4,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID4);
% % whos ulku_xrn1del_NETseq_plus_2_readin
% % 
% % disp('Finished importing raw Ulku NETseq xrn1del sample 2')
% % 
% % 
% % %% Load in the SPT4FRB Rapa 
% % %Anchor Away SPT4 (I think)
% % 
% % dataloaddir_end = 'AnchorAway/';
% % 
% % dataloaddir_SPT4AA = [dataloaddir_start dataloaddir_end];
% % 
% % 
% % 
% % fileID1 = fopen([dataloaddir_SPT4AA 'Spt4FRB_Rapa_rpt1_negativestrand.wig']);
% % ulku_SPT4AA_NETseq_minus_1_readin = textscan(fileID1,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID1);
% % whos ulku_SPT4AA_NETseq_minus_1_readin
% % 
% % 
% % 
% % fileID2 = fopen([dataloaddir_SPT4AA 'Spt4FRB_Rapa_rpt1_positivestrand.wig']);
% % ulku_SPT4AA_NETseq_plus_1_readin = textscan(fileID2,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID2);
% % whos ulku_SPT4AA_NETseq_plus_1_readin
% % 
% % disp('Finished importing raw Ulku NETseq SPT4AA sample 1')
% % 
% % 
% % 
% % fileID3 = fopen([dataloaddir_SPT4AA 'Spt4FRB_Rapa_rpt2_negativestrand.wig']);
% % ulku_SPT4AA_NETseq_minus_2_readin = textscan(fileID3,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID3);
% % whos ulku_SPT4AA_NETseq_minus_2_readin
% % 
% % 
% % 
% % fileID4 = fopen([dataloaddir_SPT4AA 'Spt4FRB_Rapa_rpt2_positivestrand.wig']);
% % ulku_SPT4AA_NETseq_plus_2_readin = textscan(fileID4,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID4);
% % whos ulku_SPT4AA_NETseq_plus_2_readin
% % 
% % disp('Finished importing raw Ulku NETseq SPT4AA sample 2')
% % 
% % 
% % %% Load in the SPT5FRB Rapa 
% % %Anchor Away SPT5 (I think)
% % 
% % dataloaddir_end = 'AnchorAway/';
% % 
% % dataloaddir_SPT5AA = [dataloaddir_start dataloaddir_end];
% % 
% % 
% % 
% % fileID1 = fopen([dataloaddir_SPT5AA 'Spt5FRB_Rapa_rpt1_negativestrand.wig']);
% % ulku_SPT5AA_NETseq_minus_1_readin = textscan(fileID1,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID1);
% % whos ulku_SPT5AA_NETseq_minus_1_readin
% % 
% % 
% % 
% % fileID2 = fopen([dataloaddir_SPT5AA 'Spt5FRB_Rapa_rpt1_positivestrand.wig']);
% % ulku_SPT5AA_NETseq_plus_1_readin = textscan(fileID2,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID2);
% % whos ulku_SPT5AA_NETseq_plus_1_readin
% % 
% % disp('Finished importing raw Ulku NETseq SPT5AA sample 1')
% % 
% % 
% % 
% % fileID3 = fopen([dataloaddir_SPT5AA 'Spt5FRB_Rapa_rpt2_negativestrand.wig']);
% % ulku_SPT5AA_NETseq_minus_2_readin = textscan(fileID3,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID3);
% % whos ulku_SPT5AA_NETseq_minus_2_readin
% % 
% % 
% % 
% % fileID4 = fopen([dataloaddir_SPT5AA 'Spt5FRB_Rapa_rpt2_positivestrand.wig']);
% % ulku_SPT5AA_NETseq_plus_2_readin = textscan(fileID4,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID4);
% % whos ulku_SPT5AA_NETseq_plus_2_readin
% % 
% % disp('Finished importing raw Ulku NETseq SPT4AA sample 2')
% % 
% % 
% % 
% % %% Load in the SPT4FRB DMSO
% % %Anchor Away SPT4 control (I think)
% % 
% % dataloaddir_end = 'AnchorAway/';
% % 
% % dataloaddir_SPT4AADMSO = [dataloaddir_start dataloaddir_end];
% % 
% % 
% % 
% % fileID1 = fopen([dataloaddir_SPT4AADMSO 'Spt4FRB_DMSO_negativestrand.wig']);
% % ulku_SPT4AADMSO_NETseq_minus_1_readin = textscan(fileID1,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID1);
% % whos ulku_SPT4AADMSO_NETseq_minus_1_readin
% % 
% % 
% % 
% % fileID2 = fopen([dataloaddir_SPT4AADMSO 'Spt4FRB_DMSO_positivestrand.wig']);
% % ulku_SPT4AADMSO_NETseq_plus_1_readin = textscan(fileID2,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID2);
% % whos ulku_SPT4AADMSO_NETseq_plus_1_readin
% % 
% % disp('Finished importing raw Ulku NETseq SPT4AA DMSO')
% % 
% % 
% % %% Load in the SPT5FRB DMSO
% % %Anchor Away SPT5 control (I think)
% % 
% % dataloaddir_end = 'AnchorAway/';
% % 
% % dataloaddir_SPT5AADMSO = [dataloaddir_start dataloaddir_end];
% % 
% % 
% % 
% % fileID1 = fopen([dataloaddir_SPT5AADMSO 'Spt5FRB_DMSO_negativestrand.wig']);
% % ulku_SPT5AADMSO_NETseq_minus_1_readin = textscan(fileID1,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID1);
% % whos ulku_SPT5AADMSO_NETseq_minus_1_readin
% % 
% % 
% % 
% % fileID2 = fopen([dataloaddir_SPT5AADMSO 'Spt5FRB_DMSO_positivestrand.wig']);
% % ulku_SPT5AADMSO_NETseq_plus_1_readin = textscan(fileID2,'%s %f %f %f', ...
% %     'CommentStyle','#');
% % fclose(fileID2);
% % whos ulku_SPT5AADMSO_NETseq_plus_1_readin
% % 
% % disp('Finished importing raw Ulku NETseq SPT5AA DMSO')
% 
% 
% 
% 
% %% NB the DMSO SPT4FRB and SPT5FRB may be combined
% %Only one repeat for each but they /it{should} be the same but also valid
% %to consider them each separately - they should identify any effect of the
% %FRB tag (I think).
% 
% 
% %% Calculate the max point for each strand
% %(might later need to extend it to the true points as could cause problems
% %if annotated regions lie outside the range of measurement)
% %Also, might need to extend it to all the different samples
% %(in case additional signal is measured at the end of chromosomes in a
% %specific sample)
% 
% %Do a separate loop for each sample so that adding additional samples is
% %less hassle
% 
% %choosing to have the same extent for each stand for simplicity
% 
% max_chromosome_extent = zeros(NoChromosomes,1);
% min_chromosome_loc = NaN*ones(NoChromosomes,1);
% 
% %WT sample 1
% 
% for chrctr=1:1:NoChromosomes
% 
%     %minus strand
%     
%     temp_chromo_locs = strfind(weissman_set1del_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [weissman_set1del_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
%         weissman_set1delNETseq_minus_1_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [weissman_set1delNETseq_minus_1_readin{2}(temp_chromo_locs), ...
%         weissman_set1delNETseq_minus_1_readin{3}(temp_chromo_locs)]));
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
%     temp_chromo_locs = strfind(weissman_set1del_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
%     
%     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
%     
%     temp_max1 = max(max( ...
%         [weissman_set1del_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
%         weissman_set1del_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
%     
%     temp_min1 = min(min( ...
%         [weissman_set1del_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
%         weissman_set1del_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
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
% % %WT sample 2
% % 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %minus strand
% %     
% %     temp_chromo_locs = strfind(weissman_set1delNETseq_minus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [weissman_set1del_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
% %         weissman_set1del_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [weissman_set1del_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
% %         weissman_set1del_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% %     
% %     %plus strand
% %         
% %     temp_chromo_locs = strfind(weissman_set1del_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [weissman_set1del_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
% %         weissman_set1del_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [weissman_set1del_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
% %         weissman_set1del_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% % end
% % 
% % 
% % 
% % 
% % %dst1del sample 1
% % 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %minus strand
% %     
% %     temp_chromo_locs = strfind(ulku_dst1del_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_dst1del_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_dst1del_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_dst1del_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_dst1del_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% %     
% %     %plus strand
% %     
% %     temp_chromo_locs = strfind(ulku_dst1del_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_dst1del_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_dst1del_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_dst1del_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_dst1del_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% % end
% % 
% % 
% % %dst1del sample 2
% % 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %minus strand
% %     
% %     temp_chromo_locs = strfind(ulku_dst1del_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_dst1del_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_dst1del_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_dst1del_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_dst1del_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% %     
% %     %plus strand
% %         
% %     temp_chromo_locs = strfind(ulku_dst1del_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_dst1del_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_dst1del_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_dst1del_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_dst1del_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% % end
% % 
% % 
% % 
% % 
% % %spt4del sample 1
% % 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %minus strand
% %     
% %     temp_chromo_locs = strfind(ulku_spt4del_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_spt4del_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_spt4del_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_spt4del_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_spt4del_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% %     
% %     %plus strand
% %     
% %     temp_chromo_locs = strfind(ulku_spt4del_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_spt4del_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_spt4del_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_spt4del_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_spt4del_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% % end
% % 
% % 
% % %spt4del sample 2
% % 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %minus strand
% %     
% %     temp_chromo_locs = strfind(ulku_spt4del_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_spt4del_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_spt4del_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_spt4del_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_spt4del_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% %     
% %     %plus strand
% %         
% %     temp_chromo_locs = strfind(ulku_spt4del_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_spt4del_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_spt4del_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_spt4del_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_spt4del_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% % end
% % 
% % 
% % 
% % 
% % %xrn1del sample 1
% % 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %minus strand
% %     
% %     temp_chromo_locs = strfind(ulku_xrn1del_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_xrn1del_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_xrn1del_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_xrn1del_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_xrn1del_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% %     
% %     %plus strand
% %     
% %     temp_chromo_locs = strfind(ulku_xrn1del_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_xrn1del_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_xrn1del_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_xrn1del_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_xrn1del_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% % end
% % 
% % 
% % %xrn1del sample 2
% % 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %minus strand
% %     
% %     temp_chromo_locs = strfind(ulku_xrn1del_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_xrn1del_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_xrn1del_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_xrn1del_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_xrn1del_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% %     
% %     %plus strand
% %         
% %     temp_chromo_locs = strfind(ulku_xrn1del_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_xrn1del_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_xrn1del_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_xrn1del_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_xrn1del_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% % end
% % 
% % 
% % 
% % 
% % 
% % %SPT4AA sample 1
% % 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %minus strand
% %     
% %     temp_chromo_locs = strfind(ulku_SPT4AA_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_SPT4AA_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT4AA_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_SPT4AA_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT4AA_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% %     
% %     %plus strand
% %     
% %     temp_chromo_locs = strfind(ulku_SPT4AA_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_SPT4AA_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT4AA_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_SPT4AA_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT4AA_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% % end
% % 
% % 
% % %SPT4AA sample 2
% % 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %minus strand
% %     
% %     temp_chromo_locs = strfind(ulku_SPT4AA_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_SPT4AA_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT4AA_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_SPT4AA_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT4AA_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% %     
% %     %plus strand
% %         
% %     temp_chromo_locs = strfind(ulku_SPT4AA_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_SPT4AA_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT4AA_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_SPT4AA_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT4AA_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% % end
% % 
% % 
% % 
% % 
% % 
% % %SPT5AA sample 1
% % 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %minus strand
% %     
% %     temp_chromo_locs = strfind(ulku_SPT5AA_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_SPT5AA_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT5AA_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_SPT5AA_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT5AA_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% %     
% %     %plus strand
% %     
% %     temp_chromo_locs = strfind(ulku_SPT5AA_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_SPT5AA_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT5AA_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_SPT5AA_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT5AA_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% % end
% % 
% % 
% % %SPT5AA sample 2
% % 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %minus strand
% %     
% %     temp_chromo_locs = strfind(ulku_SPT5AA_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_SPT5AA_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT5AA_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_SPT5AA_NETseq_minus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT5AA_NETseq_minus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% %     
% %     %plus strand
% %         
% %     temp_chromo_locs = strfind(ulku_SPT5AA_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_SPT5AA_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT5AA_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_SPT5AA_NETseq_plus_2_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT5AA_NETseq_plus_2_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% % end
% % 
% % 
% % 
% % %SPT4AA DMSO
% % 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %minus strand
% %     
% %     temp_chromo_locs = strfind(ulku_SPT4AADMSO_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_SPT4AADMSO_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT4AADMSO_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_SPT4AADMSO_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT4AADMSO_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% %     
% %     %plus strand
% %     
% %     temp_chromo_locs = strfind(ulku_SPT4AADMSO_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_SPT4AADMSO_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT4AADMSO_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_SPT4AADMSO_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT4AADMSO_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% % end
% % 
% % 
% % 
% % %SPT5AA DMSO
% % 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %minus strand
% %     
% %     temp_chromo_locs = strfind(ulku_SPT5AADMSO_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_SPT5AADMSO_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT5AADMSO_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_SPT5AADMSO_NETseq_minus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT5AADMSO_NETseq_minus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% %     
% %     %plus strand
% %     
% %     temp_chromo_locs = strfind(ulku_SPT5AADMSO_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     temp_chromo_locs = ~cellfun(@isempty,temp_chromo_locs);
% %     
% %     temp_max1 = max(max( ...
% %         [ulku_SPT5AADMSO_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT5AADMSO_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     temp_min1 = min(min( ...
% %         [ulku_SPT5AADMSO_NETseq_plus_1_readin{2}(temp_chromo_locs), ...
% %         ulku_SPT5AADMSO_NETseq_plus_1_readin{3}(temp_chromo_locs)]));
% %     
% %     
% %     max_chromosome_extent(chrctr) = ...
% %         max(max_chromosome_extent(chrctr),temp_max1);
% %     min_chromosome_loc(chrctr) = ...
% %         min(min_chromosome_loc(chrctr),temp_min1);
% %     
% % end
% 
% 
% 
% %%
% 
% weissman_set1del_NETseq_1 = cell(NoChromosomes,2);
% weissman_set1del_NETseq_2 = cell(NoChromosomes,2);
% 
% for chrctr=1:NoChromosomes
%     
%     for strandctr = 1:2
%    
%         weissman_set1del_NETseq_1{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
%         weissman_set1del_NETseq_2{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
%     
%     end
%     
% end
% 
% 
% % Ulku_dst1del_NETseq_1 = cell(NoChromosomes,2);
% % Ulku_dst1del_NETseq_2 = cell(NoChromosomes,2);
% % 
% % for chrctr=1:NoChromosomes
% %     
% %     for strandctr = 1:2
% %    
% %         Ulku_dst1del_NETseq_1{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
% %         Ulku_dst1del_NETseq_2{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
% %     
% %     end
% %     
% % end
% % 
% % 
% % Ulku_spt4del_NETseq_1 = cell(NoChromosomes,2);
% % Ulku_spt4del_NETseq_2 = cell(NoChromosomes,2);
% % 
% % for chrctr=1:NoChromosomes
% %     
% %     for strandctr = 1:2
% %    
% %         Ulku_spt4del_NETseq_1{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
% %         Ulku_spt4del_NETseq_2{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
% %     
% %     end
% %     
% % end
% % 
% % 
% % Ulku_xrn1del_NETseq_1 = cell(NoChromosomes,2);
% % Ulku_xrn1del_NETseq_2 = cell(NoChromosomes,2);
% % 
% % for chrctr=1:NoChromosomes
% %     
% %     for strandctr = 1:2
% %    
% %         Ulku_xrn1del_NETseq_1{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
% %         Ulku_xrn1del_NETseq_2{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
% %     
% %     end
% %     
% % end
% % 
% % 
% % Ulku_SPT4AA_NETseq_1 = cell(NoChromosomes,2);
% % Ulku_SPT4AA_NETseq_2 = cell(NoChromosomes,2);
% % 
% % for chrctr=1:NoChromosomes
% %     
% %     for strandctr = 1:2
% %    
% %         Ulku_SPT4AA_NETseq_1{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
% %         Ulku_SPT4AA_NETseq_2{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
% %     
% %     end
% %     
% % end
% % 
% % 
% % Ulku_SPT5AA_NETseq_1 = cell(NoChromosomes,2);
% % Ulku_SPT5AA_NETseq_2 = cell(NoChromosomes,2);
% % 
% % for chrctr=1:NoChromosomes
% %     
% %     for strandctr = 1:2
% %    
% %         Ulku_SPT5AA_NETseq_1{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
% %         Ulku_SPT5AA_NETseq_2{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
% %     
% %     end
% %     
% % end
% % 
% % 
% % Ulku_SPT4AADMSO_NETseq_1 = cell(NoChromosomes,2);
% % Ulku_SPT5AADMSO_NETseq_1 = cell(NoChromosomes,2);
% % 
% % for chrctr=1:NoChromosomes
% %     
% %     for strandctr = 1:2
% %    
% %         Ulku_SPT4AADMSO_NETseq_1{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
% %         Ulku_SPT5AADMSO_NETseq_1{chrctr,strandctr} = zeros(max_chromosome_extent(chrctr),1);
% %     
% %     end
% %     
% % end
% % 
% 
% 
% %% Process the WT data
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
%     tempmatches = strfind(weissman_set1del_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = weissman_set1del_NETseq_minus_1_readin{2}(tempmatches_processed);
%     temp_ends = weissman_set1del_NETseq_minus_1_readin{3}(tempmatches_processed);
%     temp_vals = weissman_set1del_NETseq_minus_1_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         weissman_set1del_NETseq_1{chrctr,2}(...
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
%     tempmatches = strfind(weissman_set1del_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
%     
%     tempmatches_processed = ~cellfun(@isempty,tempmatches);
%     
%     temp_starts = weissman_set1del_NETseq_plus_1_readin{2}(tempmatches_processed);
%     temp_ends = weissman_set1del_NETseq_plus_1_readin{3}(tempmatches_processed);
%     temp_vals = weissman_set1del_NETseq_plus_1_readin{4}(tempmatches_processed);
%     
%     
%     for ictr=1:1:sum(tempmatches_processed)
%        
%         weissman_set1del_NETseq_1{chrctr,1}(...
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
% %     %sample2
% %     
% %     %minus strand
% %     tempmatches = strfind(weissman_set1del_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = weissman_set1del_NETseq_minus_2_readin{2}(tempmatches_processed);
% %     temp_ends = weissman_set1del_NETseq_minus_2_readin{3}(tempmatches_processed);
% %     temp_vals = weissman_set1del_NETseq_minus_2_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         weissman_set1del_NETseq_2{chrctr,2}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 2 minus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     
% %     %plus strand
% %     tempmatches = strfind(weissman_set1del_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = weissman_set1del_NETseq_plus_2_readin{2}(tempmatches_processed);
% %     temp_ends = weissman_set1del_NETseq_plus_2_readin{3}(tempmatches_processed);
% %     temp_vals = weissman_set1del_NETseq_plus_2_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         weissman_set1del_NETseq_2{chrctr,1}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 2 plus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     
%  end
% 
% 
% %%
% 
% save([dataloaddir 'weissman_set1del_NETseq_1.mat'],'weissman_set1del_NETseq_1','-v7.3')
% 
% 
% %%
% 
% save([dataloaddir_WT 'weissman_set1del_NETseq_2.mat'],'weissman_set1del_NETseq_2','-v7.3')
% 
% 
% 
% 
% 
% 
% %% Process the dst1del data
% 
% %NB I am going from the minimum chromosome location +1 to the end location
% %might mean things are out by a basepair depending on how things were
% %supposed to be interpreted - referencing in the original file starts from
% %0 so I am assuming this is the correct way
% 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %sample1
% %     
% %     %minus strand
% %     tempmatches = strfind(ulku_dst1del_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_dst1del_NETseq_minus_1_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_dst1del_NETseq_minus_1_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_dst1del_NETseq_minus_1_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_dst1del_NETseq_1{chrctr,2}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 1 minus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     %plus strand
% %     tempmatches = strfind(ulku_dst1del_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_dst1del_NETseq_plus_1_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_dst1del_NETseq_plus_1_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_dst1del_NETseq_plus_1_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_dst1del_NETseq_1{chrctr,1}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 1 plus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     
% %     
% %     %sample2
% %     
% %     %minus strand
% %     tempmatches = strfind(ulku_dst1del_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_dst1del_NETseq_minus_2_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_dst1del_NETseq_minus_2_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_dst1del_NETseq_minus_2_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_dst1del_NETseq_2{chrctr,2}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 2 minus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     
% %     %plus strand
% %     tempmatches = strfind(ulku_dst1del_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_dst1del_NETseq_plus_2_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_dst1del_NETseq_plus_2_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_dst1del_NETseq_plus_2_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_dst1del_NETseq_2{chrctr,1}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 2 plus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     
% % end
% % 
% % 
% % %%
% % 
% % save([dataloaddir_dst1 'ulku_dst1del_NETseq_1.mat'],'Ulku_dst1del_NETseq_1','-v7.3')
% % 
% % 
% % %%
% % 
% % save([dataloaddir_dst1 'ulku_dst1del_NETseq_2.mat'],'Ulku_dst1del_NETseq_2','-v7.3')
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % %% Process the spt4del data
% % 
% % %STILL TO BE UPDATED
% % 
% % %NB I am going from the minimum chromosome location +1 to the end location
% % %might mean things are out by a basepair depending on how things were
% % %supposed to be interpreted - referencing in the original file starts from
% % %0 so I am assuming this is the correct way
% % 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %sample1
% %     
% %     %minus strand
% %     tempmatches = strfind(ulku_spt4del_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_spt4del_NETseq_minus_1_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_spt4del_NETseq_minus_1_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_spt4del_NETseq_minus_1_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_spt4del_NETseq_1{chrctr,2}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 1 minus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     %plus strand
% %     tempmatches = strfind(ulku_spt4del_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_spt4del_NETseq_plus_1_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_spt4del_NETseq_plus_1_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_spt4del_NETseq_plus_1_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_spt4del_NETseq_1{chrctr,1}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 1 plus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     
% %     
% %     %sample2
% %     
% %     %minus strand
% %     tempmatches = strfind(ulku_spt4del_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_spt4del_NETseq_minus_2_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_spt4del_NETseq_minus_2_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_spt4del_NETseq_minus_2_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_spt4del_NETseq_2{chrctr,2}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 2 minus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     
% %     %plus strand
% %     tempmatches = strfind(ulku_spt4del_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_spt4del_NETseq_plus_2_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_spt4del_NETseq_plus_2_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_spt4del_NETseq_plus_2_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_spt4del_NETseq_2{chrctr,1}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 2 plus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     
% % end
% % 
% % 
% % %%
% % 
% % save([dataloaddir_spt4 'ulku_spt4del_NETseq_1.mat'],'Ulku_spt4del_NETseq_1','-v7.3')
% % 
% % 
% % %%
% % 
% % save([dataloaddir_spt4 'ulku_spt4del_NETseq_2.mat'],'Ulku_spt4del_NETseq_2','-v7.3')
% % 
% % 
% % 
% % 
% % 
% % %% Process the xrn1del data
% % 
% % %NEED TO FINISH UPDATING THIS!
% % 
% % %NB I am going from the minimum chromosome location +1 to the end location
% % %might mean things are out by a basepair depending on how things were
% % %supposed to be interpreted - referencing in the original file starts from
% % %0 so I am assuming this is the correct way
% % 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %sample1
% %     
% %     %minus strand
% %     tempmatches = strfind(ulku_xrn1del_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_xrn1del_NETseq_minus_1_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_xrn1del_NETseq_minus_1_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_xrn1del_NETseq_minus_1_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_xrn1del_NETseq_1{chrctr,2}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 1 minus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     %plus strand
% %     tempmatches = strfind(ulku_xrn1del_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_xrn1del_NETseq_plus_1_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_xrn1del_NETseq_plus_1_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_xrn1del_NETseq_plus_1_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_xrn1del_NETseq_1{chrctr,1}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 1 plus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     
% %     
% %     %sample2
% %     
% %     %minus strand
% %     tempmatches = strfind(ulku_xrn1del_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_xrn1del_NETseq_minus_2_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_xrn1del_NETseq_minus_2_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_xrn1del_NETseq_minus_2_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_xrn1del_NETseq_2{chrctr,2}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 2 minus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     
% %     %plus strand
% %     tempmatches = strfind(ulku_xrn1del_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_xrn1del_NETseq_plus_2_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_xrn1del_NETseq_plus_2_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_xrn1del_NETseq_plus_2_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_xrn1del_NETseq_2{chrctr,1}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 2 plus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     
% % end
% % 
% % 
% % %%
% % 
% % save([dataloaddir_xrn1 'ulku_xrn1del_NETseq_1.mat'],'Ulku_xrn1del_NETseq_1','-v7.3')
% % 
% % 
% % %%
% % 
% % save([dataloaddir_xrn1 'ulku_xrn1del_NETseq_2.mat'],'Ulku_xrn1del_NETseq_2','-v7.3')
% % 
% % 
% % 
% % %% Process the SPT4AA data
% % 
% % %NB I am going from the minimum chromosome location +1 to the end location
% % %might mean things are out by a basepair depending on how things were
% % %supposed to be interpreted - referencing in the original file starts from
% % %0 so I am assuming this is the correct way
% % 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %sample1
% %     
% %     %minus strand
% %     tempmatches = strfind(ulku_SPT4AA_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_SPT4AA_NETseq_minus_1_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_SPT4AA_NETseq_minus_1_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_SPT4AA_NETseq_minus_1_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_SPT4AA_NETseq_1{chrctr,2}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 1 minus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     %plus strand
% %     tempmatches = strfind(ulku_SPT4AA_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_SPT4AA_NETseq_plus_1_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_SPT4AA_NETseq_plus_1_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_SPT4AA_NETseq_plus_1_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_SPT4AA_NETseq_1{chrctr,1}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 1 plus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     
% %     
% %     %sample2
% %     
% %     %minus strand
% %     tempmatches = strfind(ulku_SPT4AA_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_SPT4AA_NETseq_minus_2_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_SPT4AA_NETseq_minus_2_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_SPT4AA_NETseq_minus_2_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_SPT4AA_NETseq_2{chrctr,2}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 2 minus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     
% %     %plus strand
% %     tempmatches = strfind(ulku_SPT4AA_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_SPT4AA_NETseq_plus_2_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_SPT4AA_NETseq_plus_2_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_SPT4AA_NETseq_plus_2_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_SPT4AA_NETseq_2{chrctr,1}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 2 plus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     
% % end
% % 
% % 
% % %%
% % 
% % save([dataloaddir_SPT4AA 'ulku_SPT4AA_NETseq_1.mat'],'Ulku_SPT4AA_NETseq_1','-v7.3')
% % 
% % 
% % %%
% % 
% % save([dataloaddir_SPT4AA 'ulku_SPT4AA_NETseq_2.mat'],'Ulku_SPT4AA_NETseq_2','-v7.3')
% % 
% % 
% % 
% % 
% % %% Process the SPT5AA data
% % 
% % %NB I am going from the minimum chromosome location +1 to the end location
% % %might mean things are out by a basepair depending on how things were
% % %supposed to be interpreted - referencing in the original file starts from
% % %0 so I am assuming this is the correct way
% % 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %sample1
% %     
% %     %minus strand
% %     tempmatches = strfind(ulku_SPT5AA_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_SPT5AA_NETseq_minus_1_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_SPT5AA_NETseq_minus_1_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_SPT5AA_NETseq_minus_1_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_SPT5AA_NETseq_1{chrctr,2}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 1 minus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     %plus strand
% %     tempmatches = strfind(ulku_SPT5AA_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_SPT5AA_NETseq_plus_1_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_SPT5AA_NETseq_plus_1_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_SPT5AA_NETseq_plus_1_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_SPT5AA_NETseq_1{chrctr,1}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 1 plus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     
% %     
% %     %sample2
% %     
% %     %minus strand
% %     tempmatches = strfind(ulku_SPT5AA_NETseq_minus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_SPT5AA_NETseq_minus_2_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_SPT5AA_NETseq_minus_2_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_SPT5AA_NETseq_minus_2_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_SPT5AA_NETseq_2{chrctr,2}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 2 minus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     
% %     %plus strand
% %     tempmatches = strfind(ulku_SPT5AA_NETseq_plus_2_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_SPT5AA_NETseq_plus_2_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_SPT5AA_NETseq_plus_2_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_SPT5AA_NETseq_plus_2_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_SPT5AA_NETseq_2{chrctr,1}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 2 plus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     
% % end
% % 
% % 
% % %%
% % 
% % save([dataloaddir_SPT5AA 'ulku_SPT5AA_NETseq_1.mat'],'Ulku_SPT5AA_NETseq_1','-v7.3')
% % 
% % 
% % %%
% % 
% % save([dataloaddir_SPT5AA 'ulku_SPT5AA_NETseq_2.mat'],'Ulku_SPT5AA_NETseq_2','-v7.3')
% % 
% % 
% % 
% % 
% % %% Process the SPT4AADMSO data
% % 
% % %NB I am going from the minimum chromosome location +1 to the end location
% % %might mean things are out by a basepair depending on how things were
% % %supposed to be interpreted - referencing in the original file starts from
% % %0 so I am assuming this is the correct way
% % 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %sample1
% %     
% %     %minus strand
% %     tempmatches = strfind(ulku_SPT4AADMSO_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_SPT4AADMSO_NETseq_minus_1_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_SPT4AADMSO_NETseq_minus_1_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_SPT4AADMSO_NETseq_minus_1_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_SPT4AADMSO_NETseq_1{chrctr,2}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 1 minus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     %plus strand
% %     tempmatches = strfind(ulku_SPT4AADMSO_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_SPT4AADMSO_NETseq_plus_1_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_SPT4AADMSO_NETseq_plus_1_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_SPT4AADMSO_NETseq_plus_1_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_SPT4AADMSO_NETseq_1{chrctr,1}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 1 plus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     
% %     
% % end
% % 
% % 
% % %%
% % 
% % save([dataloaddir_SPT4AADMSO 'ulku_SPT4AADMSO_NETseq_1.mat'],'Ulku_SPT4AADMSO_NETseq_1','-v7.3')
% % 
% % 
% % 
% % 
% % %% Process the SPT5AADMSO data
% % 
% % %NB I am going from the minimum chromosome location +1 to the end location
% % %might mean things are out by a basepair depending on how things were
% % %supposed to be interpreted - referencing in the original file starts from
% % %0 so I am assuming this is the correct way
% % 
% % for chrctr=1:1:NoChromosomes
% % 
% %     %sample1
% %     
% %     %minus strand
% %     tempmatches = strfind(ulku_SPT5AADMSO_NETseq_minus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_SPT5AADMSO_NETseq_minus_1_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_SPT5AADMSO_NETseq_minus_1_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_SPT5AADMSO_NETseq_minus_1_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_SPT5AADMSO_NETseq_1{chrctr,2}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 1 minus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     %plus strand
% %     tempmatches = strfind(ulku_SPT5AADMSO_NETseq_plus_1_readin{1},Chromosome_Names{chrctr});
% %     
% %     tempmatches_processed = ~cellfun(@isempty,tempmatches);
% %     
% %     temp_starts = ulku_SPT5AADMSO_NETseq_plus_1_readin{2}(tempmatches_processed);
% %     temp_ends = ulku_SPT5AADMSO_NETseq_plus_1_readin{3}(tempmatches_processed);
% %     temp_vals = ulku_SPT5AADMSO_NETseq_plus_1_readin{4}(tempmatches_processed);
% %     
% %     
% %     for ictr=1:1:sum(tempmatches_processed)
% %        
% %         Ulku_SPT5AADMSO_NETseq_1{chrctr,1}(...
% %             (temp_starts(ictr)+1):temp_ends(ictr) ) = temp_vals(ictr);
% %         
% %         if mod(ictr,10000) == 0
% %             disp(['Completed ' num2str(ictr) ' of ' ...
% %                 num2str(length(temp_starts)) ...
% %                 ' 1 plus iterations for chromosome' num2str(chrctr) '.'])
% %         end
% %         
% %     end
% %     
% %     
% %     
% % end
% % 
% % 
% % %%
% % 
% % save([dataloaddir_SPT5AADMSO 'ulku_SPT5AADMSO_NETseq_1.mat'],'Ulku_SPT5AADMSO_NETseq_1','-v7.3')
% 
% 

%% NOTES

%Still needs to be run to pick out obvious bugs

%Also needs to be check for any persisting silly errors like writing the
%wrong data to the save file or processing the wrong data in a given
%section