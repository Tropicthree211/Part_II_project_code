 %shape normalisation 
 %% find mean
MEAN = mean(Total_NET_Seq_Matrix, 2) %find mean signal per gene 

%% subtract list of means from entire matrix 
MINUS_MEAN = bsxfun(@minus, Total_NET_Seq_Matrix, MEAN)   ;

%% Find standard dev
STD = std(Total_NET_Seq_Matrix, 0, 2);

%% Divide through by STD
Normalised_NETseq_Matrix = bsxfun(@rdivide, MINUS_MEAN, STD) ;

%% Image Normalised NETseq matrix 
NETseq_Heatmap = figure 
imagesc(Normalised_NETseq_Matrix)
truesize(NETseq_Heatmap)

