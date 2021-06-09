%% Histogram to visualise 
figure
histogram(Normalised_NETseq_Matrix(:))

%% Histogram of unfiltered Netseqmatrix
figure
histogram(Total_NET_Seq_Matrix(:))

%% Visualise Netseqmatrix as heatmap
%need to cutoff high values otherwise see nothing 
figure
imagesc(Total_NET_Seq_Matrix(1:500, :), [0 20])
truesize
%see high expression genes as yellow streaks consistently above 20 


%% visualisation of normalised matrix 
figure 
imagesc(Normalised_NETseq_Matrix(1:500, :), [-3 10]) %10 seemed like a valid cutoff from histogram 
truesize

%% Visualise log plot 
figure
Log_NETseq_matrix = log(NETseq_matrix_no_negative);
imagesc(Log_NETseq_matrix(1:500,:)) %fewer genes see easier change values to see all genes 
truesize

%% 