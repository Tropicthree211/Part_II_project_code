%% averaging every 10bps Bins
%make the matrix 
NETseq_Bins10bp_avg = zeros(Total_gene_number, 100);
for k = 1:Total_gene_number
    for x = 1:100
        Bin = mean(Total_NET_Seq_Matrix(k,10*x-9:10*x));
        NETseq_Bins10bp_avg(k, x) = Bin  ;
    end
end
    