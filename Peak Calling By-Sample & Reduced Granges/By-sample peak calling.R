library(Signac)
#File path
fragpath <- "/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Sample/GSM6058856_B004-A-004-D_20200817_fragments.tsv.gz"
file.exists(fragpath)

#Show the first five lines of the fragment file
head(read_tsv(fragpath, col_names = F, n_max = 5))
#The col_names function tells read_tsv() : The first line of the file is not column names but data.
#n_max = 5 means show the first five lines。

#Count the number of fragments per cell in the fragment file
require(Signac)
total_counts <- CountFragments(fragpath)

#Screening cell barcodes with more than 100 fragments
barcodes <- total_counts$CB[total_counts$frequency_count > 100]

#Create a fragment object
frags <- Signac::CreateFragmentObject(path=fragpath, cells=barcodes)

#peak calling
#peak calling for MACS3 in terminal
system('bash -c "/Users/yulongqiu/macs3_env/bin/macs3 callpeak -t /Users/yulongqiu/Desktop/Biostats/Master_Thesis/Sample/GSM6058856_B004-A-004-D_20200817_fragments.tsv.gz -g 2.7e+09 -f BED --nomodel --extsize 200 --shift -100 -n SeuratProject --outdir /Users/yulongqiu/Desktop/Biostats/Master_Thesis/peaks"')








