library(cluster)
library(factoextra)
library(tidyverse)

files = c("allFMEgenes_follow_up_matrix.csv", "non_zeroFMEgenes_follow_up_matrix.csv", "negativeFMEgenes_follow_up_matrix.csv", "mostneg_100FMEgenes_follow_up_matrix.csv", "mostneg_50FMEgenes_follow_up_matrix.csv", "mostneg_10FMEgenes_follow_up_matrix.csv", "positiveFMEgenes_follow_up_matrix.csv","mostpos_100FMEgenes_follow_up_matrix.csv", "mostpos_50FMEgenes_follow_up_matrix.csv","mostpos_10FMEgenes_follow_up_matrix.csv")

data = read.csv("/home/levinsj/Fetal_dir/AminData/FME_lasso/allFMEgenes_follow_up_matrix.csv", header=TRUE, row.names=1)

df = data.frame("ID" = rownames(data))
print(df)

for (i in files){
	data = read.csv(paste0("/home/levinsj/Fetal_dir/AminData/FME_lasso/",i), header=TRUE, row.names=1)
	data= scale(data)
	fviz_nbclust(data, kmeans, method = "silhouette")
	d1 <- dist(data, method = "euclidean") # distance matrix
	fit1 <- hclust(d1, method="ward")
	groups1 <- cutree(fit1, k=3)
	#print(names(groups1))
	#print(table(groups1))
	df[i] <- groups1
	data.total.followup = data.frame (data, groups1)
	write.csv (data, file = " data.total.followup.csv")
}
write.csv(df, file = "/home/levinsj/Fetal_dir/AminData/FME_lasso/clustering.csv")
