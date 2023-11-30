library(survival)
library(survminer)

data = read.csv("~/Fetal_dir/AminData/FME_lasso/updated_followup_with_KIM1_clusters.csv", row.names =1)

print(data)

pdf("~/Fetal_dir/AminData/FME_lasso/KIM1.pdf", height = 6, width = 6)
ggsurvplot(survfit(Surv(time = data$Time, event = data$Event_final)~data$KIM1_follow_up_matrix.csv),data=data,pval=TRUE)
dev.off()

