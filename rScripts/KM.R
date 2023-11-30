library(survival)
library(survminer)

data = read.csv("~/Fetal_dir/AminData/FME_lasso/updated_followup_with_new_clusters.csv", row.names =1)

pdf("~/Fetal_dir/AminData/FME_lasso/nonZeroFME.pdf", height = 6, width = 6)
ggsurvplot(survfit(Surv(time = data$Time, event = data$Event_final)~data$non_zeroFMEgenes_follow_up_matrix.csv),data=data,pval=TRUE)
dev.off()

pdf("~/Fetal_dir/AminData/FME_lasso/negativFME.pdf", height = 6, width = 6)
ggsurvplot(survfit(Surv(time = data$Time, event = data$Event_final)~data$negativeFMEgenes_follow_up_matrix.csv),data=data,pval=TRUE)
dev.off()

pdf("~/Fetal_dir/AminData/FME_lasso/positiveFME.pdf", height = 6, width = 6)
ggsurvplot(survfit(Surv(time = data$Time, event = data$Event_final)~data$positiveFMEgenes_follow_up_matrix.csv),data=data,pval=TRUE)
dev.off()

pdf("~/Fetal_dir/AminData/FME_lasso/100FME.pdf", height = 6, width = 6)
ggsurvplot(survfit(Surv(time = data$Time, event = data$Event_final)~data$mostneg_100FMEgenes_follow_up_matrix.csv),data=data,pval=TRUE)
dev.off()

pdf("~/Fetal_dir/AminData/FME_lasso/50FME.pdf", height = 6, width = 6)
ggsurvplot(survfit(Surv(time = data$Time, event = data$Event_final)~data$mostneg_50FMEgenes_follow_up_matrix.csv),data=data,pval=TRUE)
dev.off()

pdf("~/Fetal_dir/AminData/FME_lasso/10FME.pdf", height = 6, width = 6)
ggsurvplot(survfit(Surv(time = data$Time, event = data$Event_final)~data$mostneg_10FMEgenes_follow_up_matrix.csv),data=data,pval=TRUE)
dev.off()


pdf("~/Fetal_dir/AminData/FME_lasso/allFME.pdf", height = 6, width = 6)
ggsurvplot(survfit(Surv(time = data$Time, event = data$Event_final)~data$allFMEgenes_follow_up_matrix.csv),data=data,pval=TRUE)
dev.off()


pdf("~/Fetal_dir/AminData/FME_lasso/pos100FME.pdf", height = 6, width = 6)
ggsurvplot(survfit(Surv(time = data$Time, event = data$Event_final)~data$mostpos_100FMEgenes_follow_up_matrix.csv),data=data,pval=TRUE)
dev.off()

pdf("~/Fetal_dir/AminData/FME_lasso/pos50FME.pdf", height = 6, width = 6)
ggsurvplot(survfit(Surv(time = data$Time, event = data$Event_final)~data$mostpos_50FMEgenes_follow_up_matrix.csv),data=data,pval=TRUE)
dev.off()

pdf("~/Fetal_dir/AminData/FME_lasso/pos10FME.pdf", height = 6, width = 6)
ggsurvplot(survfit(Surv(time = data$Time, event = data$Event_final)~data$mostpos_10FMEgenes_follow_up_matrix.csv),data=data,pval=TRUE)
dev.off()

