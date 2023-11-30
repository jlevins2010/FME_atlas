library(data.table)
library(rmeta)

all_traits = c("sumstats_formatted_PASS_Type_2_Diabetes", "sumstats_formatted_PASS_Type_1_Diabetes", "sumstats_formatted_PASS_Coronary_Artery_Disease","sumstats_formatted_PASS_BMI1", "sumstats_formatted_PASS_heightZ", "sumstats_formatted_PASS_CystatinC", "sumstats_formatted_PASS_Creatinine", "sumstats_formatted_PASS_SYSTOLICadjMEDz","sumstats_formatted_PASS_DIASTOLICadjMEDz", "sumstats_formatted_PASS_HbA1c", "sumstats_formatted_PASS_HYPERTENSION_DIAGNOSED","sumstats_formatted_PASS_Urea","sumstats_formatted_PASS_Calcium","sumstats_formatted_PASS_Phosphate", "sumstats_formatted_PASS_VitaminD", "sumstats_formatted_PASS_T2D")

DM_traits = c("sumstats_formatted_PASS_Type_2_Diabetes", "sumstats_formatted_PASS_Type_1_Diabetes", "sumstats_formatted_PASS_HbA1c","sumstats_formatted_PASS_T2D")

endoPhenotypes = c("sumstats_formatted_PASS_BMI1", "sumstats_formatted_PASS_heightZ", "sumstats_formatted_PASS_CystatinC", "sumstats_formatted_PASS_Creatinine", "sumstats_formatted_PASS_SYSTOLICadjMEDz","sumstats_formatted_PASS_DIASTOLICadjMEDz", "sumstats_formatted_PASS_Urea","sumstats_formatted_PASS_Calcium","sumstats_formatted_PASS_Phosphate", "sumstats_formatted_PASS_VitaminD")

DiseasePheno = c("sumstats_formatted_PASS_Coronary_Artery_Disease", "sumstats_formatted_PASS_HYPERTENSION_DIAGNOSED")

annot_cell = "/home/levinsj/Fetal_dir/Linker/"
results_cell = "/home/levinsj/Fetal_dir/Linker/baselineLD_v2.1_annots"
base_path = "/home/levinsj/Fetal_dir/Analysis/referenceFiles/LDSC/baselineLD_v2.1_annots"
annot_names = list.files(results_cell)

run_enrichment_analysis = function(annot_cell, results_cell, base_path,annotation, traits, index_in_results=NULL, base_index = NULL){
	base <- data.frame(fread(paste0("zcat ", base_path, "/", "baselineLD.", 22, ".annot.gz")))
	if(is.null(base_index)){
		base_index = ncol(base) - 4
    }
    cell_path = paste0(annot_cell, "/", annotation)
    res = paste0(results_cell, "/", annotation, "/ABC_Road_KID/", traits[1], ".sumstats.results")
	tab2 = read.table(res,header=T)
	if(is.null(index_in_results)){index_in_results = 1:(nrow(tab2) - base_index)}
	enrich_table = matrix(0, length(index_in_results), 3)
	cat("Number of annotations together with baseline : ", nrow(tab2), "\n")
	annot_names = as.character(tab2$Category[index_in_results])
	Mref = 5961159
	for(id in 1:length(index_in_results)){
		meta_enr = NULL;
		meta_enrstat    = NULL;
		for(trait_id in 1:length(traits)){
			result.file=paste0(results_cell, "/", annotation, "/ABC_Road_KID/", traits[trait_id], ".sumstats.results")
			res=read.table(result.file,header=T)
			logfile = paste(results_cell, "/", annotation, "/ABC_Road_KID/", traits[trait_id],".sumstats.log", sep="")
			log = read.table(logfile,h=F,fill=T)
			h2g = as.numeric(as.character(log[which(log$V4=="h2:"),5]))
			myenrstat = (h2g/Mref)*((res[index_in_results[id],3]/res[index_in_results[id],2])-(1-res[index_in_results[id],3])/(1-res[index_in_results[id],2])) #step1
			myenrstat_z  = qnorm(res[index_in_results[id],7]/2) #step2
			myenrstat_sd = myenrstat/myenrstat_z #step3
			meta_enrstat = rbind(meta_enrstat, c(myenrstat, myenrstat_sd));
			meta_enr  = rbind(meta_enr, c(res[index_in_results[id],5],res[index_in_results[id],6] ));
			cat("meta_enr", meta_enr, " ", traits[trait_id])								    }
		test_eni1=meta.summaries(meta_enr[,1], meta_enr[,2],method="random")
		test_eni2=meta.summaries(meta_enrstat[,1], meta_enrstat[,2],method="random")
		cat("Printing enrichment results for annotation:", annot_names[id], "\n")
		cat(test_eni1$summary, " ",test_eni1$se.summary, " ", 2*pnorm(-abs(test_eni2$summary/test_eni2$se.summary)), "\n");
		enrich_table[id, ] = c(test_eni1$summary, test_eni1$se.summary, 2*pnorm(-abs(test_eni2$summary/test_eni2$se.summary)))
	}
	rownames(enrich_table) = annot_names
	return(enrich_table)
}

library(data.table)
library(rmeta)

out2 = c()
for(m in 1:length(annot_names)){
	out2 = rbind(out2, run_enrichment_analysis(annot_cell, results_cell,base_path = base_path, annotation = annot_names[m], traits = all_traits[6], index_in_results = 1))
}

rownames(out2) = annot_names
save(out2, file = "/home/levinsj/Fetal_dir/Linker/enrichmentScores_enrich.rda")
