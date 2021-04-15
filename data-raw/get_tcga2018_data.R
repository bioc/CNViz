library(cBioPortalData)
library(dplyr)

cbio <- cBioPortal()
studies <- getStudies(cbio)

#get studyIds
study_names <-
  data.frame(cbind(gsub(" \\(TCGA, PanCancer Atlas\\)", "",studies$name[grepl("tcga_pan_can_atlas", studies$studyId)]),
                   studies$studyId[grepl("tcga_pan_can_atlas", studies$studyId)]))
colnames(study_names) <- c("Cancer", "studyId")

#loop through studyIds to get all_tcga2018_data
all_tcga2018_data <- data.frame()
for(study in study_names$studyId){
  cbio_table <- getDataByGenePanel(cbio, study, genePanelId = "IMPACT468",
                          molecularProfileId = paste0(study, "_gistic"),
                     sampleListId = paste0(study, "_cna"))
  cbio_dat <- data.frame(cbio_table[[1]], stringsAsFactors = FALSE)
  cbio_summ <- cbio_dat %>% group_by(hugoGeneSymbol) %>%
    summarise(Gain = sum(value ==1)/n(),
              Amplification = sum(value == 2)/n(),
              ShallowDeletion = sum(value == -1)/n(),
              DeepDeletion = sum(value == -2)/n())
  cbio_summ$sample_size <- rep(paste0(" (N = ", length(unique(cbio_dat$uniquePatientKey)), ")"), nrow(cbio_summ))
  cbio_summ$studyId <- rep(study, nrow(cbio_summ))
  all_tcga2018_data <- rbind(all_tcga2018_data, cbio_summ)
}

#add sample sizes to get cbio_studies
ss <- all_tcga2018_data %>% group_by(studyId, sample_size) %>% summarise()
cbio_studies <- inner_join(study_names, ss, by = c("studyId"))
cbio_studies$Cancer <- paste0(cbio_studies$Cancer, cbio_studies$sample_size)
cbio_studies <- dplyr::select(cbio_studies, Cancer, studyId)
