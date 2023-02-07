## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages('MACP')

## ---- eval=FALSE--------------------------------------------------------------
#  
#  if(!requireNamespace("devtools", quietly = TRUE)) {
#  install.packages("devtools")
#  }
#  devtools::install_github("BabuLab-UofR/MACP")

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(MACP)
library(dplyr)
library(tidyr)

## -----------------------------------------------------------------------------
# Loading the demo data
data(exampleData)
dim(exampleData)
# Inspect the data 
glimpse(exampleData)

## -----------------------------------------------------------------------------
data_p1 = data_filtering(exampleData)
# Inspect the number of retained proteins 
dim(data_p1)

## -----------------------------------------------------------------------------
x <- data_p1
# Assign column 10 to zeros
x[,10] <- NA
data_p2 <- impute_MissingData(x)

## -----------------------------------------------------------------------------
data_p3 <- scaling(data_p1)

## -----------------------------------------------------------------------------
data_p3 <- keepMT(data_p3)
# Inspect the number of retained proteins 
dim(data_p3)

## ----message = FALSE, warning = FALSE-----------------------------------------
set.seed(100)
scored_PPI <- calculate_PPIscore(data_p3,
                                corr_removal = FALSE)

## -----------------------------------------------------------------------------
data("refcpx")

## -----------------------------------------------------------------------------
# separate the interaction pairs
PPI_pairs <-
  scored_PPI %>%
  separate(PPI, c("p1", "p2"), sep = "~") %>% select(p1,p2)

# Generate reference set of positive and negative interactions
class_labels <-
    generate_refInt(PPI_pairs[,c(1,2)],refcpx)
table(class_labels$label)

## ---- eval= FALSE-------------------------------------------------------------
#  # for example to convert mouse complexes to human complexes
#  ## load the mouse complexes
#  data("refcpx")
#  orth_mapping <- orthMappingCpx (refcpx,
#    input_species = "mouse",
#    output_species = "human",
#    input_taxid = "10090",
#    output_taxid = "9606")

## ----warning = FALSE, results="hide", message=FALSE---------------------------
set.seed(101)
predPPI_ensemble <- 
  ensemble_model(scored_PPI,
                  class_labels,
                  classifier = c("glm", "svmRadial", "ranger"),
                  cv_fold = 5,
                  plots = FALSE,
                  verboseIter = FALSE,
                  filename=file.path(tempdir(),"plots.pdf"))

# Subset predicted interactions 
pred_interactions <- predPPI_ensemble$predicted_interactions

## -----------------------------------------------------------------------------
roc_object <-
  inner_join(class_labels, pred_interactions, by ="PPI")
roc_object <-
  pROC::roc(roc_object$label,roc_object$Positive)
pROC::coords(roc_object, x="best", input="threshold", best.method="youden")

# Extract high-confidence network based on the cut-off reported from ROC curve
# Note that best-threshold can change depending on the input data
ThreshNet_PPI <- filter(pred_interactions, Positive > 0.5) 

## ----warning = FALSE, message=FALSE-------------------------------------------
ThreshNet_PPI <- 
  separate(ThreshNet_PPI, PPI, c("p1","p2"), sep = "~")
denoisedPPI <- get_DenoisedNet(ThreshNet_PPI)
dim(denoisedPPI)

## -----------------------------------------------------------------------------
finalPPI <- subcellular.mtPPI(denoisedPPI, organism = "mouse")
dim(finalPPI)

## ---- eval = FALSE------------------------------------------------------------
#  finalPPI <- finalPPI[, -4] # drop the last column
#  # Set the directory to your current directory
#  setwd("user's current directory")
#  write.table(finalPPI, file = "ppi_input_ ClusterONE.txt",
#              quote = FALSE,
#              col.names = F, row.names = F, sep = "\t")

## -----------------------------------------------------------------------------
pred_cpx <- get_clusters(csize = 2, 
                         d = 0.3, p = 2,
                         max_overlap = 0.8,
                         tpath =file.path(system.file("extdata", 
                                                      package = "MACP")))
dim(pred_cpx)

## -----------------------------------------------------------------------------
pred_cpx_mcl <- 
  MCL_clustering(finalPPI, # High-confidence interactions
                 pred_cpx, # Putative complexes produced by clusterONE
                 inflation = 9, 
                 csize =2)
dim(pred_cpx_mcl)

## -----------------------------------------------------------------------------
# first load the reference complex
data("refcpx")
Clust_tuning_result <-
  cluster_tuning(refcpx, csize = 3, 
                d = c(0.3,0.4),
                p = c(2, 2.5),
                max_overlap = c(0.6,0.7),
                tpath =
                  file.path(system.file("extdata", package = "MACP")))

## -----------------------------------------------------------------------------
pred_cpx_optimized <- get_clusters(csize = 2, 
                         d = 0.4, p = 2.5,
                         max_overlap = 0.7,
                         tpath =file.path(system.file("extdata", 
                                                      package = "MACP")))
dim(pred_cpx_optimized)

## -----------------------------------------------------------------------------
mcl_tuning_result <- 
  MCL_tuning(finalPPI,
             pred_cpx_optimized, 
             refcpx,
             inflation = c(6,8,9,10))

## -----------------------------------------------------------------------------
final_clusters <- 
  MCL_clustering(finalPPI,
             pred_cpx_optimized, 
             inflation = 10, 
             csize = 2)
dim(final_clusters)

## ---- warning=FALSE, message=FALSE, eval=FALSE--------------------------------
#  # extract the predicted complexes
#  enrich_result <-
#    enrichmentCPX(pred_cpx,
#                  threshold = 0.05,
#                  sources = "GO:BP",
#                  p.corrction.method = "bonferroni",
#                  custom_bg = NULL,
#                  org = "mmusculus")
#  head(enrich_result[, c(1,4,12)], n = 4)

## ----message = FALSE, warning = FALSE, results=FALSE--------------------------
# Load the input data
data("exampleData")
# Known reference complexes
data("refcpx")
# Perform prediction
Prediction_output <-
  predPPI_MACP(exampleData,
  refcpx,
  keepMT =TRUE, # keep mt proteins
  subcellular_mtPPI = TRUE,
  tpath = tempdir())


## ---- eval=FALSE--------------------------------------------------------------
#  ig <-
#     igraph::graph_from_data_frame(Prediction_output$filteredPPI)
#  RCy3::createNetworkFromIgraph(ig,"myIgraph")

## -----------------------------------------------------------------------------
sessionInfo()

