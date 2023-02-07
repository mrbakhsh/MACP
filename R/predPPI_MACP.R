  #' predPPI_MACP
  #' @title Predict Protein-Protein Interactions and Putative Complexes
  #' @param data A data matrix with rows including proteins and fractions
  #' along the columns. see \code{\link{exampleData}}.
  #' @param refcpx A list of known reference complexes.
  #' see \code{\link{getCPX}}.
  #' @param tpath A character string indicating the path to the project
  #' directory. If the directory is
  #' missing, it will be stored in the Temp directory.
  #' @param data_processing If TRUE, removes proteins for which peptide only
  #' detected in one fraction (i.e., "one-hit-wonders") across the
  #' co-elution table, common contaminants (e.g., keratins) only for mouse
  #' and human organisms and frequent flyers. Defaults to TRUE.
  #' See \code{\link{data_filtering}}.
  #' @param data_imputing if TRUE, imputes missing values in protein elution
  #' profile matrix via average of adjacent rows. This function is
  #' not applicable for missing values present in the first or last column.
  #' Defaults to TRUE. See \code{\link{impute_MissingData}}.
  #' @param keepMT if TRUE, removes all the non-mitochondrial proteins by
  #' mapping the co-eluted proteins from chromatography fractions to MitoCarta
  #' database. Note that this function is only applicable to
  #' mouse or human organisms.Defaults to FALSE. See \code{\link{keepMT}}.
  #' @param scaling If TRUE, performs column and row-wise normalization.
  #' Defaults to TRUE. See \code{\link{scaling}}.
  #' @param pcc If TRUE, computes pairwise protein profile similarity
  #'  using Pearson correlation metric. Defaults to TRUE.
  #'  See \code{\link{calculate_PPIscore}}.
  #' @param PCCN If TRUE, computes pairwise protein profile similarity using
  #' Pearson correlation plus noise. Defaults to TRUE.
  #' See \code{\link{calculate_PPIscore}}.
  #' @param rept Poisson iterations, defaults to 10. Defaults to TRUE.
  #' See \code{\link{calculate_PPIscore}}.
  #' @param pcc_p If TRUE, computes P-value of the Pearson correlation.
  #' Defaults to TRUE.
  #' See \code{\link{calculate_PPIscore}}.
  #' @param cosine If TRUE, computes pairwise protein profile similarity
  #' using cosine metric. Defaults to TRUE.
  #' See \code{\link{calculate_PPIscore}}.
  #' @param minfo If TRUE, computes pairwise protein profile similarity
  #' using mutual information. Defaults to TRUE.
  #' See \code{\link{calculate_PPIscore}}.
  #' @param jaccard If TRUE, computes pairwise protein profile similarity
  #' using jaccard metric. Defaults to TRUE.
  #' See \code{\link{calculate_PPIscore}}.
  #' @param apex If TRUE, computes pairwise protein profile similarity
  #' using apex. Defaults to TRUE. See \code{\link{calculate_PPIscore}}.
  #' @param bayesian If TRUE, computes pairwise protein profile similarity using
  #' Bayes correlation based on zero-count distribution. Defaults to TRUE.
  #' See \code{\link{calculate_PPIscore}}.
  #' @param wcc If TRUE, computes pairwise protein profile similarity
  #' using weighted cross correlation. Defaults to TRUE.
  #' See \code{\link{calculate_PPIscore}}.
  #' @param spearman if TRUE, computes pairwise protein profile similarity using
  #' spearman correlation. Defaults to TRUE.
  #' See \code{\link{calculate_PPIscore}}.
  #' @param kendall if TRUE, computes pairwise protein profile similarity using
  #' kendall correlation. Defaults to TRUE.
  #' See \code{\link{calculate_PPIscore}}.
  #' @param bicor if TRUE, computes pairwise protein profile similarity using
  #' biweight midcorrealtion (bicor) correlation. Defaults to TRUE.
  #' See \code{\link{calculate_PPIscore}}.
  #' @param weighted_rank if TRUE, computes pairwise protein profile similarity
  #' using weighted rank measure. Defaults to TRUE.
  #' See \code{\link{calculate_PPIscore}}.
  #' @param dice if TRUE, computes pairwise protein profile similarity
  #' using dice measure. Defaults to TRUE.
  #' See \code{\link{calculate_PPIscore}}.
  #' @param euclidean if TRUE, computes pairwise protein profile similarity
  #' using euclidean measure. Defaults to TRUE.
  #' See \code{\link{calculate_PPIscore}}.
  #' @param manhattan if TRUE, computes pairwise protein profile similarity
  #' using manhattan measure. Defaults to TRUE.
  #' See \code{\link{calculate_PPIscore}}.
  #' @param canberra if TRUE, computes pairwise protein profile similarity
  #' using canberra measure. Defaults to TRUE.
  #' See \code{\link{calculate_PPIscore}}.
  #' @param avg.distance if TRUE, computes pairwise protein profile similarity
  #' using avg.distance measure. Defaults to TRUE.
  #' See \code{\link{calculate_PPIscore}}.
  #' @param corr_removal If TRUE, removes protein pairs with
  #' correlation scores < the user defined threshold; defaults to FALSE.
  #' See \code{\link{calculate_PPIscore}}.
  #' @param corr_cutoff user defined threshold for correlation similarity
  #' scores. Defaults to 0.5.See \code{\link{calculate_PPIscore}}.
  #' @param classifier The type of classifier to use for ensemble or
  #' individual model. See \code{caret} for the available classifiers.
  #' Defaults to c("glm", "svmRadial", "ranger").
  #' See \code{\link{ensemble_model}}.
  #' @param verboseIter Logical value, indicating whether to check the status
  #' of training process;defaults to FALSE. See \code{\link{ensemble_model}}.
  #' @param cv_fold  Number of partitions for cross-validation; defaults to 5.
  #' See \code{\link{ensemble_model}}.
  #' @param plots Logical value, indicating whether to plot the performance
  #' of the learning algorithm using k-fold cross-validation;
  #' defaults to FALSE. These plots are :
  #' \itemize{ \item{pr_plot} - Precision-recall PLOT
  #'   \item{roc_plot} - ROC plot
  #'   \item{point_plot} - Point plot showing accuracy,
  #'   F1-score , positive predictive value (PPV), sensitivity (SE) and MCC.}
  #'   See \code{\link{ensemble_model}}.
  #' @param subcellular_mtPPI if TRUE, removes PPIs occurring between outer mt
  #' membrane (OMM) and matrix, between intermembrane space (IMS) and matrix,
  #' as well as between any subcellular mt compartment (except OMM) and
  #' cytosolic proteins as they deemed to be erroneous. Defaults to FALSE.
  #' See \code{\link{subcellular.mtPPI}}.
  #' @param organism Organism under study (i.e., mouse or human).
  #' Defaults to mouse. See \code{\link{subcellular.mtPPI}}.
  #' @param csize An integer, the minimum size of the predicted complexes.
  #' Defaults to 2. See \code{\link{get_clusters}}.
  #' @param d A number, density of predicted complexes. Defaults to 0.3.
  #' See \code{\link{get_clusters}}.
  #' @param p An integer, penalty value for the inclusion of each node.
  #' Defaults to 2.
  #' See \code{\link{get_clusters}}.
  #' @param max_overlap A number, specifies the maximum allowed
  #' overlap between two clusters. Defaults to 0.8.
  #' See \code{\link{get_clusters}}.
  #' @param inflation MCL inflation parameter. Defaults to 9.
  #' @return Return following data sets in the current directory including:
  #'  \itemize{
  #'  \item{unfilteredPPIs} - Unfiltered interactions
  #'  \item{filteredPPI} - High-confidence interactions
  #'  defined by ROC threshold.
  #'  \item{High_confidence interactions_with_mt_sublocalization} - if
  #'  subcellular.mtPPI is TRUE, it return high-confidene PPIs with mt
  #'  sublocalization status.
  #'  \item{predicted_cpx_clusterONE} - Putative complexes generated
  #'  by clusterONE.
  #'  \item{predicted_cpx_clusterONE_MCL} - Putative complexes generated
  #'  by clusterONE and MCL.
  #'  \item{Best_roc_curve_cutoff} - Best cutoff generated from ROC curve.}
  #' @description This function first begins by executing several
  #' pre-processing steps to improve the quality of the raw data, followed
  #' by computing similarity between protein pairs using their co-elution
  #' profiles. Computed features and and class labels generated
  #' from reference complexes are then fed into an individual or
  #' ensemble of ML classifiers.These models then generate a weighted protein
  #' interaction network in which edge weights between protein nodes represent
  #' the ML model's probability estimate for interaction.
  #' High-confidence PPIs resulted from ROC-curve cutoff analysis is then
  #' denoised and finally are partitioned via two-stage clustering,
  #' first by ClusterONE,then by MCL clustering.
  #' @importFrom utils write.table
  #' @export
  #' @examples
  #' data("exampleData")
  #' data("refcpx")
  #' output <-
  #' predPPI_MACP(exampleData,
  #' refcpx,
  #' keepMT =TRUE, # keep mt proteins
  #' subcellular_mtPPI = TRUE,
  #' tpath =tempdir())# Path for the output file




  predPPI_MACP <-
    function(data,
             refcpx,
             tpath = tempdir(),
             data_processing = TRUE,
             data_imputing = TRUE,
             scaling = TRUE,
             keepMT = FALSE,
             pcc = TRUE,
             PCCN = TRUE,
             pcc_p = TRUE,
             spearman = TRUE,
             kendall = TRUE,
             bicor = TRUE,
             weighted_rank = TRUE,
             cosine = TRUE,
             jaccard = TRUE,
             dice = TRUE,
             apex = TRUE,
             minfo = TRUE,
             bayesian = TRUE,
             wcc = TRUE,
             euclidean = TRUE,
             manhattan = TRUE,
             canberra = TRUE,
             avg.distance = TRUE,
             rept = 10,
             corr_removal = FALSE,
             corr_cutoff = 0.5,
             classifier = c("glm", "svmRadial", "ranger"),
             verboseIter = TRUE,
             cv_fold = 5,
             plots = FALSE,
             subcellular_mtPPI = FALSE,
             organism = "mouse",
             csize = 3,
             d = 0.3,
             p = 2,
             max_overlap = 0.8,
             inflation = 9){



      if (!is.matrix(data)) {
        data <- as.matrix(data)
      }

      if(is.character(data) == TRUE){
        stop("matrix must include numerical variables")
      }

      if(is.null(row.names(data))){
        stop("Please specify the row.names")
      }

      if(!is.list(refcpx)){
        stop("Reference complexes must be list")
      }

      PPI <- NULL
      Positive <- NULL
      p1 <- NULL
      p2 <- NULL
      pred_dat <- NULL



      data_n <- data

      if(data_processing) {
        data_n <-
          data_filtering(data_n)
        message("Number of retained proteins:", nrow(data_n))
      }
      if(data_imputing){
        data_n <-
          impute_MissingData(data_n)
        message("Imputation completes ...")
      }
      if(scaling){
        data_n <-
          scaling(data_n)
        message("Scaling completes ...")
      }
      if(keepMT){
        data_n <-
          keepMT(data_n)
        message("Number of mitochondrial (mt) proteins:",  nrow(data_n))
      }


      # compute similarity
      scored_PPI <-
        calculate_PPIscore(data_n,
                           pcc = pcc,
                           PCCN = PCCN,
                           pcc_p = pcc_p,
                           spearman = spearman,
                           kendall = kendall,
                           bicor = bicor,
                           weighted_rank = weighted_rank,
                           cosine = cosine,
                           jaccard = jaccard,
                           dice = dice,
                           apex = apex,
                           minfo = minfo,
                           bayesian = bayesian,
                           wcc = wcc,
                           euclidean = euclidean,
                           manhattan = manhattan,
                           canberra = canberra,
                           avg.distance = avg.distance,
                           rept = rept,
                           corr_removal = corr_removal,
                           corr_cutoff = corr_cutoff)
      message("Computing similarity completes ...")

      # separate the interaction pairs
      PPI_pairs <-
        scored_PPI %>%
        separate(PPI, c("p1", "p2"), sep = "~") %>% select(p1,p2)

      # now generate reference interactions
      class_labels <-
        generate_refInt(PPI_pairs[,c(1,2)],refcpx)

      message("Building class labels completes ...")


      output <- list()

      predPPI <- ensemble_model(scored_PPI,
                                  class_labels,
                                  cv_fold = cv_fold,
                                  classifier = classifier,
                                  verboseIter = verboseIter,
                                  plots = plots,
                                  filename = file.path(tpath, "plots.pdf"))
      message("Prediction completes ...")

      # cut-off selection based on ROC-curve
      pred_interactions <- predPPI$predicted_interactions
      roc_object <-
        inner_join(class_labels, pred_interactions, by ="PPI")
      roc_object <-
        pROC::roc(roc_object$label,roc_object$Positive)
      s <- pROC::coords(roc_object, x="best", input="threshold",
                        best.method="youden")
      message("The ROC-curve cutoff is:", s$threshold)

      # Extract high-confidence network based on the cut-off reported from ROC curve
      pred_dat_filt <-
        filter(pred_interactions, Positive >= s$threshold) %>%
        separate(PPI, c("p1", "p2"), sep = "~")
      message("Network denoising...")
      pred_dat_filt <- get_DenoisedNet(pred_dat_filt)

      if(subcellular_mtPPI) {

        message("filtering non-relevant mtPPIs...")
        mt_filt <-
          subcellular.mtPPI(pred_dat_filt, organism = organism)
        pred_dat_filt <- mt_filt[, -4]
      }


      # save the output files
      fname <- file.path(tpath, "unfilteredPPIs.txt")
      write.table(predPPI$predicted_interactions, file = fname,
                  row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

      fname <- file.path(tpath, "ppi_input_ClusterONE.txt")
      write.table(pred_dat_filt,
                  file = fname,
                  row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)


      # generate clusters by clusterONE
      predcpx <-
        get_clusters(csize =csize,
                     d= d,
                     p = p,
                     max_overlap = max_overlap,
                     tpath = tpath)
      message("ClusterONE clustering completes...")


      fname <- file.path(tpath, "predicted_complexes_clusterONE.txt")
      write.table(predcpx, file = fname,
                  row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

      # generate clusters by MCL
      final_clusters <-
        MCL_clustering(pred_dat_filt,
                       predcpx,
                       inflation = inflation,
                       csize = csize)
      message("MCL clustering completes...")
      fname <- file.path(tpath, "predicted_complexes_clusterONE_MCL.txt")
      write.table(final_clusters, file = fname,
                  row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)



      output[["unfilteredPPIs"]] <- predPPI$predicted_interactions
      output[["filteredPPI"]] <- pred_dat_filt
      output[["High_confidence interactions_with_mt_sublocalization"]] <- mt_filt
      output[["predicted_cpx_clusterONE"]] <- predcpx
      output[["predicted_cpx_clusterONE_MCL"]] <- final_clusters
      output[["Best_roc_curve_cutoff"]] <- s$threshold




      return(output)

    }

