  .m_df <- function(x, sname = "WCC"){
    rowname <- NULL
    name <- NULL
    datPPI <-
      x %>% as.data.frame %>% rownames_to_column() %>%
      pivot_longer(-rowname) %>% na.omit()%>%
      mutate(`PPI` =paste(`rowname`, `name`, sep = "~")) %>%
      select(4,3)
    colnames(datPPI)[2] <- sname
    return(datPPI)

  }
  .wtd_rank = function(mat) {
    ranks = apply(mat, 2, rank, ties = "average")
    # weight the ranks
    # calculate the savage scores
    n = nrow(mat)
    reciprocals = 1 / seq_len(n)
    savage = sapply(seq_len(n), function(i) sum(reciprocals[i:n]))
    # replace each rank with the savage score
    savages = ranks
    savages[] = savage[ranks]
    # calculate pearson correlation
    cor = cor(savages, method = "p")
    return(cor)
  }

  .int = function(mat) {
    present = !is.na(mat) & mat > 0
    intersect = crossprod(present)
    J = intersect
    return(J)
  }

  .wcc_f <- function(x,y)
    ptw::wcc(x, y, trwdth = 1)

  .jaccard = function(mat) {
    present = !is.na(mat) & mat > 0
    intersect = crossprod(present)
    union = nrow(mat) - crossprod(!present)
    J = (2 * intersect) / union
    return(J)
  }

  .dice = function(mat) {
    present = !is.na(mat) & mat > 0
    J = 1- as.matrix(arules::dissimilarity(t(present), method = 'dice'))
    return(J)
  }


  .coapex_f <- function(x,y)
    ifelse((max(x) == max(y)), 1, 0)

  # Bayes correlation
  .Bayes_Corr <- function(alpha0, beta0, X){
    nrowsX <- nrow(X)
    k <- ncol(X)
    cs <- colSums(X)
    alphas <- matrix(rep(alpha0,nrowsX), nrow=nrowsX, byrow=TRUE) + X
    betas  <- matrix(rep(beta0,nrowsX), nrow=nrowsX, byrow=TRUE) +
      matrix(rep(cs,nrowsX), nrow=nrowsX, byrow=TRUE) - X
    alphasPLUSbetas <- alphas + betas

    # First BIG product term for covariance formula
    Psi <- alphas/alphasPLUSbetas -
      matrix(rep(rowSums(alphas/alphasPLUSbetas)/k, k), ncol=k,
             byrow=FALSE)

    # Covariance matrix
    cov_mtrx <- Psi %*% t(Psi) / k

    # Variances (this is a column vector of length = nrowsX)
    var_vec <-
      as.matrix( ( rowSums( (alphas*betas)/
                              ( (alphasPLUSbetas^2)*
                                  (alphasPLUSbetas+1) ) )
                   + rowSums(Psi^2) )/k )

    Bcorrvals <- cov_mtrx / sqrt( var_vec %*% t(var_vec) )
    diag(Bcorrvals) <- 1
    return(Bcorrvals)
  }

  # Computing the Bayesian correlations assuming third (zero count-motivated)
  # prior
  .Bayes_Prior3 <- function(X){
    d <- dim(X)
    cs <- colSums(X)
    alpha0 <- (cs+1)/(max(cs)+1)
    beta0 <- rep(1,d[2])
    Bcorrvals <- .Bayes_Corr(alpha0,beta0,X)
    return(Bcorrvals)
  }


  .PCCN <- function(x, rept = 10){



    A <- t(x)
    M <-
      ncol(A)
    N <-
      nrow(A)
    i <- 0
    # message("Compute PCC ...")
    pb <-
      txtProgressBar(min = 1, max = rept, style = 3)
    repeat{
      A.rpoisson <-
        apply(A, c(1, 2), function(x){rpois(1, lambda = x)})
      C.rpoisson <-
        A.rpoisson + 1/M
      B.rpoisson <-
        C.rpoisson/rowSums(C.rpoisson)
      i <- i + 1
      setTxtProgressBar(pb, i)
      B.cor <-
        cor(B.rpoisson, use = "pairwise.complete.obs")
      B.cor[is.na(B.cor)] <-0
      if(i == 1){
        PCC.mat <- B.cor
      }else{
        PCC.mat <- PCC.mat + B.cor
      }
      if(i == rept){
        break
      }
    }
    close(pb)
    PCC.mat.avg <-
      PCC.mat/rept
    return(PCC.mat.avg)

  }

  .merge.all <- function(x, y) {
    merge(x, y, all=TRUE, by="PPI")
  }



  #' calculate_PPIscore
  #' @title Calculate Pairwise Protein Profile Similarity using Different
  #' Metrics
  #' @param x A co-elution data matrix with proteins in rows and
  #' fractions in columns.
  #' @param pcc If TRUE, computes pairwise protein profile similarity
  #'  using Pearson correlation metric.
  #' @param PCCN If TRUE, computes pairwise protein profile similarity using
  #' Pearson correlation plus noise.This function is adapted from the
  #' PCCN function in the SMED package.
  #' @param rept Poisson iterations, defaults to 10.
  #' @param pcc_p If TRUE, computes P-value of the Pearson correlation.
  #' @param cosine If TRUE, computes pairwise protein profile similarity
  #' using cosine metric.
  #' @param minfo If TRUE, computes pairwise protein profile similarity
  #' using mutual information.
  #' @param jaccard If TRUE, computes pairwise protein profile similarity
  #' using jaccard metric.
  #' @param apex If TRUE, computes pairwise protein profile similarity
  #' using apex.
  #' @param bayesian If TRUE, computes pairwise protein profile similarity using
  #' Bayes correlation based on zero-count distribution.
  #' @param wcc If TRUE, computes pairwise protein profile similarity
  #' using weighted cross correlation.
  #' @param spearman if TRUE, computes pairwise protein profile similarity using
  #' spearman correlation.
  #' @param kendall if TRUE, computes pairwise protein profile similarity using
  #' kendall correlation.
  #' @param bicor if TRUE, computes pairwise protein profile similarity using
  #' biweight midcorrealtion (bicor) correlation.
  #' @param weighted_rank if TRUE, computes pairwise protein profile similarity
  #' using weighted rank measure.
  #' @param dice if TRUE, computes pairwise protein profile similarity
  #' using dice measure.
  #' @param euclidean if TRUE, computes pairwise protein profile similarity
  #' using euclidean measure.
  #' @param manhattan if TRUE, computes pairwise protein profile similarity
  #' using manhattan measure.
  #' @param canberra if TRUE, computes pairwise protein profile similarity
  #' using canberra measure.
  #' @param avg.distance if TRUE, computes pairwise protein profile similarity
  #' using avg.distance measure.
  #' @param corr_removal If TRUE, removes protein pairs with
  #' correlation scores < the user defined threshold ; defaults to FALSE.
  #' @param corr_cutoff user defined threshold for correlation similarity
  #' scores. Defaults to 0.5.
  #' @return A data frame containing the calculated features for all possible
  #' protein pairs.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @importFrom Hmisc rcorr
  #' @importFrom stats cor
  #' @importFrom utils setTxtProgressBar
  #' @importFrom stats rpois
  #' @importFrom utils txtProgressBar
  #' @importFrom dplyr mutate
  #' @importFrom dplyr select
  #' @importFrom dplyr %>%
  #' @importFrom lsa cosine
  #' @importFrom WGCNA mutualInfoAdjacency
  #' @importFrom tidyr pivot_longer
  #' @importFrom tibble rownames_to_column
  #' @importFrom stats na.omit
  #' @importFrom stats dist
  #' @description This function first removes proteins pairs for which two
  #' proteins never occurred in the same fractions, then computes
  #' pairwise protein similarity using up to 18 metrics (by default all the 18
  #' measures are activated). This function also
  #' provides users with an option to choose an appropriate co-fractionation
  #' correlation score cut-off using the `corr_cutoff` argument,
  #' if argument `corr_removal` is set to TRUE.
  #' @export
  #' @examples
  #' M1<-matrix(rnorm(36),nrow=6)
  #' M1 <- abs(M1)
  #' rownames(M1) <- c("A","B","C","D","E","F")
  #' scored_Data <- calculate_PPIscore(M1)




  calculate_PPIscore <-
    function(x,
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
             corr_cutoff = 0.5){

      if (!is.matrix(x)) {
        x <- as.matrix(x)
      }

      if(is.character(x) == TRUE){
        stop("matrix must include numerical variables")
      }

      if(is.null(row.names(x))){
        stop("Please specify the row.names")
      }

      intersection <- NA
      PPI <- NA

      x <- x[sort(rownames(x)), ]


      f_matrices <- list()


      if(pcc) {
        m <- cor(t(x), method = "pearson",
                 use = "pairwise.complete.obs")
        row.names(m) <- colnames(m) <- rownames(x)
        m[lower.tri(m, diag=TRUE)] <- NA
        pcc_d <-.m_df(m, sname = "PCC")
        f_matrices[["pcc"]] <- pcc_d
        message("Pearson scoring completes ...")

      }
      if(pcc_p) {

        m <-
          rcorr(t(x))$P
        row.names(m) <- colnames(m) <- rownames(x)
        m[lower.tri(m, diag=TRUE)] <- NA
        pcc_p_d <- .m_df(m, sname = "PCC_P")
        f_matrices[["pcc_p"]] <- pcc_p_d

        message("Pearson-P value scoring completes ...")


      }

      if(spearman) {
        m <- cor(t(x), method = "spearman",
                 use = "pairwise.complete.obs")
        row.names(m) <- colnames(m) <- rownames(x)
        m[lower.tri(m, diag=TRUE)] <- NA
        spear_d <- .m_df(m, sname = "spearman")
        f_matrices[["spear_d"]] <- spear_d

        message("Spearman scoring completes ...")


      }

      if(kendall) {
        m <- cor(t(x), method = "kendall",
                 use = "pairwise.complete.obs")
        row.names(m) <- colnames(m) <- rownames(x)
        m[lower.tri(m, diag=TRUE)] <- NA
        kendall_d <- .m_df(m, sname = "kendall")
        f_matrices[["kendall_d"]] <- kendall_d

        message("Kendall scoring completes ...")


      }

      if (bicor) {
        m = WGCNA::bicor(t(x), use = 'pairwise.complete.obs')
        m[lower.tri(m, diag=TRUE)] <- NA
        bicor_d <- .m_df(m, sname = "bicor")
        f_matrices[["bicor_d"]] <- bicor_d

        message("Biweight midcorrelation (bicor) scoring completes ...")


      }

      if (weighted_rank) {
        m = .wtd_rank(t(x))
        m[lower.tri(m, diag=TRUE)] <- NA
        weighted_rank_d <- .m_df(m, sname = "weighted_rank")
        f_matrices[["weighted_rank_d"]] <- weighted_rank_d

        message("Weighted rank correlation scoring completes ...")

      }

      if(wcc){
        m <-
          as.matrix(proxy::dist(x,.wcc_f))
        m[lower.tri(m, diag=TRUE)] <- NA
        wcc_d <- .m_df(m, sname = "WCC")
        f_matrices[["wcc"]] <- wcc_d

        message("Weighted cross correlation scoring completes ...")

      }

      if(cosine){
        m = cosine(t(x))
        m[lower.tri(m, diag=TRUE)] <- NA
        cos_d <- .m_df(m, sname = "Cosine")
        f_matrices[["cosine"]] <- cos_d
        message("Cosine scoring completes ...")

      }
      if(jaccard) {
        m = .jaccard(t(x))
        m[lower.tri(m, diag=TRUE)] <- NA
        jacc_d <- .m_df(m, sname = "Jaccard")
        f_matrices[["jaccard_d"]] <- jacc_d

        message("Jaccard scoring completes ...")

      }
      if(dice) {
        m = .dice(t(x))
        m[lower.tri(m, diag=TRUE)] <- NA
        dice_d <- .m_df(m, sname = "dice")
        f_matrices[["dice_d"]] <- dice_d

        message("Dice scoring completes ...")

      }


      if(apex) {
        m <-
          as.matrix(proxy::dist(x, .coapex_f))
        m[lower.tri(m, diag=TRUE)] <- NA
        cppx <- .m_df(m, sname = "co_apex")
        f_matrices[["coapex"]] <- cppx
        message("Apex scoring completes ...")

      }

      if(bayesian) {
        m <- .Bayes_Prior3(x)
        m[lower.tri(m, diag=TRUE)] <- NA
        bayes_d <- .m_df(m, sname = "Bayes")
        f_matrices[["bayesian"]] <- bayes_d

        message("Bayesian scoring completes ...")


      }
      if(minfo) {
        m = mutualInfoAdjacency(t(x))$AdjacencyUniversalVersion1
        m[lower.tri(m, diag=TRUE)] <- NA
        MI_d <- .m_df(m, sname = "MI")
        f_matrices[["minfo"]] <- MI_d

        message("Mutual information scoring completes ...")

      }
      if (manhattan) {
        m =  as.matrix(dist(x, method = "manhattan"))
        m[lower.tri(m, diag=TRUE)] <- NA
        manhattan_d <- .m_df(m, sname = "manhattan")
        manhattan_d$manhattan =
          manhattan_d$manhattan/max(manhattan_d$manhattan)
        manhattan_d$manhattan =
          1- manhattan_d$manhattan
        f_matrices[["manhattan_d"]] <- manhattan_d

        message("Manhattan distance scoring completes ...")

      }

      if (canberra) {
        m =  as.matrix(dist(x, method = "canberra"))
        m[lower.tri(m, diag=TRUE)] <- NA
        canberra_d <- .m_df(m, sname = "canberra")
        canberra_d$canberra =
          canberra_d$canberra/max(canberra_d$canberra)
        canberra_d$canberra =
          1- canberra_d$canberra
        f_matrices[["canberra_d"]] <- canberra_d
        message("Canberra distance scoring completes ...")

      }

      if(avg.distance) {
        m =  as.matrix(philentropy::distance(x, method = "avg"))
        row.names(m) <- colnames(m) <- rownames(x)
        m[lower.tri(m, diag=TRUE)] <- NA
        avg.distance_d <- .m_df(m, sname = "avg.distance")
        avg.distance_d$avg.distance =
          avg.distance_d$avg.distance/max(avg.distance_d$avg.distance)
        avg.distance_d$avg.distance =
          1- avg.distance_d$avg.distance

        f_matrices[["avg.distance_d"]] <- avg.distance_d

        message("Average distance scoring completes ...")
      }

      if (euclidean) {
        m = as.matrix(dist(x, method = 'euclidean'))
        m[lower.tri(m, diag=TRUE)] <- NA
        euclidean_d <- .m_df(m, sname = "euclidean")
        euclidean_d$euclidean =
          euclidean_d$euclidean/max(euclidean_d$euclidean)
        euclidean_d$euclidean =
          1- euclidean_d$euclidean
        f_matrices[["euclidean_d"]] <- euclidean_d

        message("Euclidean scoring completes ...")

      }

      if(PCCN){
        m <- .PCCN(x, rept = rept)
        m[lower.tri(m, diag=TRUE)] <- NA
        PCCN_d <- .m_df(m, sname = "PCCN")
        f_matrices[["PCCN"]] <- PCCN_d
        message("Pearson with noise scoring completes ...")

      }

      datPPI <-
        Reduce(.merge.all, f_matrices)

      # remove pairs never observed in the same fraction
        m = .int(t(x))
        m[lower.tri(m, diag=TRUE)] <- NA
        int_d <-
          .m_df(m, sname = "intersection") %>%
          filter(intersection >= 1)
        datPPI <-
          datPPI %>%
          filter(PPI %in% int_d$PPI)

      if(corr_removal){
      # remove PPI with at least 0.5 cut-off across all the scores
      c <- colnames(datPPI)
      dnames <- c("PPI", "PCC_P","co_apex")
      c <- setdiff(c, dnames)
      datPPI <- datPPI[which(rowSums(datPPI[,c] >= corr_cutoff) >= 1), ]
      }

      return(datPPI)
    }


