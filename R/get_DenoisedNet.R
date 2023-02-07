
  #' get_DenoisedNet
  #' @title Denoising Predicted Protein-Protein Interactions
  #' @param ppi Interactions data containing id1-id2-weight triplets.
  #' @return  A data.frame containing denoised network.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @importFrom dplyr rowwise
  #' @importFrom stats median
  #' @importFrom dplyr left_join
  #' @description This function removes the noise in the form of false
  #' positive edges in the predicted networks using network topology.
  #' @export
  #' @examples
  #' # high-confidence network as input
  #' ppi <-
  #' read.table(system.file("extdata/ppi_input_ClusterONE.txt",
  #' package="MACP"),
  #' quote="\"", comment.char="")
  #' # Perform network denoising
  #' denoisetNet <- get_DenoisedNet(ppi)

  get_DenoisedNet <-
    function(ppi) {

      . <- NA
      Int1.Edges <- NA
      Int2.Edges <- NA
      edgeweight <- NA
      hscore <- NA
      intersection <- NA
      ints <- NA
      lengthp1 <- NA
      lengthp2 <- NA
      p1 <- NA
      p2 <- NA

      colnames(ppi) <-
        c("p1", "p2", "edgeweight")

      InteractorAEdges <- aggregate(ppi[,2] ~ ppi[,1],
                                    ppi,
                                    FUN = function(x)
                                      paste(unique(x), collapse = ";"))
      colnames(InteractorAEdges) <- c("p1","Int1.Edges")
      InteractorAEdges$Int1.Edges <-
        vapply(lapply(strsplit(InteractorAEdges$Int1.Edges, ";"), unique),
               paste, character(1L), collapse = ";")

      InteractorBEdges <- aggregate(ppi[,1] ~ ppi[,2],
                                    ppi,
                                    FUN = function(x)
                                      paste(unique(x), collapse = ";"))
      colnames(InteractorBEdges) <- c("p2","Int2.Edges")
      InteractorBEdges$Int2.Edges <-
        vapply(lapply(strsplit(InteractorBEdges$Int2.Edges, ";"), unique),
               paste, character(1L), collapse = ";")


      denoisedNet <-
        left_join(ppi, InteractorAEdges , by = "p1") %>%
        left_join(., InteractorBEdges, "p2") %>%
        mutate(ints = mapply(function(x,y)
          paste(length(intersect(x, y)), collapse = ";"),
          strsplit(.$Int1.Edges, ";"), strsplit(.$Int2.Edges, ";"))) %>%
        mutate(ints = as.numeric(ints)) %>%
        mutate(lengthp1 = as.numeric(lengths(strsplit(Int1.Edges, ";")))) %>%
        mutate(lengthp2 = as.numeric(lengths(strsplit(Int2.Edges, ";")))) %>%
        rowwise() %>%
        mutate(hscore = min(ints/lengthp1,ints/lengthp2)) %>%
        filter(hscore > 0) %>%
        select(`p1`,`p2`,`edgeweight`)

      denoisedNet$edgeweight <- as.numeric(denoisedNet$edgeweight)


      return(denoisedNet)

    }
