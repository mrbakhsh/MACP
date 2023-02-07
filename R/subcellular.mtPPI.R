  #' subcellular.mtPPI
  #' @title Keep Mitochondrial (mt) Proteins
  #' @param ppi Interactions data containing id1-id2-weight triplets.
  #' @param organism Organism under study (i.e., mouse or human).
  #' Defaults to mouse.
  #' @return Filtered PPI netwrok.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @description This function removes PPIs occurring between outer mt
  #' membrane (OMM) and matrix, between intermembrane space (IMS) and matrix,
  #' as well as between any subcellular mt compartment (except OMM) and
  #' cytosolic proteins as they deemed to be erroneous
  #' @export
  #' @examples
  #' ppi <-
  #' read.table(system.file("extdata/ppi_input_ClusterONE.txt",
  #' package="MACP"),
  #' quote="\"", comment.char="")
  #' filtered_mtEdges <- subcellular.mtPPI(ppi)


  subcellular.mtPPI <- function(ppi, organism = "mouse") {

    mito.mito <- NULL
    sub_mito1 <- NULL
    sub_mito2 <- NULL


    colnames(ppi) <-
      c("p1", "p2", "edgeweight")

    # first open mitochondrial proteins
    mito_prot <- read.csv(
      system.file("extdata/MitoCarta_protein.csv", package = "MACP"),
      header = TRUE
    )

    back_p <- union(ppi$p1, ppi$p2)
    overlap_p <- intersect(back_p, mito_prot$p)
    if(length(overlap_p) == 0){
      stop("There is no overlap between data input and MitoCarta database")
    }



    if(organism == "mouse"){

    # open mouse sub.compartmetn
    m_sub.comp <- read.csv(
      system.file("extdata/m_mito_subcellcompartment.csv", package = "MACP"),
      header = TRUE
    )
    overlap_p <- intersect(back_p, m_sub.comp$UniProt)
    if(length(overlap_p) == 0){
      stop("There is no overlap between data input and MitoCarta database")
    }
    colnames(m_sub.comp)[1] <- "p1"
    ppi <- inner_join(ppi, m_sub.comp, by = "p1")
    colnames(m_sub.comp)[1] <- "p2"
    ppi <- inner_join(ppi, m_sub.comp, by = "p2")
    colnames(ppi)[4:5] <- c("sub_mito1", "sub_mito2")
    ppi <-
      unite(ppi, mito.mito, c(sub_mito1, sub_mito2), sep = "~")
    ppi_f <-
      ppi %>%
      filter(mito.mito != "Matrix~IMS") %>%
      filter(mito.mito != "IMS~Matrix") %>%
      filter(mito.mito != "MOM~Matrix") %>%
      filter(mito.mito != "Matrix~MOM")
    }

    if(organism == "human"){

      # open mouse sub.compartmetn
      h_sub.comp <- read.csv(
        system.file("extdata/h_mito_subcellcompartment.csv", package = "MACP"),
        header = TRUE
      )
      overlap_p <- intersect(back_p, h_sub.comp$p1)
      if(length(overlap_p) == 0){
        stop("There is no overlap between data input and MitoCarta database")
      }

      colnames(h_sub.comp)[1] <- "p1"
      ppi <- inner_join(ppi, h_sub.comp, by = "p1")
      colnames(h_sub.comp)[1] <- "p2"
      ppi <- inner_join(ppi, h_sub.comp, by = "p2")
      colnames(ppi)[4:5] <- c("sub_mito1", "sub_mito2")
      ppi <-
        unite(ppi, mito.mito, c(sub_mito1, sub_mito2), sep = "~")
      ppi_f <-
        ppi %>%
        filter(mito.mito != "Matrix~IMS") %>%
        filter(mito.mito != "IMS~Matrix") %>%
        filter(mito.mito != "MOM~Matrix") %>%
        filter(mito.mito != "Matrix~MOM")
    }


    return(ppi_f)

  }


