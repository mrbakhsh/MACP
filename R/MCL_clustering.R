  .mappedPPI_cpx <- # this function map each complexes to edge weight
    function(m, complexes) {

      output <- list()

      #p_value <- numeric(0)
      for (i in seq_along(complexes)) {
        complex_name <- names(complexes)[i]
        complex <- complexes[[i]]

        # subset complex to proteins present in this network
        nodes <- colnames(m)
        overlap <- intersect(complex, nodes)


        # calculate median PCC for intra-complex interactions
        idxing_mat <- t(combn(overlap, 2))
        edge_weights <- na.omit(m[idxing_mat])
        edge_weights <- as.data.frame(edge_weights)
        d <- cbind(idxing_mat,edge_weights )
        d <- filter(d, edge_weights != 0)

        output[[i]] <- d

      }

      return(output)
    }

  .m_list <- # convert df to adjacency matrix
    function(x) {

      output <- list()

      #p_value <- numeric(0)
      for (i in seq_along(x)) {

        m <-
          igraph::get.adjacency(
            igraph::graph_from_data_frame(x[[i]], directed = FALSE),
            attr = "edge_weights",
            sparse = FALSE
          )

        output[[i]] <- m

      }

      return(output)
    }

  .cpx_duplciate_remove <- # remove subset MCL clustering
    function(MCLclust) {

      . <- NA
      size <- NA
      ID <- NA
      ID1 <- NA
      ClustID <- NA

    df <-
      MCLclust %>%
      mutate(size = lengths(strsplit(.$Members, "\\s+"))) %>%
      arrange(size) %>% select(1,2)

    list.cpx <- #split tge list
      strsplit(df$Members, "\\s+")
    names(list.cpx) <- #chnage the name of the complexes (can be skiped)
      df$ClustID
    nms <-  ### all the name combination
      combn( names(list.cpx), 2,
             FUN = paste0, collapse = "~",
             simplify = FALSE)
    complex.list <-  #all the cpx combination
      combn( list.cpx , 2 , simplify = FALSE )
    intersect <- lapply( #get the intersection
      complex.list ,
      function(x)
        length( intersect( x[[1]] , x[[2]] ) ) )
    length.complex1 <- lapply(#get the length of complex1 1
      complex.list ,
      function(x) length(( x[[1]]) ) )
    Names <-
      data.frame(matrix(unlist(nms),
                        nrow = length(nms), byrow = T))
    Intersect <-
      data.frame(matrix(unlist(intersect),
                        nrow = length(intersect), byrow = T))
    Length.cpx1 <-
      data.frame(matrix(unlist(length.complex1),
                        nrow = length(length.complex1), byrow = T))
    dfcomb <- cbind(Names,Intersect,Length.cpx1)
    colnames(dfcomb) <- c("ID", "Intersect","Length.cpx1")
    dfcomb <-
      dfcomb %>%
      filter(intersect > 0) %>%
      separate( ID, c("ID1","ID2"),
                sep = "~", remove = FALSE)

    complete.subset1 <-
      dfcomb %>%
      group_by(ID1) %>%
      arrange(Length.cpx1) %>%
      filter(Intersect == Length.cpx1)
    s <- unique(complete.subset1$ID1)
    dff <- df[!df$ClustID %in% s,]
    dff$ClustID <-
      as.numeric(dff$ClustID)
    dff <-
      dff %>%
      arrange(ClustID) %>%
      select(-1) %>% rowid_to_column("ClustID")

    return(dff)
  }


  #' MCL_clustering
  #' @title MCL clustering
  #' @param hc_ppi High-confidence interactions data containing
  #' id1-id2-weight triplets.
  #' @param predcpx A data.frame containing predicted complexes resulted from
  #' \code{\link[MACP]{get_clusters}}.
  #' @param inflation MCL inflation parameter. Defaults to 9.
  #' @param csize  An integer, the minimum size of the predicted complexes.
  #' Defaults to 2.
  #' @return List of refined complexes.
  #' @description This function applies MCL clustering to further refine
  #' the predicted subnetworks produced by ClusterONE.
  #' @importFrom tidyr unnest
  #' @author Matineh Rahmatbakhsh
  #' @export
  #' @importFrom tibble rowid_to_column
  #' @examples
  #' # open high-confidence network
  #' hc_ppi <-
  #' read.delim(
  #' system.file("extdata/ppi_input_ClusterONE.txt", package = "MACP"),
  #' header = FALSE)
  #' # predict complexes by ClusterONE
  #' predcpx <-
  #' get_clusters(csize = 3, d = 0.3, p = 2,
  #' max_overlap = 0.8,
  #' tpath = file.path(system.file("extdata", package = "MACP")))
  #' # Break down big complexes by MCL
  #' MCL_clusters <- MCL_clustering(hc_ppi, predcpx, inflation = 4, csize = 2)




  MCL_clustering <-
    function(hc_ppi,
             predcpx,
             inflation = 9,
             csize = 2){

      clustid <- NULL
      id <- NULL
      . <- NULL
      V1 <- NULL
      Members <- NA
      . <- NA

      colInput <-
        c("ClustID", "Members")
      if(!all(colInput %in% colnames(predcpx))){
        missingCol <-
          setdiff(colInput,
                  colnames(predcpx)[match(colInput, colnames(predcpx))])
        stop("Input data missing: ", paste(missingCol, collapse = ", "))
      }

      if(!is.data.frame(predcpx)){
        predcpx <- as.data.frame(predcpx)
      }
      if(!is.data.frame(hc_ppi)){
        hc_ppi <- as.data.frame(hc_ppi)
      }


      colnames(hc_ppi)[3] <- "V3"

      # break down big complexes
      predcpx_breakdown <-
        predcpx %>%
        filter(lengths(strsplit(Members, " ")) >= 15)



      # convert to matrix
      m <-
        adjacency_matrix <- #create adjacency matrix
        igraph::get.adjacency(
          igraph::graph_from_data_frame(hc_ppi, directed = FALSE),
          attr = "V3",
          sparse = FALSE
        )


      cpx_l <- strsplit(predcpx_breakdown$Members, "\\s+")
      names(cpx_l) <-
        predcpx_breakdown$ClustID

      #map to edge weight
      ss <-
        .mappedPPI_cpx(m, cpx_l)
      d <-
        .m_list(ss)

      # retrive column names
      cols <-
        lapply(d, function(x) colnames(x))
      cols <- lapply(cols, function(x) paste(x, collapse = ";"))

      df.cols <-
        do.call(rbind, cols) %>%
        as.data.frame(.) %>%
        rowid_to_column("id") %>%
        mutate(V1 = strsplit(as.character(V1), ";"))  %>%
        unnest(V1) %>%
        dplyr::rename(uniprot=V1) %>%
        arrange(id) %>%
        na.omit(.)



      ## apply mcl clustering
      mcl.clust <-
        lapply(d, function(m)
          MCL::mcl(m, addLoops = T,
                   inflation = inflation,
                   allow1 = FALSE,
                   max.iter = 100, ESM = FALSE))
      membership2 <-
        lapply(mcl.clust, function(x) x[["Cluster"]])
      mem <- lapply(membership2, function(x) paste(x, collapse = ";"))


      df <-
        do.call(rbind, mem) %>%
        as.data.frame(.) %>%
        rowid_to_column("id") %>%
        mutate(V1 = strsplit(as.character(V1), ";"))  %>%
        unnest(V1) %>%
        dplyr::rename(member=V1) %>%
        arrange(id) %>%
        na.omit(.)


      df <-
        cbind(df, df.cols[,2])
      colnames(df)[2:3] <-
        c("clustid","Member")

      df <-
        df %>%
        filter(clustid > 0)

      MCLclust <- aggregate(list(df[,3]),
                       by = list(df$id, df$clustid),
                       FUN = function(x) paste(unique(x), collapse = " "))
      MCLclust <-
        MCLclust %>% rownames_to_column("ClustID") %>%
        dplyr::select(1,4)
      colnames(MCLclust)[2] <- "Members"
      MCLclust$Members <-
        vapply(lapply(strsplit(MCLclust$Members, " "), unique),
               paste, character(1L), collapse = " ")

      predcpx_notfilt <-
        predcpx %>%
        filter(lengths(strsplit(Members, " ")) < 15)

      # final clusters
      MCLclust<-
        rbind(predcpx_notfilt,MCLclust)

      #remove redundancy
      MCLclust <- .cpx_duplciate_remove(MCLclust)

      #size filtering
      MCLclust <-
        MCLclust %>%
        filter(lengths(strsplit(MCLclust$Members, " ")) >= csize)




      return(MCLclust)
    }
