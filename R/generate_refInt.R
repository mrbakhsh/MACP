
  # create adjacency matrix from data frame
  .adj_m_d <- function(x){

    m <-
      adjacency_matrix <- #create adjacency matrix
      get.adjacency(
        graph_from_data_frame(x, directed = FALSE),
        attr = "class",
        sparse = FALSE
      )
    diag(m) <- 1

    return(m)
  }


  # Create refint_interactions
  .refint <- function(x) {
    values <- NA

    ppiCpx <-
      lapply(x, function(cpx){
        m <-
          apply(combn(cpx, 2),2, sort)
        s <-
          paste(m[1, ], m[2, ], sep = "~")
        return(s)
      })

    dat.refINT <-
      aggregate(data = stack(ppiCpx),
                `ind` ~ `values`,
                FUN = paste,
                collapse = ";")

    dat.refINT <- dat.refINT %>%
      separate(values, c("Prot_1", "Prot_2"), sep = "~") %>%
      select(1,2) %>%
      mutate(class = 1)

    return(dat.refINT)
  }


  #' generate_refInt
  #' @title Generate Class Labels for Data Input Based on Gold Reference Set
  #' @param x A data frame with interacting proteins in the first two
  #' columns.
  #' @param refcpx A list containing gold standard protein complexes.
  #' @return A Data frame containing class labels for protein pairs in the
  #' data input. If protein pairs involve in same protein complexes are
  #' assigned to Positive, otherwise Negative.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @importFrom stats aggregate
  #' @importFrom utils stack
  #' @importFrom tidyr separate
  #' @importFrom igraph graph_from_data_frame
  #' @importFrom igraph get.adjacency
  #' @importFrom stats na.omit
  #' @importFrom utils combn
  #' @importFrom dplyr sample_n
  #' @description This function creates class labels for protein pairs in the
  #' same order in the data input based on gold reference set.
  #' @export


  generate_refInt <- function(x,
                              refcpx){
    . <- NA
    if(!is.list(refcpx))
      stop("Refcpx Must Be List")

    node_columns = c(1, 2)
    col1 <- node_columns[1]
    col2 <- node_columns[2]
    prot_1 <- as.character(x[[col1]])
    prot_2 <- as.character(x[[col2]])
    w_prot <- union(prot_1,prot_2)


    # Keep cpx having same members as our data
    refcpx <-
      lapply(refcpx, function (x) x[x %in% w_prot])

    refcpx <-
      refcpx[lapply(refcpx,length)>=2]

    if(length(refcpx) == 0)
      stop("No Overlap Between Data Input & gold reference set ...")

    # Create interactions
    refintDF <- .refint(refcpx)

    # Create adjacency matrix from data.frame
    adj_m <-  .adj_m_d(refintDF)


    class.lab_index <-
      prot_1 %in% rownames(adj_m) &
      prot_2 %in% rownames(adj_m)

    idxing_dat <- cbind(prot_1[class.lab_index], prot_2[class.lab_index])
    label <- rep(NA, nrow(x))
    label[class.lab_index] <- adj_m[idxing_dat]

    data_labeled <- cbind(x, label)
    data_labeled <-
      data_labeled %>%
      mutate(PPI = paste(.[,1], .[,2], sep = "~")) %>%
      select(4,3) %>%
      mutate(label = ifelse(label == 1, "Positive",
                            "Negative"))
    data_labeled <- na.omit(data_labeled)

    return(data_labeled)
  }


