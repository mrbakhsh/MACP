

  #' orthMappingCpx
  #' @title Protein Complex Ortholog Mapping
  #' @param datInput A list containing reference complexes (i.e.,
  #' CORUM complexes). Note that the members of each complexes must be
  #' represented by UniProt accession identifier.
  #' @param input_species Name of the input species (e.g., "mouse","fly").
  #' See \code{\link[orthogene]{map_species}} to return a full list of
  #' available species.
  #' @param output_species Name of the output species (e.g., "human").
  #' See \code{\link[orthogene]{map_species}} to return a full list of
  #' available species.
  #' @param input_taxid A numeric value that specifies the NCBI
  #' taxonomy identifier (TaxId) for input organism (e.g., 10090).
  #' @param output_taxid A numeric value that specifies the NCBI
  #' taxonomy identifier (TaxId) for output organism.
  #' @return A list containing complexes, whose members
  #' converted to output_species.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @description This function uses \code{\link[orthogene]{convert_orthologs}}
  #' function to support ortholog mapping of protein complexes between any pair
  #' of 700+ species.
  #' @export
  #' @examples
  #' data("refcpx")
  #' orth_mapping <- orthMappingCpx (refcpx,
  #' input_species = "mouse",
  #' output_species = "human",
  #' input_taxid = "10090",
  #' output_taxid = "9606")





  orthMappingCpx <-
    function(datInput,
             input_species,
             output_species,
             input_taxid,
             output_taxid){

      if(!is.list(datInput)){
        stop("Reference complexes must be list")
      }

      . <- NULL
      V1 <- NULL
      id <- NULL
      input_gene <- NULL
      conv_gene <- NULL

      # convert list to data frame
      agg_mem <-
        lapply(datInput, function(x) paste(x, collapse = ";"))


      df_cpx <-
        do.call(rbind, agg_mem) %>%
        as.data.frame(.) %>%
        rownames_to_column("id") %>%
        mutate(V1 = strsplit(as.character(V1), ";"))  %>%
        unnest(V1) %>%
        dplyr::rename(accession=V1) %>%
        arrange(id) %>%
        na.omit(.)

      # fetch proteome of input species
      prot_input <- protti::fetch_uniprot_proteome(input_taxid,
                                             columns = c("accession",
                                                         "gene_primary"),
                                             reviewed = TRUE)

      df_cpx_g <-
        inner_join(df_cpx, prot_input, by = "accession")

      # now orthologous mapping
      orth_Df <- orthogene::convert_orthologs(unique(df_cpx_g$gene_primary),
                                   input_species = input_species,
                                   output_species = output_species)

      orth_Df <-
        orth_Df %>%
        tibble::rownames_to_column("conv_gene") %>%
        rename(gene_primary = input_gene) %>%
        select(1,3)



      df_cpx_g2 <-
        inner_join(df_cpx_g, orth_Df, by = "gene_primary") %>%
        select(1,4) %>%
        rename(gene_primary=conv_gene)


      # fetch proteome of output species
      prot_output <- protti::fetch_uniprot_proteome(output_taxid,
                                                   columns = c("accession",
                                                               "gene_primary"),
                                                   reviewed = TRUE)

      fout <-
        inner_join(df_cpx_g2,prot_output,by = "gene_primary") %>%
        select(1,3)


      #conver to list of cpx
      fout_cpx <- aggregate(list(fout$accession),
                            by = list(fout$id),
                            FUN = function(x) paste(unique(x), collapse = " "))
      colnames(fout_cpx) <- c("name", "member")

      fout_cpx$member <-
        vapply(lapply(strsplit(fout_cpx$member, " "), unique),
               paste, character(1L), collapse = " ")
      fout_cpx <-
        fout_cpx %>%
        filter(lengths(strsplit(fout_cpx$member, " ")) >= 3)

      fout_cpx_list <-
        strsplit(fout_cpx$member, " ")
      names(fout_cpx_list) <- fout_cpx$name

      # remove redudancy
      fout_cpx_list <-
        EliminateCpxRedundance(fout_cpx_list)

      return(fout_cpx_list)



    }


