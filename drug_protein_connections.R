#' Title
#'
#' @param filePath a character vector, Full path where all the input files are located and the output to be saved
#' @param diseaseName a character vector, indicating the disease name
#' @param dr_Disease a vector of strings (character vectors), indicating all the drug names (generic names) involved in the disease
#'
#' @returns an dataframe as an .rda object, indicating all the drug names and involved gene symbols (targets, enzymes, carriers and transporters)
#' @export
#'
#' @examples
#' drugNames <- c("Aclidinium","Albuterol","Arformoterol","Budesonide","Carmoterol","Cefpodoxime")
#' drug_protein_mapping(filePath = "path/to/your/files/,
#' diseaseName = "COPD",
#' dr_Disease = drugNames)
drug_protein_mapping <- function(filePath, diseaseName = "COPD", dr_Disease) {
  library(dplyr)
  library(data.table)

  # filePath = "C:\\Users\\Hp\\OneDrive - Imam university\\Documents\\Az projects\\DrugBank Projects\\drugbank 5.1.13\\"

  # read drug targets
  dt = fread(paste0(filePath, "uniprot links_target.csv")) %>%
    dplyr::mutate("Protein_role" = "Target")

  # read drug enzymes
  denz = fread(paste0(filePath, "uniprot links_enzymes.csv")) %>%
    dplyr::mutate("Protein_role" = "Enzyme")

  # read drug carrier
  dcar = fread(paste0(filePath, "uniprot links_carrier.csv")) %>%
    dplyr::mutate("Protein_role" = "Carrier")

  # read drug transporter
  dtrans = fread(paste0(filePath, "uniprot links_transporter.csv")) %>%
    dplyr::mutate("Protein_role" = "Transporter")

  # read all the human protein-gs mapping file
  prot.df = fread(paste0(filePath, "HUMAN_9606_idmapping_Gene_Names.txt"))


  drug_all_Protein <- rbind(dt, denz, dcar, dtrans)

  drug_all_Protein %>% fwrite(paste0(filePath, "drug_related_proteins_all.csv"))



  # save all COPD drug-related proteins (Targets, Enzymes, Carrier, and Transporter)
  result <- drug_all_Protein %>%
    dplyr::filter(Name %in% dr_Disease) %>% # filter disease-related drug proteins
    dplyr::inner_join(x = .,
                      y = prot.df,
                      by = c("UniProt ID" = "UNIPROT"))
  # print statistics
  cat(
    paste0(
      "Total disease drugs: ",
      dr_Disease %>% length(),
      "\n",
      "Total drugs-related protein: ",
      drug_all_Protein %>%
        dplyr::filter(Name %in% dr_Disease) %>%
        dplyr::pull(`UniProt ID`) %>% unique() %>%
        length(),
      "\n",
      "Total drugs-related geneSymbol: ",
      result$SYMBOL %>%
        unique() %>%
        length()
    )
  )
  result %>%
    save(file = paste0(filePath, "DrugBank_", diseaseName, "_geneSymbols.rda"))
}

# # read Disease drugs only
# dr_Disease <- fread(paste0(filePath, "COPD_Drugs_Drugbank.csv")) %>%
#   dplyr::filter(Drug != "") %>%
#   dplyr::pull(Drug) %>%
#   unique()
