library(dplyr)
library(data.table)

filePath = "C:\\Users\\Hp\\OneDrive - Imam university\\Documents\\Az projects\\DrugBank Projects\\drugbank 5.1.13\\"

# read drug targets
dt = fread(paste0(filePath, "uniprot links_target.csv")) %>%
  dplyr::mutate("Protein_role" = "Target")

# read drug enzymes
denz = fread(paste0(filePath, "uniprot links_enzymes.csv"))%>%
  dplyr::mutate("Protein_role" = "Enzyme")

# read drug carrier
dcar = fread(paste0(filePath, "uniprot links_carrier.csv"))%>%
  dplyr::mutate("Protein_role" = "Carrier")

# read drug transporter
dtrans = fread(paste0(filePath, "uniprot links_transporter.csv"))%>%
  dplyr::mutate("Protein_role" = "Transporter")



rbind(
  dt, denz, dcar, dtrans
) %>% fwrite(paste0(filePath, "drug_related_proteins_all.csv"))

