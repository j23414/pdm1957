#! /usr/bin/env Rscript

library(magrittr)

data <- readr::read_delim("influenza_na.dat", delim = "\t", 
                          col_names = c("gb", "host", "segnum", "subtype",
                                        "country", "date", "len", "strain",
                                        "x1", "x2", "x3"))

cdata <- data %>%
  subset(date < 1966) %>%
  subset(!grepl("Influenza B", strain)) %>%
  subset(!grepl("Influenza C", strain)) %>%
  subset(!grepl("unidentified influenza virus", strain)) %>%
  dplyr::arrange(desc(len)) %>%
  dplyr::mutate(
    x1 = NULL,
    x2 = NULL,
    strain = gsub("Influenza A virus \\(","",strain), 
    strain = gsub("\\(.*","", strain) %>% gsub("\\)", "", .)
  ) %>% 
  dplyr::mutate(
    HH = dplyr::case_when(segnum==4 ~ gb),
    NN = dplyr::case_when(segnum==6 ~ gb),
    PB2 = dplyr::case_when(segnum==1 ~ gb),
    PB1 = dplyr::case_when(segnum==2 ~ gb),
    PA = dplyr::case_when(segnum==3 ~ gb),
    NP = dplyr::case_when(segnum==5 ~ gb),
    M = dplyr::case_when(segnum==7 ~ gb),
    NS = dplyr::case_when(segnum==8 ~ gb)
  ) %>%
  dplyr::group_by(strain, host, subtype, country, date) %>%
  dplyr::summarise(
    HH  = toString(HH)  %>% gsub("NA, ", "", .) %>% gsub(", NA", "", .) %>% gsub("NA", "", .),
    NN  = toString(NN)  %>% gsub("NA, ", "", .) %>% gsub(", NA", "", .) %>% gsub("NA", "", .),
    PB2 = toString(PB2) %>% gsub("NA, ", "", .) %>% gsub(", NA", "", .) %>% gsub("NA", "", .),
    PB1 = toString(PB1) %>% gsub("NA, ", "", .) %>% gsub(", NA", "", .) %>% gsub("NA", "", .),
    PA  = toString(PA)  %>% gsub("NA, ", "", .) %>% gsub(", NA", "", .) %>% gsub("NA", "", .),
    NP  = toString(NP)  %>% gsub("NA, ", "", .) %>% gsub(", NA", "", .) %>% gsub("NA", "", .),
    M   = toString(M)   %>% gsub("NA, ", "", .) %>% gsub(", NA", "", .) %>% gsub("NA", "", .),
    NS  = toString(NS)  %>% gsub("NA, ", "", .) %>% gsub(", NA", "", .) %>% gsub("NA", "", .)
  ) %>%
  subset(HH!="" & NN!="" & PB1!="" & PB2!="" & PA!="" & NP!="" & M!="" & NS!="")

writexl::write_xlsx(cdata,"pdm1957.xlsx")

# === Generate the genbank id list for analysis

for (seg in c("HH", "NN", "PB2", "PB1", "PA", "NP", "M", "NS")){
  cat(paste(seg, " is processing...\n"))
  out <- data.frame(seg=cdata[[seg]]) %>%
    dplyr::mutate(
      seg = gsub(",.*", "", seg)
    )
  readr::write_delim(out, path = paste(seg, "ids", sep = "."), col_names = F)
}


