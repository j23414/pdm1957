#! /usr/bin/env Rscript

library(magrittr)

data <- readr::read_delim("influenza_na.dat", delim="\t", 
                          col_names = c("gb", "host", "segnum", "subtype","country", "date", "len","strain","x1","x2","x3"))

cdata <- data %>%
  subset(date<1965) %>%
  subset(!grepl("Influenza B", strain)) %>%
  subset(!grepl("Influenza C", strain)) %>%
  subset(!grepl("unidentified influenza virus", strain)) %>%
  dplyr::mutate(
    #len=NULL,
    x1=NULL,
    x2=NULL,
    strain=gsub("Influenza A virus \\(","",strain), 
    #%>% 
      #gsub("Influenza B virus \\(","",.) %>%
      #gsub("Influenza C virus \\(","",.),
    strain=gsub("\\(.*","", strain) %>% gsub("\\)", "", .)
  ) %>% 
  dplyr::mutate(
    HH=dplyr::case_when(segnum==4~gb),
    NN=dplyr::case_when(segnum==6~gb),
    PB2=dplyr::case_when(segnum==1~gb),
    PB1=dplyr::case_when(segnum==2~gb),
    PA=dplyr::case_when(segnum==3~gb),
    NP=dplyr::case_when(segnum==5~gb),
    M=dplyr::case_when(segnum==7~gb),
    NS=dplyr::case_when(segnum==8~gb)
  )%>%
  dplyr::group_by(strain, host, subtype, country, date) %>%
  dplyr::summarise(
    HH=toString(HH) %>% gsub("NA, ","",.) %>% gsub(", NA","", .),
    NN=toString(NN) %>% gsub("NA, ","",.) %>% gsub(", NA","", .),
    HH=toString(HH) %>% gsub("NA, ","",.) %>% gsub(", NA","", .),
    HH=toString(HH) %>% gsub("NA, ","",.) %>% gsub(", NA","", .),
    HH=toString(HH) %>% gsub("NA, ","",.) %>% gsub(", NA","", .),
    HH=toString(HH) %>% gsub("NA, ","",.) %>% gsub(", NA","", .),
    HH=toString(HH) %>% gsub("NA, ","",.) %>% gsub(", NA","", .),
    HH=toString(HH) %>% gsub("NA, ","",.) %>% gsub(", NA","", .),
  )

writexl::write_xlsx(cdata,"pdm1957.xlsx")
