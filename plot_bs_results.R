
#### INITIALIZE ####
list.of.packages <- c("vcfR", 
                      "tidyverse",
                      "janitor",
                      "ggrepel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(tidyverse)
library(ggrepel)
library(janitor)
library(vcfR)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("plot_R.r")

#### FUNCTION ####
get_newest_filename <-
  function(DataPath,
           FilePattern){
    list.files(path = DataPath,
               pattern = FilePattern,
               full.names = TRUE) %>%
      file.info() %>%
      rownames_to_column(var="file_name") %>%
      arrange(desc(mtime)) %>%
      head(1) %>%
      pull(var=file_name)
  }

#### USER DEFINED VARIABLES ####
fdrCUTOFF <- 0.05
fstFILE <- get_newest_filename(DataPath="./",
                               FilePattern="*fst.txt")
selFILE <- get_newest_filename(DataPath="./",
                               FilePattern="*sel")
# strFILE <- "../opihiSK2017_WholeIsland.A4.5.5.Fltr21.1.FilterOurWacky_noNA_noOutliers_Str"
idFILE <- "../sample_name_decode.tsv"
vcfPATH = "../data_files/Pfalcifer.3.7.Fltr20.2.randSNPperLoc.vcf"

#### FOLLOWING BAYESCAN MAN ####
# plot
results <- plot_bayescan(fstFILE, 
                         FDR=fdrCUTOFF)

#print # outliers and SNP ids
results$nb_outliers
results$outliers

mydata=read.table(selFILE,
                  colClasses="numeric")

parameter="Fst1"
plot(density(mydata[[parameter]]),xlab=parameter,main=paste(parameter,"posterior distribution"))

bs_plot <- plot_bayescan(fstFILE)

#### READ IN DATA ####

chromData <-
  vcfR2tidy(read.vcfR(vcfPATH))$fix %>%
  clean_names() %>%
  mutate(locus_id = row_number()) %>%
  select(locus_id,
         chrom,
         pos)

sampleMetaData <-
  read_tsv(idFILE)

posteriorDATA <- 
  read.table(selFILE,
             colClasses="numeric") %>%
  mutate(step = as.numeric(row.names(.))) %>%
  rename(Fst01=Fst1,
         Fst02=Fst2,
         Fst03=Fst3,
         Fst04=Fst4,
         # Fst05=Fst5,
         # Fst06=Fst6,
         # Fst07=Fst7,
         # Fst08=Fst8,
         # Fst09=Fst9
         ) %>%
  pivot_longer(cols = contains("Fst"),
               names_to = "sample",
               values_to = "fst") %>%
  mutate(sample = str_remove(sample, "Fst")) %>%
  left_join(sampleMetaData,
            by = c("sample" = "sample_id_bscan")) %>%
  select(step,
         logL,
         sample:date_collection)

fstData <- 
  read_delim(fstFILE,
             delim = " ",
             skip=1,
             col_names=FALSE) %>%
  rename(locus_id = X1,
         prob = X3,
         log10_po = X4,
         qval = X5,
         alpha = X6,
         fst = X7) %>%
  select(-X2,
         -X8) %>%
  #error need to convert fdrCUTOFF
  mutate(outlier = case_when(qval <= fdrCUTOFF ~ TRUE,
                             TRUE ~ FALSE),
         outlier_locus_id = case_when(qval <= fdrCUTOFF ~ locus_id,
                                       TRUE ~ NA_real_)) %>%
  left_join(chromData) %>%
  select(chrom,
         pos,
         locus_id,
         outlier_locus_id,
         prob:outlier)

#### PLOT FST POSTERIORS by SAMPLE ####
posteriorDATA %>%
  ggplot(aes(y=fst,
             x=step,
             color = sample_id_detailed)) +
  geom_line() +
  labs(title = "Trace Plot",
       subtitle = "There should be no trends.") +
  facet_grid(sample_id_detailed ~ .,
             scales = "free_y") 

posteriorDATA %>%
  ggplot(aes(x=fst,
             fill = sample_id_detailed)) +
  geom_histogram(bins=200,
                 alpha = 0.7) 

posteriorDATA %>%
  ggplot(aes(x=fst,
             fill = sample_id_detailed)) +
  geom_histogram(bins=200)+
  facet_wrap(sample_id_detailed ~ .)

#### PLOT OUTLIERS ####
fstData %>%
  ggplot(aes(x=qval,
             y=fst,
             color = outlier)) +
  geom_point(size = 5) + 
  geom_vline(xintercept = 0.05,
             linetype="dashed") + 
  geom_vline(xintercept = 0.01,
             linetype="dashed") + 
  geom_vline(xintercept = 0.001,
             linetype="dashed") +
  scale_x_continuous(trans='log10') +
  geom_label_repel(aes(label = outlier_locus_id))

fstData %>%
  filter(outlier == TRUE)

#### OUTPUT FST FILE W OUTLIER INFO ####
write_tsv(fstData,
          "snps_w_outliers.tsv")

