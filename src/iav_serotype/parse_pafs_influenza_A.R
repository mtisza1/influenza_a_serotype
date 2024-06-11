#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=5) {
  stop("5 arguments must be supplied \
       flu info db, paf file, sample name, \
       out directory, score threshold.", call.=FALSE)
} else if (length(args)==5) {
  sprintf("arguments found. Running.")
}

Flu_A_DB <- args[1]

paf_file <- args[2]

samp <- args[3]

out_dir <- args[4]

score_thresh <- args[5]

paf_dt <- fread(paf_file, select = 1:12, sep = "\t", header = F, 
                fill = T, 
                col.names = 
                  c("qname", "qlength", "qstart", "qend", "strand",
                    "tname", "tlength", "tstart", "tend", 
                    "num_matches", "align_length", "mapq"))

flu_info_dt <- fread(Flu_A_DB, sep = "\t", header = T)

merge_dt <- merge(paf_dt, flu_info_dt,
                  by.x = "tname", by.y = "accession")

assigment_dt <- merge_dt %>%
  group_by(qname, tname, serotype, segment, strand) %>%
  summarize(tot_read_length = sum(qlength),
            tot_align = sum(align_length),
            tot_match = sum(num_matches),
            ANI = tot_match / tot_align,
            AF = tot_align / tot_read_length,
            align_score = ANI * AF) %>%
  arrange(qname, desc(align_score))


sum_dt <- assigment_dt %>%
  ungroup() %>%
  group_by(qname, serotype, segment) %>%
  summarize(n = n(),
            top_score = max(align_score),
            avg_score = mean(align_score)) %>%
  arrange(qname, desc(top_score)) %>%
  ungroup() %>%
  group_by(qname) %>%
  arrange(desc(top_score)) %>%
  slice_max(top_score, n = 2) %>%
  mutate(read_assignment = case_when(
    (max(top_score - 0.003)) >= min(top_score) |
      n_distinct(serotype) == 1 ~ first(serotype),
    TRUE ~ "ambiguous")) %>%
  filter(max(top_score) >= score_thresh)

write.table(sum_dt,
            file = sprintf("%s/%s_read_summary.tsv", 
                           out_dir, samp),
            quote = F, row.names = F, col.names = T,
            sep = "\t")

assignp <- sum_dt %>%
  slice_max(top_score, with_ties = F) %>%
  ggplot(aes(x = read_assignment)) +
  geom_bar() +
  theme(axis.text.x = 
          element_text(angle = 90, 
                       vjust = 0.5, 
                       hjust=1))

ggsave(assignp,
       file = sprintf("%s/%s_read_serotype_assignment.pdf", 
                      out_dir, samp)
)

sum_dt %>% 
  ungroup() %>%
  group_by(qname) %>%
  slice_max(top_score, with_ties = F) %>% 
  group_by(read_assignment) %>%
  select(read_assignment, qname) %>%
  group_walk(~ write.table(.x, paste0(out_dir, "/", samp, 
                                      "_", .y$read_assignment, ".txt"),
                           quote = F, row.names = F, col.names = F,
                           sep = "\t"))


