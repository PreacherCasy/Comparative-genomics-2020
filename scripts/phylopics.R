library(dplyr)
library(ggplot2)
library(ggtree)
library(magrittr)
library(pheatmap)

#### Basic statistics
## Per genome gene distribution
gene_num <- read.delim(file.path('pivots', 'gene_number.tsv'),
                       header = FALSE,
                       sep = '\t')


ggplot(gene_num,
        aes(x = reorder(V1, V2),
            y = V2)) +
  geom_bar(stat = 'identity',
           fill = "#5C88DAFF") +
  theme_bw()  +
  theme(axis.text.y = element_text(color = 'black', 
                                   size = 13),
        axis.title.y = element_text(face = "bold", color = "black", 
                                    size = 16),
        panel.background = element_rect(fill = "#DCDCDC"),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(size = 0.1, linetype = 'dashed',
                                        colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_text(face = "bold",
                                    color = "black", 
                                    size = 16),
        legend.title = element_text(face = "bold",
                                    color = "black", 
                                    size = 16),
        legend.text = element_text(color = "black", 
                                   size = 14)) +
  scale_y_continuous(breaks = seq(0, 600, by = 100)) +
  ggsci::scale_fill_startrek() +
  labs(x = 'Genomes',
       y = 'No of genes')

## Per genome gene distribution - chromosome genes only
contig_num <- read.delim('gene_content_table.tsv',
                         sep = '\t',
                         header = TRUE,
                         stringsAsFactors = FALSE)
contig_num[is.na(contig_num$ContigType), 'AssemblyType'] <- 'Chromosome'
contig_num[contig_num$ContigType %in% c('chromosome', 'plasmid'), 'AssemblyType'] <- 'Complete'
contig_num[is.na(contig_num)] <- 'Chromosome'


contig_num$Assembly %<>% as.factor()

ggplot(contig_num %>%
         filter(ContigType != 'plasmid'),
       aes(x = reorder(Assembly, NoOfGenes),
           y = NoOfGenes,
           fill = AssemblyType)) +
  geom_bar(stat = 'identity') +
  theme_bw()  +
  theme(axis.text.y = element_text(color = 'black', 
                                   size = 13),
        axis.title.y = element_text(face = "bold", color = "black", 
                                    size = 16),
        panel.background = element_rect(fill = "#DCDCDC"),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(size = 0.1, linetype = 'dashed',
                                        colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_text(face = "bold",
                                    color = "black", 
                                    size = 16),
        legend.title = element_text(face = "bold",
                                    color = "black", 
                                    size = 16),
        legend.text = element_text(color = "black", 
                                   size = 14)) +
  scale_y_continuous(breaks = seq(0, 600, by = 100)) +
  ggsci::scale_fill_startrek() +
  labs(x = 'Chromosome contigs',
       y = 'No of genes',
       fill = 'Assembly type')

## Genome length distribution
genome_lengths <- read.delim(
  file.path('pivots', 'sibeliaz_header.tsv'),
  sep = '\t',
  header = TRUE,
  stringsAsFactors = FALSE
)

ggplot(genome_lengths,
       aes(x = reorder(Description, Size),
           y = Size / 1000)) +
  geom_bar(stat = 'identity',
           fill = "#5C88DAFF") +
  theme_bw()  +
  theme(axis.text.y = element_text(color = 'black', 
                                   size = 13),
        axis.title.y = element_text(face = "bold", color = "black", 
                                    size = 16),
        panel.background = element_rect(fill = "#DCDCDC"),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(size = 0.1, linetype = 'dashed',
                                        colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_text(face = "bold",
                                    color = "black", 
                                    size = 16),
        legend.title = element_text(face = "bold",
                                    color = "black", 
                                    size = 16),
        legend.text = element_text(color = "black", 
                                   size = 14)) +
  ggsci::scale_fill_startrek() +
  labs(x = 'Genomes',
       y = 'Length, kbp')


## Shared block coverage

ggplot(genome_lengths,
       aes(x = reorder(Description, -Size),
           y = 4400 / Size)) +
  geom_bar(stat = 'identity',
           fill = "#5C88DAFF") +
  theme_bw()  +
  theme(axis.text.y = element_text(color = 'black', 
                                   size = 13),
        axis.title.y = element_text(face = "bold", color = "black", 
                                    size = 16),
        panel.background = element_rect(fill = "#DCDCDC"),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(size = 0.1, linetype = 'dashed',
                                        colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_text(face = "bold",
                                    color = "black", 
                                    size = 16),
        legend.title = element_text(face = "bold",
                                    color = "black", 
                                    size = 16),
        legend.text = element_text(color = "black", 
                                   size = 14)) +
  ggsci::scale_fill_startrek() +
  labs(x = 'Genomes',
       y = 'Length, kbp')


#### OrthoFinder
## visualize the gene presence/absence matrix

of_pamatr <- read.delim(file.path('pivots',
                                  'Orthogroups.GeneCount.tsv'),
                        sep = '\t', 
                        header = TRUE,
                        row.names = 1)

of_pamatr %>% nrow()

of_pamatr %<>% as.matrix() %>% 
  t() %>%
  as.data.frame() %>%
  slice(-c(nrow(.))) %>%
  `rownames<-`(sapply(rownames(.), function(x) strsplit(x, '\\.') %>%
                        unlist() %>%
                        dplyr::first()))

of_pamatr %>% pheatmap(cluster_cols = FALSE,
                       show_colnames = FALSE, 
                       color = c('white', 'yellow2', 'orange2', 'orangered'),
                       legend_breaks = 0:3)

of_pamatr %>% View()
colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, 'RdYlBu')))(100)
?pheatmap

#### Roary
## visualize the gene presence/absence matrix

r_pamatr <- read.delim(file.path('pivots', 
                                 'gene_presence_absence.csv'),
                       sep = ',',
                       header = TRUE)

r_pamatr_fplot <- r_pamatr %>% dplyr::select(colnames(.)[grepl('GCF', colnames(.))])
r_pamatr_fplot %<>% sapply(function(x) ifelse(x == '', 0, 1)) %>%
  as.matrix() %>% 
  t() %>% 
  as.data.frame() %>% 
  `rownames<-`(sapply(rownames(.), function(x) x %>% 
                        strsplit('\\.') %>% 
                        unlist() %>% 
                        dplyr::first()))

r_pamatr_fplot %>% 
  pheatmap(cluster_cols = FALSE,
           show_colnames = FALSE,
           color = c('yellow2','orangered'),
           legend_breaks = 0:1)

## Roary pangenome profile
summary_dir <- file.path(getwd(), 'pivots', 'summaries')
files <- summary_dir %>% list.files()

sapply(1:length(files),
       function(x) assign(files[x] %>% 
                            strsplit('\\.') %>% 
                            unlist() %>% 
                            dplyr::first(),
                          read.delim(file.path(summary_dir, files[x]),
                                     header = FALSE,
                                     sep = '\t'),
                          envir = .GlobalEnv))

summary_to_plot <- tibble('I=1.00' = summary_100$V3,
                          'I=0.90' = summary_90$V3,
                          'I=0.85' = summary_85$V3,
                          'I=0.80' = summary_80$V3,
                          'I=0.75' = summary_75$V3,
                          'I=0.70' = summary_70$V3,
                          'I=0.60' = summary_60$V3,
                          'I=0.50' = summary_50$V3,
                          'Pangenome component' = summary_100$V1) %>% 
                   tibble::column_to_rownames('Pangenome component') %>% 
                   as.matrix() %>% 
                   t() %>% 
                   as.data.frame()

summary_to_plot[, 1:4] %<>% sapply(function(x) x / summary_to_plot$`Total genes`)
summary_to_plot %<>% tibble::rownames_to_column('Identity')
summary_to_plot <- data.frame(`Identity` = rep(summary_to_plot[,'Identity'], 4),
            `Pangenome component` = c(rep('Core genes', 8),
                                      rep('Soft core genes', 8),
                                      rep('Shell genes', 8),
                                      rep('Cloud genes', 8)),
            `Number` = c(summary_to_plot[, 'Core genes'],
                  summary_to_plot[, 'Soft core genes'],
                  summary_to_plot[, 'Shell genes'],
                  summary_to_plot[, 'Cloud genes']))

ggplot(summary_to_plot) +
  geom_bar(aes(x = forcats::fct_rev(Identity),
               y = Number, 
               fill = `Pangenome.component`),
           stat = 'identity') +
  theme_bw()  +
  theme(axis.text.y = element_text(color = 'black', 
                                    size = 13),
         axis.title.y = element_blank(),
         panel.background = element_rect(fill = "#DCDCDC"),
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_line(size = 0.1, linetype = 'dashed',
                                         colour = "black"),
         axis.text.x = element_text(color = 'black', 
                                    angle = 60, vjust = 1, 
                                    size = 14, hjust = 1),
         axis.title.x = element_text(face = "bold",
                                     color = "black", 
                                     size = 16),
         legend.title = element_text(face = "bold",
                                     color = "black", 
                                     size = 16),
         legend.text = element_text(color = "black", 
                                    size = 14)) +
  labs(fill = 'Pangenome component',
       x = 'Identity threshold') +
  ggsci::scale_fill_startrek()

