library(cluster)
library(dplyr)
library(hash)
library(magrittr)
library(Quartet)

roary_matrix <- read.delim(file.path('trees',
                                     'gene_presence_absence_local_backup.csv'),
                           header = TRUE,
                           sep = ',',
                           stringsAsFactors = FALSE)

roary_matrix %>% View()

for_clustering <- roary_matrix %>% dplyr::select(c(1, grep('GCF', colnames(roary_matrix)))) %>% tibble::column_to_rownames(var = 'Gene')
for_clustering[for_clustering != ''] <- 1
for_clustering[for_clustering == ''] <- 0
for_clustering %<>% sapply(as.numeric)

jaccard <- function(M, user1, user2) {
  sums = rowSums(M[,c(user1, user2)])
  
  similarity = length(sums[sums==2])
  total = length(sums[sums==1]) + similarity
  
  return(similarity/total)
}

jaccard_matrix <- sapply(1:ncol(for_clustering),
                         function(x) sapply(1:ncol(for_clustering),
                                            function(y) jaccard(for_clustering, x, y)))

jaccard_matrix %<>% `colnames<-`(colnames(for_clustering)) %>% `rownames<-`(colnames(for_clustering))

wpgma_roary <- phangorn::wpgma(jaccard_matrix)
ape::write.tree(wpgma_roary, file.path('trees', 'wpgma_roary.newick'))


orthof_matrix <- read.delim(file.path(
  'trees', 'Orthogroups.GeneCount.tsv'
  ),
  header = TRUE,
  sep = '\t',
  stringsAsFactors = FALSE
  )

orthof_matrix %<>% tibble::column_to_rownames(
  var = 'Orthogroup'
  ) %>% 
  dplyr::select(!Total)

orthof_matrix %<>% sapply(
  function(x) 
    ifelse(
      x > 0, 1, 0
      )
  )

jaccard_matrix <- sapply(1:ncol(orthof_matrix),
                         function(x) sapply(1:ncol(orthof_matrix),
                                            function(y) jaccard(orthof_matrix, x, y)))

jaccard_matrix %<>% `colnames<-`(colnames(orthof_matrix)) %>% `rownames<-`(colnames(orthof_matrix))
wpgma_orthof <- phangorn::wpgma(jaccard_matrix)

ape::write.tree(
  wpgma_orthof,
  file.path(
  'trees', 'wpgma_orthofinder.newick'
),
)

tree_files <- list.files('trees') %>% .[!grepl('csv', .) & !grepl('tsv', .)]

trees <- sapply(tree_files, function(x)
                ape::read.tree(file.path('trees', x)))

tree_comparison <- function(tree1, tree2) {
  statuses <- QuartetStatus(tree1, tree2)
  return(SimilarityMetrics(statuses, similarity = TRUE))
    
}



all_comparisons <- sapply(
  combn(5, 2) %>%
    as.data.frame(),
  function(x) tree_comparison(trees[[x[1]]], trees[[x[2]]])
  ) %>% 
  as.matrix() %>% 
  t() %>% 
  as.data.frame()

buildNames <- function(x) {
  tree_names <- hash(seq_len(5), c('Roary_binary',
                                   'OrthoFinder_RAxML',
                                   'Roary_RAxML',
                                   'Roary_WPGMA',
                                   'OrthoFinder_WPGMA'))
  linenames <- sapply(
      as.data.frame(),
    combn(5,2) %>% 
    function(x) paste(
      values(tree_names,
             key = x[1]),
      'Vs',
      values(tree_names,
             key = x[2]),
      sep = ''
      ) 
      )
  return(linenames)
}


all_comparisons %<>% `colnames<-`(colnames(tree_comparison(tree1, tree2))) %>% `rownames<-`(buildNames())

write.table(all_comparisons, file.path('trees', 'quartet_report.tsv'), sep = '\t', col.names = TRUE, row.names = FALSE)


