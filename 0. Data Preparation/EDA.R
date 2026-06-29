#===============================================================================
# TESTING

library(plotly)
library(stats)
library(ggplot2)

cell_by_count <- TF_array_filtered %>% 
  group_by(cell_ontology_class) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count))

main_type <- TF_array_filtered %>% 
  filter(cell_ontology_class == "neutrophil") %>% 
  select(all_of(tf_names)) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(across(everything(),
                if (normalisation_method == 1) {
                  ~ . / sum(., na.rm = TRUE)} else {
                    ~ . / sqrt(sum(. ^ 2, na.rm = TRUE))})) %>% 
  t()

prin_comp <- prcomp(main_type, rank. = 3)

components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3

tot_explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
tot_explained_variance_ratio <- 100 * sum(tot_explained_variance_ratio)

tit = paste("Total Explained Variance =", tot_explained_variance_ratio)

other_type <- TF_array_filtered %>% 
  filter(cell_ontology_class == "neuron") %>% 
  select(all_of(tf_names)) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(across(everything(),
                if (normalisation_method == 1) {
                  ~ . / sum(., na.rm = TRUE)} else {
                    ~ . / sqrt(sum(. ^ 2, na.rm = TRUE))})) %>% 
  t()

main_type <- as.data.frame(as.matrix(main_type) %*% prin_comp$rotation)
main_type$group <- "Main"

other_type <- as.data.frame(as.matrix(other_type) %*% prin_comp$rotation)
other_type$group <- "Other"

df <- bind_rows(main_type, other_type)

fig <- plot_ly(
  df,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~group,
  colors = c("Main" = "red", "Other" = "blue")) %>%
  add_markers(size = 12)

fig <- fig %>%
  layout(
    title = tit,
    scene = list(bgcolor = "#e5ecf6")
  )

fig


#===============================================================================
