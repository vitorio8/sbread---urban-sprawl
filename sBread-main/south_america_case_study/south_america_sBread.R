library(tidyverse)
library(sf)
library(spdep)
library(igraph)
library(gridExtra)
library(kableExtra)
library(here)

here::i_am("south_america_case_study/south_america_sBread.R")

import::here(barabasi.game, .from = igraph)
import::here(sourceCpp, .from = Rcpp)
import::here(mclapply, .from = parallel)

setwd("..")
import::here(sBread.from.igraph, .from = "Rcpp_code/sBread_from_igraph.R")
import::here(similarity.sequence, .from = "Rcpp_code/convergence_diagnostics.R")
import::here(discrete.ess, .from = "Rcpp_code/discrete_ess.R")
sourceCpp(here("Rcpp_code", "main.cpp"))
sourceCpp(here("Rcpp_code", "analyze_posterior.cpp"))
setwd("south_america_case_study")

# cell size in metres
cellsize <- 200000
# Time interval between snapshots
delta_t <- 10

#Baselines
baseline.cellsize <- 200000
baseline.delta_t <- 10
baseline.transmission_rate <- 0.1

#rate of transmission between neighbors
transmission_rate <- baseline.transmission_rate * 
  (delta_t / baseline.delta_t) / ((cellsize / baseline.cellsize) ^ 2)


seed <- 4
set.seed(seed)
K <- 2
mut.rate <- 0.001 * delta_t
num.chains <- 1

length.burnin <- 3 * 10 ^ 7
sample.size.by.chain <- 10 ^ 3
sample.interval <- 3 * 10 ^ 4



data_path <- here("south_america_case_study", "data")


# Load polygon, crop and reproject to an equal area CRS
crs_eq_area <- st_crs("+proj=tcea +lon_0=-57.3046875 +datum=WGS84 +units=m +no_defs")

america <- st_read(here(data_path, "ne_10m_land", "ne_10m_land.shp")) %>%
  st_crop(xmin = -100 , ymin = -60, xmax = -20, ymax = 15) %>%
  st_cast("POLYGON") %>%
  filter(!is.null(geometry)) %>%
  mutate(area = st_area(geometry) %>% as.vector()) %>% 
  arrange(desc(area)) %>%
  slice(1:2) %>%
  st_difference(st_bbox(c(xmin = -100, ymin = 6, 
                          xmax = -81, ymax = 20),
                        crs = 4326) %>% st_as_sfc()) %>%
  st_transform(crs_eq_area)

ggplot(america) +
  geom_sf() +
  theme_minimal()


# Generate the hexagonal grid
bbox <- america %>% st_bbox() %>% st_as_sfc()
grid <- st_make_grid(bbox, what = "polygons", 
                     cellsize = cellsize, square = FALSE, flat_topped = FALSE) %>%
  st_as_sf()

st_geometry(grid) = "geom"

# Filter grid cells intersecting with South America 
america_grid <- grid %>%
  st_filter(america, .predicates = st_intersects) %>%
  mutate(id = row_number())

current_ie <- st_read(here(data_path, "language_polygons_current.json")) %>%
  filter(family == "Indo-European") %>%
  st_transform(crs_eq_area)


historical_ie <- read_csv(here(data_path, "historical_seeds_ie.csv")) %>%
  st_as_sf(coords = c("lon", "lat"), crs = st_crs(4326)) %>%
  st_transform(crs_eq_area)

# Arrival of first languages
t_initial <- historical_ie %>% 
  filter(year == min(year)) %>%
  pull(year)

# Current snapshot
t_final <- 1990

t_snapshots <- seq(t_initial, t_final, by = delta_t)

# For each new arrival find the closest year in t_snapshots
historical_ie <- historical_ie %>% 
  rowwise() %>%
  mutate(match_year = min(t_snapshots[t_snapshots >= year]))

history <- sapply(t_snapshots, function (t) {
  if (t == t_initial) {
    ifelse(st_intersects(america_grid, historical_ie %>%
                           filter(match_year == t), sparse = F) %>%
             rowSums() > 0L, 1, 0)
  } else if (t == t_final) {
    ifelse(st_intersects(america_grid, current_ie, sparse = F) %>% 
             rowSums() > 0L, 1, 0)
  } else {
    ifelse(st_intersects(america_grid, historical_ie %>%
                           filter(match_year == t), sparse = F) %>% 
             rowSums() > 0L, 1, -1)}}, simplify = T) %>% t 

# Join the first and last row of the matrix with the grid 
america_grid <- america_grid %>% left_join(
  data.frame(initial = history[1, ],
             current = history[nrow(history), ]) %>%
    mutate(id = row_number()), by = "id")

america_initial_map <- ggplot(america_grid) +
  geom_sf(mapping = aes(fill = factor(initial))) +
  ggtitle("1509") +
  labs(fill = NULL) +
  scale_fill_discrete(labels = c("indigenous", "Indo-European")) +
  theme_minimal() 

america_current_map <- ggplot(america_grid) +
  geom_sf(mapping = aes(fill = factor(current))) +
  ggtitle("1990") +
  labs(fill = NULL) +
  scale_fill_discrete(labels = c("indigenous", "Indo-European")) +
  theme_minimal() 

grid.arrange(america_initial_map, america_current_map, ncol = 1, nrow = 2)

adj_list <- poly2nb(america_grid) %>% 
  nb2mat(style = "B") %>%
  graph_from_adjacency_matrix(mode = "undirected") %>%
  as_adj_list()


###modify the adj_list for sBread
num.node <- length(adj_list)

adj_list_modified <- list()
for(i in 1:num.node){
  adj_list_modified[[i]] <- c(as.vector(adj_list[[i]] - 1), i - 1)
}

weight_list <- list()
for(i in 1:num.node){
  num.neighbors.include.self <- length(adj_list_modified[[i]])
  weight_list[[i]] <- rep(transmission_rate, num.neighbors.include.self)
  weight_list[[i]][num.neighbors.include.self] <- 1 - transmission_rate * (num.neighbors.include.self - 1)
}

T <- nrow(history)
N <- ncol(history)

#####################sBread algorithm

sample.histories <- sBread(seed, K, mut.rate, num.chains, length.burnin, sample.size.by.chain,
                           sample.interval, adj_list_modified, weight_list, history, TRUE)

posterior.ie.language <- compute_posterior(K, T, N, sample.histories)


#################for plot

america_grid <- america_grid %>% left_join(
  data.frame(mid_state1 = posterior.ie.language[[2]][T/4,]) %>%
    mutate(id = row_number()), by = "id")

america_mid_map1 <- ggplot(america_grid) +
  geom_sf(mapping = aes(fill = mid_state1)) +
  scale_fill_gradientn(colours = c("#F8766D", "white", "#00BFC4"),
                       values = c(0, 0.5, 1),
                       limits=c(0, 1.001)) +
  ggtitle("mid1") +
  labs(fill = NULL) +
  theme_minimal()

america_grid <- america_grid %>% left_join(
  data.frame(mid_state2 = posterior.ie.language[[2]][T/2,]) %>%
    mutate(id = row_number()), by = "id")

america_mid_map2 <- ggplot(america_grid) +
  geom_sf(mapping = aes(fill = mid_state2)) +
  scale_fill_gradientn(colours = c("#F8766D", "white", "#00BFC4"),
                       values = c(0, 0.5, 1),
                       limits=c(0, 1.001)) +
  ggtitle("mid2") +
  labs(fill = NULL) +
  theme_minimal() 

america_grid <- america_grid %>% left_join(
  data.frame(mid_state3 = posterior.ie.language[[2]][3*T/4,]) %>%
    mutate(id = row_number()), by = "id")

america_mid_map3 <- ggplot(america_grid) +
  geom_sf(mapping = aes(fill = mid_state3)) +
  scale_fill_gradientn(colours = c("#F8766D", "white", "#00BFC4"),
                       values = c(0, 0.5, 1),
                       limits=c(0, 1.001)) +
  ggtitle("mid3") +
  labs(fill = NULL) +
  theme_minimal()

map_grid <- grid.arrange(america_initial_map, america_mid_map1, america_mid_map2,
                         america_mid_map3, america_current_map, ncol = 5, nrow = 1)
map_grid
# ggsave(file="map_grid.pdf", plot=map_grid, width=20, height=8, dpi=4000)
