#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(clue)
library(colorspace)
library(cowplot)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {
  stop("No directory path provided", call. = FALSE)
}

# Set up paths
coordinates_path <- file.path(args[1], "coordinates.tsv")
output_path <- file.path(args[1], "coordinates_annotated.tsv")

# Constants
FONT_SIZE = 20
DISTANCE_CUTOFF = 5
levels = c("emd3228", "sphere", "emoji", "Sphere [Theory]", "Rectangle [Theory]")
labels = c("Ribosome (EMD-3228)", "Sphere", "Emoji", "Sphere [Theory]", "Rectangle [Theory]")

colorMapping = c(
  `Ribosome (EMD-3228)` = "#E41A1C",
  Sphere = "#377EB8",
  Emoji = "#4DAF4A",
  `Sphere [Theory]` = "#984EA3"
)

ground_truth = fread("../particle_lists/TS_037_cyto_ribosomes.csv")
fas_particles = fread("../particle_lists/TS_037_fas.csv")
colnames(ground_truth) = c("x", "y", "z")
colnames(fas_particles) = c("x", "y", "z")
total_particles = rbind(ground_truth, fas_particles)
colnames(total_particles) = c("x", "y", "z")
classes = c(rep("Ribosome", nrow(ground_truth)), rep("FAS", nrow(fas_particles)))

estimate = fread(coordinates_path)

if(nrow(estimate) == 0){
  fwrite(data.table(), output_path)
}

if(grepl(pattern = "pytme", coordinates_path)){
  colnames(estimate) = c("z", "y", "x", "score", "name")
}else{
  colnames(estimate) = c("x", "y", "z", "score", "name")
}

name = basename(unique(estimate$name))
score = estimate$score
estimate = estimate[, c("x", "y", "z")]

print(estimate)
print(ground_truth)
distances = as.matrix(dist(rbind(estimate, ground_truth)))
distancesTotal = as.matrix(dist(rbind(estimate, total_particles)))

# Ribosome assignment
n1 = nrow(estimate)
n2 = nrow(ground_truth)
distances = distances[(n1+1):(n1+n2), 1:n1]
assignment = solve_LSAP(distances, maximum = FALSE)
nearest_neighbors = as.vector(assignment)
nearest_neighbor_distances = sapply(seq_along(nearest_neighbors), function(i) {
  distances[i, nearest_neighbors[i]]
})
estimate[nearest_neighbors, assignment := 1:nrow(ground_truth)]
estimate[nearest_neighbors, distance := nearest_neighbor_distances]
nearest_neighbors = apply(distances, 2, which.min)
nearest_neighbor_distances = sapply(seq_along(nearest_neighbors), function(i) {
  distances[nearest_neighbors[i], i]
})
estimate[, assignmentNotUnique := nearest_neighbors]
estimate[, distanceNotUnique := nearest_neighbor_distances]

# Ribosome + FAS assignment
n2 = nrow(total_particles)
distancesTotal = distancesTotal[(n1+1):(n1+n2), 1:n1]
assignment = clue::solve_LSAP(distancesTotal, maximum = FALSE)
nearest_neighbors = as.vector(assignment)
nearest_neighbor_distances = sapply(seq_along(nearest_neighbors), function(i) {
  distancesTotal[i, nearest_neighbors[i]]
})
estimate[nearest_neighbors, assignmentClass := classes]
estimate[nearest_neighbors, distanceBoth := nearest_neighbor_distances]

estimate[, score := score]
estimate$group = gsub(name, pattern = "(.*)_(\\d+)", replacement = "\\1")
estimate$name = gsub(name, pattern = "(.*)_(\\d+)", replacement = "\\2")
estimate

fwrite(estimate, output_path)
