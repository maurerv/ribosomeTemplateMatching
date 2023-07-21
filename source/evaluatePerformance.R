library(data.table)
library(ggplot2)
library(clue)
library(parallel)
library(colorspace)
library(cowplot)

setwd("/Users/vmaurer/src/ribosomeSpheres")

FONT_SIZE = 20
DISTANCE_CUTOFF = 20

ground_truth = fread("particle_lists/TS_037_cyto_ribosomes.csv")
fas_particles = fread("particle_lists/TS_037_fas.csv")

colnames(ground_truth) = c("x", "y", "z")
colnames(fas_particles) = c("x", "y", "z")
total_particles = rbind(ground_truth, fas_particles)
colnames(total_particles) = c("x", "y", "z")
classes = c(rep("Ribosome", nrow(ground_truth)), rep("FAS", nrow(fas_particles)))

estimate_paths = list.files("templates",recursive = T, full.names = T, pattern = "coordinates.tsv")

estimates = mclapply(estimate_paths, function(estimate_path){
  estimate = fread(estimate_path)
  if(nrow(estimate) == 0){
    return(data.table())
  }
  colnames(estimate) = c("x", "y", "z", "score", "name")
  name = basename(unique(estimate$name))
  score = estimate$score
  estimate = estimate[, c("x", "y", "z")]
  
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
}, mc.cores = 8)
estimates = rbindlist(estimates)
# fwrite(estimates, "/Users/vmaurer/src/ribosomeSpheres/mappings.csv")

estimates = fread("/Users/vmaurer/src/ribosomeSpheres/mappings.csv")

estimates[, name := as.numeric(name)]
estimates[, Class := factor(group, 
                            levels = c("7p6z", "sphere", "emoji"), 
                            labels = c("Ribosome (7p6z)", "Sphere", "Emoji"))
]

estimates[is.na(distance), distance := DISTANCE_CUTOFF * 10]
estimates[is.na(distanceBoth), distanceBoth := DISTANCE_CUTOFF * 10]
estimates[is.na(assignmentClass), assignmentClass := ""]

estimates[, true_positives := cumsum(distance < DISTANCE_CUTOFF), by = .(group, name)]
estimates[, false_positives := cumsum(distance >= DISTANCE_CUTOFF), by = .(group, name)]

estimates[, true_positivesFAS := cumsum(distanceBoth < DISTANCE_CUTOFF & assignmentClass == "FAS"), 
          by = .(group, name)]
estimates[, true_positivesRibosome := cumsum(distanceBoth < DISTANCE_CUTOFF & assignmentClass == "Ribosome"), 
          by = .(group, name)]
estimates[, false_positivesTotal := cumsum(distanceBoth >= DISTANCE_CUTOFF), 
          by = .(group, name)]

p1 = ggplot(estimates, aes(x = false_positives + true_positives, 
                      y = true_positives / nrow(ground_truth), 
                      color = name, group = name))+
  geom_line(linewidth = 1.5)+
  theme_bw(base_size = FONT_SIZE)+
  scale_color_continuous_divergingx(name = "Radius [Voxel]", 
                                    palette = "RdYlBu", 
                                    mid = 10,
                                    limits = c(min(estimates$name), max(estimates$name)),
                                    breaks = c(min(estimates$name), 5, 10, 15, max(estimates$name))
  )+
  scale_y_continuous(labels = scales::percent)+
  facet_grid(~Class)+
  xlab("Number of picked particles")+
  ylab("Recall")+
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "cm"))


plot_data = estimates[
  , .(tp = sum(distance < DISTANCE_CUTOFF), 
      fp = sum(distance >= DISTANCE_CUTOFF)
  ), by = .(Class, name)
]
plot_data[, recall := tp / nrow(ground_truth), by = .(Class, name)]
plot_data[, precision := tp / (tp + fp), by = .(Class, name)]
p2 = ggplot(plot_data, aes(x = name, y = precision, group = Class, color = Class))+
  geom_point()+
  geom_line()+
  ylab("true postives")+
  scale_color_brewer(name = "Template", palette = "Set1")+
  theme_bw(base_size = FONT_SIZE)+
  xlab("Template radius [Voxel]")+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(breaks = c(min(plot_data$name), 5, 10, 15, max(plot_data$name)))+
  ylab("Precision")+
  theme(legend.position = "bottom")

p3 = ggplot(plot_data[tp > 0, mean(tp), by = .(Class, name)], 
            aes(x = name, y = V1, group = Class, color = Class))+
  geom_point()+
  geom_line()+
  ylab("true postives")+
  scale_color_brewer(name = "Template", palette = "Set1")+
  theme_bw(base_size = FONT_SIZE)+
  xlab("Template radius")+
  ylab("Average number of picks per true particle")+
  theme(legend.position = "bottom")
p3

p4 = ggplot(estimates, aes(x = true_positivesRibosome / nrow(ground_truth), 
                           y = true_positivesFAS / nrow(fas_particles), 
                           color = name, group = name))+
  geom_line(linewidth = 1.5)+
  theme_bw(base_size = FONT_SIZE)+
  scale_color_continuous_divergingx(name = "Radius [Voxel]", 
                                    palette = "RdYlBu", 
                                    mid = 10,
                                    limits = c(min(estimates$name), max(estimates$name)),
                                    breaks = c(min(estimates$name), 5, 10, 15, max(estimates$name))
  )+
  scale_y_continuous(labels = scales::percent, limits = c(0, .6))+
  scale_x_continuous(labels = scales::percent)+
  facet_grid(~Class)+
  xlab("Sampled Recall FAS")+
  ylab("Sampled Recall Ribosome")+
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "cm"))
p4

upper = cowplot::plot_grid(plotlist = list(p1), labels = "A", label_size = 20)
middle = cowplot::plot_grid(plotlist = list(p2, p3), labels = c("B", "C"), label_size = 20)
lower = cowplot::plot_grid(plotlist = list(p4), labels = c("D"), label_size = 20)
total = cowplot::plot_grid(plotlist = list(upper, middle), nrow = 2, labels = NA)

ggsave("/Users/vmaurer/src/ribosomeSpheres/plots/templateMatchingRibosome.pdf",
       total, width = 16, height = 12)
total = cowplot::plot_grid(plotlist = list(upper, middle, lower), nrow = 3, labels = NA)
ggsave("/Users/vmaurer/src/ribosomeSpheres/plots/templateMatchingRibosomeTotalPicking.pdf",
       total, width = 16, height = 18)


radialAverages = fread("/Users/vmaurer/Desktop/radialAverages.tsv")
bins = sort(unique(radialAverages$bin))
vals = seq(min(bins), max(bins), .1)
sphereTheory = data.table(
  class = "Sphere [Theory]", bin = vals, value = abs(besselJ(vals, nu = 1)/vals)
)
sphereTheory[is.na(value), value := 1]
# rectangleTheory = data.table(
#   class = "Rectangle [Theory]", bin = vals, value = abs(besselJ(vals, nu = 0)**3)
# )
# rectangleTheory[is.na(value), value := 1]

radialAverages = rbindlist(list(radialAverages, sphereTheory))
radialAverages[, value := (value - min(value))/(max(value) - min(value)), by = class]

radialAverages[, Class := factor(class, 
                                 levels = c("ribosome", "sphere", "emoji", "Sphere [Theory]", "Rectangle [Theory]"), 
                                 labels = c("Ribosome (7p6z)", "Sphere", "Emoji", "Sphere [Theory]", "Rectangle [Theory]"))
]
radialAverages[, Linetype := factor(class, 
                                 levels = c("ribosome", "sphere", "emoji", "Sphere [Theory]", "Rectangle [Theory]"), 
                                 labels = c("Practical", "Practical", "Practical", "Theory", "Theory"))
]
p6 = ggplot(radialAverages[bin < 10], aes(x = bin, y = value, color = Class))+
  geom_line(mapping = aes(linetype = Linetype))+
  # geom_line()+
  geom_point()+
  scale_color_brewer(name = "Template", palette = "Set1")+
  ylab("Power spectrum average")+
  xlab("Radius [Voxel]")+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(breaks = seq(0, 10, 2))+
  theme_bw(base_size = FONT_SIZE)+
  theme(legend.position = "bottom")+
  guides(linetype = FALSE)
p6
ggsave("/Users/vmaurer/src/ribosomeSpheres/plots/radialAverages.pdf",
       p6, width = 8, height = 6)

lgd = cowplot::get_legend(p6)
middle = cowplot::plot_grid(plotlist = list(
  p2 + theme(legend.position = "None"), 
  p6+ theme(legend.position = "None")), labels = c("B", "C"), label_size = 20, ncol = 2
)
middleTotal = cowplot::plot_grid(plotlist = list(middle, lgd), nrow = 2, rel_heights = c(.9, .1))
middle = cowplot::plot_grid(plotlist = list(p2, p6), labels = c("B", "C"), label_size = 20)
total = cowplot::plot_grid(plotlist = list(upper, middleTotal, lower), nrow = 3, labels = NA)
cowplot::ggsave2("/Users/vmaurer/src/ribosomeSpheres/plots/templateMatchingRibosomeTotal.pdf",
       total, width = 16, height = 18)
