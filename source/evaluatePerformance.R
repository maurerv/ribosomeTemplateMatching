library(data.table)
library(ggplot2)
library(clue)
library(parallel)
library(colorspace)
library(cowplot)
library(dplyr)

setwd("/Users/vmaurer/src/ribosomeSpheres")

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

# See process_coordinates.R to generate data loaded by this function
load_data = function(path){
  estimate_paths = list.files(
    path = path,
    recursive = T, 
    full.names = T, 
    pattern = "coordinates_annotated.tsv"
  )
  estimates = rbindlist(lapply(estimate_paths, function(x){
    temp = fread(x)
    if(!nrow(temp)){
      return(data.table())
    }
    return(temp)
  }))
  estimates[, name := as.numeric(name)]
  estimates[, Class := factor(group, levels = levels, labels = labels)]
  # estimates = estimates[!is.na(Class)]
  
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
  return(estimates)
}

ground_truth = fread("particle_lists/TS_037_cyto_ribosomes.csv")
fas_particles = fread("particle_lists/TS_037_fas.csv")

colnames(ground_truth) = c("x", "y", "z")
colnames(fas_particles) = c("x", "y", "z")
total_particles = rbind(ground_truth, fas_particles)
colnames(total_particles) = c("x", "y", "z")
classes = c(rep("Ribosome", nrow(ground_truth)), rep("FAS", nrow(fas_particles)))

# Data used for generating main test figures
estimates = load_data("templates_inverted_2k")

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
  scale_y_continuous(labels = scales::percent, limits = c(0, .5))+
  facet_grid(~Class)+
  xlab("Number of picked particles")+
  ylab("Recall")+
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "cm"))

p1_precision = ggplot(estimates, 
                      aes(x = false_positives + true_positives,
                          y = true_positives / (false_positives + true_positives),
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
  ylab("Precision")+
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "cm"))

p1_precision_recall = ggplot(estimates, 
                      aes(x = true_positives / nrow(ground_truth),
                          y = true_positives / (false_positives + true_positives),
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
  xlab("Recall")+
  ylab("Precision")+
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
  geom_point(size = 3)+
  geom_line(linewidth = 1.5)+
  ylab("true postives")+
  scale_color_manual(name = "Template",
                     values = c(`Ribosome (EMD-3228)` = "#A7D9CB", `Sphere` = "#2188c9", `Emoji` = "#fca349"))+
  theme_bw(base_size = FONT_SIZE)+
  xlab("Template radius [Voxel]")+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(breaks = c(min(plot_data$name), 5, 10, 15, max(plot_data$name)))+
  ylab("Precision")+
  theme(legend.position = "bottom")
ggsave("/Users/vmaurer/Desktop/comparison_template.pdf", p2, width = 8, height = 6)

p2 = ggplot(plot_data, aes(x = name, y = precision, group = Class, color = Class))+
  geom_point(size = 3)+
  geom_line(linewidth = 1.5)+
  ylab("true postives")+
  scale_color_manual(name = "Template", values = colorMapping)+
  theme_bw(base_size = FONT_SIZE)+
  xlab("Template radius [Voxel]")+
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05 * max(plot_data$precision)))+
  scale_x_continuous(breaks = c(min(plot_data$name), 5, 10, 15, max(plot_data$name)))+
  ylab("Precision")+
  theme(legend.position = "bottom")

p3 = ggplot(plot_data[tp > 0, mean(tp), by = .(Class, name)],
            aes(x = name, y = V1, group = Class, color = Class))+
  geom_point(size = 3)+
  geom_line(linewidth = 1.5)+
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
  scale_y_continuous(labels = scales::percent, limits = c(0, .5))+
  scale_x_continuous(labels = scales::percent, limits = c(0, .5))+
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


radialAverages = fread("source/radialAverages.tsv")
bins = sort(unique(radialAverages$bin))
vals = seq(min(bins), max(bins), .1)

sphereTheory = data.table(
  class = "Sphere [Theory]", bin = vals, value = abs(besselJ(vals, nu = 1)/vals)
)

a = 10
vals = seq(0, 0.5, .001) * 2 * pi * a
sphereTheory = data.table(
  class = "Sphere [Theory]", bin = vals,
  value = 4 * pi * a ** 3 * abs(besselJ(vals, nu = 1)/vals )
)
max_value = max(sphereTheory$value, na.rm = T)
sphereTheory[is.na(value), value := max_value]

radialAverages = rbindlist(list(radialAverages, sphereTheory))
radialAverages[, value := (value - min(value))/(max(value) - min(value)), by = class]

radialAverages[, Class := factor(class, levels, labels)]
radialAverages[, Linetype := factor(class, levels, labels)]
radialAverages[grepl(pattern = "theory", ignore.case = T, class), Linewidth := 1.5]
radialAverages[!grepl(pattern = "theory", ignore.case = T, class), Linewidth := 1]
plot_data = radialAverages[bin < 10]
p6 = ggplot(plot_data, aes(x = bin, y = value, color = Class))+
  geom_line(data = plot_data[!grepl(pattern = "theory", ignore.case = T, class)],
            linewidth = 1)+
  geom_line(data = plot_data[grepl(pattern = "theory", ignore.case = T, class)],
            linetype = "twodash", linewidth = 2)+
  geom_point(data = plot_data[!grepl(pattern = "theory", ignore.case = T, class)], size = 3)+
  scale_color_manual(name = element_blank(), values = colorMapping)+
  ylab("Average Fourier Magnitude")+
  xlab("Distance [Voxel]")+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(breaks = seq(0, 10, 2))+
  theme_bw(base_size = FONT_SIZE - 2)+
  theme(legend.position = "bottom")+
  guides(linetype = FALSE)
p6
ggsave("plots/radialAverages.pdf", p6, width = 8, height = 6)


lgd = cowplot::get_legend(p6)
middle = cowplot::plot_grid(plotlist = list(
  p2 + theme(legend.position = "None"),
  p6+ theme(legend.position = "None")), labels = c("B", "C"), label_size = 20, ncol = 2
)
middleTotal = cowplot::plot_grid(plotlist = list(middle, lgd), nrow = 2, rel_heights = c(.9, .1))
middle = cowplot::plot_grid(plotlist = list(p2, p6), labels = c("B", "C"), label_size = 20)
total = cowplot::plot_grid(plotlist = list(upper, middleTotal, lower), nrow = 3, labels = NA)
cowplot::ggsave2("plots/templateMatchingRibosomeTotal.pdf",
       total, width = 16, height = 18)

cowplot::ggsave2("plots/precisionByTemplate.pdf",
                 p1_precision, width = 16, height = 6)

ggsave("plots/fas_precisionrecall.pdf", p4, width = 16, height = 6)

lgd = cowplot::get_legend(p1)
total = cowplot::plot_grid(plotlist = list(
  p1 + theme(legend.position = "None"),
  p1_precision + theme(legend.position = "None"),
  lgd
  ), nrow = 3, labels = c("A", "B", ""), rel_heights = c(.9, .9, .1),
  label_size = 20
)
total
cowplot::ggsave2("plots/templateMatchingRibosomeTotal.pdf",
       total, width = 16, height = 12
)


angles_980 = load_data("templates_inverted_980")
angles_2k = load_data("templates_inverted_2k")
angles_15k = load_data("templates_inverted_15k")
angles_980[, angles := 980]
angles_2k[, angles := 1944]
angles_15k[, angles := 15192]
pytom_data = rbindlist(list(angles_980, angles_2k, angles_15k))

angles_980 = load_data("templates_pytme_980")
angles_2k = load_data("templates_pytme_2k")
angles_15k = load_data("templates_pytme_15k")
angles_980[, angles := 980]
angles_2k[, angles := 1944]
angles_15k[, angles := 15192]
pytme_data = rbindlist(list(angles_980, angles_2k, angles_15k))

pytom_data[, tool := "PyTom"]
pytme_data[, tool := "pyTME"]
plot_data = rbindlist(list(pytom_data, pytme_data))

# Curve shape is consistent across number of picks
# plot_data = plot_data[true_positives + false_positives < 2000]
# plot_data = plot_data[true_positives + false_positives < 750]

plot_data = plot_data[!is.na(Class)]
plot_data[is.na(distance), distance := DISTANCE_CUTOFF * 10]
plot_data[is.na(distanceBoth), distanceBoth := DISTANCE_CUTOFF * 10]
plot_data[is.na(assignmentClass), assignmentClass := ""]
plot_data = plot_data[
  , .(tp = sum(distance < DISTANCE_CUTOFF),
      fp = sum(distance >= DISTANCE_CUTOFF)
  ), by = .(Class, name, angles, tool)
]
plot_data[, recall := tp / nrow(ground_truth), by = .(Class, name)]
plot_data[, precision := tp / (tp + fp), by = .(Class, name)]
plot_data[, tool := factor(tool, levels = c("PyTom", "pyTME"), labels = c("PyTom", "pyTME"))]
si_p1 = ggplot(plot_data, aes(x = name, y = precision, group = Class, color = Class))+
  geom_point(size = 3)+
  geom_line(linewidth = 1.5)+
  scale_color_manual(name = "Template", values = colorMapping)+
  theme_bw(base_size = FONT_SIZE)+
  xlab("Template radius [Voxel]")+
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.1 * max(plot_data$precision)))+
  scale_x_continuous(breaks = c(min(plot_data$name), 5, 10, 15, max(plot_data$name)))+
  ylab("Precision")+
  facet_grid(angles~tool)+
  theme(legend.position = "bottom")
cowplot::ggsave2("plots/templateMatchingPrecisionComparison.pdf",
                 si_p1, width = 14, height = 14)

temp = ggplot(plot_data[tool == "PyTom"], aes(x = name, y = precision, group = Class, color = Class))+
  geom_point(size = 3)+
  geom_line(linewidth = 1.5)+
  scale_color_manual(name = "Template", values = colorMapping)+
  theme_bw(base_size = FONT_SIZE)+
  xlab("Template radius [Voxel]")+
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.1 * max(plot_data$precision)))+
  scale_x_continuous(breaks = c(min(plot_data$name), 5, 10, 15, max(plot_data$name)))+
  ylab("Precision")+
  facet_grid(~angles)+
  theme(legend.position = "bottom")
cowplot::ggsave2("plots/templateMatchingPrecision_pyTOM.pdf",
                 temp, width = 14, height = 7)


angles_temp = load_data("templates_pytme_2k")
angles_temp[is.na(distance), distance := DISTANCE_CUTOFF * 10]
angles_temp[is.na(distanceBoth), distanceBoth := DISTANCE_CUTOFF * 10]
angles_temp[is.na(assignmentClass), assignmentClass := ""]
angles_temp = angles_temp[
  , .(tp = sum(distance < DISTANCE_CUTOFF),
      fp = sum(distance >= DISTANCE_CUTOFF)
  ), by = .(group, name)
]
angles_temp[, recall := tp / nrow(ground_truth), by = .(group, name)]
angles_temp[, precision := tp / (tp + fp), by = .(group, name)]
p7 = ggplot(angles_temp[group == "ha"], aes(x = name, y = precision, group = group))+
  geom_point(size = 3)+
  geom_line(linewidth = 1.5)+
  theme_bw(base_size = FONT_SIZE - 2)+
  xlab("Template radius [Voxel]")+
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.1 * max(plot_data$precision)))+
  scale_x_continuous(breaks = c(min(plot_data$name), 5, 10, 15, max(plot_data$name)))+
  ylab("Precision")+
  theme(legend.position = "bottom")
ggsave("plots/ha_precision.pdf", p7, width = 8, height = 6)

