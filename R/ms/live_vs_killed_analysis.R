# A little setup
point.size <- 0.5
line.size <- 0.25
old_theme <- theme_get() 
theme_set(old_theme +
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank()))
  

# Create v0 column with minimum values
min.v0s <- samp_slopes %>%  
  group_by(treatment) %>%
  summarise(min.v0s = min(v0[v0 - v0.se > 0], na.rm=TRUE))
min.v0.killed <- min.v0s$min.v0s[min.v0s$treatment == "killed"]
min.v0.live <- min.v0s$min.v0s[min.v0s$treatment == "live"]


# Pivot to a wider format with columns for live and killed activities, and then
#   adjust activities so that 
#   a.) activities below detection limit are set to detection limit, and 
#   b.) activities that have been adjusted like that are marked
live_killed <- samp_slopes %>%
  select(depth.mbsf, substrate, treatment, class, enzyme, v0, v0.se) %>%
  pivot_wider(id_cols = c(depth.mbsf, substrate, class, enzyme), # note id_cols defaults to everything not in names_from and values_from
              names_from = treatment,
              values_from = c(v0, v0.se)) %>%
  mutate(v0.adj_killed = case_when(v0_killed < min.v0.killed ~ min.v0.killed,
                                   v0_killed >= min.v0.killed ~ v0_killed),
         v0.adj_live = case_when(v0_live < min.v0.live ~ min.v0.live,
                                 v0_live >= min.v0.live ~ v0_live)) %>%
  mutate(v0.min_live = case_when(v0.adj_live - v0.se_live <= 1e-2 ~ 1e-2,
                                 TRUE ~ v0.adj_live - v0.se_live),
         v0.min_killed = case_when(v0.adj_killed - v0.se_killed <= 1e-2 ~ 1e-2,
                                   TRUE ~ v0.adj_killed - v0.se_killed),
         live.se.is.truncated = case_when((v0.min_live == 1e-2)  ~ TRUE,
                                     TRUE ~ FALSE),
         killed.se.is.truncated = case_when((v0.min_killed == 1e-2)  ~ TRUE,
                                            TRUE ~ FALSE))

# Plot of activity retained after autoclaving vs activity before autoclaving
p_live_killed <- ggplot(live_killed, aes(x=v0.adj_live, y=v0.adj_killed)) + 
  geom_smooth(aes(group=1), method = "lm", color="black") +
  geom_crossbar(aes(xmin=v0.min_live, xmax=v0.adj_live, colour = live.se.is.truncated), size = line.size) +
  geom_crossbar(aes(xmin=v0.adj_live, xmax=v0.adj_live+v0.se_live), size = line.size) +
  geom_crossbar(aes(ymin=v0.min_killed, ymax=v0.adj_killed, colour=killed.se.is.truncated), size = line.size) +
  geom_crossbar(aes(ymin=v0.adj_killed, ymax=v0.adj_killed+v0.se_killed), size = line.size) +
  geom_point(size=point.size) +
  scale_x_log10(name = expression(paste("live ", v[0], ", ", "nmol ", g^{-1}, " sed ", hr^{-1}))) + 
  scale_y_log10(name = expression(paste("autoclaved ", v[0], ", ", "nmol ", g^{-1}, " sed ", hr^{-1}))) + 
  scale_colour_manual(values = c("black", "gray75"), guide=FALSE) +
  theme_bw() + 
  theme(text = element_text(size=9),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


if(print.plots) {
    print(p_live_killed)
}
# if(save.plots) {
#   ggsave("plots/live_vs_killed.png", p_live_killed, height=2.5, width=3, units = "in", dpi=300)
# }

######
# Calculate average fraction of enzyme activity lost
#####
frac_killed <- live_killed %>%
  filter(depth.mbsf < 50 & v0_live > v0.min_live) %>% 
  mutate(frac.killed = v0.adj_killed / v0.adj_live) %>% 
  mutate(class.fac = fct_reorder(factor(class), frac.killed, median, .desc = TRUE))

total_frac_killed <- frac_killed %>%
  summarise(median.frac = mean(frac.killed, na.rm=TRUE),
           iqr.25 = quantile(frac.killed, 0.25),
           iqr.75 = quantile(frac.killed, 0.75),
           n = n())
print(total_frac_killed)
  
# For all classes
frac_killed_by_class <- frac_killed %>%
  group_by(class) %>%
  summarise(median.frac = mean(frac.killed, na.rm=TRUE),
            iqr.25 = quantile(frac.killed, 0.25),
            iqr.75 = quantile(frac.killed, 0.75),
            n = n())
print(frac_killed_by_class)
conover.test::conover.test(frac_killed$frac.killed, frac_killed$class, kw=TRUE) # See below; nonparametric tests are necessary

class_lines <- data.frame(x=c(0.5, 1.5), xend=c(2.5, 3.5), y=c(2.65, 2.5), yend=c(2.65, 2.5))

# Plot all classes
p_classes <- ggplot(frac_killed, aes(x=class.fac, y=frac.killed)) + 
  geom_boxplot(size = line.size, outlier.shape = NA) + 
  geom_point(position = position_dodge2(width = 0.3), size = point.size) + 
  scale_y_continuous(name = "activity after autoclaving", labels = scales::percent) + 
  geom_segment(data=  class_lines, aes(x=x, xend=xend, y=y, yend=yend)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=-45, hjust =0))
if(print.plots) {
  print(p_classes)
}

# For just the peptidases
frac_killed_peptidase_summ <- frac_killed %>%
  filter(class == "peptidase") %>%
  group_by(enzyme) %>%
  summarise(median.frac = mean(frac.killed, na.rm=TRUE),
            iqr.25 = quantile(frac.killed, 0.25),
            iqr.75 = quantile(frac.killed, 0.75),
            n = n()) %>%
  arrange(desc(median.frac)) 
print(frac_killed_peptidase_summ)
conover.test::conover.test(frac_killed$frac.killed[frac_killed$class == "peptidase"], frac_killed$enzyme[frac_killed$class == "peptidase"], kw=TRUE) # See below; nonparametric tests are necessary

# Show differences with lines
peptidase_lines <- data.frame(x=c(1.5, 3.5), xend = c(4.5, 6.5), y=c(1.6, 1.5), yend = c(1.6, 1.5))

frac_killed_peptidase <- frac_killed %>%
  filter(class=="peptidase") %>%
  mutate(enzyme.fac = factor(enzyme),
         enzyme.fac = fct_reorder(enzyme.fac, frac.killed, median, .desc = TRUE))
  

p_frac_peptidases_killed <- ggplot(frac_killed_peptidase, aes(x=enzyme.fac, y=frac.killed)) +
  geom_boxplot(size = line.size, outlier.shape = NA) +
  geom_point(position=position_dodge2(width=0.3), size = point.size) + 
  geom_segment(data = peptidase_lines, aes(x=x, xend=xend, y=y, yend=yend)) + 
  #geom_text(data = peptidase_labels, aes(x=enzyme.fac, y=median.frac, label = label)) + 
  scale_y_continuous(name = "activity after autoclaving", labels = scales::percent) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=-45, hjust =0))
if(print.plots) {
  print(p_frac_peptidases_killed)
}


spacing <- theme(plot.margin = unit(c(0.5, 0.3, 0 ,0.5), "cm")) # Create custom spacing around the plots to make the panel labels look right

tiff("plots/activity_differences_combined.tiff", height = 2.5, width = 7.3, units = "in", res= 300, compression = "lzw")
cowplot::plot_grid(p_live_killed + spacing,
                   p_classes + spacing, 
                   p_frac_peptidases_killed + spacing, 
                   labels=c("A", "B", "C"), 
                   rel_widths = c(1, 0.7, 1), 
                   label_size = 10, 
                   nrow = 1)
dev.off()




















#######
# Some linear models
#######

frac_killed_lm_all <- lm(frac.killed ~ class, data= frac_killed)
#plot(frac_killed_lm_all) # These are fairly skewed and have some weird leverage things going on
kruskal.test(frac.killed ~ class, data= frac_killed) # p=1.739e-5
summary(lm(rank(frac.killed) ~ rank(class), data= frac_killed))
conover.test::conover.test(frac_killed$frac.killed, frac_killed$class, kw=TRUE)

live_killed_mod <- lm(v0.adj_killed ~ v0.adj_live, data = live_killed)
summary(live_killed_mod)

live_killed_mod_subs <- lm(v0.adj_killed ~ v0.adj_live + substrate, data = live_killed)
summary(live_killed_mod_subs)

#  Test the idea that glycosylases have a different relationship between live and killed than peptidases
gp <- live_killed %>%
  filter(class != "phosphatase")

ggplot(gp, aes(x=v0.adj_live, y=v0.adj_killed, shape = class)) +
  geom_point() +
  geom_smooth(method = "lm", aes(linetype = class), colour = "black") + 
  scale_x_log10(name = expression(paste("live ", v[0], ", nmol g ", sed^{-1}, ", ", hr^{-1}))) +
  scale_y_log10(name = expression(paste("autoclaved ", v[0], ", nmol g ", sed^{-1}, ", ", hr^{-1}))) + 
  scale_shape_manual(values = c(19, 1))
ggsave("plots/live_vs_killed_by_class.png", height = 3, width = 4, units = "in", dpi = 300)

frac_surviving <- live_killed %>%
  mutate(frac.surviving = v0.adj_killed / v0.adj_live) # remember that v0.adj_killed actually represents what is left after killing - thus, the fraction surviving

ggplot(frac_surviving, aes(x=class, y=frac.surviving, colour = enzyme)) + 
  geom_boxplot() +
  #geom_point(position = position_dodge(width = 0.1)) + 
  scale_y_log10()

# Linear models of fraction surviving
class_mod <- lm(frac.surviving ~ class, data = frac_surviving %>% filter (depth.mbsf < 50))
summary(class_mod)
print(TukeyHSD(aov(class_mod)))
plot(TukeyHSD(aov(class_mod)))

peptidase_mod <- lm(frac.surviving ~ enzyme, data = frac_surviving %>% filter(class == "peptidase" & depth.mbsf < 50))
summary(peptidase_mod)
plot(TukeyHSD(aov(peptidase_mod)))
print(TukeyHSD(aov(peptidase_mod)))



median_frac_surviving <- frac_surviving %>%
  subset(depth.mbsf < 50) %>%
  group_by(class) %>%
  summarise(median.frac.surviving = quantile(frac.surviving, 0.5),
            low.iqr = quantile(frac.surviving, 0.25),
            hi.iqr = quantile(frac.surviving, 0.75),
            mean = mean(frac.surviving, na.rm = TRUE),
            sd = sd(frac.surviving, na.rm = TRUE))

print(median_frac_surviving)


summary_peptidases_surviving <- frac_surviving %>%
  filter(depth.mbsf < 50 & class == "peptidase") %>%
  group_by(enzyme) %>%
  summarise(median.frac.surviving = quantile(frac.surviving, 0.5),
            low.iqr = quantile(frac.surviving, 0.25),
            hi.iqr = quantile(frac.surviving, 0.75),
            mean = mean(frac.surviving, na.rm = TRUE),
            sd = sd(frac.surviving, na.rm = TRUE))
print(summary_peptidases_surviving)

frac_surviving %>%
  filter(depth.mbsf < 50 & class == "peptidase" & enzyme != "ornithyl AP") %>%
  summarise(median.frac.surviving = quantile(frac.surviving, 0.5),
            low.iqr = quantile(frac.surviving, 0.25),
            hi.iqr = quantile(frac.surviving, 0.75),
            mean = mean(frac.surviving, na.rm = TRUE),
            sd = sd(frac.surviving, na.rm = TRUE))

# Create linear models of fraction surviving by class
class_mod <- lm(frac.surviving ~ class, data = median_frac_surviving)


ggplot()