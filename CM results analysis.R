rm(list=ls())

#remotes::install_github("jvparidon/lmerMultiMember")
library(lmerMultiMember)
library(lme4)
library(car)
library(jtools)
library(parameters)
library(ggplot2)
library(emmeans)
library(dplyr)

analysis_data_path <- "path/CM_analysis_df.csv"
plot_save_path <- "path/"

analysis_data <- read.csv(analysis_data_path)



analysis_data$add_train_relevel <- relevel(as.factor(analysis_data$add_train), "TRUE")
analysis_data$treatment_relevel <- relevel(as.factor(analysis_data$treatment), "TR")
analysis_data$release_cond_relevel <- relevel(as.factor(analysis_data$release_cond), "solo")

### Sample sizes

analysis_data[which(!is.na(analysis_data$Distance_closest_analysis)),] %>%
  group_by(release_cond, treatment, add_train) %>%
  summarize(count_dist_closest = sum(!is.na(Distance_closest_analysis)),
            count_sub_site = n_distinct(bird, site))

analysis_data[which(!is.na(analysis_data$Distance_average_analysis)),] %>%
  group_by(release_cond, treatment, add_train) %>%
  summarize(count_dist_average = sum(!is.na(Distance_average_analysis)),
            count_sub_site = n_distinct(bird, site))

analysis_data[which(!is.na(analysis_data$HEI_analysis)),] %>%
  group_by(release_cond, treatment, add_train) %>%
  summarize(count_hei = sum(!is.na(HEI_analysis)),
            count_sub_site = n_distinct(bird, site))

############ Main analysis plots


p <- ggplot(data=analysis_data, 
            aes(x=treatment, color=release_cond, y=Distance_average_analysis))+
  facet_wrap(vars(add_train), 
             labeller = labeller(add_train = c("FALSE" = "Forgetting\ngroup", 
                                               "TRUE" = "Extra training\ngroup")))+
  geom_boxplot()+
  scale_x_discrete(limits=c("TR", "CM"), 
                   labels=c("End of\ntraining", "Memory\ntest"))+
  geom_point(position = position_dodge(width = 0.75))+
  scale_y_log10()+
  theme_nice()+
  scale_color_manual(values=c("blue", "red"))+
  labs(x="", y="Average distance to baseline tracks (m)", color="")
p
#ggsave(paste(plot_save_path, "Average_distance_baseline.png", sep=""), plot = p, width = 6, height = 4.5, dpi = 300)

p2 <- ggplot(data=analysis_data, 
             aes(x=treatment, color=release_cond, y=Distance_closest_analysis))+
  facet_wrap(vars(add_train), 
             labeller = labeller(add_train = c("FALSE" = "Forgetting\ngroup", 
                                               "TRUE" = "Extra training\ngroup")))+
  geom_boxplot()+
  scale_x_discrete(limits=c("TR", "CM"), 
                   labels=c("End of\ntraining", "Memory\ntest"))+
  geom_point(position = position_dodge(width = 0.75))+
  scale_y_log10()+
  theme_nice()+
  scale_color_manual(values=c("blue", "red"))+
  labs(x="", y="Mean distance to nearest baseline track (m)", color="")
p2
#ggsave(paste(plot_save_path, "Closest_distance_baseline.png", sep=""), , plot = p2, width = 6, height = 4.5, dpi = 300)

p2a <- ggplot(data=analysis_data[analysis_data$add_train==FALSE,], 
             aes(x=treatment, color=release_cond, y=Distance_closest_analysis))+
  geom_boxplot()+
  scale_x_discrete(limits=c("TR", "CM"), 
                   labels=c("End of\ntraining", "Memory\ntest"))+
  geom_point(position = position_dodge(width = 0.75))+
  scale_y_log10()+
  theme_nice()+
  scale_color_manual(values=c("blue", "red"))+
  labs(x="", y="Distance to baseline routes (m)", color="")
p2a
#ggsave(paste(plot_save_path, "Fig1c.png", sep=""), plot = p2a, width = 4, height = 6.5, dpi = 300)

p2b <- ggplot(data=analysis_data[analysis_data$add_train==T,], 
              aes(x=treatment, color=release_cond, y=Distance_closest_analysis))+
  geom_boxplot()+
  scale_x_discrete(limits=c("TR", "CM"), 
                   labels=c("End of\ntraining", "Memory\ntest"))+
  geom_point(position = position_dodge(width = 0.75))+
  scale_y_log10()+
  theme_nice()+
  scale_color_manual(values=c("blue", "red"))+
  labs(x="", y="Distance to baseline routes (m)", color="")
p2b
#ggsave(paste(plot_save_path, "Fig2.png", sep=""), plot = p2b, width = 3, height = 4, dpi = 300)


p3 <- ggplot(data=analysis_data, 
             aes(x=treatment, color=release_cond, y=HEI_analysis))+
  facet_wrap(vars(add_train), 
             labeller = labeller(add_train = c("FALSE" = "Forgetting\ntreatment", 
                                               "TRUE" = "Extra training\ntreatment")))+
  geom_boxplot()+
  scale_x_discrete(limits=c("TR", "CM"), 
                   labels=c("End of\ntraining", "Memory\ntest"))+
  scale_color_manual(values=c("blue", "red"))+
  theme_nice()+
  ylim(c(0.27, 1))+
  geom_point(position = position_dodge(width = 0.75))+
  labs(x="", y="Homing Efficiency Index", color="")
p3
#ggsave(paste(plot_save_path, "Fig3.png", sep=""), 
#plot = p3, width = 5, height = 4, dpi = 300)


W_merge <- weights_from_vector(c(paste(analysis_data$bird, 
                                       analysis_data$partner, sep = ",")))
W_merge <- W_merge[-1,]


mod_dist_closest <- lmerMultiMember::lmer(data=analysis_data, 
                                              log(Distance_closest_analysis)~ 
                                                release_cond*treatment*add_train
                                              + site +(1|names)
                                              +(1|Pair_ID),
                                              memberships = list(names = W_merge))
anova(mod_dist_closest)
summary(mod_dist_closest)
parameters::p_value(mod_dist_closest, method="wald")



mod_dist_average <- lmerMultiMember::lmer(data=analysis_data, 
                                          log(Distance_average_analysis)~ 
                                            release_cond*treatment*add_train_relevel
                                          + site +(1|names)
                                          +(1|Pair_ID),
                                          memberships = list(names = W_merge))
anova(mod_dist_average)
summary(mod_dist_average)
parameters::p_value(mod_dist_average, method="wald")


analysis_data$transformed_hei <- log((1-analysis_data$HEI_analysis)/analysis_data$HEI_analysis)
mod_hei <- lmerMultiMember::lmer(data=analysis_data, 
                                 transformed_hei~ 
                                   release_cond_relevel*treatment*add_train
                                 + site +(1|names)
                                 +(1|Pair_ID),
                                 memberships = list(names = W_merge))
anova(mod_hei)
summary(mod_hei)
parameters::p_value(mod_hei, method="wald")

############ ASSUMPTIONS CHECK
model_to_check <- mod_dist_average
#MODS TO CHECK
#mod_dist_closest
#mod_dist_average
#mod_hei


# Check residuals
residuals <- resid(model_to_check)
fitted_values <- fitted(model_to_check)

# 1. QQ-plot for normality
qqnorm(residuals)
qqline(residuals)

# 2. Kolmogorov-Smirnov test for normality
ks.test(residuals, "pnorm", mean(residuals), sd(residuals))

# 3. Residuals vs. Fitted values plot
plot(fitted_values, residuals)
abline(h = 0, col = "red")
title("Residuals vs Fitted Values")

# 4. Residuals independence plot
plot(residuals)
title("Residuals Plot")



############ BEST INDS ANALYSIS

df_pairs <- analysis_data[which(analysis_data$release_cond=="pair"),]

df_best_inds_closest <- analysis_data[which(!is.na(analysis_data$Distance_closest_analysis)),] %>%
  group_by(Pair_ID, date, site) %>%
  filter(n_distinct(bird) == 2) %>%  # Ensure both individuals have a record on that day
  slice_min(Distance_closest_analysis, n = 1) %>%    # Keep the row with the highest performance value
  ungroup()

df_pair_bestinds_closest <- bind_rows(df_pairs[which(!is.na(df_pairs$Distance_closest_analysis)),], 
                                      df_best_inds_closest)


df_pair_bestinds_closest %>%
  group_by(release_cond, treatment, add_train) %>%
  summarize(count_dist_closest = sum(!is.na(Distance_closest_analysis)),
            count_sub_site = n_distinct(bird, site))



p4 <- ggplot(data=df_pair_bestinds_closest, 
             aes(x=treatment_addtrain, color=release_cond, y=Distance_closest_analysis))+
  geom_boxplot()+
  scale_x_discrete(limits=c("TR", "CM_N", "ATR", "CM_A"), 
                   labels=c("End of \ntraining", "Memory test", "End of \nextra training", "Memory test; \nextra training \ngroup"))+  
  geom_point(position = position_dodge(width = 0.75))+
  scale_y_log10()+
  theme_nice()+
  scale_color_manual(values=c("black", "red"), labels=c("pair", "best solo"))+
  labs(x="", y="Average distance to nearest baseline track (m)", color="")
p4


W_merge_2 <- weights_from_vector(c(paste(df_pair_bestinds_closest$bird, 
                                         df_pair_bestinds_closest$partner, sep = ",")))
W_merge_2 <- W_merge_2[-1,]

mod_dist_closest_best_inds <- lmerMultiMember::lmer(data=df_pair_bestinds_closest, 
                                                    log(Distance_closest_analysis)~ 
                                                      release_cond*treatment*add_train_relevel
                                                    + site +(1|names)
                                                    +(1|Pair_ID),
                                                    memberships = list(names = W_merge_2))
anova(mod_dist_closest_best_inds)
summary(mod_dist_closest_best_inds)
parameters::p_value(mod_dist_closest_best_inds)


df_best_inds_average <- analysis_data[which(!is.na(analysis_data$Distance_average_analysis)),] %>%
  group_by(Pair_ID, date, site) %>%
  filter(n_distinct(bird) == 2) %>%  # Ensure both individuals have a record on that day
  slice_min(Distance_average_analysis, n = 1) %>%    # Keep the row with the highest performance value
  ungroup()


df_pair_bestinds_average <- bind_rows(df_pairs[which(!is.na(df_pairs$Distance_average_analysis)),], 
                                      df_best_inds_average)


df_pair_bestinds_average %>%
  group_by(release_cond, treatment, add_train) %>%
  summarize(count_dist_av = sum(!is.na(Distance_average_analysis)),
            count_sub_site = n_distinct(bird, site))


p5 <- ggplot(data=df_pair_bestinds_average, 
             aes(x=treatment_addtrain, color=release_cond, y=Distance_average_analysis))+
  geom_boxplot()+
  scale_x_discrete(limits=c("TR", "CM_N", "ATR", "CM_A"), 
                   labels=c("End of \ntraining", "Memory test", "End of \nextra training", "Memory test; \nextra training \ngroup"))+  
  geom_point(position = position_dodge(width = 0.75))+
  scale_y_log10()+
  theme_nice()+
  scale_color_manual(values=c("black", "red"), labels=c("pair", "best solo"))+
  labs(x="", y="Average distance to baseline tracks (m)", color="")
p5


W_merge_2b <- weights_from_vector(c(paste(df_pair_bestinds_average$bird, 
                                         df_pair_bestinds_average$partner, sep = ",")))
W_merge_2b <- W_merge_2b[-1,]

mod_dist_average_best_inds <- lmerMultiMember::lmer(data=df_pair_bestinds_average, 
                                                    log(Distance_average_analysis)~ 
                                                      release_cond*treatment*add_train_relevel
                                                    + site +(1|names)
                                                    +(1|Pair_ID),
                                                    memberships = list(names = W_merge_2b))
anova(mod_dist_average_best_inds)
summary(mod_dist_average_best_inds)
parameters::p_value(mod_dist_average_best_inds)

############ ASSUMPTIONS CHECK
model_to_check <- mod_dist_average_best_inds
#MODS TO CHECK
#mod_dist_closest_best_inds
#mod_dist_average_best_inds


# Check residuals
residuals <- resid(model_to_check)
fitted_values <- fitted(model_to_check)

# 1. QQ-plot for normality
qqnorm(residuals)
qqline(residuals)

# 2. Kolmogorov-Smirnov test for normality
ks.test(residuals, "pnorm", mean(residuals), sd(residuals))

# 3. Residuals vs. Fitted values plot
plot(fitted_values, residuals)
abline(h = 0, col = "red")
title("Residuals vs Fitted Values")

# 4. Residuals independence plot
plot(residuals)
title("Residuals Plot")


