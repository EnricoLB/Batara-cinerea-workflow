setwd("/Users/en_gl/Documents/Unesp/Mestrado/B.cinerea/")

library(MASS)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(emmeans)
library(nlme)
library(tidyr)
library(multcompView)
library(multcomp)
library(ggExtra)
library(ggrepel)
library(cowplot)

# ============================================================
# READ & PREPARE MORPHOLOGICAL DATA
# ============================================================
morfo_raw <- read.csv(
  "Análises/Morfológica/morfo_batara.csv",
  header = TRUE,
  sep    = ",",
  na.strings = c("", "NA"),
  check.names = FALSE
)

# Rename columns to internal names
morfo <- morfo_raw
names(morfo)[names(morfo) == "Subespécie"]                  <- "ssp"
names(morfo)[names(morfo) == "Sexo"]                        <- "sex"
names(morfo)[names(morfo) == "Comprimento da Asa (cm)"]     <- "wing"
names(morfo)[names(morfo) == "Comprimento da Cauda (cm)"]   <- "tail"
names(morfo)[names(morfo) == "Comprimento da Crista (cm)"]  <- "crest"
names(morfo)[names(morfo) == "Comprimento do Bico (mm)"]    <- "beak"
names(morfo)[names(morfo) == "Largura do Bico (mm)"]        <- "beak_w"
names(morfo)[names(morfo) == "Altura do Bico (mm)"]         <- "beak_h"
names(morfo)[names(morfo) == "Curvatura do Bico (mm)"]      <- "beak_c"
names(morfo)[names(morfo) == "Tarso - Metatarso (mm)"]      <- "tarsus"
names(morfo)[names(morfo) == "Museu"]                       <- "researcher"

# Convert measurements to numeric
morfo[c("wing","tail","beak","tarsus")] <- lapply(
  morfo[c("wing","tail","beak","tarsus")], as.numeric
)

morfo$ssp  <- factor(morfo$ssp)
morfo$sex  <- factor(morfo$sex)

morfo <- na.omit(morfo)

head(morfo)
summary(morfo)
str(morfo)


# ============================================================
# MALES vs FEMALES — MORPHO (t-test / Mann-Whitney)
# ============================================================
variaveis_m <- morfo[morfo$sex == "M", ]
variaveis_f <- morfo[morfo$sex == "F", ]

variaveis_a_comparar <- c("wing", "tail", "beak", "tarsus")

# Shapiro-Wilk normality test
resultados_shapiro <- lapply(variaveis_a_comparar, function(var) {
  list(
    M = shapiro.test(variaveis_m[[var]]),
    F = shapiro.test(variaveis_f[[var]])
  )
})
names(resultados_shapiro) <- variaveis_a_comparar
resultados_shapiro

# t-tests
resultados_t <- lapply(variaveis_a_comparar, function(var) {
  t.test(variaveis_m[[var]], variaveis_f[[var]], var.equal = FALSE)
})
names(resultados_t) <- variaveis_a_comparar
resultados_t

# Mann-Whitney
resultados_mann_whitney <- lapply(variaveis_a_comparar, function(var) {
  wilcox.test(variaveis_m[[var]], variaveis_f[[var]])
})
names(resultados_mann_whitney) <- variaveis_a_comparar
resultados_mann_whitney


# ============================================================
# MANOVA — MORPHO (by subspecies, split by sex)
# ============================================================
morfo_m <- morfo[morfo$sex == "M", ]
morfo_f <- morfo[morfo$sex == "F", ]

# Drop crest if present (as in original)
morfo_m$crest <- NULL
morfo_f$crest <- NULL

manova_model_m <- manova(cbind(wing, tail, beak, tarsus) ~ ssp, data = morfo_m)
summary(manova_model_m)
summary.aov(manova_model_m)

manova_model_f <- manova(cbind(wing, tail, beak, tarsus) ~ ssp, data = morfo_f)
summary(manova_model_f)
summary.aov(manova_model_f)

# Post-hoc Tukey (example: tarsus in males)
TukeyHSD(aov(tarsus ~ ssp, data = morfo_m))


# ============================================================
# TABLE 3 — Descriptive stats per subspecies (wing, females)
# ============================================================
resultado_wing_f <- morfo_f %>%
  dplyr::filter(!is.na(ssp) & !is.na(wing)) %>%
  dplyr::group_by(ssp) %>%
  dplyr::summarise(
    n          = n(),
    Media      = mean(wing, na.rm = TRUE),
    Minimo     = min(wing,  na.rm = TRUE),
    Maximo     = max(wing,  na.rm = TRUE),
    DesvioPadrao = sd(wing, na.rm = TRUE)
  )
print(resultado_wing_f)


# ============================================================
# HIERARCHICAL MODELS — MALES (researcher as random effect)
# ============================================================
data_m <- morfo[morfo$sex == "M", ]
data_m$ssp <- factor(data_m$ssp)
data_m$ssp <- relevel(data_m$ssp, ref = levels(data_m$ssp)[1])  # adjust ref if needed

# Wing
data_m$wing <- as.numeric(data_m$wing)
data_m <- na.omit(data_m)
model_w <- lme(wing ~ ssp, random = ~1 | researcher, data = data_m)
summary(model_w)
emm <- emmeans(model_w, ~ ssp)
posthoc_w <- pairs(emm, adjust = "tukey")
print(posthoc_w)

# Tail
data_m$tail <- as.numeric(data_m$tail)
data_m <- na.omit(data_m)
model_tl <- lme(tail ~ ssp, random = ~1 | researcher, data = data_m)
summary(model_tl)
emm <- emmeans(model_tl, ~ ssp)
posthoc_tl <- pairs(emm, adjust = "tukey")
print(posthoc_tl)

# Beak
data_m$beak <- as.numeric(data_m$beak)
data_m <- na.omit(data_m)
model_bk <- lme(beak ~ ssp, random = ~1 | researcher, data = data_m)
summary(model_bk)
emm <- emmeans(model_bk, ~ ssp)
posthoc_bk <- pairs(emm, adjust = "tukey")
print(posthoc_bk)

# Tarsus
data_m$tarsus <- as.numeric(data_m$tarsus)
data_m <- na.omit(data_m)
model_tr <- lme(tarsus ~ ssp, random = ~1 | researcher, data = data_m)
summary(model_tr)
emm <- emmeans(model_tr, ~ ssp)
posthoc_tr <- pairs(emm, adjust = "tukey")
print(posthoc_tr)


# ============================================================
# SUMMARY TABLE — p-values and estimates from LME models
# ============================================================
models_list  <- list(model_w, model_tl, model_bk, model_tr)
measurements <- c("wing", "tail", "beak", "tarsus")

summary_results <- lapply(models_list, function(model) {
  if (!"tTable" %in% names(summary(model))) return(NA)
  coef_table <- summary(model)$tTable
  ssp_rows   <- grep("^ssp", rownames(coef_table))
  data.frame(
    term      = rownames(coef_table)[ssp_rows],
    estimate  = coef_table[ssp_rows, "Value"],
    p_value   = coef_table[ssp_rows, "p-value"]
  )
})

res_list <- Map(function(df, meas) {
  if (is.data.frame(df)) {
    df$Measurement <- meas
    df
  } else {
    data.frame(term = NA, estimate = NA, p_value = NA, Measurement = meas)
  }
}, summary_results, measurements)

summary_table <- do.call(rbind, res_list)
summary_table$term <- sub("^ssp", "", summary_table$term)
summary_table <- summary_table[, c("Measurement", "term", "estimate", "p_value")]
print(summary_table)

# Plot p-values
ggplot(summary_table, aes(x = Measurement, y = p_value)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "Comparison of p-values for each Measurement",
       y = "p-value", x = "Measurement") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ============================================================
# LDA — MORPHO
# ============================================================
my_data <- morfo
my_data[c("wing","tail","beak","tarsus")] <- lapply(
  my_data[c("wing","tail","beak","tarsus")], as.numeric
)
my_data <- na.omit(my_data)

my_data_f <- my_data[my_data$sex == "F", ]
my_data_m <- my_data[my_data$sex == "M", ]

# Perform LDA
lda_result <- lda(ssp ~ wing + tail + beak + tarsus, data = my_data)
print(lda_result)

lda_pred <- predict(lda_result)
lda_data <- data.frame(lda_pred$x, ssp = my_data$ssp, sex = my_data$sex)

# Centroids
centroids_ssp <- aggregate(cbind(LD1, LD2) ~ ssp, data = lda_data, FUN = mean)

# Vector scaling
lda_scaling <- as.data.frame(lda_result$scaling)
lda_scaling$variable <- rownames(lda_scaling)
max_loading_magnitude <- max(sqrt(lda_scaling$LD1^2 + lda_scaling$LD2^2))
vec_scale <- 1.5 / max_loading_magnitude

# Plotting config
pad_range <- function(v, p = 0.08) {
  r <- range(v, na.rm = TRUE); d <- diff(r)
  c(r[1] - p*d, r[2] + p*d)
}
xlims <- pad_range(lda_data$LD1)
ylims <- pad_range(lda_data$LD2)

# Define colors dynamically from subspecies levels
ssp_levels <- levels(my_data$ssp)
# Assign colors — adjust palette as needed
base_cols <- c("black","darkorange","darkgreen","steelblue","purple","brown")
cols <- setNames(base_cols[seq_along(ssp_levels)], ssp_levels)

# Main plot

# Make sure nothing above ends with a dangling +

p_main <- ggplot(lda_data, aes(x = LD1, y = LD2)) +
  stat_ellipse(aes(fill = ssp, color = ssp),
               geom = "polygon", alpha = 0.12,
               level = 0.75, show.legend = FALSE) +
 # stat_ellipse(aes(color = ssp),
  #             level = 0.95, linewidth = 0.5,
   #            linetype = "dashed", show.legend = FALSE) +
  geom_point(aes(fill = ssp, shape = ssp),
             colour = "white", stroke = 0.4, size = 2.8, alpha = 0.85) +
  geom_point(data = centroids_ssp, aes(x = LD1, y = LD2, fill = ssp),
             shape = 21, size = 5, color = "white", stroke = 0.8,
             show.legend = FALSE) +
  geom_label_repel(data = centroids_ssp, aes(x = LD1, y = LD2, label = ssp),
                   size = 3, label.size = 0, seed = 1, min.segment.length = 0,
                   fill = "white", alpha = 0.9, box.padding = 0.25,
                   point.padding = 0.2, segment.color = "grey50",
                   show.legend = FALSE) +
  geom_segment(data = lda_scaling,
               aes(x = 0, y = 0,
                   xend = LD1 * vec_scale,
                   yend = LD2 * vec_scale),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "blue", linewidth = 0.8, inherit.aes = FALSE) +
  geom_text_repel(data = lda_scaling,
                  aes(x = LD1 * vec_scale,
                      y = LD2 * vec_scale,
                      label = variable),
                  color = "blue", fontface = "bold", size = 3.5,
                  nudge_x = 0.1, nudge_y = 0.1, inherit.aes = FALSE) +
  scale_fill_manual(values = cols, breaks = names(cols), name = "Subspecies") +
  scale_colour_manual(values = cols, breaks = names(cols), guide = "none") +
  scale_shape_manual(values = c(21, 22, 24, 25, 23)[seq_along(ssp_levels)],
                     name = "Subspecies") +
  scale_x_continuous(limits = xlims) +
  scale_y_continuous(limits = ylims) +
  labs(x = "LDA 1", y = "LDA 2") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 10) +
  theme(
    legend.position    = "right",
    legend.title       = element_blank(),
    legend.text        = element_text(face = "italic"),
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(linewidth = 0.2, colour = "grey92"),
    plot.background    = element_rect(fill = "transparent", colour = NA),
    panel.background   = element_rect(fill = "transparent", colour = "grey70"),
    axis.title         = element_text(colour = "grey20")
  )
p_top <- ggplot(lda_data, aes(x = LD1, fill = ssp)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(limits = xlims, expand = c(0, 0)) +
  theme_void() + theme(legend.position = "none")

p_right <- ggplot(lda_data, aes(x = LD2, fill = ssp)) +
  geom_density(alpha = 0.3) + coord_flip() +
  scale_fill_manual(values = cols) +
  scale_x_continuous(limits = ylims, expand = c(0, 0)) +
  theme_void() + theme(legend.position = "none")

legend           <- get_legend(p_main)
p_main_no_legend <- p_main + theme(legend.position = "none")
placeholder      <- ggplot() + theme_void()

bottom_row <- plot_grid(p_main_no_legend, p_right, ncol = 2,
                        align = "h", axis = "tb", rel_widths = c(4, 1))
top_row    <- plot_grid(p_top, placeholder, ncol = 2, rel_widths = c(4, 1))

final_plot_morpho <- plot_grid(
  plot_grid(top_row, bottom_row, nrow = 2, rel_heights = c(1, 4)),
  legend, ncol = 2, rel_widths = c(1, 0.15)
)
final_plot_morpho



# Save
output_file_morpho <- "/Users/en_gl/Documents/Unesp/Mestrado/B.cinerea/Plots/lda_morpho_plot_hd.png"
save_plot(output_file_morpho, final_plot_morpho,
          base_height = 6, base_width = 8, dpi = 600)


# ============================================================
# READ & PREPARE VOCAL DATA
# ============================================================
vocal_raw <- read.csv(
  "Vocal/Combined data/batara.csv",
  header    = TRUE,
  sep       = ",",
  na.strings = c("", "NA"),
  check.names = TRUE   # keeps R-safe names like Low.Freq..Hz.
)

vocal <- vocal_raw

# Rename to internal names
names(vocal)[names(vocal) == "Low.Freq..Hz."]   <- "low.fq"
names(vocal)[names(vocal) == "High.Freq..Hz."]  <- "high.fq"
names(vocal)[names(vocal) == "Peak.Freq..Hz."]  <- "peak.fq"
names(vocal)[names(vocal) == "Delta.Time..s."]  <- "time"
names(vocal)[names(vocal) == "Center.Time..s."] <- "center.time"
names(vocal)[names(vocal) == "Dur.90...s."]     <- "dur.90"
names(vocal)[names(vocal) == "BW.90...Hz."]     <- "bw.90"
names(vocal)[names(vocal) == "source"]          <- "Taxon"

vocal$peak.fq <- as.numeric(vocal$peak.fq)
vocal$Taxon   <- as.factor(vocal$Taxon)

head(vocal)
summary(vocal)
str(vocal)


# ============================================================
# MALES vs FEMALES — VOCAL (ANOVA per variable)
# ============================================================
variaveis_a_comparar_v <- c("low.fq","high.fq","peak.fq","time","center.time","dur.90","bw.90")

resultados_aov <- lapply(variaveis_a_comparar_v, function(var) {
  formula <- as.formula(paste(var, "~ Taxon"))
  aov(formula, data = vocal)
})
names(resultados_aov) <- variaveis_a_comparar_v
lapply(resultados_aov, summary)

# Tukey example — peak.fq
TukeyHSD(resultados_aov$peak.fq)


# ============================================================
# MANOVA — VOCAL
# ============================================================
manova_vocal <- manova(
  cbind(low.fq, high.fq, peak.fq, time, center.time, dur.90, bw.90) ~ Taxon,
  data = vocal
)
summary(manova_vocal)
summary.aov(manova_vocal)

TukeyHSD(aov(peak.fq ~ Taxon, data = vocal))


# ============================================================
# TUKEY TABLE — all significant variables
# ============================================================
get_p <- function(v) {
  fit <- aov(as.formula(paste(v, "~ Taxon")), data = vocal)
  summary(fit)[[1]][["Pr(>F)"]][1]
}
pvals    <- sapply(variaveis_a_comparar_v, get_p)
alpha    <- 0.05
sig_vars <- names(pvals)[pvals < alpha]
sig_vars

tukey_list <- lapply(sig_vars, function(v) {
  fit  <- aov(as.formula(paste(v, "~ Taxon")), data = vocal)
  out  <- TukeyHSD(fit)$Taxon
  list(var = v, fit = fit, tukey = out)
})
names(tukey_list) <- sig_vars

tukey_df <- do.call(rbind, lapply(tukey_list, function(x) {
  df             <- as.data.frame(x$tukey)
  df$comparison  <- rownames(df)
  df$variable    <- x$var
  rownames(df)   <- NULL
  df
}))
tukey_df <- tukey_df %>% select(variable, comparison, diff, lwr, upr, `p adj`)
print(tukey_df)


# ============================================================
# TABLE 3 — Descriptive stats vocal (peak.fq by Taxon)
# ============================================================
resultado_vocal <- vocal %>%
  dplyr::filter(!is.na(Taxon)) %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarise(
    n            = n(),
    Media        = mean(peak.fq, na.rm = TRUE),
    Minimo       = min(peak.fq,  na.rm = TRUE),
    Maximo       = max(peak.fq,  na.rm = TRUE),
    DesvioPadrao = sd(peak.fq,   na.rm = TRUE)
  )
print(resultado_vocal)


# Significance letters for peak.fq boxplot
modelo_pk  <- aov(peak.fq ~ Taxon, data = vocal)
tukey_pk   <- TukeyHSD(modelo_pk)
letras     <- multcompLetters4(modelo_pk, tukey_pk)
letras_df  <- as.data.frame.list(letras$Taxon)
letras_df$Taxon <- rownames(letras_df)
medias     <- aggregate(peak.fq ~ Taxon, data = vocal, mean)
medias     <- merge(medias, letras_df, by = "Taxon")

ggplot(vocal, aes(x = Taxon, y = peak.fq, fill = Taxon)) +
  geom_boxplot(alpha = 0.7) +
  geom_text(data = medias, aes(x = Taxon, y = peak.fq + 100, label = Letters),
            size = 5, vjust = 0) +
  labs(y = "Peak frequency (Hz)", x = "Taxon") +
  theme_minimal() +
  theme(legend.position = "none")


# ============================================================
# LDA — VOCAL
# ============================================================

sapply(my_data_v[, c("low.fq","high.fq","peak.fq","time","center.time","dur.90","bw.90")], class)



my_data_v <- vocal

# Rebuild Taxon factor cleanly
my_data_v$Taxon <- factor(as.character(my_data_v$Taxon))

# Define only the columns the LDA actually needs
lda_cols <- c("Taxon", "low.fq", "high.fq", "peak.fq",
              "time", "center.time", "dur.90", "bw.90")

# Subset to only those columns BEFORE na.omit
my_data_v <- my_data_v[, lda_cols]

# Now na.omit only removes rows where these specific columns have NAs
my_data_v <- na.omit(my_data_v)

# Rebuild Taxon after omit
my_data_v$Taxon <- droplevels(factor(as.character(my_data_v$Taxon)))

# Confirm — should show 158 / 224 / 87
print(table(my_data_v$Taxon))
print(nrow(my_data_v))

# Now run LDA
lda_result_v <- lda(
  Taxon ~ low.fq + high.fq + peak.fq + time + center.time + dur.90 + bw.90,
  data = my_data_v
)
print(lda_result_v)

coeficientes    <- lda_result_v$scaling
coef_ld1        <- coeficientes[, 1]
coef_ld1_sorted <- sort(coef_ld1, decreasing = TRUE)
print(coef_ld1_sorted)

lda_pred_v  <- predict(lda_result_v)
lda_data_v  <- data.frame(lda_pred_v$x, Taxon = my_data_v$Taxon)
lda_data_v$ssp <- lda_data_v$Taxon

centroids_v <- aggregate(cbind(LD1, LD2) ~ ssp, data = lda_data_v, FUN = mean)

lda_scaling_v <- as.data.frame(lda_result_v$scaling)
lda_scaling_v$variable <- rownames(lda_scaling_v)
max_load_v  <- max(sqrt(lda_scaling_v$LD1^2 + lda_scaling_v$LD2^2))
vec_scale_v <- 2.5 / max_load_v

xlims_v <- pad_range(lda_data_v$LD1)
ylims_v <- pad_range(lda_data_v$LD2)

taxon_levels <- levels(my_data_v$Taxon)
cols_v <- setNames(base_cols[seq_along(taxon_levels)], taxon_levels)

# Main vocal LDA plot
p_main_v <- ggplot(lda_data_v, aes(x = LD1, y = LD2)) +
  stat_density_2d(aes(fill = ssp, color = ssp), geom = "polygon", alpha = 0.10, bins = 3) +
  stat_density_2d(aes(colour = ssp), bins = 6, linewidth = 0.4, alpha = 0.8, show.legend = FALSE) +
  geom_point(aes(fill = ssp, shape = ssp),
             colour = "white", stroke = 0.4, size = 2.8, alpha = 0.85) +
  geom_point(data = centroids_v, aes(x = LD1, y = LD2, fill = ssp),
             shape = 21, size = 5, color = "white", stroke = 0.8, show.legend = FALSE) +
  geom_label_repel(data = centroids_v, aes(x = LD1, y = LD2, label = ssp),
                   size = 3, label.size = 0, seed = 1, min.segment.length = 0,
                   fill = "white", alpha = 0.9, box.padding = 0.25,
                   segment.color = "grey50", show.legend = FALSE) +
  geom_segment(data = lda_scaling_v,
               aes(x = 0, y = 0, xend = LD1 * vec_scale_v, yend = LD2 * vec_scale_v),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "blue", linewidth = 0.8, inherit.aes = FALSE) +
  geom_text_repel(data = lda_scaling_v,
                  aes(x = LD1 * vec_scale_v, y = LD2 * vec_scale_v, label = variable),
                  color = "blue", fontface = "bold", size = 3.5,
                  max.overlaps = Inf, inherit.aes = FALSE) +
  scale_fill_manual(values = cols_v, breaks = names(cols_v), name = "Taxon") +
  scale_colour_manual(values = cols_v, breaks = names(cols_v), guide = "none") +
  scale_shape_manual(values = c(21, 22, 24, 25, 23)[seq_along(taxon_levels)],
                     name = "Taxon") +
  scale_x_continuous(limits = xlims_v) +
  scale_y_continuous(limits = ylims_v) +
  labs(x = "LDA 1", y = "LDA 2") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 10) +
  theme(
    legend.position  = "right",
    legend.title     = element_blank(),
    legend.text      = element_text(face = "italic"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, colour = "grey92"),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = "grey70"),
    axis.title       = element_text(colour = "grey20")
  )

p_top_v <- ggplot(lda_data_v, aes(x = LD1, fill = ssp)) +
  geom_density(alpha = 0.3) + scale_fill_manual(values = cols_v) +
  scale_x_continuous(limits = xlims_v, expand = c(0,0)) +
  theme_void() + theme(legend.position = "none")

p_right_v <- ggplot(lda_data_v, aes(x = LD2, fill = ssp)) +
  geom_density(alpha = 0.3) + coord_flip() +
  scale_fill_manual(values = cols_v) +
  scale_x_continuous(limits = ylims_v, expand = c(0,0)) +
  theme_void() + theme(legend.position = "none")

legend_v           <- get_legend(p_main_v)
p_main_v_no_legend <- p_main_v + theme(legend.position = "none")
placeholder_v      <- ggplot() + theme_void()

bottom_row_v <- plot_grid(p_main_v_no_legend, p_right_v, ncol = 2,
                          align = 'h', axis = 'tb', rel_widths = c(4, 1))
top_row_v    <- plot_grid(p_top_v, placeholder_v, ncol = 2, rel_widths = c(4, 1))

final_plot_vocal <- plot_grid(
  plot_grid(top_row_v, bottom_row_v, nrow = 2, rel_heights = c(1, 4)),
  legend_v, ncol = 2, rel_widths = c(1, 0.15)
)
final_plot_vocal

# Save
output_file_vocal <- "/Users/en_gl/Documents/Unesp/Mestrado/B.cinerea/Plots/lda_vocal_plot_hd.png"
save_plot(output_file_vocal, final_plot_vocal,
          base_height = 6, base_width = 8, dpi = 600)

# ============================================================
# LDA DIAGNOSABILITY & CONFUSION MATRICES
# ============================================================

# --- 1. Morphological Data ---
# Extract posterior probabilities and predicted classes
lda_pred_morpho <- predict(lda_result)

# Build the confusion matrix
conf_matrix_morpho <- table(True = my_data$ssp, Predicted = lda_pred_morpho$class)

# Calculate overall misclassification rate
misclass_rate_morpho <- sum(my_data$ssp != lda_pred_morpho$class) / nrow(my_data)

cat("\n--- MORPHOLOGICAL LDA RESULTS ---\n")
print(conf_matrix_morpho)
cat("Overall Misclassification Rate:", round(misclass_rate_morpho * 100, 1), "%\n")


# --- 2. Vocal Data ---
# Extract posterior probabilities and predicted classes
lda_pred_vocal <- predict(lda_result_v)

# Build the confusion matrix
conf_matrix_vocal <- table(True = my_data_v$Taxon, Predicted = lda_pred_vocal$class)

# Calculate overall misclassification rate
misclass_rate_vocal <- sum(my_data_v$Taxon != lda_pred_vocal$class) / nrow(my_data_v)

cat("\n--- VOCAL LDA RESULTS ---\n")
print(conf_matrix_vocal)
cat("Overall Misclassification Rate:", round(misclass_rate_vocal * 100, 1), "%\n")
