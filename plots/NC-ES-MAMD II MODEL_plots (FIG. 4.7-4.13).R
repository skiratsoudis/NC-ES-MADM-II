# ==========================================================
# NES-MADM II — EXTRA Case Study (Figures 7–13)
# Standalone plotting script (NO file I/O, NO ggsave)
# ==========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

# -------------------------
# Criteria legend (X1..X6)
# -------------------------
crit_map <- tibble::tibble(
  criterion = paste0("X", 1:6),
  description = c(
    "LPI overall score (2018) — Benefit",
    "Trade openness: Trade (% of GDP) (2018) — Benefit",
    "Internet users (% of population) (2018) — Benefit",
    "GDP per capita, PPP (constant 2021 intl$) (2018) — Benefit",
    "Inflation, consumer prices (annual %) (2018) — Cost",
    "Tariff rate, applied, weighted mean, all products (%) (2018) — Cost"
  )
)

legend_plot <- function() {
  leg <- crit_map %>%
    mutate(label = paste0(criterion, " = ", description),
           y = rev(row_number()))
  ggplot(leg, aes(x = 0, y = y, label = label)) +
    geom_text(hjust = 0, size = 3.15) +
    xlim(0, 1) +
    theme_void() +
    theme(plot.margin = margin(0, 6, 0, 6))
}

attach_legend <- function(p) {
  p / legend_plot() + plot_layout(heights = c(1, 0.33))
}

# -------------------------
# Embedded data (from: NES-MADM II_ EXTRA CASE STUDY results.xlsx)
# -------------------------
weights <- tibble::tibble(
  criterion = c("X1","X2","X3","X4","X5","X6"),
  x_OBJ = c(0.069287251577547, 0.217044829118105, 0.0847127471163314,
            0.264286558671642, 0.167003466386275, 0.197665147130099),
  x_SBJ_base = c(0.102037099216594, 0.14635777908206, 0.118879303030948,
                 0.279308021577778, 0.179251144103471, 0.174166653),
  x_INT_base = c(0.0377894692028733, 0.169794244671432, 0.0538291515371841,
                 0.394562822652501, 0.160008993399678, 0.184015318536332),
  x_SBJ_trained = c(0.108788375045292, 0.22908554709659, 0.113531792824589,
                    0.196093246371197, 0.166764892924952, 0.18573614573738),
  x_INT_trained = c(0.0411294329649066, 0.271310743431061, 0.0524790293299528,
                    0.282785062218058, 0.151966514899204, 0.200329217156817)
) %>% mutate(delta_x_INT = x_INT_trained - x_INT_base)

alts <- tibble::tibble(
  alternative = c("Hubonia","Richovia","Balancedia","Openland",
                  "Digitopia","Tariffstan","Inflatia","Fragilia"),
  P_Y_base = c(0.178041029168492, 0.192238224971801, 0.149076101189709, 0.162494587654604,
               0.167166129590907, 0.0591106707730405, 0.0465658329514923, 0.0453074237009544),
  P_Y_trained = c(0.196543708161364, 0.163382717950519, 0.142485669397523, 0.179771594473729,
                  0.158011258357173, 0.0599823698740461, 0.0494694257025312, 0.0503532560831152),
  P_target = c(0.65, 0.03, 0.05, 0.18, 0.07, 0.015, 0.003, 0.002)
) %>%
  mutate(
    delta_P = P_Y_trained - P_Y_base,
    rank_base = rank(-P_Y_base, ties.method = "min"),
    rank_trained = rank(-P_Y_trained, ties.method = "min")
  )

elasticities <- tibble::tibble(
  Alternative = c("Hubonia","Richovia","Balancedia","Openland","Digitopia","Tariffstan","Inflatia","Fragilia"),
  X1=c(0.0472527157305562,0.0220410317042731,0.0305950842271099,0.0442811722525906,0.0299880126458259,0.063197,0.080461,0.067754),
  X2=c(0.404918819451715,0.110705589253177,0.177717488413812,0.402452389182611,0.171703779751547,0.301545,0.329071,0.287366),
  X3=c(0.0408367283312708,0.0598322286841575,0.0505516720130304,0.0372061402590813,0.0599131303745472,0.077198,0.104006,0.030653),
  X4=c(0.178068409821054,0.556945006289356,0.314399208395185,0.171318617242097,0.336668140042121,0.210050,0.226394,0.083406),
  X5=c(0.123201488090255,0.137985112804489,0.181661771680366,0.130049940238329,0.169098095343723,0.348010,0.000000,0.165823),
  X6=c(0.205721838575149,0.112490,0.245074775270497,0.214692,0.232628,0.000000,0.260068,0.364997)
)

sens <- tibble::tibble(
  Alternative = c("Hubonia","Richovia","Balancedia","Openland","Digitopia","Tariffstan","Inflatia","Fragilia"),
  X1 = c(0.0110632240025346,-0.0286674067481689,-0.0137983742133697,0.00520828648989914,-0.0161828576696213,0.013712,0.014331,0.014334),
  X2 = c(0.114628581049036,-0.114542326739225,-0.0582128007312173,0.102911149875395,-0.0687030342401225,0.033991,-0.005455,-0.006416),
  X3 = c(-0.0201553852697342,0.0105822230721985,-0.00241821095776752,-0.0241848213829534,0.0103459192927994,-0.004882,0.015224,0.015386),
  X4 = c(-0.104957659503593,0.228427038105416,0.0229724380192513,-0.102187918621758,0.0434193750077266,-0.043578,-0.018377,-0.025527),
  X5 = c(-0.033901991081304,-0.013697,0.025372,-0.023625,0.016233,0.010114,-0.004112,-0.001384),
  X6 = c(0.005707,-0.077268,0.034326,0.013901,0.027477,-0.009357,-0.001611,0.003825)
)

# -------------------------
# Figure 7: x_INT Base vs Trained
# -------------------------
w_long <- weights %>%
  select(criterion, Base = x_INT_base, Trained = x_INT_trained) %>%
  pivot_longer(cols = c(Base, Trained), names_to = "Phase", values_to = "Weight") %>%
  mutate(criterion = factor(criterion, levels = paste0("X",1:6)))

fig7 <- ggplot(w_long, aes(x = criterion, y = Weight, fill = Phase)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  scale_y_continuous(labels = number_format(accuracy = 0.001)) +
  labs(title = "Figure 7. Integrated weights x_INT (Base vs Trained) — EXTRA scenario",
       x = NULL, y = "Weight") +
  theme_minimal(base_size = 12)

print(attach_legend(fig7))

# -------------------------
# Figure 8: P(Y) Base vs Trained
# -------------------------
p_long <- alts %>%
  select(alternative, Base = P_Y_base, Trained = P_Y_trained) %>%
  pivot_longer(cols = c(Base, Trained), names_to = "Phase", values_to = "P") %>%
  mutate(alternative = factor(alternative, levels = alts %>% arrange(desc(P_Y_trained)) %>% pull(alternative)))

fig8 <- ggplot(p_long, aes(x = alternative, y = P, fill = Phase)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(title = "Figure 8. Alternative probabilities P(Y): Base vs Trained — EXTRA scenario",
       x = NULL, y = "Probability") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

print(attach_legend(fig8))

# -------------------------
# Figure 9: ΔP = Trained − Base
# -------------------------
fig9 <- ggplot(alts %>% mutate(alternative = factor(alternative, levels = alternative[order(delta_P)])),
               aes(x = alternative, y = delta_P)) +
  geom_hline(yintercept = 0) +
  geom_col(width = 0.70) +
  scale_y_continuous(labels = number_format(accuracy = 0.001)) +
  labs(title = "Figure 9. Probability redistribution ΔP(Y) = Trained − Base — EXTRA scenario",
       x = NULL, y = "ΔP") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

print(attach_legend(fig9))

# -------------------------
# Figure 10: Trained vs Target
# -------------------------
pt_long <- alts %>%
  select(alternative, Trained = P_Y_trained, Target = P_target) %>%
  pivot_longer(cols = c(Trained, Target), names_to = "Series", values_to = "P") %>%
  mutate(alternative = factor(alternative, levels = alts %>% arrange(desc(P_target)) %>% pull(alternative)))

fig10 <- ggplot(pt_long, aes(x = alternative, y = P, fill = Series)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(title = "Figure 10. Trained P(Y) vs supervision target P_target — EXTRA scenario",
       x = NULL, y = "Probability") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

print(attach_legend(fig10))

# -------------------------
# Figure 11: Elasticities heatmap
# -------------------------
elas_long <- elasticities %>%
  pivot_longer(cols = starts_with("X"), names_to = "criterion", values_to = "elasticity") %>%
  mutate(
    criterion = factor(criterion, levels = paste0("X",1:6)),
    Alternative = factor(Alternative, levels = rev(elasticities$Alternative))
  )

fig11 <- ggplot(elas_long, aes(x = criterion, y = Alternative, fill = elasticity)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.3f", elasticity)), size = 3, color="white") +
  labs(title = "Figure 11. Elasticities (final): Alternatives × Criteria — EXTRA scenario",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12)

print(attach_legend(fig11))

# -------------------------
# Figure 12: Sensitivity heatmap (dP/dx)
# -------------------------
sens_long <- sens %>%
  pivot_longer(cols = starts_with("X"), names_to = "criterion", values_to = "dpdx") %>%
  mutate(
    criterion = factor(criterion, levels = paste0("X",1:6)),
    Alternative = factor(Alternative, levels = rev(sens$Alternative))
  )

fig12 <- ggplot(sens_long, aes(x = criterion, y = Alternative, fill = dpdx)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.3f", dpdx)), size = 3, color="yellow") +
  labs(title = "Figure 12. Local sensitivities dP/dx: Alternatives × Criteria — EXTRA scenario",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12)

print(attach_legend(fig12))

# -------------------------
# Figure 13: Global sensitivity Σ|dP/dx|
# -------------------------
glob <- sens_long %>%
  group_by(criterion) %>%
  summarise(GlobalSensitivity = sum(abs(dpdx), na.rm = TRUE), .groups = "drop")

fig13 <- ggplot(glob, aes(x = criterion, y = GlobalSensitivity)) +
  geom_col(width = 0.70) +
  scale_y_continuous(labels = number_format(accuracy = 0.001)) +
  labs(title = "Figure 13. Global sensitivity: Σ_y |∂P(y)/∂x_μ| by criterion — EXTRA scenario",
       x = NULL, y = "Global sensitivity") +
  theme_minimal(base_size = 12)

print(attach_legend(fig13))
