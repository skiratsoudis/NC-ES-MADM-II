# Figure 14 — Training effect on ADI and NMGI across scenarios
# No file I/O. Produces the figure in the active plotting device.

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
})

df <- data.frame(
  Scenario = c("S0","S0","S1","S1","S2","S2","EXTRA","EXTRA"),
  Phase    = c("Base","Trained","Base","Trained","Base","Trained","Base","Trained"),
  ADI      = c(0.089163,0.102706, 0.078801,0.092041, 0.069755,0.080178, 0.061057,0.056760),
  NMGI     = c(0.035508,0.035960, 0.041477,0.035319, 0.042949,0.040465, 0.046103,0.047258)
)

df_long <- df %>%
  pivot_longer(cols = c("ADI","NMGI"), names_to = "Metric", values_to = "Value") %>%
  mutate(
    Phase = factor(Phase, levels = c("Base","Trained")),
    Scenario = factor(Scenario, levels = c("S0","S1","S2","EXTRA")),
    Metric = factor(Metric, levels = c("ADI","NMGI"))
  )

p <- ggplot(df_long, aes(x = Phase, y = Value, group = Scenario)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.4) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  labs(
    x = NULL,
    y = NULL,
    title = "Figure 14. Training effect on ADI and NMGI across scenarios",
    subtitle = "S0–S2: empirical TradeHub-2018 regimes; EXTRA: mechanism demo under supervision"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

print(p)
