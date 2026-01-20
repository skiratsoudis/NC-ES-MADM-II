# Figures generated: Figures 4.1-4.6 (as referenced in the manuscript)
# No file I/O — figures rendered in active graphics device


suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)  # for plot stacking (plot / legend)
})

# --------------------
# Embedded results data
# --------------------
csv_data <- list(
  S0 = list(
    Weights = "criterion,x_OBJ,x_SBJ_base,x_INT_base,x_SBJ_trained,x_INT_trained\nLPI_Overall_2018,0.0260249495923246,0.213527182002479,0.0334464826708094,0.207660428931826,0.029900613916119\nTrade_Openness_%GDP_2018,0.131238845867537,0.153343840215721,0.121124105404443,0.160717534330803,0.116700019748099\nInternet_Users_%Pop_2018,0.140721847225536,0.153343840215721,0.102033621718552,0.169146477721765,0.103866764458129\nGDPpc_PPP_const2021$_2018,0.363042591208651,0.171543683634539,0.486967614911565,0.156102280190371,0.488727292739085\nInflation_CPI_%_2018,0.102093017222221,0.121708846935467,0.054019065666313,0.123334693406271,0.0554999084382001\nTariff_WM_AllProducts_%_MRV,0.23687874888373,0.186532606996073,0.202409109628317,0.183038585418965,0.205305400700369\n",
    Alternatives = "alternative,P_Y_base,P_Y_trained\nAfghanistan,0.0524703318354977,0.0440646312867626\nAngola,0.0413098072796611,0.0398937230035963\nArgentina,0.0822702231414602,0.0838829692752737\nAustralia,0.175451085341727,0.183204069173761\nAustria,0.190679625192888,0.200921039879047\nBangladesh,0.0323964780150571,0.026289967508393\nBelarus,0.109115764739274,0.117447145812153\nBelgium,0.195388741009365,0.20835560597423\nBenin,0.0348680511056126,0.0303063701030648\nBolivia,0.0860508923394585,0.065634477983718\n",
    System_Base = "S_Y,S_Y_given_X,I_Y_given_X,S_X,S_XY,J_XY,NMI,CES,CSF,ADI,NMGI\n3.02573463224701,2.89079468040967,0.134939951837338,2.23387570690442,5.25961033915143,0.134939951837338,0.0513122610847847,0.0003416500230804,0.044596956297486,0.0891631167743238,0.0355078050951722\n",
    System_Trained = "S_Y,S_Y_given_X,I_Y_given_X,S_X,S_XY,J_XY,NMI,CES,CSF,ADI,NMGI,DMI\n2.9807454155678,2.85390591434741,0.126839501220394,2.05769462838988,5.03844004395768,0.126839501220394,0.0503493317393895,-0.0147462026600183,0.042552978948,0.102705624785241,0.0359602196942132,0.940740740740741\n",
    Elasticities_Final = "Alternative,LPI_Overall_2018,Trade_Openness_%GDP_2018,Internet_Users_%Pop_2018,GDPpc_PPP_const2021$_2018,Inflation_CPI_%_2018,Tariff_WM_AllProducts_%_MRV\nAfghanistan,0.044458,0.183508,0.07278,0.10311,0.158258,0.437886\nAngola,0.064443,0.303725,0.082624,0.085534,0.115855,0.347819\nArgentina,0.053281,0.26426,0.093049,0.339553,0.091878,0.158978\nAustralia,0.014505,0.101445,0.168686,0.5627,0.039971,0.112692\nAustria,0.01519,0.107623,0.157812,0.569833,0.044087,0.105454\nBangladesh,0.11452,0.469965,0.060798,0.021675,0.074216,0.258826\nBelarus,0.052985,0.274227,0.115194,0.28294,0.064955,0.2097\nBelgium,0.014434,0.105267,0.16138,0.552588,0.042013,0.124318\nBenin,0.10791,0.44442,0.052907,0.025608,0.087035,0.28212\nBolivia,0.054458,0.271474,0.12267,0.248631,0.067477,0.235291\n",
    dP_dxSBJ_Final = "Alternative,dP/dLPI_Overall_2018,dP/dTrade_Openness_%GDP_2018,dP/dInternet_Users_%Pop_2018,dP/dGDPpc_PPP_const2021$_2018,dP/dInflation_CPI_%_2018,dP/dTariff_WM_AllProducts_%_MRV\nAfghanistan,0.003089,0.018317,-0.010268,-0.069837,0.046085,0.065424\nAngola,0.004186,0.033612,0.00042,-0.022708,0.008425,-0.010405\nArgentina,0.005027,0.060829,0.025543,0.122361,-0.037774,-0.099838\nAustralia,0.001117,0.009048,0.01653,0.10252,-0.0093,-0.02068\nAustria,0.00091,0.00746,0.012975,0.091345,-0.008803,-0.02038\nBangladesh,-0.00057,-0.00412,-0.01087,-0.01458,0.00619,0.02395\nBelarus,0.004832,0.06606,0.023561,0.053292,-0.01814,-0.04782\nBelgium,0.00086,0.00692,0.01201,0.08587,-0.00878,-0.01998\nBenin,-0.00032,-0.00266,-0.00788,-0.01218,0.00624,0.01699\nBolivia,0.002788,0.058318,0.022682,0.043114,-0.016286,-0.040946\n"
  ),
  S1 = list(
    Weights = "criterion,x_OBJ,x_SBJ_base,x_INT_base,x_SBJ_trained,x_INT_trained\nLPI_Overall_2018,0.0260249495923246,0.282762824180078,0.0506162676610129,0.273895399864027,0.0449221399893874\nTrade_Openness_%GDP_2018,0.131238845867537,0.20639497196336,0.18630894901052,0.218594414239789,0.180795944429309\nInternet_Users_%Pop_2018,0.140721847225536,0.121838975673184,0.0907163901938774,0.135659127851689,0.0926948118082019\nGDPpc_PPP_const2021$_2018,0.363042591208651,0.121838975673184,0.429792403195361,0.10300306633306,0.425606048305088\nInflation_CPI_%_2018,0.102093017222221,0.101803459373975,0.0508596932889824,0.0998517845471592,0.0514462850324135\nTariff_WM_AllProducts_%_MRV,0.23687874888373,0.165360793136219,0.191706296650246,0.169,0.2045347704356\n",
    Alternatives = "alternative,P_Y_base,P_Y_trained\nAfghanistan,0.056998625153174,0.0479767441860465\nAngola,0.0455406018052335,0.0439509101703295\nArgentina,0.0788050383426513,0.080630779392338\nAustralia,0.164325175193267,0.172688024711696\nAustria,0.182411434677097,0.193363844753663\nBangladesh,0.0358246158159851,0.0284685671916526\nBelarus,0.11265158869991,0.12314309113739\nBelgium,0.198902931216161,0.207903608864387\nBenin,0.0385414557608166,0.0343490304709141\nBolivia,0.085998533336704,0.0675233981225831\n",
    System_Base = "S_Y,S_Y_given_X,I_Y_given_X,S_X,S_XY,J_XY,NMI,CES,CSF,ADI,NMGI\n3.06015834865638,2.92257233995171,0.137586008704668,2.33736124665627,5.39751959531265,0.137586008704668,0.0509810828314853,0.0115869476310556,0.0449596162106908,0.0788005455728754,0.0414767618159557\n",
    System_Trained = "S_Y,S_Y_given_X,I_Y_given_X,S_X,S_XY,J_XY,NMI,CES,CSF,ADI,NMGI,DMI\n3.01617343119434,2.88294143200228,0.133231999192057,2.17833823461811,5.19451166581245,0.133231999192057,0.0512969914960764,-0.0028272266326824,0.0441717884257016,0.0920413315276462,0.035318625044901,0.866666666666667\n",
    Elasticities_Final = "Alternative,LPI_Overall_2018,Trade_Openness_%GDP_2018,Internet_Users_%Pop_2018,GDPpc_PPP_const2021$_2018,Inflation_CPI_%_2018,Tariff_WM_AllProducts_%_MRV\nAfghanistan,0.061346,0.261116,0.059656,0.082471,0.134738,0.400672\nAngola,0.081388,0.356929,0.076192,0.07357,0.091693,0.320228\nArgentina,0.075917,0.342561,0.084238,0.278293,0.076132,0.142859\nAustralia,0.024085,0.160668,0.155146,0.523923,0.036221,0.099957\nAustria,0.025687,0.175354,0.147379,0.531803,0.039993,0.079784\nBangladesh,0.137966,0.552898,0.047332,0.019022,0.061237,0.182545\nBelarus,0.076166,0.352728,0.101633,0.239574,0.055411,0.174488\nBelgium,0.024414,0.174773,0.14862,0.513814,0.040554,0.097825\nBenin,0.127571,0.527017,0.041302,0.021547,0.071307,0.211256\nBolivia,0.075264,0.345289,0.111205,0.213174,0.057541,0.197528\n",
    dP_dxSBJ_Final = "Alternative,dP/dLPI_Overall_2018,dP/dTrade_Openness_%GDP_2018,dP/dInternet_Users_%Pop_2018,dP/dGDPpc_PPP_const2021$_2018,dP/dInflation_CPI_%_2018,dP/dTariff_WM_AllProducts_%_MRV\nAfghanistan,0.002877,0.017629,-0.013269,-0.079031,0.049977,0.068682\nAngola,0.004103,0.034931,-0.000675,-0.029057,0.006986,-0.017909\nArgentina,0.004955,0.058696,0.022917,0.103517,-0.032408,-0.089484\nAustralia,0.001098,0.009135,0.015693,0.089907,-0.008707,-0.018148\nAustria,0.000886,0.007325,0.012101,0.08017,-0.008303,-0.017904\nBangladesh,-0.00057,-0.004105,-0.011207,-0.016053,0.006507,0.025428\nBelarus,0.004675,0.063629,0.022013,0.046487,-0.015777,-0.042809\nBelgium,0.000819,0.006744,0.011047,0.075296,-0.007737,-0.016488\nBenin,-0.00031,-0.002542,-0.008178,-0.013267,0.006716,0.017763\nBolivia,0.00269,0.055258,0.020676,0.041079,-0.016253,-0.038284\n"
  ),
  S2 = list(
    Weights = "criterion,x_OBJ,x_SBJ_base,x_INT_base,x_SBJ_trained,x_INT_trained\nLPI_Overall_2018,0.0260249495923246,0.361807676084703,0.0746143076716728,0.348250675548829,0.0663860769310841\nTrade_Openness_%GDP_2018,0.131238845867537,0.261657606990498,0.272110721264777,0.281670296918306,0.270769628411272\nInternet_Users_%Pop_2018,0.140721847225536,0.0738983294311263,0.0635937503005046,0.0909337180517308,0.064342058586058\nGDPpc_PPP_const2021$_2018,0.363042591208651,0.0709623138182195,0.314802044429682,0.0604827181556787,0.323529677243483\nInflation_CPI_%_2018,0.102093017222221,0.0637899083737141,0.0364726790493884,0.062414410520141,0.0361633868028022\nTariff_WM_AllProducts_%_MRV,0.23687874888373,0.167883865302,0.238406497284,0.156248180805315,0.238809172025301\n",
    Alternatives = "alternative,P_Y_base,P_Y_trained\nAfghanistan,0.064287385010605,0.0550935550935551\nAngola,0.0506639313523396,0.0493793793793794\nArgentina,0.0740614623387504,0.0757717717717718\nAustralia,0.150625264732633,0.158118118118118\nAustria,0.172058820414513,0.182433933933934\nBangladesh,0.0389030073636341,0.0297922922922923\nBelarus,0.120816123715006,0.13271046046046\nBelgium,0.195865394972681,0.207039039039039\nBenin,0.0437948330941491,0.0384894894894895\nBolivia,0.0889237760056892,0.071172972972973\n",
    System_Base = "S_Y,S_Y_given_X,I_Y_given_X,S_X,S_XY,J_XY,NMI,CES,CSF,ADI,NMGI\n3.09020769137778,2.95941194309974,0.130795748278044,2.30073208347918,5.39093977485696,0.130795748278044,0.0485235891248768,0.0211982064874585,0.0423258068729506,0.0697547920637451,0.0429492342629792\n",
    System_Trained = "S_Y,S_Y_given_X,I_Y_given_X,S_X,S_XY,J_XY,NMI,CES,CSF,ADI,NMGI,DMI\n3.05558272112216,2.92239354606592,0.133189175056246,2.21811231857205,5.27369503969421,0.133189175056246,0.0505112094704705,0.0101073619643453,0.0435890023099866,0.0801775988182533,0.0404651507540377,0.866666666666667\n",
    Elasticities_Final = "Alternative,LPI_Overall_2018,Trade_Openness_%GDP_2018,Internet_Users_%Pop_2018,GDPpc_PPP_const2021$_2018,Inflation_CPI_%_2018,Tariff_WM_AllProducts_%_MRV\nAfghanistan,0.078946,0.340544,0.03606,0.054593,0.082475,0.407383\nAngola,0.096202,0.411959,0.042032,0.054158,0.062081,0.333568\nArgentina,0.101333,0.431079,0.064042,0.220756,0.053384,0.129406\nAustralia,0.034428,0.211075,0.115934,0.469238,0.032663,0.136662\nAustria,0.03664,0.22575,0.110081,0.476214,0.035503,0.115812\nBangladesh,0.171901,0.620006,0.029023,0.014635,0.040074,0.124361\nBelarus,0.106147,0.43646,0.08025,0.184202,0.043659,0.149282\nBelgium,0.035024,0.218771,0.110897,0.459343,0.035609,0.140356\nBenin,0.162584,0.592458,0.025645,0.01657,0.049368,0.153375\nBolivia,0.102191,0.431003,0.092198,0.157751,0.043316,0.17354\n",
    dP_dxSBJ_Final = "Alternative,dP/dLPI_Overall_2018,dP/dTrade_Openness_%GDP_2018,dP/dInternet_Users_%Pop_2018,dP/dGDPpc_PPP_const2021$_2018,dP/dInflation_CPI_%_2018,dP/dTariff_WM_AllProducts_%_MRV\nAfghanistan,0.001987,0.013648,-0.024963,-0.121786,0.052764,0.067478\nAngola,0.003734,0.035381,-0.00922,-0.055284,0.004061,-0.030232\nArgentina,0.004767,0.055834,0.014167,0.061274,-0.02097,-0.069456\nAustralia,0.001194,0.010182,0.012927,0.055975,-0.006152,-0.012749\nAustria,0.000946,0.007972,0.010076,0.05047,-0.006028,-0.011229\nBangladesh,-0.00048,-0.003633,-0.011262,-0.018846,0.006152,0.028648\nBelarus,0.004423,0.060576,0.016432,0.03262,-0.011661,-0.03306\nBelgium,0.000833,0.006876,0.008869,0.047786,-0.005329,-0.010771\nBenin,-0.00025,-0.002099,-0.007962,-0.015603,0.006642,0.019934\nBolivia,0.002655,0.053971,0.018604,0.025432,-0.010182,-0.022732\n"
  )
)

read_embedded <- function(txt) read.csv(text = txt, stringsAsFactors = FALSE, check.names = FALSE)
S <- lapply(csv_data, function(x) lapply(x, read_embedded))

# -------------------------
# X1..X6 mapping (legend)
# -------------------------
crit_map <- tibble(
  CriterionFull = c(
    "LPI_Overall_2018",
    "Trade_Openness_%GDP_2018",
    "Internet_Users_%Pop_2018",
    "GDPpc_PPP_const2021$_2018",
    "Inflation_CPI_%_2018",
    "Tariff_WM_AllProducts_%_MRV"
  ),
  X = paste0("X", 1:6),
  Description = c(
    "LPI overall score (2018)",
    "Trade openness: Trade (% of GDP) (2018)",
    "Internet users (% of population) (2018)",
    "GDP per capita, PPP (constant 2021 intl$) (2018)",
    "Inflation, consumer prices (annual %) (2018)",
    "Applied tariff rate, weighted mean (all products, %) (2018)"
  )
)

# Legend plot WITHOUT row_number()/n() inside aes()
make_x_legend_plot <- function(title = "Criterion legend (X-codes)") {
  leg_df <- crit_map %>%
    transmute(Line = paste0(X, " = ", Description))
  # Force order top-to-bottom by factor levels (no data-masking functions)
  leg_df$Line <- factor(leg_df$Line, levels = rev(leg_df$Line))
  
  ggplot(leg_df, aes(x = 0, y = Line)) +
    geom_text(aes(label = Line), hjust = 0, size = 3.6) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0),
      plot.margin = margin(2, 5, 2, 5)
    )
}

# Rename helper (X = old_full_name)
rename_to_X <- function(df, first_col = "Alternative") {
  names(df)[1] <- first_col
  full   <- crit_map$CriterionFull
  xnames <- crit_map$X
  
  if (all(full %in% names(df))) {
    df <- df %>% rename(!!!setNames(full, xnames))  # new=X, old=full
  } else {
    # fallback: just enforce X1..X6 on columns 2..7
    names(df)[2:7] <- xnames
  }
  df
}

theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

# --------------------
# Data preparation
# --------------------
alt_long <- bind_rows(lapply(names(S), function(sc){
  df <- S[[sc]]$Alternatives
  df$Scenario <- sc
  df
}))

w_long <- bind_rows(lapply(names(S), function(sc){
  df <- S[[sc]]$Weights
  df$Scenario <- sc
  df
})) %>%
  left_join(crit_map, by = c("criterion" = "CriterionFull")) %>%
  mutate(CriterionX = factor(X, levels = crit_map$X))

sys_long <- bind_rows(lapply(names(S), function(sc){
  b <- S[[sc]]$System_Base;    b$Regime <- "Base"
  t <- S[[sc]]$System_Trained; t$Regime <- "Trained"
  b$Scenario <- sc; t$Scenario <- sc
  bind_rows(b, t)
}))

diag_vars <- c("NMI","CES","CSF","ADI","NMGI")
sys_diag_long <- sys_long %>%
  select(Scenario, Regime, all_of(diag_vars)) %>%
  pivot_longer(cols = all_of(diag_vars), names_to = "Index", values_to = "Value")

elas_S0_long <- rename_to_X(S[["S0"]]$Elasticities_Final, "Alternative") %>%
  pivot_longer(cols = starts_with("X"), names_to = "X", values_to = "Elasticity") %>%
  mutate(X = factor(X, levels = crit_map$X))

sens_S2_long <- rename_to_X(S[["S2"]]$dP_dxSBJ_Final, "Alternative") %>%
  pivot_longer(cols = starts_with("X"), names_to = "X", values_to = "dPdx") %>%
  mutate(X = factor(X, levels = crit_map$X))

glob_sens <- bind_rows(lapply(names(S), function(sc){
  df <- rename_to_X(S[[sc]]$dP_dxSBJ_Final, "Alternative")
  df %>%
    select(-Alternative) %>%
    summarise(across(everything(), ~sum(abs(.), na.rm = TRUE))) %>%
    pivot_longer(cols = everything(), names_to = "X", values_to = "Impact") %>%
    mutate(Scenario = sc, X = factor(X, levels = crit_map$X))
}))

# --------------------
# Figures (printed; no saving)
# --------------------

# CS1
p1 <- ggplot(alt_long, aes(x = reorder(alternative, P_Y_trained), y = P_Y_trained)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~Scenario, ncol = 3) +
  labs(title = "Figure CS1. Alternative-level probability mass P(Y) (trained)",
       x = NULL, y = "P(Y)") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_pub
print(p1 / make_x_legend_plot() + plot_layout(heights = c(4, 1.2)))

# CS2
p2 <- ggplot(w_long, aes(x = CriterionX, y = x_INT_trained)) +
  geom_col() +
  facet_wrap(~Scenario, ncol = 3) +
  labs(title = "Figure CS2. Integrated weights x_INT (trained)",
       x = NULL, y = "Weight") +
  theme_pub
print(p2 / make_x_legend_plot() + plot_layout(heights = c(4, 1.2)))

# CS3
p3 <- ggplot(sys_diag_long, aes(x = Index, y = Value, fill = Regime)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~Scenario, ncol = 3) +
  labs(title = "Figure CS3. System diagnostics (Base vs Trained)",
       x = NULL, y = "Value") +
  theme_pub
print(p3 / make_x_legend_plot() + plot_layout(heights = c(4, 1.2)))

# CS4
p4 <- ggplot(elas_S0_long, aes(x = X, y = Alternative, fill = Elasticity)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.3f", Elasticity)), size = 3,color="white") +
  labs(title = "Figure CS4. Elasticity allocation heatmap (S0, final)",
       x = NULL, y = NULL) +
  theme_pub
print(p4 / make_x_legend_plot() + plot_layout(heights = c(4.5, 1.4)))

# CS5
p5 <- ggplot(sens_S2_long, aes(x = X, y = Alternative, fill = dPdx)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.3f", dPdx)), size = 3) +
  labs(title = "Figure CS5. Sensitivity heatmap dP/dx_SBJ (S2, final)",
       x = NULL, y = NULL) +
  scale_fill_gradient2(labels = number_format(accuracy = 0.001)) +
  theme_pub
print(p5 / make_x_legend_plot() + plot_layout(heights = c(4.5, 1.4)))

# CS6
p6 <- ggplot(glob_sens, aes(x = X, y = Impact)) +
  geom_col() +
  facet_wrap(~Scenario, ncol = 3) +
  labs(title = "Figure CS6. Global sensitivity by criterion (sum |dP/dx_SBJ|)",
       x = NULL, y = "Impact") +
  theme_pub
print(p6 / make_x_legend_plot() + plot_layout(heights = c(4, 1.2)))

