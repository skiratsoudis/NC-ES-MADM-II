################################################################################
# app.R — NC-ES–MADM II II Pro (Office look, compact, full deliverables, cancelable compute)
## Model produced by LT COL Dr. Kiratsoudis Sideris
## All rights reserved by the author.
################################################################################

# ============ 0) Packages ============
req_pkgs <- c(
  "shiny","shinyWidgets","bslib","shinyjs","DT","ggplot2","plotly",
  "openxlsx","readxl","dplyr","tidyr","purrr","scales","shinycssloaders","tibble"
)
ensure_packages <- function(pkgs=req_pkgs){
  missing <- pkgs[!suppressWarnings(sapply(pkgs, requireNamespace, quietly=TRUE))]
  if(length(missing)) tryCatch(install.packages(missing, repos="https://cloud.r-project.org"), error=function(e) NULL)
  invisible(lapply(pkgs, function(p) suppressMessages(require(p, character.only=TRUE))))
}
ensure_packages()

suppressWarnings({
  library(shiny); library(shinyWidgets); library(bslib); library(shinyjs)
  library(DT); library(ggplot2); library(plotly)
  library(openxlsx); library(readxl)
  library(dplyr); library(tidyr); library(purrr); library(scales)
  library(shinycssloaders); library(tibble)
})

try({
  Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
}, silent=TRUE)

`%||%` <- function(a,b) if(!is.null(a)) a else b
to_num <- function(v){ v <- suppressWarnings(as.numeric(v)); v[!is.finite(v)] <- NA_real_; v }
parse_num <- function(x, default=NA_real_){
  if(is.numeric(x)) return(x)
  x2 <- suppressWarnings(as.numeric(x)); if(is.finite(x2)) return(x2)
  x3 <- suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
  if(is.finite(x3)) x3 else default
}

# ============ 0.A) App Metadata (Intro Section Only) ============
APP_VERSION <- "v1.0.0"
APP_RELEASE_DATE <- "2026-01-03"
APP_PRODUCER <- "Model produced by LT COL Dr. Kiratsoudis Sideris"
APP_FUNCTION <- paste0(
  "NC-ES–MADM II Pro platform for entropy–synergy decision support, ",
  "including objective/subjective/integrated weighting, probability structure, ",
  "indices (NMI, CES, CSF, ADI, NMGI), and full sensitivity diagnostics, ",
  "with cancelable live computation."
)

# ============ 1) Sample + Template ============
generate_sample_data <- function(M=6,N=8,K=3,seed=7){
  set.seed(seed)
  X <- data.frame(criterion=paste0("X",1:M),
                  type=c("Benefit","Cost","Benefit","Benefit","Cost","Benefit"))
  Y <- data.frame(alternative=paste0("Y",1:N))
  Xi <- matrix(round(runif(M*N,10,100),2), nrow=M, ncol=N)
  Xi <- data.frame(criterion=X$criterion, Xi, check.names=FALSE)
  names(Xi)[-1] <- Y$alternative
  W <- matrix(runif(M*K), nrow=M, ncol=K); W <- apply(W,2,function(z) z/sum(z))
  SBJ <- data.frame(criterion=X$criterion, W, check.names=FALSE); names(SBJ)[-1] <- paste0("DM",1:K)
  DM_Rho <- data.frame(DM=paste0("DM",1:K), rho=c(0.9,0.8,1.0)[1:K])
  rnk <- sample(1:N,N); tau <- 0.25; lg <- -rnk/tau; lg <- lg - max(lg)
  P_target <- exp(lg)/sum(exp(lg))
  WIN <- data.frame(alternative=Y$alternative, rank=rnk, P_target=round(P_target,6))
  list(X=X,Y=Y,Xi=Xi,SBJ=SBJ,DM_Rho=DM_Rho,WIN=WIN)
}
build_template_workbook <- function(d){
  wb <- createWorkbook()
  addWorksheet(wb,"X");        writeData(wb,"X",d$X);        freezePane(wb,"X", firstRow=TRUE)
  addWorksheet(wb,"Y");        writeData(wb,"Y",d$Y);        freezePane(wb,"Y", firstRow=TRUE)
  addWorksheet(wb,"Xi");       writeData(wb,"Xi",d$Xi);      freezePane(wb,"Xi", firstRow=TRUE)
  addWorksheet(wb,"SBJ");      writeData(wb,"SBJ",d$SBJ);    freezePane(wb,"SBJ", firstRow=TRUE)
  addWorksheet(wb,"DM_Rho");   writeData(wb,"DM_Rho",d$DM_Rho); freezePane(wb,"DM_Rho", firstRow=TRUE)
  addWorksheet(wb,"WIN");      writeData(wb,"WIN",d$WIN);    freezePane(wb,"WIN", firstRow=TRUE)
  wb
}

# ============ 2) Model ============
normalize_simplex <- function(x, eps=1e-15){ x[is.na(x)] <- 0; s <- sum(x); if(!is.finite(s)||s<=eps) rep(1/length(x),length(x)) else x/s }
linear_inversion_cost <- function(x_row, eps=1e-9){ m <- max(x_row, na.rm=TRUE); (m+eps) - x_row }
orient_and_normalize <- function(Xi, is_cost, eps=1e-9){
  Xi_p <- Xi
  for(i in seq_len(nrow(Xi))) if(is_cost[i]) Xi_p[i,] <- linear_inversion_cost(Xi[i,], eps)
  Xi_p[Xi_p<0] <- 0; rs <- rowSums(Xi_p)
  rho <- sweep(Xi_p, 1, ifelse(rs>0,rs,1), "/")
  rho[!is.finite(rho)] <- 0
  if(any(rowSums(rho)==0)) rho[rowSums(rho)==0,] <- 1/ncol(Xi)
  list(Xi_prime=Xi_p, rho=rho)
}
compute_entropy_div <- function(rho, log2=TRUE){
  h <- apply(rho,1,function(r){ pr <- pmax(r,1e-15); N <- length(pr)
  if(log2) -sum(pr*(log(pr)/log(2)))/(log(N)/log(2)) else -sum(pr*log(pr))/log(N) })
  list(h=h, d=1-h)
}
objective_weights <- function(d) normalize_simplex(d)
integrated_weights <- function(x_OBJ,x_SBJ) normalize_simplex(x_OBJ*x_SBJ)
fallback_consensus_from_panel <- function(W, rho_dm=NULL){
  M <- nrow(W); K <- ncol(W); if(is.null(rho_dm)||length(rho_dm)!=K) rho_dm <- rep(1,K)
  rho_dm <- pmax(rho_dm,1e-15); rho_dm <- rho_dm/sum(rho_dm)
  lg <- rep(0,M); for(k in 1:K) lg <- lg + rho_dm[k]*log(pmax(W[,k],1e-15))
  normalize_simplex(exp(lg))
}
neural_consensus_fast <- function(W, rho_dm=NULL, x_OBJ=NULL, lambda_ent=0, lambda_align=0){
  M <- nrow(W); K <- ncol(W); if(is.null(rho_dm)||length(rho_dm)!=K) rho_dm <- rep(1,K)
  rho_dm <- pmax(as.numeric(rho_dm),1e-12); rho_dm <- rho_dm/sum(rho_dm)
  logmix <- rep(0,M); for(k in 1:K) logmix <- logmix + rho_dm[k]*log(pmax(W[,k],1e-12))
  if(!is.null(x_OBJ) && lambda_align>0) logmix <- logmix + lambda_align*log(pmax(x_OBJ,1e-12))
  temp <- 1 + max(lambda_ent,0); normalize_simplex(exp(logmix/temp))
}
conditional_entropy_row <- function(r, log2=TRUE){ pr <- pmax(r,1e-15); if(log2) -sum(pr*(log(pr)/log(2))) else -sum(pr*log(pr)) }
overall_scores <- function(rho, x_INT) as.numeric(t(rho)%*%x_INT) |> normalize_simplex()
overall_entropy_Y <- function(PY, log2=TRUE){ pr <- pmax(PY,1e-15); if(log2) -sum(pr*(log(pr)/log(2))) else -sum(pr*log(pr)) }
total_cond_entropy <- function(S_rows,x_INT) sum(x_INT*S_rows)
normalized_total_cond_entropy <- function(SYx,SY) SYx/pmax(SY,1e-15)
joint_prob <- function(rho,x_INT) rho * x_INT
joint_entropy <- function(PXY, log2=TRUE){ pr <- pmax(as.numeric(PXY),1e-15); if(log2) -sum(pr*(log(pr)/log(2))) else -sum(pr*log(pr)) }
entropy_SX <- function(x_INT, log2=TRUE){ pr <- pmax(x_INT,1e-15); if(log2) -sum(pr*(log(pr)/log(2))) else -sum(pr*log(pr)) }
mutual_information <- function(SX,SY,SXY) SX + SY - SXY
NMI <- function(J, SX, SY) 2*J/pmax(SX+SY,1e-15)
CES <- function(S_rows,SY) mean((SY-S_rows)/pmax(SY,1e-15))
CSF <- function(SYx,SY) 1 - (SYx/pmax(SY,1e-15))
ADI <- function(SY,N,log2=TRUE){ denom <- if(log2) (log(N)/log(2)) else log(N); 1 - (SY/pmax(denom,1e-15)) }
NMGI <- function(NMI_val,CES_val,CSF_val,ADI_val){
  X <- pmax(c(NMI=NMI_val, CES=CES_val, CSF=CSF_val, ADI=ADI_val),1e-15)
  p <- X/sum(X); d <- 1 + ( - (1/log(4))*sum(p*log(p)) - 1 ) # 1 - entropy
  as.numeric(d*mean(X))
}
J_xint_wrt_xsbj <- function(x_OBJ,x_SBJ){
  z <- x_OBJ*x_SBJ; s <- sum(z); M <- length(z)
  J <- matrix(0,M,M)
  for(mu in 1:M) for(i in 1:M) J[mu,i] <- if(mu==i) x_OBJ[mu]*(s - z[mu])/(s^2) else - z[mu]*x_OBJ[i]/(s^2)
  J
}

# ============ 3) UI ============
light_theme <- bs_theme(
  version=5, bootswatch="flatly",
  base_font=font_google("Inter", TRUE),
  heading_font=font_google("Inter", TRUE),
  code_font=font_google("JetBrains Mono", TRUE),
  primary="#2563eb"
)
dark_theme <- bs_theme(
  version=5, bootswatch="cyborg",
  base_font=font_google("Inter", TRUE),
  heading_font=font_google("Inter", TRUE),
  code_font=font_google("JetBrains Mono", TRUE),
  primary="#1dd1a1"
)

tab_icon <- function(name, ic) tagList(icon(ic), HTML(paste("&nbsp;", name)))

ui <- navbarPage(
  title="NC-ES–MADM II",
  theme=light_theme, id="nav",
  header = tagList(
    useShinyjs(),
    tags$style(HTML("
      /* Office-like tightness & no mysterious top space */
      .navbar{font-weight:600;}
      .navbar + .container-fluid { padding-top: 0 !important; }
      .container-fluid { padding-top: 0 !important; }
      .tab-content { margin-top: 0 !important; }
      .bslib-value-box{min-width:160px;margin-bottom:8px}
      .tight-row{margin-left:-6px;margin-right:-6px}
      .tight-row > [class^='col-']{padding-left:6px;padding-right:6px;margin-bottom:12px;}
      .ribbon{background:var(--bs-body-bg);border-bottom:1px solid rgba(0,0,0,.08);
              padding:6px 12px;display:flex;gap:12px;align-items:center;}
      .status-dot{width:10px;height:10px;border-radius:50%;display:inline-block;margin-right:6px;background:#9CA3AF}
      .status-dot.running{background:#22c55e}.status-dot.stopped{background:#ef4444}

      /* Intro (PES-like) */
      @keyframes softGlowIntro {
        0%   { box-shadow:0 0 18px 4px rgba(255,255,0,0.35); }
        50%  { box-shadow:0 0 28px 8px rgba(255,255,0,0.70); }
        100% { box-shadow:0 0 18px 4px rgba(255,255,0,0.35); }
      }
      .intro-wrap{
        max-width: 1040px;
        margin: 18px auto 0 auto;
        padding: 18px 18px;
      }
      .intro-hero{
        background: rgba(0,0,0,0.05);
        border: 1px solid rgba(0,0,0,0.08);
        border-radius: 18px;
        padding: 18px 18px;
      }
      .intro-img{
        width: 180px; height: 180px;
        border-radius: 50%;
        overflow: hidden;
        border: 4px solid rgba(255,255,0,0.55);
        background: rgba(0,0,0,0.06);
        animation: softGlowIntro 3s ease-in-out infinite;
        margin: 0 auto;
      }
      .intro-pill{
        display:inline-flex;
        align-items:center;
        gap:18px;
        padding:10px 14px;
        border-radius: 14px;
        background: rgba(0,0,0,0.04);
        border: 1px solid rgba(0,0,0,0.08);
        font-weight: 700;
      }
    "))
  ),
  
  # ---------- NEW: Intro / Instructions (PES-MADM II style) ----------
  tabPanel(tab_icon("Instructions","info-circle"),
           div(class="intro-wrap",
               div(class="intro-hero",
                   div(class="intro-img",
                       tags$img(
                         src = "https://cdn.thecollector.com/wp-content/uploads/2024/07/thinker-auguste-rodin-what-so-special.jpg?width=1400&quality=70",
                         style = "width:100%; height:100%; object-fit:cover;"
                       )
                   ),
                   tags$div(style="text-align:center; margin-top:14px;",
                            tags$h2("NC-ES–MADM II Decision Support Platform",
                                    style="margin:0; font-weight:800; letter-spacing:0.3px;"),
                            tags$div(APP_PRODUCER,
                                     style="margin-top:6px; font-weight:800; color:#b8860b; font-size:16px;"),
                            tags$div(style="margin-top:12px;",
                                     tags$span(class="intro-pill",
                                               tags$span(paste0("Version: ", APP_VERSION)),
                                               tags$span(paste0("Release date: ", APP_RELEASE_DATE))
                                     )
                            ),
                            tags$div(style="margin:12px auto 0 auto; max-width:880px; font-size:13.5px; line-height:1.45;",
                                     APP_FUNCTION
                            ),
                            tags$hr(style="margin:14px 0 10px 0; opacity:0.25;"),
                            tags$div(style="text-align:left; max-width:920px; margin:0 auto;",
                                     tags$h5("Operational Steps (Quick Guide)"),
                                     tags$ul(
                                       tags$li("Prepare or download the provided template .xlsx (tab: Inputs → Download Template)."),
                                       tags$li("Populate X, Y, Xi, SBJ, DM_Rho (and optionally WIN) sheets exactly as specified."),
                                       tags$li("Upload the Excel file, set computation options, then select Compute / Refresh."),
                                       tags$li("Use Cancel (Stop) to interrupt long computations at any point."),
                                       tags$li("Export the full results workbook from Summary & Export.")
                                     ),
                                     tags$h5("Required Excel Sheets"),
                                     tags$ul(
                                       tags$li("X: criterion, type (Benefit/Cost)"),
                                       tags$li("Y: alternative"),
                                       tags$li("Xi: criterion + alternative columns (raw scores)"),
                                       tags$li("SBJ: criterion + DM columns (panel subjective weights per DM)"),
                                       tags$li("DM_Rho: DM, rho (DM reliabilities)"),
                                       tags$li("WIN (optional): alternative with rank and/or P_target")
                                     )
                            )
                   )
               )
           )
  ),
  
  # ---------- Inputs
  tabPanel(tab_icon("Inputs","folder-open"),
           sidebarLayout(
             sidebarPanel(width=4,
                          fileInput("file_xlsx","Upload Excel (.xlsx)", accept=".xlsx"),
                          helpText("Sheets required: X, Y, Xi, SBJ, DM_Rho. Optional WIN: rank and/or P_target."),
                          awesomeCheckbox("use_sample","Use built-in sample", value=TRUE),
                          switchInput("dark_mode","Dark mode", value=FALSE),
                          hr(),
                          textInput("epsilon", "ε (linear inversion safeguard)", "1e-9"),
                          switchInput("log2", "Entropies in log2", value=TRUE),
                          hr(),
                          switchInput("use_nn", "Neural consensus", value=TRUE),
                          numericInput("lambda_ent", "λ_ent (smoothing)", 0.0, min=0, step=0.01),
                          numericInput("lambda_align", "λ_align (align with x^OBJ)", 0.0, min=0, step=0.01),
                          hr(),
                          checkboxInput("use_win", "Use supervision if WIN present", value=TRUE),
                          numericInput("win_tau", "τ for softmax (rank→P_target)", value=0.20, min=0.01, step=0.01),
                          div(class="mb-2",
                              actionButton("run","Compute / Refresh", class="btn btn-primary"),
                              actionButton("cancel","Stop", class="btn btn-outline-danger"),
                              shinyWidgets::progressBar(id="comp_pb", value=0, total=100, display_pct=TRUE, size="sm")
                          ),
                          downloadButton("dl_template","Download Template (.xlsx)")
             ),
             mainPanel(
               h5("Preview"),
               tabsetPanel(
                 tabPanel("Criteria",     withSpinner(DTOutput("tbl_criteria"), type=4)),
                 tabPanel("Alternatives", withSpinner(DTOutput("tbl_alts"), type=4)),
                 tabPanel("Xi (Raw)",     withSpinner(DTOutput("tbl_xi"), type=4)),
                 tabPanel("Subjective Weights", withSpinner(DTOutput("tbl_wsbj"), type=4)),
                 tabPanel("DM Reliabilities",   withSpinner(DTOutput("tbl_rho"), type=4)),
                 tabPanel("WIN (optional)",     withSpinner(DTOutput("tbl_win"), type=4)),
                 tabPanel("Diagnostics",        tableOutput("tbl_checks"))
               )
             )
           )
  ),
  
  # ---------- Dashboard (compact, no gap)
  tabPanel(tab_icon("Dashboard","dashboard"),
           div(class="ribbon",
               span(id="status_dot", class="status-dot stopped"),
               strong(textOutput("run_status", inline=TRUE)),
               span(" — Live compute with Cancel")
           ),
           div(class="tight-row",
               fluidRow(
                 column(3, bslib::value_box(title="NMGI (Base)", value=textOutput("kpi_nmgi0"))),
                 column(3, bslib::value_box(title="NMGI (Trained)", value=textOutput("kpi_nmgi1"))),
                 column(3, bslib::value_box(title="CSF (Trained)", value=textOutput("kpi_csf1"))),
                 column(3, bslib::value_box(title="Panel DMI", value=textOutput("kpi_dmi")))
               ),
               fluidRow(
                 column(6, bslib::card(title="Alternative Probabilities — Baseline vs Trained",
                                       withSpinner(plotlyOutput("plt_py_compare", height="330px"), type=4))),
                 column(6, bslib::card(title="System Indices — Baseline vs Trained",
                                       withSpinner(plotlyOutput("plt_indices_compare", height="330px"), type=4)))
               ),
               fluidRow(
                 column(6, bslib::card(title="Integrated Weights — Baseline vs Trained",
                                       withSpinner(plotlyOutput("plt_intw_compare", height="330px"), type=4))),
                 column(6, bslib::card(title="Subjective Consensus — Baseline vs Trained",
                                       withSpinner(plotlyOutput("plt_sbj_compare", height="330px"), type=4)))
               ),
               fluidRow(
                 column(6, bslib::card(title="Change in Alternative Probabilities (Δ P(Y))",
                                       withSpinner(plotlyOutput("plt_delta_py", height="330px"), type=4))),
                 column(6, bslib::card(title="Top Influential Criteria (|∑ν ∂Pν/∂x^SBJμ|)",
                                       withSpinner(plotlyOutput("plt_top_influence", height="330px"), type=4)))
               ),
               fluidRow(
                 column(4,
                        bslib::card(
                          bslib::card_header("Radar Controls"),
                          selectInput("radar_alt","Single alternative", choices=NULL),
                          checkboxInput("radar_all","Show all alternatives (grid)", TRUE),
                          sliderInput("radar_cols","Grid columns", min=2,max=4,value=3,step=1)
                        )
                 ),
                 column(8, bslib::card(title="Criteria Profiles — Baseline vs Trained",
                                       withSpinner(plotlyOutput("plt_radar_wall", height="520px"), type=4)))
               )
           )
  ),
  
  # ---------- Neural Consensus
  tabPanel(tab_icon("Neural Consensus","project-diagram"),
           fluidRow(
             column(6, bslib::card(title="Final Subjective Consensus",
                                   withSpinner(plotlyOutput("plt_sbj_final", height="360px"), type=4))),
             column(6, bslib::card(title="Status / Engine", verbatimTextOutput("txt_nn_status")))
           ),
           bslib::card(withSpinner(DTOutput("tbl_sbj"), type=4))
  ),
  
  # ---------- Objective & Integrated
  tabPanel(tab_icon("Objective & Integrated","sliders-h"),
           fluidRow(
             column(6, bslib::card(title="Per-Criterion Entropy & Diversification",
                                   withSpinner(plotlyOutput("plt_entropy", height="320px"), type=4))),
             column(6, bslib::card(title="Objective vs Integrated Weights",
                                   withSpinner(plotlyOutput("plt_compare_w_all", height="320px"), type=4)))
           ),
           bslib::card(withSpinner(DTOutput("tbl_h_d_xobj"), type=4)), br(),
           bslib::card(withSpinner(DTOutput("tbl_intw"), type=4))
  ),
  
  # ---------- Probability Structure
  tabPanel(tab_icon("Probability Structure","project-diagram"),
           fluidRow(
             column(6, bslib::card(title="System Entropies / Info",
                                   tableOutput("tbl_sys_entropy"))),
             column(6, bslib::card(title="Cumulative Mass of Trained P(Y)",
                                   withSpinner(plotlyOutput("plt_py_cum", height="300px"), type=4)))
           ),
           bslib::card(title="Joint Probability Heatmap (Trained)",
                       withSpinner(plotlyOutput("plt_heat_pxy", height="520px"), type=4))
  ),
  
  # ---------- Sensitivity (FULL)
  tabPanel(tab_icon("Sensitivity","tachometer-alt"),
           fluidRow(
             column(6, bslib::card(title="Elasticity η (Alternative × Criterion)",
                                   withSpinner(DTOutput("tbl_elasticity"), type=4))),
             column(6, bslib::card(title="Jacobian ∂P/∂x^SBJ",
                                   withSpinner(DTOutput("tbl_dP_dxsbj"), type=4)))
           ),
           fluidRow(
             column(6,
                    bslib::card(
                      bslib::card_header("Worst-case ℓ1 budget"),
                      sliderInput("wc_eps", "ε for x^INT", min=0, max=0.5, value=0.05, step=0.01),
                      tableOutput("tbl_wc_bounds")
                    )
             ),
             column(6,
                    bslib::card(
                      bslib::card_header("Stochastic Uncertainty"),
                      sliderInput("sigma_int", "σ for x^INT (indep.)", min=0, max=0.2, value=0.05, step=0.005),
                      withSpinner(plotlyOutput("plt_uncertainty", height="330px"), type=4)
                    )
             )
           )
  ),
  
  # ---------- Supervision & Tuning
  tabPanel(tab_icon("Supervision & Tuning","wrench"),
           fluidRow(
             column(6, bslib::card(title="Fit to Target (if WIN provided)",
                                   withSpinner(plotlyOutput("plt_py_vs_ptarget", height="320px"), type=4))),
             column(6, bslib::card(title="WIN Metrics", tableOutput("tbl_win_metrics")))
           ),
           hr(),
           fluidRow(
             column(4, numericInput("grid_lambda_align_max", "Grid λ_align max", 2.0, min=0, step=0.1)),
             column(4, numericInput("grid_lambda_ent_max",   "Grid λ_ent max",   0.5, min=0, step=0.05)),
             column(4, actionButton("btn_grid", "Auto-tune λ (FAST)", class="btn btn-warning"))
           ),
           bslib::card(withSpinner(plotlyOutput("plt_grid", height="420px"), type=4))
  ),
  
  # ---------- Notation & Help
  tabPanel(tab_icon("Notation & Help","question-circle"),
           p("Model equations, symbol meanings, and code mapping."),
           downloadButton("dl_notation","Download Notation (.xlsx)"),
           br(), br(),
           bslib::card(withSpinner(DTOutput("tbl_notation"), type=4))
  ),
  
  # ---------- Summary & Export (merged & compact)
  tabPanel(tab_icon("Summary & Export","file-export"),
           fluidRow(
             column(8, bslib::card(title="Executive Summary", verbatimTextOutput("txt_overall"))),
             column(4, bslib::card(title="Download", downloadButton("dl_results","Download Results (.xlsx)")))
           ),
           bslib::card(withSpinner(DTOutput("tbl_insight"), type=4)),
           verbatimTextOutput("txt_summary")
  )
)

# ============ 4) Server ============
server <- function(input, output, session){
  
  observe({ session$setCurrentTheme(if (isTRUE(input$dark_mode)) dark_theme else light_theme) })
  
  cancel_flag <- reactiveVal(FALSE)
  running <- reactiveVal(FALSE)
  observeEvent(input$cancel, {
    cancel_flag(TRUE); running(FALSE); updateProgressBar(session,"comp_pb",value=0)
    shinyjs::html("run_status","Stopped"); shinyjs::runjs("$('#status_dot').removeClass('running').addClass('stopped');")
  })
  output$run_status <- renderText(if (isTRUE(running())) "Running…" else "Ready")
  
  output$dl_template <- downloadHandler(
    filename=function() "template_esmadm2.xlsx",
    content=function(file){ d <- generate_sample_data(); wb <- build_template_workbook(d); saveWorkbook(wb,file,overwrite=TRUE) }
  )
  
  data_input <- reactive({
    if (isTRUE(input$use_sample) || is.null(input$file_xlsx)) {
      generate_sample_data()
    } else {
      path <- input$file_xlsx$datapath
      out <- list(
        X = read_excel(path,"X"),
        Y = read_excel(path,"Y"),
        Xi = read_excel(path,"Xi"),
        SBJ = read_excel(path,"SBJ"),
        DM_Rho = read_excel(path,"DM_Rho")
      )
      out$WIN <- tryCatch(read_excel(path,"WIN"), error=function(e) NULL)
      out
    }
  })
  
  checks <- reactive({
    d <- data_input(); sbj <- d$SBJ; if ("criterion" %in% names(sbj)) sbj <- dplyr::select(sbj, -criterion)
    col_sums <- if (ncol(sbj)>0) colSums(sbj) else numeric(0)
    list(
      "SBJ: column sums" = paste(names(col_sums), sprintf("%.4f", col_sums)),
      "DM_Rho length matches #DMs" = length(d$DM_Rho$rho) == ncol(sbj),
      "WIN present" = !is.null(d$WIN),
      "Cost flags parsed" = all(!is.na(match(tolower(d$X$type), tolower(c("benefit","cost","c","κ","kostos","κοστος")))))
    )
  })
  output$tbl_checks    <- renderTable({ x <- checks(); data.frame(Check=names(x), Status=sapply(x, function(v) paste(v, collapse=" | ")), check.names=FALSE) }, striped=TRUE, bordered=TRUE, spacing="m")
  output$tbl_criteria  <- renderDT(datatable(data_input()$X,  options=list(pageLength=6, scrollX=TRUE)))
  output$tbl_alts      <- renderDT(datatable(data_input()$Y,  options=list(pageLength=6, scrollX=TRUE)))
  output$tbl_xi        <- renderDT(datatable(data_input()$Xi, options=list(pageLength=6, scrollX=TRUE)))
  output$tbl_wsbj      <- renderDT(datatable(data_input()$SBJ, options=list(pageLength=6, scrollX=TRUE)))
  output$tbl_rho       <- renderDT(datatable(data_input()$DM_Rho, options=list(pageLength=6, scrollX=TRUE)))
  output$tbl_win       <- renderDT({ d <- data_input(); if (is.null(d$WIN)) return(NULL); datatable(d$WIN, options=list(pageLength=10, scrollX=TRUE)) })
  output$txt_nn_status <- renderText({
    if (!isTRUE(input$use_nn)) return("Neural consensus: OFF — geometric pooling in effect.")
    "Neural engine: FAST (closed-form + supervised; no torch)."
  })
  
  # ---- Main compute (cancel-aware)
  results <- eventReactive(input$run, {
    cancel_flag(FALSE); running(TRUE)
    shinyjs::html("run_status","Running…")
    shinyjs::runjs("$('#status_dot').removeClass('stopped').addClass('running');")
    
    pb <- function(v){ updateProgressBar(session,"comp_pb", value=max(0,min(100,round(v)))) }
    
    withProgress(message = "Computing NC–ES-MADM II ...", value = 0, {
      on.exit({ running(FALSE); shinyjs::html("run_status","Ready"); shinyjs::runjs("$('#status_dot').removeClass('running').addClass('stopped');") }, add=TRUE)
      if (cancel_flag()) stop("Cancelled by user")
      
      d <- data_input(); incProgress(0.05); pb(6)
      criteria <- d$X$criterion; alts <- d$Y$alternative
      updateSelectInput(session,"radar_alt", choices=alts, selected = alts[1])
      
      Xi_df <- d$Xi; if ("criterion" %in% names(Xi_df)) Xi_df <- dplyr::select(Xi_df, -criterion)
      Xi <- as.matrix(Xi_df); if (!all(colnames(Xi) == alts)) Xi <- as.matrix(d$Xi %>% dplyr::select(any_of(c("criterion", alts))) %>% dplyr::select(-criterion))
      is_cost <- tolower(d$X$type) %in% c("cost","c","κ","kostos","κοστος")
      SBJ <- d$SBJ; if ("criterion" %in% names(SBJ)) SBJ <- dplyr::select(SBJ, -criterion)
      W_SBJ <- as.matrix(SBJ); rho_dm <- to_num(d$DM_Rho$rho)
      eps_val <- parse_num(input$epsilon, 1e-9)
      
      on <- orient_and_normalize(Xi, is_cost, eps=eps_val); incProgress(0.18); pb(22)
      Xi_prime <- on$Xi_prime; rho <- on$rho
      
      ent <- compute_entropy_div(rho, log2 = isTRUE(input$log2))
      h <- ent$h; dvec <- ent$d; x_OBJ <- objective_weights(dvec); incProgress(0.30); pb(35)
      if (cancel_flag()) stop("Cancelled by user")
      
      # helper
      metrics_for <- function(x_INT){
        S_rows <- apply(rho,1,conditional_entropy_row, log2=isTRUE(input$log2))
        P_Y <- overall_scores(rho, x_INT)
        S_Y <- overall_entropy_Y(P_Y, log2=isTRUE(input$log2))
        S_Y_given_X <- total_cond_entropy(S_rows, x_INT)
        I_Y_given_X <- normalized_total_cond_entropy(S_Y_given_X, S_Y)
        P_XY <- joint_prob(rho, x_INT)
        S_XY <- joint_entropy(P_XY, log2=isTRUE(input$log2))
        S_X <- entropy_SX(x_INT, log2=isTRUE(input$log2))
        J_XY <- mutual_information(S_X, S_Y, S_XY)
        nmi <- NMI(J_XY, S_X, S_Y); ces <- CES(S_rows,S_Y); csf <- CSF(S_Y_given_X,S_Y); adi <- ADI(S_Y, ncol(rho), log2=isTRUE(input$log2))
        nmgi <- NMGI(nmi, ces, csf, adi)
        list(P_Y=P_Y,S_Y=S_Y,S_Y_given_X=S_Y_given_X,I_Y_given_X=I_Y_given_X,
             P_XY=P_XY,S_XY=S_XY,S_X=S_X,J_XY=J_XY,NMI=nmi,CES=ces,CSF=csf,ADI=adi,NMGI=nmgi,S_rows=S_rows)
      }
      
      # Targets (if any)
      WIN <- d$WIN; P_target <- NULL
      if (!is.null(WIN) && isTRUE(input$use_win)) {
        WIN2 <- WIN %>% right_join(data.frame(alternative = alts), by = "alternative")
        if ("P_target" %in% names(WIN2)) {
          P_target <- normalize_simplex(to_num(WIN2$P_target))
        } else if ("rank" %in% names(WIN2)) {
          rnk <- to_num(WIN2$rank); rnk[is.na(rnk)] <- max(rnk, na.rm=TRUE)+5
          logits <- -rnk / max(0.01, input$win_tau); logits <- logits - max(logits, na.rm=TRUE)
          P_target <- exp(logits); P_target[!is.finite(P_target)] <- 0; P_target <- P_target/sum(P_target)
        }
      }
      
      # Baseline
      x_SBJ0 <- if (isTRUE(input$use_nn)) neural_consensus_fast(W_SBJ, rho_dm, x_OBJ,
                                                                lambda_ent = input$lambda_ent %||% 0, lambda_align = input$lambda_align %||% 0)
      else fallback_consensus_from_panel(W_SBJ, rho_dm=rho_dm)
      x_INT0 <- integrated_weights(x_OBJ, x_SBJ0)
      met0 <- metrics_for(x_INT0); incProgress(0.20); pb(55)
      if (cancel_flag()) stop("Cancelled by user")
      
      # Supervised nudged (if targets)
      x_SBJ1 <- tryCatch({
        if (!is.null(P_target) && isTRUE(input$use_nn)) {
          iters <- 260; lr <- 0.5; fd_eps <- 1e-4
          v <- log(pmax(x_SBJ0,1e-12))
          loss <- function(vv){
            x <- exp(vv - max(vv)); x <- x/sum(x)
            panel <- 0; for(k in seq_len(ncol(W_SBJ)))
              panel <- panel + (rho_dm[k]/sum(rho_dm))*sum(x*(log(pmax(x,1e-12)) - log(pmax(W_SBJ[,k],1e-12))))
            align <- if (input$lambda_align>0) input$lambda_align * sum(x*(log(pmax(x,1e-12))-log(pmax(x_OBJ,1e-12)))) else 0
            ent   <- if (input$lambda_ent>0)   input$lambda_ent   * sum(x*log(pmax(x,1e-12))) else 0
            x_INT <- normalize_simplex(x_OBJ*x); PY <- as.numeric(t(rho)%*%x_INT)
            ce <- -sum(P_target*log(pmax(PY,1e-12)))
            panel + align + ent + ce
          }
          for(t in seq_len(iters)){
            if (cancel_flag()) stop("Cancelled by user")
            g <- numeric(length(v))
            for (j in seq_along(v)) { v1<-v; v1[j]<-v1[j]+fd_eps; v2<-v; v2[j]<-v2[j]-fd_eps; g[j]<-(loss(v1)-loss(v2))/(2*fd_eps) }
            v <- pmax(pmin(v - lr*g, 30), -30)
          }
          x <- exp(v - max(v)); x/sum(x)
        } else x_SBJ0
      }, error=function(e) x_SBJ0)
      x_INT1 <- integrated_weights(x_OBJ, x_SBJ1)
      met1 <- metrics_for(x_INT1); pb(85)
      
      # Panel agreement
      P_Y_by_DM <- NULL; tau_by_DM <- NULL; DMI <- NA_real_
      if (!is.null(W_SBJ) && ncol(W_SBJ)>=1){
        K <- ncol(W_SBJ); N <- ncol(rho); P_Y_by_DM <- matrix(NA_real_, nrow=N, ncol=K)
        for(k in 1:K){
          x_int_k <- integrated_weights(x_OBJ, normalize_simplex(W_SBJ[,k]))
          P_Y_by_DM[,k] <- as.numeric(t(rho)%*%x_int_k)
        }
        rk_cons <- rank(-met1$P_Y, ties.method="average")
        tau_by_DM <- sapply(1:K, function(k){ rk <- rank(-P_Y_by_DM[,k], ties.method="average"); suppressWarnings(cor(rk,rk_cons,method="kendall")) })
        DMI <- mean(tau_by_DM, na.rm=TRUE)
      }
      
      # Sensitivities (final)
      dP_dxint <- t(rho)
      eta1 <- matrix(NA_real_, nrow=length(met1$P_Y), ncol=length(x_INT1))
      for(nu in seq_along(met1$P_Y)) for(mu in seq_along(x_INT1))
        eta1[nu,mu] <- (x_INT1[mu]/pmax(met1$P_Y[nu],1e-15))*rho[mu,nu]
      J_x <- J_xint_wrt_xsbj(x_OBJ,x_SBJ1)
      dP_dxsbj <- dP_dxint %*% J_x
      
      # WIN metrics
      win_metrics <- NULL
      if (!is.null(P_target)) {
        eps <- 1e-12; ce <- -sum(P_target*log(pmax(met1$P_Y,eps)))
        kl <- sum(P_target*log(pmax(P_target,eps)/pmax(met1$P_Y,eps)))
        brier <- mean((met1$P_Y - P_target)^2)
        top1_pred <- which.max(met1$P_Y); top1_true <- which.max(P_target)
        win_metrics <- list(CE=ce, KL=kl, Brier=brier, Top1=as.numeric(top1_pred==top1_true),
                            Top3=as.numeric(top1_true %in% order(met1$P_Y, decreasing=TRUE)[1:min(3,length(met1$P_Y))]),
                            P_target=P_target)
      }
      incProgress(0.12); pb(100)
      
      list(criteria_names=criteria, alternative_names=alts,
           Xi=Xi, Xi_prime=Xi_prime, rho=rho, h=h, d=dvec, x_OBJ=x_OBJ,
           x_SBJ0=x_SBJ0, x_INT0=x_INT0, met0=met0,
           x_SBJ1=x_SBJ1, x_INT1=x_INT1, met1=met1,
           P_Y_by_DM=P_Y_by_DM, tau_by_DM=tau_by_DM, DMI=DMI,
           dP_dxint=dP_dxint, eta1=eta1, J_x=J_x, dP_dxsbj=dP_dxsbj,
           WIN=d$WIN, WIN_metrics=win_metrics)
    })
  }, ignoreInit=FALSE)
  
  # ---- KPI
  output$kpi_nmgi0 <- renderText({ req(results()); sprintf("%.3f", results()$met0$NMGI) })
  output$kpi_nmgi1 <- renderText({ req(results()); sprintf("%.3f", results()$met1$NMGI) })
  output$kpi_csf1  <- renderText({ req(results()); sprintf("%.3f", results()$met1$CSF) })
  output$kpi_dmi   <- renderText({ req(results()); sprintf("%.3f", results()$DMI) })
  
  # ---- Dashboard visuals
  output$plt_py_compare <- renderPlotly({
    req(results()); r <- results()
    df <- tibble(alternative=rep(r$alternative_names,2),
                 value=c(r$met0$P_Y, r$met1$P_Y),
                 kind=rep(c("Baseline","Trained"), each=length(r$alternative_names)))
    ggplotly(ggplot(df, aes(y=reorder(alternative,value), x=value, fill=kind)) +
               geom_col(position="dodge") + labs(y="Alternative", x="P(Y)") +
               scale_x_continuous(labels=percent_format(accuracy=0.1)) + guides(fill="none"))
  })
  output$plt_indices_compare <- renderPlotly({
    req(results()); r <- results()
    base <- c(NMI=r$met0$NMI, CSF=r$met0$CSF, ADI=r$met0$ADI, CES=r$met0$CES, NMGI=r$met0$NMGI)
    fin  <- c(NMI=r$met1$NMI, CSF=r$met1$CSF, ADI=r$met1$ADI, CES=r$met1$CES, NMGI=r$met1$NMGI)
    df <- bind_rows(tibble(Index=names(base), Value=as.numeric(base), Kind="Baseline"),
                    tibble(Index=names(fin),  Value=as.numeric(fin),  Kind="Trained"))
    ggplotly(ggplot(df, aes(x=Index, y=Value, fill=Kind)) + geom_col(position="dodge") + labs(x=NULL,y="Index"))
  })
  output$plt_intw_compare <- renderPlotly({
    req(results()); r <- results()
    df <- tibble(criterion=rep(r$criteria_names,2),
                 value=c(r$x_INT0, r$x_INT1),
                 kind=rep(c("Baseline","Trained"), each=length(r$criteria_names)))
    ggplotly(ggplot(df, aes(y=reorder(criterion,value), x=value, fill=kind)) +
               geom_col(position="dodge") + labs(y="Criterion", x="x^INT") +
               scale_x_continuous(labels=percent_format(accuracy=0.1)) + guides(fill="none"))
  })
  output$plt_sbj_compare <- renderPlotly({
    req(results()); r <- results()
    df <- tibble(criterion=rep(r$criteria_names,2),
                 value=c(r$x_SBJ0, r$x_SBJ1),
                 kind=rep(c("Baseline","Trained"), each=length(r$criteria_names)))
    ggplotly(ggplot(df, aes(y=reorder(criterion,value), x=value, fill=kind)) +
               geom_col(position="dodge") + labs(y="Criterion", x="x^SBJ") +
               scale_x_continuous(labels=percent_format(accuracy=0.1)) + guides(fill="none"))
  })
  output$plt_delta_py <- renderPlotly({
    req(results()); r <- results()
    df <- tibble(alternative=r$alternative_names, Delta=r$met1$P_Y - r$met0$P_Y)
    ggplotly(ggplot(df, aes(y=reorder(alternative,Delta), x=Delta)) + geom_col() +
               labs(y="Alternative", x="Δ P(Y) (Trained − Baseline)") +
               scale_x_continuous(labels=percent_format(accuracy=0.1)))
  })
  output$plt_top_influence <- renderPlotly({
    req(results()); r <- results()
    infl <- colSums(abs(r$dP_dxsbj))
    df <- tibble(criterion=r$criteria_names, impact=infl)
    ggplotly(ggplot(df, aes(y=reorder(criterion,impact), x=impact)) + geom_col() +
               labs(y="Criterion", x="Total |∑ν ∂Pν/∂x^SBJμ|"))
  })
  
  # Radar wall
  make_radar <- function(cat, base, new, title){
    pad <- function(v) c(v, v[1]); padcats <- c(cat, cat[1])
    plot_ly(type='scatterpolar', mode='lines+markers') |>
      add_trace(r=pad(as.numeric(base)), theta=padcats, name="Baseline") |>
      add_trace(r=pad(as.numeric(new)),  theta=padcats, name="Trained") |>
      layout(title=list(text=title), polar=list(radialaxis=list(visible=TRUE)), showlegend=FALSE)
  }
  output$plt_radar_wall <- renderPlotly({
    req(results()); r <- results(); cats <- r$criteria_names
    if (isTRUE(input$radar_all)) {
      plots <- lapply(seq_along(r$alternative_names), function(j){
        b <- r$met0$P_XY[,j]; b <- b/sum(b); t <- r$met1$P_XY[,j]; t <- t/sum(t)
        make_radar(cats,b,t,r$alternative_names[j])
      })
      subplot(plots, nrows=ceiling(length(plots)/input$radar_cols), shareX=TRUE, shareY=TRUE, margin=0.03)
    } else {
      j <- match(input$radar_alt, r$alternative_names)
      b <- r$met0$P_XY[,j]; b <- b/sum(b); t <- r$met1$P_XY[,j]; t <- t/sum(t)
      make_radar(cats,b,t,paste0("Profile — ",r$alternative_names[j]))
    }
  })
  
  # Neural tab
  output$tbl_sbj <- renderDT({ req(results()); tibble(criterion=results()$criteria_names, x_SBJ=results()$x_SBJ1) |>
      datatable(options=list(pageLength=10, scrollX=TRUE), rownames=FALSE) })
  output$plt_sbj_final <- renderPlotly({
    req(results()); df <- tibble(criterion=results()$criteria_names, x_SBJ=results()$x_SBJ1)
    ggplotly(ggplot(df, aes(y=reorder(criterion,x_SBJ), x=x_SBJ)) + geom_col() +
               labs(y="Criterion", x="x^SBJ (Trained)") + scale_x_continuous(labels=percent_format(accuracy=0.1)))
  })
  
  # Objective & Integrated
  output$tbl_h_d_xobj <- renderDT({ req(results()); tibble(criterion=results()$criteria_names, h_mu=results()$h, d_mu=results()$d, x_OBJ=results()$x_OBJ) |>
      datatable(options=list(pageLength=10, scrollX=TRUE), rownames=FALSE) })
  output$plt_entropy <- renderPlotly({
    req(results()); df <- tibble(criterion=results()$criteria_names, h_mu=results()$h, d_mu=results()$d) |>
      pivot_longer(cols=c(h_mu,d_mu), names_to="metric", values_to="value")
    ggplotly(ggplot(df, aes(y=reorder(criterion,value), x=value, fill=metric)) + geom_col(position="dodge") +
               labs(y="Criterion", x="Entropy / Diversification") + guides(fill=guide_legend(title=NULL)))
  })
  output$plt_compare_w_all <- renderPlotly({
    req(results()); r <- results()
    df <- tibble(criterion=r$criteria_names, `x^OBJ`=r$x_OBJ, `x^INT (Base)`=r$x_INT0, `x^INT (Trained)`=r$x_INT1) |>
      pivot_longer(-criterion, names_to="scheme", values_to="value")
    ggplotly(ggplot(df, aes(y=reorder(criterion,value), x=value, fill=scheme)) + geom_col(position="dodge") +
               labs(y="Criterion", x="Weight") + scale_x_continuous(labels=percent_format(accuracy=0.1)))
  })
  output$tbl_intw <- renderDT({ req(results()); tibble(criterion=results()$criteria_names, x_OBJ=results()$x_OBJ, x_INT_base=results()$x_INT0, x_INT_trained=results()$x_INT1) |>
      datatable(options=list(pageLength=10, scrollX=TRUE), rownames=FALSE) })
  
  # Probability structure
  output$plt_py_cum <- renderPlotly({
    req(results()); df <- tibble(alternative=results()$alternative_names, P_Y=results()$met1$P_Y) |>
      arrange(desc(P_Y)) |> mutate(cum=cumsum(P_Y))
    ggplotly(ggplot(df, aes(x=seq_along(P_Y), y=cum)) + geom_line() + geom_point() +
               labs(x="Rank position", y="Cumulative P(Y)"))
  })
  output$tbl_sys_entropy <- renderTable({
    req(results()); r <- results(); with(r$met1, data.frame(
      `S(Y|X)`=S_Y_given_X, `I(Y|X)`=I_Y_given_X, `S(X)`=S_X, `S(Y)`=S_Y, `S(X,Y)`=S_XY, `J(X;Y)`=J_XY, check.names=FALSE))
  }, striped=TRUE, bordered=TRUE, spacing="m")
  output$plt_heat_pxy <- renderPlotly({
    req(results()); r <- results()
    df <- as.data.frame(r$met1$P_XY); colnames(df) <- r$alternative_names; df$criterion <- r$criteria_names
    df_long <- df %>% pivot_longer(-criterion, names_to="alternative", values_to="prob")
    base <- ggplot(df_long, aes(x=alternative, y=criterion, fill=prob)) +
      geom_tile() + scale_fill_viridis_c() + labs(x="Alternative", y="Criterion", fill="P(Xμ,Yν)") +
      theme(axis.text.x=element_text(angle=45, hjust=1))
    ggplotly(base) %>% layout(margin=list(l=80,r=20,b=80,t=40))
  })
  
  # Sensitivity
  output$tbl_elasticity <- renderDT({
    req(results()); r <- results()
    E <- r$eta1; colnames(E) <- r$criteria_names; rownames(E) <- r$alternative_names
    datatable(round(E, 4), options=list(pageLength=10, scrollX=TRUE))
  })
  output$tbl_dP_dxsbj <- renderDT({
    req(results()); r <- results()
    G <- r$dP_dxsbj; colnames(G) <- paste0("∂/∂", r$criteria_names); rownames(G) <- r$alternative_names
    datatable(round(G, 6), options=list(pageLength=10, scrollX=TRUE))
  })
  output$tbl_wc_bounds <- renderTable({
    req(results()); r <- results()
    eps <- input$wc_eps %||% 0.05
    max_grad <- apply(r$dP_dxint, 1, function(row) max(abs(row)))
    data.frame(alternative = r$alternative_names,
               max_abs_delta = round(eps * max_grad, 6),
               check.names = FALSE)
  }, striped=TRUE, bordered=TRUE, spacing="m")
  output$plt_uncertainty <- renderPlotly({
    req(results()); r <- results()
    sigma <- input$sigma_int %||% 0.05
    M <- length(r$x_INT1); var_x <- rep(sigma^2, M)
    CovP <- r$dP_dxint %*% diag(var_x, nrow=M, ncol=M) %*% t(r$dP_dxint)
    se <- sqrt(pmax(diag(CovP), 0))
    df <- tibble(alternative = r$alternative_names, mean = r$met1$P_Y,
                 lo = pmax(r$met1$P_Y - 1.96*se, 0), hi = pmin(r$met1$P_Y + 1.96*se, 1))
    ggplotly(ggplot(df, aes(y=reorder(alternative, mean), x=mean)) +
               geom_col(alpha=0.85) + geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.25) +
               labs(y="Alternative", x="P(Y) ± 1.96·SE"))
  })
  
  # Supervision & tuning
  output$plt_py_vs_ptarget <- renderPlotly({
    req(results()); r <- results(); wm <- r$WIN_metrics
    if (is.null(wm) || is.null(wm$P_target)) return(NULL)
    df <- tibble(alternative=r$alternative_names, `Trained P(Y)`=r$met1$P_Y, `Target`=wm$P_target) |>
      pivot_longer(-alternative, names_to="kind", values_to="P")
    ggplotly(ggplot(df, aes(y=reorder(alternative,P), x=P, fill=kind)) + geom_col(position="dodge") +
               labs(y="Alternative", x="Probability"))
  })
  output$tbl_win_metrics <- renderTable({
    req(results()); wm <- results()$WIN_metrics
    if (is.null(wm)) return(data.frame(Message="No WIN metrics/training (WIN missing or disabled)."))
    data.frame(CrossEntropy=round(wm$CE,6), KL_Ptarget_PY=round(wm$KL,6),
               Brier=round(wm$Brier,6), Top1_match=wm$Top1, Top3_contains=wm$Top3, check.names=FALSE)
  }, striped=TRUE, bordered=TRUE, spacing="m")
  observeEvent(input$btn_grid, {
    req(results()); r0 <- results(); d <- data_input(); req(d$WIN); req(isTRUE(input$use_win))
    WIN2 <- d$WIN %>% right_join(data.frame(alternative=r0$alternative_names), by="alternative")
    if (!("P_target" %in% names(WIN2))) {
      if ("rank" %in% names(WIN2)) {
        rnk <- to_num(WIN2$rank); rnk[is.na(rnk)] <- max(rnk, na.rm=TRUE)+5
        logits <- -rnk / max(0.01, input$win_tau); logits <- logits - max(logits, na.rm=TRUE)
        P_target <- exp(logits); P_target[!is.finite(P_target)] <- 0; P_target <- P_target/sum(P_target)
      } else return(NULL)
    } else P_target <- normalize_simplex(to_num(WIN2$P_target))
    SBJ <- d$SBJ; if ("criterion" %in% names(SBJ)) SBJ <- dplyr::select(SBJ,-criterion); W_SBJ <- as.matrix(SBJ)
    rho_dm <- to_num(d$DM_Rho$rho)
    
    la <- seq(0, max(0,input$grid_lambda_align_max), length.out=15)
    le <- seq(0, max(0,input$grid_lambda_ent_max),   length.out=9)
    Xg <- expand.grid(lambda_align=la, lambda_ent=le); CE <- numeric(nrow(Xg))
    for(i in seq_len(nrow(Xg))){
      xa <- Xg$lambda_align[i]; xe <- Xg$lambda_ent[i]
      x_SBJ_try <- neural_consensus_fast(W_SBJ, rho_dm=rho_dm, x_OBJ=r0$x_OBJ, lambda_ent=xe, lambda_align=xa)
      x_INT_try <- integrated_weights(r0$x_OBJ, x_SBJ_try)
      P_Y_try <- as.numeric(t(r0$rho)%*%x_INT_try)
      CE[i] <- -sum(P_target*log(pmax(P_Y_try,1e-12)))
    }
    Xg$CE <- CE; best <- Xg[which.min(Xg$CE), , drop=FALSE]
    output$plt_grid <- renderPlotly({
      ggplotly(ggplot(Xg, aes(x=lambda_align, y=lambda_ent, fill=CE)) + geom_tile() + scale_fill_viridis_c() +
                 geom_point(data=best, aes(x=lambda_align, y=lambda_ent), size=3, shape=21, stroke=1.2, fill="white") +
                 labs(x="λ_align", y="λ_ent", fill="Cross-Entropy"))
    })
    showNotification(sprintf("Best λ: align=%.3f, ent=%.3f (CE=%.5f)", best$lambda_align, best$lambda_ent, best$CE),
                     type="message", duration=6)
  })
  
  # Notation
  notation_df <- reactive({
    data.frame(
      Symbol=c("Xi","type","ρ (rho)","h_μ","d_μ","x^OBJ","W_SBJ","ρ_DM","x^SBJ","x^INT","P(Y)",
               "S(Y), S(Y|X)","S(X), S(X,Y)","J(X;Y), NMI","CES","CSF","ADI","NMGI","DMI",
               "∂P/∂x^INT","Elasticity η"),
      Meaning=c(
        "Raw scores (criteria × alternatives)","Benefit/Cost per criterion","Row-normalized intensities",
        "Normalized row entropy","Diversification (1 − h_μ)","Objective weights","Panel subjective weights",
        "DM reliabilities","Consensus subjective weights","Integrated weights (normalize x^OBJ ⊙ x^SBJ)",
        "Alternative probabilities","Alt. entropy / conditional entropy","Criterion entropy / joint entropy",
        "Mutual information & normalized","Complementary Entropy Share","Criterion Support Factor",
        "Alternative Distinction Index","Composite index","Mean Kendall τ (panel agreement)",
        "Jacobian to x^INT (= ρ^T)","Log–log elasticity of P(Y) w.r.t x^INT"
      ),
      stringsAsFactors=FALSE
    )
  })
  output$tbl_notation <- renderDT(datatable(notation_df(), options=list(pageLength=25, scrollX=TRUE), rownames=FALSE))
  output$dl_notation <- downloadHandler(
    filename=function() "notation_esmadm2.xlsx",
    content=function(file){ wb <- createWorkbook(); addWorksheet(wb,"Notation"); writeData(wb,"Notation", notation_df()); saveWorkbook(wb,file,overwrite=TRUE) }
  )
  
  # Summary & Export
  output$txt_overall <- renderText({
    req(results()); r <- results()
    ord <- order(r$met1$P_Y, decreasing=TRUE)
    lead <- r$alternative_names[ord][1]
    gap  <- if (length(ord)>=2) r$met1$P_Y[ord][1] - r$met1$P_Y[ord][2] else r$met1$P_Y[ord][1]
    infl <- r$criteria_names[which.max(colSums(abs(r$dP_dxsbj)))]
    paste0("Leader: ", lead, " (margin ", sprintf("%.2f%%", 100*gap),
           "). NMGI ", sprintf("%.3f", r$met0$NMGI)," → ",sprintf("%.3f", r$met1$NMGI),
           "; CSF=",sprintf("%.3f", r$met1$CSF), "; DMI=",sprintf("%.3f", r$DMI),
           ". Most influential criterion: ", infl, ".")
  })
  output$tbl_insight <- renderDT({
    req(results()); r <- results()
    ord <- order(r$met1$P_Y, decreasing=TRUE); N <- length(r$alternative_names)
    tibble::tibble(
      Item=c("Top-1","Top-3 cumulative P(Y)","NMGI","CSF","ADI","NMI"),
      Baseline=c(r$alternative_names[order(r$met0$P_Y, decreasing=TRUE)][1],
                 sprintf("%.2f%%", 100*sum(sort(r$met0$P_Y, decreasing=TRUE)[1:min(3,N)])),
                 sprintf("%.3f", r$met0$NMGI), sprintf("%.3f", r$met0$CSF),
                 sprintf("%.3f", r$met0$ADI), sprintf("%.3f", r$met0$NMI)),
      Trained=c(r$alternative_names[ord][1],
                sprintf("%.2f%%", 100*sum(r$met1$P_Y[ord][1:min(3,N)])),
                sprintf("%.3f", r$met1$NMGI), sprintf("%.3f", r$met1$CSF),
                sprintf("%.3f", r$met1$ADI), sprintf("%.3f", r$met1$NMI))
    ) |> datatable(options=list(dom='t', pageLength=15, scrollX=TRUE), rownames=FALSE)
  })
  output$txt_summary <- renderText({
    req(results()); r <- results(); ord <- order(r$met1$P_Y, decreasing=TRUE)
    top <- paste0(r$alternative_names[ord][1:min(3,length(ord))], collapse=", ")
    paste0("Top-ranked: ", top,
           "\nNMGI(base→trained) = ", round(r$met0$NMGI,4)," → ",round(r$met1$NMGI,4),
           " | CSF=", round(r$met1$CSF,4),
           " | ADI=", round(r$met1$ADI,4),
           " | NMI=", round(r$met1$NMI,4),
           " | DMI=", round(r$DMI,4))
  })
  output$dl_results <- downloadHandler(
    filename=function() "esmadm2_results.xlsx",
    content=function(file){
      req(results()); r <- results()
      wb <- createWorkbook()
      addWorksheet(wb,"Weights")
      writeData(wb,"Weights", data.frame(criterion=r$criteria_names,
                                         x_OBJ=r$x_OBJ, x_SBJ_base=r$x_SBJ0, x_INT_base=r$x_INT0,
                                         x_SBJ_trained=r$x_SBJ1, x_INT_trained=r$x_INT1))
      addWorksheet(wb,"Alternatives")
      writeData(wb,"Alternatives", data.frame(alternative=r$alternative_names,
                                              P_Y_base=r$met0$P_Y, P_Y_trained=r$met1$P_Y))
      addWorksheet(wb,"System_Base")
      writeData(wb,"System_Base", data.frame(S_Y=r$met0$S_Y, S_Y_given_X=r$met0$S_Y_given_X, I_Y_given_X=r$met0$I_Y_given_X,
                                             S_X=r$met0$S_X, S_XY=r$met0$S_XY, J_XY=r$met0$J_XY,
                                             NMI=r$met0$NMI, CES=r$met0$CES, CSF=r$met0$CSF, ADI=r$met0$ADI, NMGI=r$met0$NMGI))
      addWorksheet(wb,"System_Trained")
      writeData(wb,"System_Trained", data.frame(S_Y=r$met1$S_Y, S_Y_given_X=r$met1$S_Y_given_X, I_Y_given_X=r$met1$I_Y_given_X,
                                                S_X=r$met1$S_X, S_XY=r$met1$S_XY, J_XY=r$met1$J_XY,
                                                NMI=r$met1$NMI, CES=r$met1$CES, CSF=r$met1$CSF, ADI=r$met1$ADI, NMGI=r$met1$NMGI, DMI=r$DMI))
      addWorksheet(wb,"Elasticities_Final")
      E <- r$eta1; colnames(E) <- r$criteria_names; rownames(E) <- r$alternative_names
      writeData(wb,"Elasticities_Final", cbind(Alternative=rownames(E), round(E,6)))
      addWorksheet(wb,"dP_dxSBJ_Final")
      G <- r$dP_dxsbj; colnames(G) <- paste0("dP/d", r$criteria_names); rownames(G) <- r$alternative_names
      writeData(wb,"dP_dxSBJ_Final", cbind(Alternative=rownames(G), round(G,8)))
      if (!is.null(r$WIN_metrics)) {
        addWorksheet(wb,"WIN_metrics")
        wm <- r$WIN_metrics
        writeData(wb,"WIN_metrics", data.frame(CrossEntropy=wm$CE, KL_Ptarget_PY=wm$KL, Brier=wm$Brier, Top1_match=wm$Top1, Top3_contains=wm$Top3))
        addWorksheet(wb,"P_vs_Target")
        writeData(wb,"P_vs_Target", data.frame(alternative=r$alternative_names, P_Y_trained=r$met1$P_Y, P_target=r$WIN_metrics$P_target))
      }
      saveWorkbook(wb,file,overwrite=TRUE)
    }
  )
}

shinyApp(ui, server)
