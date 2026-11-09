# Copyright (c) 2026 Tom Goodale. SPDX-License-Identifier: MIT
# Survival explorer: Kaplan-Meier curves and Cox models on real datasets
#
# Run from RStudio with the "Run App" button, or from the R console:
#   shiny::runApp("apps/survival-explorer")
#
# Packages: shiny, survival, ggsurvfit, broom, dplyr, ggplot2

library(shiny)
library(survival)
library(ggsurvfit)
library(broom)
library(dplyr)
library(ggplot2)

# ---- Datasets (all ship with the survival package) -------------------------
datasets <- list(
  "Lung cancer (NCCTG)" = survival::lung |>
    transmute(time, status = status - 1,
              sex = factor(sex, 1:2, c("Male", "Female")),
              age, age_group = cut(age, c(0, 60, 70, Inf), c("<60", "60-69", "70+")),
              ecog = factor(pmin(ph.ecog, 2), 0:2, c("0", "1", "2+")),
              wt.loss),
  "Veterans' lung cancer trial" = survival::veteran |>
    transmute(time, status,
              treatment = factor(trt, 1:2, c("Standard", "Test")),
              celltype, karno, age, prior = factor(prior, c(0, 10), c("No", "Yes"))),
  "Colon cancer adjuvant trial (death)" = survival::colon |>
    filter(etype == 2) |>
    transmute(time, status, rx, sex = factor(sex, 0:1, c("Female", "Male")),
              age, nodes, obstruct = factor(obstruct, 0:1, c("No", "Yes")),
              differ = factor(differ, 1:3, c("Well", "Moderate", "Poor")))
)

ui <- fluidPage(
  titlePanel("Survival explorer: Kaplan–Meier & Cox regression"),
  sidebarLayout(
    sidebarPanel(width = 3,
      selectInput("data", "Dataset", names(datasets)),
      helpText("Time is in days. status: 1 = event, 0 = censored."),
      hr(),
      h4("Kaplan–Meier"),
      selectInput("group", "Compare groups by", choices = NULL),
      checkboxInput("ci", "Show 95% confidence intervals", TRUE),
      checkboxInput("risk", "Show risk table", TRUE),
      checkboxInput("median", "Mark median survival", FALSE),
      hr(),
      h4("Cox model"),
      selectizeInput("covs", "Covariates", choices = NULL, multiple = TRUE)
    ),
    mainPanel(width = 9,
      tabsetPanel(
        tabPanel("Kaplan–Meier",
          plotOutput("km_plot", height = "520px"),
          h4("Median survival"), tableOutput("km_table"),
          h4("Log-rank test"), verbatimTextOutput("logrank")),
        tabPanel("Cox model",
          h4("Hazard ratios"), tableOutput("cox_table"),
          plotOutput("forest", height = "300px"),
          verbatimTextOutput("cox_fit")),
        tabPanel("Proportional hazards check",
          p("Scaled Schoenfeld residuals: a flat line is consistent with ",
            "proportional hazards. A small p-value suggests the HR changes over time."),
          tableOutput("zph_table"),
          plotOutput("zph_plot", height = "420px")),
        tabPanel("Data", tableOutput("data_head"))
      )
    )
  )
)

server <- function(input, output, session) {

  dat <- reactive(datasets[[input$data]])

  observeEvent(dat(), {
    vars <- setdiff(names(dat()), c("time", "status"))
    cat_vars <- vars[sapply(dat()[vars], \(x) is.factor(x) || is.character(x))]
    updateSelectInput(session, "group", choices = c("(none)" = "1", cat_vars),
                      selected = cat_vars[1])
    updateSelectizeInput(session, "covs", choices = vars, selected = vars[1:2])
  })

  km_formula <- reactive({
    req(input$group)
    as.formula(paste("Surv(time, status) ~", input$group))
  })

  output$km_plot <- renderPlot({
    fit <- survfit2(km_formula(), data = dat())
    p <- ggsurvfit(fit, linewidth = 1) +
      add_censor_mark(alpha = 0.5) +
      labs(x = "Days", y = "Survival probability") +
      theme(text = element_text(size = 15))
    if (input$ci)     p <- p + add_confidence_interval()
    if (input$median) p <- p + add_quantile(y_value = 0.5, linetype = "dashed")
    if (input$risk)   p <- p + add_risktable(risktable_stats = "n.risk", size = 4.5)
    if (input$group != "1") p <- p + add_pvalue(location = "annotation")
    p
  })

  output$km_table <- renderTable({
    fit <- survfit(km_formula(), data = dat())
    tab <- as.data.frame(summary(fit)$table)
    tab$group <- if (nrow(tab) > 1 || input$group != "1") rownames(tab) else "All"
    tab[, c("group", "records", "events", "median", "0.95LCL", "0.95UCL")]
  }, digits = 0)

  output$logrank <- renderPrint({
    if (input$group == "1") return(cat("Choose a grouping variable to compare curves."))
    survdiff(km_formula(), data = dat())
  })

  cox_fit <- reactive({
    req(length(input$covs) > 0)
    coxph(as.formula(paste("Surv(time, status) ~", paste(input$covs, collapse = " + "))),
          data = dat())
  })

  output$cox_table <- renderTable({
    tidy(cox_fit(), exponentiate = TRUE, conf.int = TRUE) |>
      transmute(Term = term, HR = estimate, `Lower 95%` = conf.low,
                `Upper 95%` = conf.high, `p-value` = format.pval(p.value, digits = 2))
  }, digits = 3)

  output$forest <- renderPlot({
    tidy(cox_fit(), exponentiate = TRUE, conf.int = TRUE) |>
      mutate(term = factor(term, levels = rev(term))) |>
      ggplot(aes(estimate, term)) +
      geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
      geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.2, orientation = "y") +
      geom_point(size = 3, colour = "#0072B2") +
      scale_x_log10() +
      labs(x = "Hazard ratio (log scale)", y = NULL) +
      theme_minimal(base_size = 15)
  })

  output$cox_fit <- renderPrint({
    f <- cox_fit()
    cat("n =", f$n, " events =", f$nevent,
        "\nConcordance =", round(f$concordance["concordance"], 3), "\n")
  })

  zph <- reactive(cox.zph(cox_fit()))

  output$zph_table <- renderTable({
    z <- zph()$table
    data.frame(Term = rownames(z), `Chi-square` = z[, "chisq"], df = z[, "df"],
               `p-value` = format.pval(z[, "p"], digits = 2), check.names = FALSE)
  }, digits = 2)

  output$zph_plot <- renderPlot({
    z <- cox.zph(cox_fit(), terms = FALSE)   # one panel per coefficient
    k <- ncol(z$y)
    par(mfrow = c(ceiling(k / 4), min(k, 4)), mar = c(4, 4, 2, 1))
    for (i in seq_len(k)) {
      plot(z[i], resid = TRUE, se = TRUE, lwd = 2, col = "#0072B2")
      abline(h = coef(cox_fit())[i], lty = 3, col = "grey40")
    }
  })

  output$data_head <- renderTable(head(dat(), 15))
}

shinyApp(ui, server)
