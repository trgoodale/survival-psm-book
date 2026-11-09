# Copyright (c) 2026 Tom Goodale. SPDX-License-Identifier: MIT
# Propensity score explorer: matching / weighting + survival outcome
#
# Run from RStudio with the "Run App" button, or from the R console:
#   shiny::runApp("apps/psm-explorer")
#
# Packages: shiny, survival, ggsurvfit, MatchIt, cobalt, broom, dplyr, ggplot2

library(shiny)
library(survival)
library(ggsurvfit)
library(MatchIt)
library(cobalt)
library(broom)
library(dplyr)
library(ggplot2)

# Same simulation as R/simulate_cohort.R, with the strength of confounding
# exposed as a parameter (1 = as in the book).
simulate_cohort <- function(n = 2000, hr_treat = 0.70, conf = 1, seed = 2026) {
  set.seed(seed)
  age      <- round(pmin(pmax(rnorm(n, 65, 10), 35), 95))
  female   <- rbinom(n, 1, 0.45)
  diabetes <- rbinom(n, 1, plogis(-1.4 + 0.04 * (age - 65)))
  bmi      <- round(rnorm(n, 28, 4.5), 1)
  smoker   <- rbinom(n, 1, 0.25)
  lp_treat <- -0.6 + conf * (0.06 * (age - 65) + 1.0 * diabetes + 0.08 * (bmi - 28) -
    0.3 * female + 0.4 * smoker)
  treated <- rbinom(n, 1, plogis(lp_treat))
  lp_death <- log(hr_treat) * treated + 0.06 * (age - 65) + 0.6 * diabetes +
    0.03 * (bmi - 28) - 0.25 * female + 0.55 * smoker
  t_event <- (-log(runif(n)) / (0.05 * exp(lp_death)))^(1 / 1.2)
  t_cens  <- pmin(runif(n, 3, 8), rexp(n, 0.04))
  data.frame(id = seq_len(n), age, female, diabetes, bmi, smoker, treated,
             time = pmin(t_event, t_cens), status = as.integer(t_event <= t_cens))
}

all_covs <- c("age", "female", "diabetes", "bmi", "smoker")

ui <- fluidPage(
  titlePanel("Propensity score explorer"),
  sidebarLayout(
    sidebarPanel(width = 3,
      h4("1. Simulate a cohort"),
      sliderInput("n", "Patients", 500, 10000, 2000, step = 500),
      sliderInput("hr", "True conditional HR for treatment", 0.3, 1.5, 0.7, step = 0.05),
      sliderInput("conf", "Strength of confounding (1 = book)", 0, 3, 1, step = 0.25),
      numericInput("seed", "Seed", 2026),
      hr(),
      h4("2. Propensity score model"),
      checkboxGroupInput("covs", "Covariates included", all_covs, selected = all_covs),
      helpText("All five are true confounders. Leave one out to see the effect of an unmeasured confounder."),
      hr(),
      h4("3. Method"),
      radioButtons("method", NULL,
                   c("Nearest-neighbour matching (ATT)" = "nearest",
                     "Inverse probability weighting (ATE)" = "iptw")),
      conditionalPanel("input.method == 'nearest'",
        sliderInput("ratio", "Controls per treated (ratio)", 1, 3, 1),
        checkboxInput("replace", "Match with replacement", FALSE),
        checkboxInput("usecal", "Use a caliper", TRUE),
        sliderInput("caliper", "Caliper (SD of logit PS)", 0.01, 1, 0.2, step = 0.01)
      )
    ),
    mainPanel(width = 9,
      tabsetPanel(
        tabPanel("Balance",
          verbatimTextOutput("sample_sizes"),
          plotOutput("love", height = "330px"),
          plotOutput("psdist", height = "330px")),
        tabPanel("Survival",
          plotOutput("km", height = "480px")),
        tabPanel("Treatment effect",
          p("Hazard ratio for treatment from each approach. The red dashed line is the ",
            "true conditional HR used in the simulation. The PS-based estimates are ",
            em("marginal"), " HRs, which are usually slightly closer to 1 (non-collapsibility)."),
          plotOutput("effects", height = "300px"),
          tableOutput("effects_table"))
      )
    )
  )
)

server <- function(input, output, session) {

  cohort <- reactive(simulate_cohort(input$n, input$hr, input$conf, input$seed))

  ps_formula <- reactive({
    covs <- if (length(input$covs)) input$covs else "1"
    as.formula(paste("treated ~", paste(covs, collapse = " + ")))
  })

  # Returns a list: a balance object for cobalt, the analysis data with weights,
  # and a description.
  design <- reactive({
    d <- cohort()
    if (input$method == "nearest") {
      m <- matchit(ps_formula(), data = d, method = "nearest", link = "linear.logit",
                   ratio = input$ratio, replace = input$replace,
                   caliper = if (input$usecal) input$caliper else NULL)
      # With replacement a control can appear in several pairs, so cluster
      # on the patient rather than the pair.
      md <- if (input$replace) get_matches(m, id = "unit_id") else match.data(m)
      if (input$replace) md$subclass <- md$unit_id
      list(bal = m, data = md, cluster = "subclass", m = m)
    } else {
      ps <- fitted(glm(ps_formula(), data = d, family = binomial))
      d$ps <- ps
      d$weights <- ifelse(d$treated == 1, 1 / ps, 1 / (1 - ps))
      d$subclass <- d$id
      list(bal = NULL, data = d, cluster = "id", ps = ps)
    }
  })

  output$sample_sizes <- renderPrint({
    des <- design()
    if (input$method == "nearest") {
      print(summary(des$m)$nn)
    } else {
      w <- des$data$weights
      t <- des$data$treated
      ess <- function(w) sum(w)^2 / sum(w^2)
      cat("Effective sample size after weighting:\n",
          " Treated:", round(ess(w[t == 1])), "of", sum(t == 1), "\n",
          " Control:", round(ess(w[t == 0])), "of", sum(t == 0), "\n",
          " Largest weight:", round(max(w), 1), "\n")
    }
  })

  # Balance is always checked on ALL five confounders, even ones left out of
  # the PS model, so you can see what an "unmeasured" confounder looks like.
  bal_obj <- reactive({
    des <- design()
    if (input$method == "nearest") {
      bal.tab(des$m, addl = as.formula(paste("~", paste(all_covs, collapse = "+"))),
              binary = "std", un = TRUE)
    } else {
      bal.tab(as.formula(paste("treated ~", paste(all_covs, collapse = "+"))),
              data = des$data, weights = des$data$weights, distance = des$data$ps,
              binary = "std", s.d.denom = "pooled", estimand = "ATE", un = TRUE)
    }
  })

  output$love <- renderPlot({
    love.plot(bal_obj(), abs = TRUE, thresholds = c(m = 0.1), var.order = "unadjusted",
              colors = c("#D55E00", "#0072B2"), shapes = c("circle", "triangle"), size = 4,
              title = "Covariate balance (all five true confounders)")
  })

  output$psdist <- renderPlot({
    des <- design()
    if (input$method == "nearest") {
      bal.plot(des$m, var.name = "distance", which = "both", type = "histogram",
               mirror = TRUE, colors = c("#E69F00", "#0072B2"))
    } else {
      bal.plot(treated ~ ps, data = des$data, weights = des$data$weights,
               var.name = "ps", which = "both", type = "histogram", mirror = TRUE,
               colors = c("#E69F00", "#0072B2"))
    }
  })

  output$km <- renderPlot({
    des <- design()
    d <- des$data
    survfit2(Surv(time, status) ~ treated, data = d, weights = weights) |>
      ggsurvfit(linewidth = 1) +
      add_confidence_interval() +
      labs(x = "Years", y = "Survival probability",
           title = if (input$method == "nearest") "Matched cohort" else "IPT-weighted cohort") +
      theme(text = element_text(size = 15))
  })

  effects <- reactive({
    d <- cohort()
    des <- design()
    ad <- des$data
    adj_f <- as.formula(paste("Surv(time, status) ~ treated +", paste(all_covs, collapse = " + ")))
    fits <- list(
      "Naive" = coxph(Surv(time, status) ~ treated, data = d),
      "Fully adjusted Cox" = coxph(adj_f, data = d),
      "PS analysis" = coxph(Surv(time, status) ~ treated, data = ad, weights = weights,
                            robust = TRUE, cluster = ad$subclass)
    )
    names(fits)[3] <- if (input$method == "nearest") "PS-matched Cox" else "IPTW Cox"
    bind_rows(lapply(fits, \(f) tidy(f, exponentiate = TRUE, conf.int = TRUE)[1, ]), .id = "Method") |>
      transmute(Method = factor(Method, levels = rev(names(fits))),
                HR = estimate, Lower = conf.low, Upper = conf.high)
  })

  output$effects <- renderPlot({
    ggplot(effects(), aes(HR, Method)) +
      geom_vline(xintercept = 1, colour = "grey60") +
      geom_vline(xintercept = input$hr, linetype = "dashed", colour = "#D55E00") +
      geom_errorbar(aes(xmin = Lower, xmax = Upper), width = 0.15, orientation = "y") +
      geom_point(size = 4, colour = "#0072B2") +
      scale_x_log10() +
      labs(x = "Hazard ratio (log scale)", y = NULL) +
      theme_minimal(base_size = 16)
  })

  output$effects_table <- renderTable(effects() |> arrange(desc(Method)), digits = 3)
}

shinyApp(ui, server)
