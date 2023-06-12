library(shiny)
library(pedtools)
library(pedprobr)
library(forrel)
library(tibble)
library(ggplot2)
library(patchwork)
library(tidyr)

CASES = c("Duo : unrel",
          "Sibs : unrel",
          "Halfsibs : unrel",
          "Grandp : unrel",
          "Sibs : halfsibs",
          "Grandp : uncle",
          "Grandp : halfsibs",
          "Uncle : halfsibs")
names(CASES) = CASES

IDS = c("A", "B")

PEDS = lapply(CASES, function(case) {
  rels = strsplit(case, " : ")[[1]]
  setNames(builtinPeds[tolower(rels)], c("Ped 1", "Ped 2"))
})

VARIABLES = c(
  "Distance (centiMorgan)" = "dist",
  "Recombination rate" = "rho",
  "Mutation rate" = "mutrate",
  "Frequency of '1'" = "freq")

debug = FALSE



ui = fluidPage(

  tags$head(tags$style(HTML("
  .well {margin-bottom: 6px; padding: 15px}
  .inline label{ display: table-cell; padding-right:3px }
  .inline .form-group { display: table-row;}
  .form-control {padding: 0px 1px 0px 8px}
  #variable .form-group {margin-bottom: 0px}
  #lastrow .form-group {margin-bottom: 0px}
  "))),

  # Application title
  titlePanel(
    title = HTML("<b>Linkage Lab:</b> Exploring kinship LR with linked markers"),
    windowTitle = "LinkageLab"),


  sidebarLayout(
    sidebarPanel(width = 3, style = "min-width:200px;",
      selectInput("comp", label = "Comparison", choices = CASES, selected = "Sibs : halfsibs"),

      HR,

      fluidRow(
        column(8, h4("Genotypes", style = "margin-top: 4px")),
        column(4, actionButton("sim", "Simulate", style = "float:right; padding:2px 10px",
                               class = "btn-sm btn-warning")),
      ),
      tags$div(class = "inline",
      fluidRow(
        column(6, style = "padding-right: 5px", wellPanel(
          style = "padding: 5px",
          p("A", style = "font-size: 18px; text-align: center; margin-bottom:0px"),
          textInput("A1", "M1:", width = "100%"),
          textInput("A2", "M2:", width = "100%")
        )),
        column(6, style = "padding-left: 5px", wellPanel(
          style = "padding: 5px",
          p("B", style = "font-size: 18px; text-align: center; margin-bottom:0px"),
          textInput("B1", "M1:", width = "100%"),
          textInput("B2", "M2:", width = "100%")
        )),
      )),

      actionButton("updatePeds", "Update pedigrees", style = "padding: 2px; margin-top:5px",
                   width = "100%", class = "btn-success"),


      HR,

      h4("Frequencies"),
      tags$div(class = "inline",
        fluidRow(
          column(6, numericInput("p1", label = "1:", value = 0.1, min = 0, max = 1, step = 0.01, width = "100%")),
          column(6, numericInput("p2", label = "2:", value = 0.2, min = 0, max = 1, step = 0.01, width = "100%")),
        ),
        fluidRow(
          column(6, numericInput("p3", label = "3:", value = 0.3, min = 0, max = 1, step = 0.01, width = "100%")),
          column(6, numericInput("p4", label = "4:", value = 0.4, min = 0, max = 1, step = 0.01, width = "100%")),
        ),
      ),

      HR,

      tags$div(class = "inline",
        numericInput("cm",    "Distance (cm):", value = 0, step = 1, min = 0, width = "100%"),
        numericInput("rho",   "Recomb rate:", value = 0, step = 0.01, min = 0, max = 1, width = "100%"),
        numericInput("mrate", "Mutation rate:", value = 0.001, step = 0.001, min = 0, max = 1, width = "100%")
      ),

      tags$div(id = "variable", selectInput("variable", "Variable", choices = VARIABLES)),
      fluidRow(id = "lastrow",
        column(6, numericInput("npoints", "Points", value = 10, min = 5, max = 50, step = 5)),
        column(6, actionButton("compute", "Plot!", class = "btn-primary", width = "100%",
                               style = "margin-top: 15px; font-size: 150%")),
      )
),
    # Plots
    mainPanel(width = 9, align = "left",
      fluidRow(
        column(width = 5, align = "left", style = "max-width:400px",
          plotOutput("pedplot1", height = "327px"),
          tags$div(style = "text-align: center", htmlOutput("ibs")),
          plotOutput("pedplot2", height = "327px")
        ),
        column(width = 7, align = "left", style = "max-width:600px",
          plotOutput("graphs", height = "674px"),
        )
      )
    )
  )
)


server = function(input, output, session) {

  genodat = reactive({
    if(debug) message("genodat")
    g = cbind(M1 = c(A = input$A1, B = input$B1),
              M2 = c(A = input$A2, B = input$B2))
    g[g == ""] = NA
    g
  })

  emptydat = cbind(M1 = c(A = NA, B = NA),
                   M2 = c(A = NA, B = NA))

  afreq = reactive(c("1" = input$p1, "2" = input$p2,
                     "3" = input$p3, "4" = 1-input$p1-input$p2-input$p3))

  observeEvent(afreq(), {
    if(debug) message("afreq-4")
    updateNumericInput(session, "p4", value = afreq()[["4"]])
  })

  peds = reactiveVal(NULL)
  peds(isolate(PEDS[[input$comp]] |> setGenos(emptydat)))

  observeEvent(input$comp, {
    if(debug) message("newcomp")
    newpeds = tryCatch(
      PEDS[[input$comp]] |> setGenos(genodat()),
      error = function(e) {
        errModal(e)
        PEDS[[input$comp]] |> setGenos(emptydat)
    })
    peds(newpeds)
  })

  observeEvent(input$updatePeds, {
    if(debug) message("update")
    newpeds = tryCatch(peds() |> setGenos(genodat()), error = errModal)
    peds(req(newpeds))
  })

  observeEvent(input$sim, {
    if(debug) message("sim")

    gt = peds()[[sample(1:2, size = 1)]] |>
      markerSim(N = 2, ids = IDS, alleles = 1:4, verbose = FALSE) |>
      getGenotypes(IDS)
    colnames(gt) = c("M1", "M2")

    updateTextInput(session, "A1", value = gt[["A", "M1"]])
    updateTextInput(session, "A2", value = gt[["A", "M2"]])
    updateTextInput(session, "B1", value = gt[["B", "M1"]])
    updateTextInput(session, "B2", value = gt[["B", "M2"]])
    peds(peds() |> setGenos(gt))
  })

  observeEvent(req(input$rho), {
    cm = req(haldane(rho = input$rho))
    if(abs(input$cm - cm) > sqrt(.Machine$double.eps))
      updateNumericInput(session, "cm", value = cm)
  })
  observeEvent(req(input$cm), {
    rho = req(haldane(cM = input$cm))
    if(abs(input$rho - rho) > sqrt(.Machine$double.eps))
      updateNumericInput(session, "rho", value = rho)
  })

  LRdat = reactiveVal(NULL)

  observeEvent({input$compute; input$variable; input$npoints}, {
    if(debug) message("compute")

    N = req(input$npoints)
    afr = afreq()
    mrate = input$mrate
    rho = input$rho

    # Checks
    err = NULL
    if(N < 3 || N > 100)
      err = "The number of points must be between 3 and 100"
    if(any(is.na(afr) | afr < 0 | afr > 1))
      err = "Allele frequencies must be between 0 and 1"
    if(is.na(mrate) || mrate < 0 || mrate > 1)
      err = "The mutation rate must be between 0 and 1"
    if(is.na(rho) || rho < 0 || rho > 0.5)
      err = "The recombination rate must be between 0 and 0.5"

    if(!is.null(err)) {
      errModal(err)
      return()
    }

    peds = tryCatch(error = errModal, warning = errModal,
      lapply(peds(), function(p)
        p |> setAfreq12(afr) |> setMutmod(model = "equal", rate = mrate)))
    req(peds)

    lik2 = function(ped, r = rho) {
      likelihood2(ped, 1, 2, rho = r)
    }

    setFr1 = function(a)
      c("1" = a, afr[2:3], "4" = 1 - a - sum(afr[2:3]))

    dat = switch(input$variable,
      dist = tibble(
        x = seq(0, 200, length = N),
        lik1 = sapply(x, function(cm) lik2(peds[[1]], r = haldane(cM = cm))),
        lik2 = sapply(x, function(cm) lik2(peds[[2]], r = haldane(cM = cm))),
      ),
      rho = tibble(
        x = seq(0, 0.5, length = N),
        lik1 = sapply(x, function(r) lik2(peds[[1]], r)),
        lik2 = sapply(x, function(r) lik2(peds[[2]], r)),
      ),
      mutrate = tibble(
        x = seq(0, 1, length = N),
        lik1 = sapply(x, function(mutr) peds[[1]] |> setMutmod(rate = mutr, update = T) |> lik2()),
        lik2 = sapply(x, function(mutr) peds[[2]] |> setMutmod(rate = mutr, update = T) |> lik2()),
      ),
      freq = tibble(
        x = seq(0.001, 1-0.001-sum(afr[2:3]), length = N),
        lik1 = sapply(x, function(a) peds[[1]] |> setAfreq12(setFr1(a)) |> lik2()),
        lik2 = sapply(x, function(a) peds[[2]] |> setAfreq12(setFr1(a)) |> lik2()),
      )
    )

    dat$LR = dat$lik1/dat$lik2
    LRdat(dat)
  })

  output$pedplot1 = renderPlot({
    if(debug) message("pedplot1")
    niceplot(peds()[[1]], title = "Ped 1", fillcol = "lightblue")
  }, execOnResize = TRUE)

  output$pedplot2 = renderPlot({
    if(debug) message("pedplot2")
    niceplot(peds()[[2]], title = "Ped 2", fillcol = "pink")
  }, execOnResize = TRUE)

  output$ibs = renderText({
    ibs = getIBS(peds()[[1]])
    HTML(sprintf("# shared alleles: <b>%d</b> and <b>%d</b>", ibs[1], ibs[2]))
  })


  output$graphs = renderPlot({
    if(debug) message("graph")

    dat = req(LRdat())
    xlab = names(VARIABLES)[match(input$variable, VARIABLES)]

    p1 = ggplot(dat, aes(x, LR)) +
      theme_classic(base_size = 16) +
      geom_line() +
      labs(x = xlab, y = "LR", title = "LR vs. variable") +
      lims(y = c(0, max(1, dat$LR))) +
      theme(plot.margin = margin(b = 10))


    likdat = rbind(data.frame(x = dat$x, lik = dat$lik1, Ped = "Ped 1"),
                   data.frame(x = dat$x, lik = dat$lik2, Ped = "Ped 2"))

    p2 = ggplot(likdat, aes(x, lik, col = as.factor(Ped))) +
      theme_classic(base_size = 16) +
      geom_line() +
      labs(x = xlab, y = "Likelihood", col = NULL,  title = "Likelihood vs. variable") +
      scale_colour_manual(values = c("blue", "red")) +
      lims(y = c(0, NA)) +
      theme(legend.key.width = unit(1, "cm"),
            plot.margin = margin(t = 10))

    p1 / p2
  }, execOnResize = TRUE)

}


# Run the application
shinyApp(ui = ui, server = server)
