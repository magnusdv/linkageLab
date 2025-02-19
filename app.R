suppressMessages(suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(pedtools)
  library(pedprobr)
  library(forrel)
  library(tibble)
  library(plotly)
}))

CASES = c("Duo : unrelated",
          "Sibs : unrelated",
          "Half-sibs : unrelated",
          "Grandparent : unrelated",
          "Sibs : half-sibs",
          "Grandparent : uncle",
          "Grandparent : half-sibs",
          "Uncle : half-sibs")
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
  "Frequency of '1' allele" = "freq")

debug = FALSE



ui = fluidPage(
  useShinyjs(),

  tags$head(tags$style(HTML("
  .well {margin-bottom: 6px; padding: 15px}
  .inline label{ display: table-cell; padding-right:3px; white-space: nowrap;}
  .inline .form-group { display: table-row}
  .form-control {padding: 3px 1px 3px 8px; height: auto; margin-top:1px; margin-bottom:1px;}
  #variable .form-group {margin-bottom: 0px; white-space: nowrap;}
  #lastrow .form-group {margin-bottom: 7px}
  @media (min-width: 1200px) { /* Adjust for large screens */
      .sidebar { max-width: 280px; }
  }
  "))),

  # Application title
  titlePanel(
    title = HTML("<b>Linkage Lab:</b> Kinship LR with linked markers"),
    windowTitle = "Linkage Lab"),

  div(style = "margin-top:-5px; margin-bottom: 8px;",
      HTML('<b>A pedagogical tool for exploring the effect of linkage in kinship testing.
           Built on the <a href="https://magnusdv.github.io/pedsuite/" target="_blank">pedsuite</a>.
           Source code: <a href="https://github.com/magnusdv/linkageLab" target="_blank">GitHub</a>.</b>')),

  sidebarLayout(
    sidebarPanel(width = 2, style = "min-width:200px;", class = "sidebar",
      h4("Comparison", style = "margin-top: 2px"),
                 selectInput("comp", label = NULL, choices = CASES, selected = "Sibs : halfsibs"),

      fluidRow(
        column(8, h4("Genotypes", style = "margin-top: 6px")),
        column(4, actionButton("sim", "Simulate", style = "float:right; padding:2px 10px; margin-top: 2px",
                               class = "btn-sm btn-warning")),
      ),
      tags$div(class = "inline",
      fluidRow(
        column(6, style = "padding-right: 5px", wellPanel(
          style = "padding: 5px; border-color: silver;",
          p("A", style = "font-size: 18px; text-align: center; margin-bottom:0px"),
          textInput("A1", "M1:", width = "100%"),
          textInput("A2", "M2:", width = "100%")
        )),
        column(6, style = "padding-left: 5px", wellPanel(
          style = "padding: 5px; border-color: silver;",
          p("B", style = "font-size: 18px; text-align: center; margin-bottom:0px"),
          textInput("B1", "M1:", width = "100%"),
          textInput("B2", "M2:", width = "100%")
        )),
      )),

      actionButton("updatePeds", "Update", style = "padding: 2px; margin-top:5px",
                   width = "100%", class = "btn-success"),
      HR,

      h5("Allele frequencies", style = "font-weight:bold;"),
      tags$div(class = "inline",
        fluidRow(
          column(6, numericInput("p1", label = "1:", value = 0.1, min = 0, max = 1, step = 0.05, width = "100%")),
          column(6, numericInput("p2", label = "2:", value = 0.2, min = 0, max = 1, step = 0.05, width = "100%")),
        ),
        fluidRow(
          column(6, numericInput("p3", label = "3:", value = 0.3, min = 0, max = 1, step = 0.05, width = "100%")),
          column(6, numericInput("p4", label = "4:", value = 0.4, min = 0, max = 1, step = 0.05, width = "100%")),
        ),
      ),

      HR,

      tags$div(class = "inline",
        numericInput("cm",    "Dist (cm):", value = 0, step = 1, min = 0, width = "100%"),
        numericInput("rho",   "Rec. rate:", value = 0, step = 0.01, min = 0, max = 1, width = "100%"),
        numericInput("mrate", "Mut. rate:", value = 0, step = 0.01, min = 0, max = 1, width = "100%")
      ),

      HR,
      tags$div(id = "variable", selectInput("variable", "Plot variable", choices = VARIABLES),
               style = "margin-bottom: 10px"),
      fluidRow(id = "lastrow",
        column(6, numericInput("npoints", "Points", value = 11, min = 5, max = 50, step = 1)),
      ),

),
    # Plots
    mainPanel(width = 10,
      fluidRow(
        column(width = 5, align = "left", style = "max-width:400px",
          plotOutput("pedplot1", height = "325px"),
          tags$div(style = "text-align: center", htmlOutput("ibs")),
          plotOutput("pedplot2", height = "325px")
        ),
        column(width = 7, align = "left", style = "max-width:600px",
          plotlyOutput("graph1", height = "325px"),
          tags$div(style = "height:20px"),
          plotlyOutput("graph2", height = "325px"),
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
    gg = fix2(g)
    if(!identical(g, gg)) {
      g = gg
      isolate({
        updateTextInput(session, "A1", value = g[["A", "M1"]])
        updateTextInput(session, "A2", value = g[["A", "M2"]])
        updateTextInput(session, "B1", value = g[["B", "M1"]])
        updateTextInput(session, "B2", value = g[["B", "M2"]])
      })
    }

    g[g == ""] = NA
    g
  })

  emptydat = cbind(M1 = c(A = NA, B = NA),
                   M2 = c(A = NA, B = NA))

  afreq = reactive(c("1" = input$p1, "2" = input$p2,
                     "3" = input$p3, "4" = 1-input$p1-input$p2-input$p3))

  observeEvent(afreq(), {
    if(debug) message("update p4 field")
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
    gt = genodat()
    newpeds = tryCatch(peds() |> setGenos(gt), error = errModal)
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

  observeEvent(input$variable, {
    if(debug) message("disable/enable")

    disab = switch(input$variable, dist =, rho = c("cm", "rho"), mutrate = "mrate", freq = "p1")
    enab = setdiff(c("cm", "rho", "mrate", "p1"), disab)

    for(w in disab) shinyjs::disable(w)
    for(w in enab)  shinyjs::enable(w)
  })

  # Main data reactive
  LRdat = reactive({
    if(debug) message("LRdat")

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

    if(input$variable != "mrate" && (is.na(mrate) || mrate < 0 || mrate > 1))
      err = "The mutation rate must be between 0 and 1"

    if(!input$variable %in% c("dist", "rho") && (is.na(rho) || rho < 0 || rho > 0.5))
      err = "The recombination rate must be between 0 and 0.5"

    if(!is.null(err)) {
      errModal(err)
      req(FALSE)
    }

    peds = tryCatch(error = errModal, warning = errModal,
      lapply(peds(), function(p)
        p |> setAfreq12(afr) |> setMutmod(model = "equal", rate = mrate)))
    req(peds)

    # Wrappers
    lik2 = function(ped, r = rho) likelihood2(ped, 1, 2, rho = r)
    setFr1 = function(a)  c("1" = a, afr[2:3], "4" = 1 - a - sum(afr[2:3]))

    # Compute likelihoods
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
        x = seq(0.0001, 1-0.0001-sum(afr[2:3]), length = N),
        lik1 = sapply(x, function(a) peds[[1]] |> setAfreq12(setFr1(a)) |> lik2()),
        lik2 = sapply(x, function(a) peds[[2]] |> setAfreq12(setFr1(a)) |> lik2()),
      )
    )

    dat$LR = dat$lik1/dat$lik2
    dat
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


  output$graph1 = renderPlotly({
    if(debug) message("graph1")
    dat = req(LRdat())
    xlab = names(VARIABLES)[match(input$variable, VARIABLES)]

    # Hack to avoid max tick at 0.9.
    # (When max LR is around 1.05, default ticks are 0,.3,.6,.9,1.2; 1.2 not shown)
    brks = function(lims) {
      b = scales::extended_breaks()(lims)
      if(max(b) < 1.21)
        b = c(0, 0.25, 0.5, 0.75, 1)
      b
    }

    p1 = ggplot(dat, aes(x, LR)) +
      theme_classic(base_size = 14) +
      geom_line() + geom_point(size = 0.8) +
      labs(x = xlab, y = "LR") +
      scale_y_continuous(limits = c(0, max(1, max(dat$LR))),
                         breaks = brks, expand = expansion(mult = c(0, 0.08)))

    ggplotly(p1, tooltip = "y") |> plotly::config(displayModeBar = FALSE)
  })

  output$graph2 = renderPlotly({
    if(debug) message("graph1")
    dat = req(LRdat())
    likdat = rbind(data.frame(x = dat$x, lik = dat$lik1, Ped = "Ped 1"),
                   data.frame(x = dat$x, lik = dat$lik2, Ped = "Ped 2"))
    likdat$text = paste("Likelihood:", likdat$x)

    xlab = names(VARIABLES)[match(input$variable, VARIABLES)]

    p2 = ggplot(likdat, aes(x, y = lik, col = as.factor(Ped), group = Ped,
                            text = paste0("x: ", sprintf("%.3g", x), "<br>y: ", sprintf("%.4g", lik)))) +
      theme_classic(base_size = 14) +
      geom_line() + geom_point(size = 0.8) +
      labs(x = xlab, y = "Likelihood", col = NULL) +
      scale_y_continuous(limits = c(0, NA),
                         expand = expansion(mult = c(0, 0.08))) +
      scale_colour_manual(values = c("blue", "red"))

    ggplotly(p2, tooltip = "text") |>
      plotly::config(displayModeBar = FALSE) |>
      plotly::layout(
        legend = list(orientation = "h", x = 1, xanchor = "right", y = 1.09,
                      yanchor = "center", bgcolor = "transparent"),
        showlegend = TRUE)
  })

}


# Run the application
shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
