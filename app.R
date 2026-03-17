suppressMessages(suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(pedtools)
  library(pedprobr)
  library(forrel)
  library(tibble)
  library(plotly)
}))

# Deployment commands:
# remotes::install_github("magnusdv/pedtools")

debug = F
prefill = F

# Prefill values for paper example
prefill_vals = list(
  amat = rbind(A = c(1,1,1,2), B = c(1,1,1,3)),
  fr1 = c(norSTR::norwayDB$vWA["16"], 0, 0),
  fr2 = norSTR::norwayDB$D12S391[c("18", "21", "22")]
)


IDS = c("A", "B")

builtinPeds = list(
  paternity = nuclearPed(father = IDS[1], children = IDS[2]),
  unrelated = pedtools::singletons(IDS),
  sibs = nuclearPed(children = IDS),
  `half-sibs` = halfSibPed() |> relabel(old = 4:5, new = IDS),
  avuncular = avuncularPed() |> relabel(old = c(3,6), new = IDS),
  grandparent = linearPed(2) |> relabel(old = c(1,5), new = IDS)
)

CASES = c("Paternity : unrelated",
          "Paternity : sibs",
          "Sibs : unrelated",
          "Sibs : half-sibs",
          "Half-sibs : unrelated",
          "Half-sibs : grandparent",
          "Half-sibs : avuncular",
          "Grandparent : unrelated",
          "Grandparent : avuncular",
          "Avuncular : unrelated")

names(CASES) = CASES

PEDS = lapply(CASES, function(case) {
  rels = strsplit(case, " : ")[[1]]
  setNames(builtinPeds[tolower(rels)], c("Ped 1", "Ped 2"))
})

VARIABLES = c("Distance (cM)" = "dist",
              "Recombination rate" = "rho",
              "Mutation rate" = "mutrate")


# UI ----------------------------------------------------------------------


ui = fluidPage(
  useShinyjs(),
  tags$head(includeHTML("GA.html")),
  tags$head(tags$style(HTML("
  .well{margin-bottom:6px;padding:15px}

  .selectize-dropdown .option {padding:2px 10px; line-height:1.2; white-space:nowrap}
  .selectize-dropdown-content {max-height: none;}

  .inline .form-group{display:table-row;margin-bottom:0}
  .inline label{display:table-cell;padding-right:3px;white-space:nowrap}
  .inline .form-control{padding:3px 1px 3px 8px;height:auto;margin:1px 0}

  #variable .form-group{margin-bottom:0}
  #lastrow .form-group{margin-bottom:7px}

  @media (min-width:1200px){.sidebar{max-width:280px}}

  #updatePeds{position:relative;transition:box-shadow .2s ease,filter .2s ease,opacity .2s ease}
  #updatePeds:disabled{opacity:.55;filter:saturate(.6);box-shadow:none}
  #updatePeds.dirty{animation:updatePulse 1.3s ease-in-out infinite;box-shadow:0 0 0 .25rem rgba(91,192,222,.30)}
  #updatePeds.dirty::after{content:'';position:absolute;top:8px;right:10px;width:10px;height:10px;border-radius:50%;background:rgb(220,53,69);box-shadow:0 0 0 4px rgba(220,53,69,.18)}
  @keyframes updatePulse{0%{box-shadow:0 0 0 .10rem rgba(91,192,222,.12)}50%{box-shadow:0 0 0 .35rem rgba(91,192,222,.30)}100%{box-shadow:0 0 0 .10rem rgba(91,192,222,.12)}}

  .cellgrid{display:grid;grid-template-columns:20px repeat(4,1fr);column-gap:1px;align-items:center}
  .cellgrid .colname{font-size:13px;font-weight:600;text-align:center;margin:0;padding:0}
  .cellgrid .rowname{font-size:13px;font-weight:600;text-align:left;margin:0;padding:0 2px 0 0}
  .cellgrid .form-group{margin:0}
  .cellgrid .form-control{padding:2px 3px;height:auto;margin:0;text-align:center}
  .afreq-cell .form-control{font-style:italic}

  .cellgrid input[type=number]{-moz-appearance:textfield}
  .cellgrid input[type=number]::-webkit-outer-spin-button,
  .cellgrid input[type=number]::-webkit-inner-spin-button{-webkit-appearance:none;margin:0}
  "))),

  tags$script(HTML("
  $(document).on('shown.bs.modal', '#shiny-modal', function() {
    $(this).find('.modal-footer button').trigger('focus');
  });
  ")),

  # Application title
  titlePanel(
    title = HTML("<b>Linkage Lab:</b> Kinship testing with linked markers"),
    windowTitle = "Linkage Lab"),

  div(style = "margin-top:-5px; margin-bottom: 8px;",
      HTML('<b>A tool for exploring how linkage affects the LR in kinship testing.
           Built on the <a href="https://magnusdv.github.io/pedsuite/" target="_blank">pedsuite</a>.
           Source code: <a href="https://github.com/magnusdv/linkageLab" target="_blank">GitHub</a>.</b>')),

  sidebarLayout(
    sidebarPanel(width = 2, style = "min-width:200px;", class = "sidebar",
      h4("Comparison", style = "margin-top: 2px"),
      selectInput("comp", label = NULL, choices = CASES, selected = "Sibs : half-sibs"),

      fluidRow(
        column(8, h4("Genotypes", style = "margin-top: 6px")),
        column(4, actionButton("sim", "Simulate", style = "float:right; padding:2px 10px; margin-top: 2px",
                               class = "btn-sm btn-warning")),
      ),

      tags$div(class = "cellgrid",
        tags$div(), tags$div(class = "colname", "A1"), tags$div(class = "colname", "A2"),
                    tags$div(class = "colname", "B1"), tags$div(class = "colname", "B2"),
        alleleRow("M1", ids = paste0("m1_", c("A1", "A2", "B1", "B2")), vals = rep("", 4)),
        alleleRow("M2", ids = paste0("m2_", c("A1", "A2", "B1", "B2")), vals = rep("", 4))
      ),

      h4("Allele frequencies", style = "margin-top:14px; margin-bottom:6px;"),
      tags$div(class = "cellgrid",
        tags$div(), tags$div(class = "colname", "1"), tags$div(class = "colname", "2"),
                    tags$div(class = "colname", "3"), tags$div(class = "colname", "4"),
        afRow("M1", ids = paste0("m1_p", 1:4), vals = c(0.1, 0.2, 0.3, 0.4)),
        afRow("M2", ids = paste0("m2_p", 1:4), vals = c(0.1, 0.2, 0.3, 0.4))
      ),

      actionButton("updatePeds", "Update", style = "padding: 2px; margin-top:10px",
                   width = "100%", class = "btn-default"),
      HR,


      tags$div(class = "inline", style = "margin-bottom: 5px;",
        numericInput("cm",    "Dist (cM):", value = 0, step = 1, min = 0, width = "100%"),
        numericInput("rho",   "Rec. rate:", value = 0, step = 0.01, min = 0, max = 0.5, width = "100%"),
        numericInput("mrate", "Mut. rate:", value = 0, step = 0.01, min = 0, max = 1, width = "100%"),
      ),

      tags$div(style = "margin-left: 5px",
        radioButtons("mapfun", NULL, choices = c("Haldane", "Kosambi"), inline = TRUE, width = "100%",
                     selected = "Kosambi")
      ),

      HR,
      tags$div(id = "variable", selectInput("variable", "Plot variable", choices = VARIABLES),
               style = "margin-bottom: 10px"),
      fluidRow(id = "lastrow",
        column(6, numericInput("npoints", "Points", value = 11, min = 5, max = 50, step = 1)),
        column(6, numericInput("xmax",    "Max", value = 200))
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



# Server function ---------------------------------------------------------


server = function(input, output, session) {

  session$onSessionEnded(stopApp)

  session$onFlushed(function() {
    if(!prefill) return()
    if(debug) message("prefill")
    amat = prefill_vals$amat
    for(i in 1:2) for(j in 1:2) {
      updateTextInput(session, sprintf("m%d_A%d", i, j), value = amat[[1, 2*(i-1) + j]])
      updateTextInput(session, sprintf("m%d_B%d", i, j), value = amat[[2, 2*(i-1) + j]])
    }
    for(k in 1:3) {
      updateNumericInput(session, paste0("m1_p", k), value = as.numeric(prefill_vals$fr1[k]))
      updateNumericInput(session, paste0("m2_p", k), value = as.numeric(prefill_vals$fr2[k]))
    }
    shinyjs::delay(200, shinyjs::click("updatePeds"))
  }, once = TRUE)

  amat = reactive({
    if(debug) message("allele matrix")

    cbind(M1.1 = c(A = input$m1_A1, B = input$m1_B1),
          M1.2 = c(A = input$m1_A2, B = input$m1_B2),
          M2.1 = c(A = input$m2_A1, B = input$m2_B1),
          M2.2 = c(A = input$m2_A2, B = input$m2_B2))
  })
  amatEmpty = cbind(M1.1 = c(A = NA, B = NA),
                    M1.2 = c(A = NA, B = NA),
                    M2.1 = c(A = NA, B = NA),
                    M2.2 = c(A = NA, B = NA))

  # Allele frequencies
  afreq = reactiveValues(afr1 = c(0.1, 0.2, 0.3, 0.4) |> setNames(1:4),
                         afr2 = c(0.1, 0.2, 0.3, 0.4) |> setNames(1:4))

  # Disable and auto-update p4 fields
  observe({
    shinyjs::disable("m1_p4")
    shinyjs::disable("m2_p4")
  })

  observeEvent(list(input$m1_p1, input$m1_p2, input$m1_p3), {
    if(debug) message("update m1-p4")
    s = sum(c(input$m1_p1, input$m1_p2, input$m1_p3), na.rm = TRUE)
    updateNumericInput(session, "m1_p4", value = 1 - s) #max(0, min(1, 1 - s)))
  }, ignoreInit = FALSE)

  observeEvent(list(input$m2_p1, input$m2_p2, input$m2_p3), {
    if(debug) message("update m2-p4")
    s = sum(c(input$m2_p1, input$m2_p2, input$m2_p3), na.rm = TRUE)
    updateNumericInput(session, "m2_p4", value = 1 - s) #max(0, min(1, 1 - s)))
  }, ignoreInit = FALSE)

  peds = reactiveVal(NULL)
  peds(isolate(PEDS[[input$comp]] |> setAmat(amatEmpty)))

  observeEvent(input$comp, {
    if(debug) message("newcomp")
    newpeds = tryCatch(
      PEDS[[input$comp]] |> setAmat(amat()),
      error = function(e) {
        errModal(e)
        PEDS[[input$comp]] |> setAmat(amatEmpty)
    })
    peds(newpeds)
  })


  # Update pedigrees ------------------------------------------------------

  dirty = reactiveValues(alleles = FALSE, afreqs = FALSE)

  observeEvent(input$updatePeds, {

    # Update genotypes
    if(dirty$alleles) {
      if(debug) message("update genotypes")
      newpeds = tryCatch(peds() |> setAmat(amat()), error = errModal)
      peds(req(newpeds))
    }

    # Update freqs
    if(dirty$afreqs) {
      if(debug) message("update afreqs")
      afreq$afr1 = c("1" = input$m1_p1, "2" = input$m1_p2, "3" = input$m1_p3,
             "4" = 1 - input$m1_p1 - input$m1_p2 - input$m1_p3)
      afreq$afr2 = c("1" = input$m2_p1, "2" = input$m2_p2, "3" = input$m2_p3,
             "4" = 1 - input$m2_p1 - input$m2_p2 - input$m2_p3)
    }
    dirty$alleles = dirty$afreqs = FALSE
  })


  alleleFields = c("m1_A1","m1_B1","m1_A2","m1_B2","m2_A1","m2_B1","m2_A2","m2_B2")
  afreqFields = c("m1_p1","m1_p2","m1_p3","m2_p1","m2_p2","m2_p3")

  observeEvent(lapply(alleleFields, function(id) input[[id]]), {dirty$alleles = TRUE}, ignoreInit = TRUE)
  observeEvent(lapply(afreqFields, function(id) input[[id]]), {dirty$afreqs = TRUE}, ignoreInit = TRUE)

  observeEvent(list(dirty$alleles, dirty$afreqs), {
    dirt = dirty$alleles || dirty$afreqs
    if(debug) message("dirty:", dirt)
    if(dirt) {
      shinyjs::enable("updatePeds")
      shinyjs::removeClass("updatePeds", "btn-default")
      shinyjs::addClass("updatePeds", "btn-info dirty")
    } else {
      shinyjs::disable("updatePeds")
      shinyjs::removeClass("updatePeds", "btn-info dirty")
      shinyjs::addClass("updatePeds", "btn-default")
    }
  }, ignoreInit = FALSE)

  observeEvent(input$sim, {
    if(debug) message("sim")

    amatSim = peds()[[sample(1:2, size = 1)]] |>
      markerSim(N = 2, ids = IDS, alleles = 1:4, verbose = FALSE) |>
      getAlleles(IDS)

    colnames(amatSim) = paste0("M", colnames(amatSim))

    for(i in 1:2) for(j in 1:2) {
      updateTextInput(session, sprintf("m%d_A%d", i, j), value = amatSim[[1, 2*(i-1) + j]])
      updateTextInput(session, sprintf("m%d_B%d", i, j), value = amatSim[[2, 2*(i-1) + j]])
    }
    shinyjs::delay(200, shinyjs::click("updatePeds"))
  })

  # Map function reactive
  mapfun = reactive(switch(input$mapfun, Haldane = haldane, Kosambi = kosambi))

  src = reactiveVal(NULL)

  # Update rho when map function or cm changes
  observeEvent(list(mapfun(), input$cm), {
    if(identical(src(), "rho")) { src(NULL); return() }
    rho = suppressWarnings(mapfun()(cM = input$cm))
    if(!isTRUE(all.equal(input$rho, rho, tolerance = 1e-5))) {
      if(debug) message("cm -> rho = ", rho)
      src("cm")
      updateNumericInput(session, "rho", value = rho)
    }
  })

  # Update cm when rho changes
  observeEvent(input$rho, {
    if(identical(src(), "cm")) { src(NULL); return() }
    cm = suppressWarnings(mapfun()(rho = input$rho))
    if(!isTRUE(all.equal(input$cm, cm, tolerance = 1e-5))) {
      if(debug) message("rho -> cm = ", cm)
      src("rho")
      updateNumericInput(session, "cm", value = cm)
    }
  }, ignoreInit = TRUE)

  # Max when variable changes
  observeEvent(input$variable, {
    maxval = switch(input$variable,
      dist = 200,
      rho = 0.5,
      mutrate = 1
    )
    freezeReactiveValue(input, "xmax")
    updateNumericInput(session, "xmax", value = maxval)
  })

  # Main data reactive
  LRdat = reactive({
    if(debug) message("LRdat")

    variable = input$variable
    N = req(input$npoints)
    xmax = req(input$xmax)
    afr1 = afreq$afr1
    afr2 = afreq$afr2
    mrate0 = input$mrate
    cm0 = input$cm
    mapfun = mapfun()
    rho0 = mapfun(cM = cm0) #isolate(input$rho)
    cm0 = min(cm0, 200)

    # Include user point?
    x0 = switch(variable, dist = cm0, rho = rho0, mutrate = mrate0)
    hasUser = !is.na(x0)
    addUser = function(v) if(hasUser) unique.default(c(x0, v)) else v

    # Checks
    err = NULL
    if(N < 3 || N > 100)
      err = "The number of points must be between 3 and 100"

    if(xmax <= 0)
      err = "The max value must be positive"

    if((variable == "rho" && xmax > 0.5) || (variable == "mutrate" && xmax > 1))
      err = "Max value is too large"

    if(!is.na(cm0) && cm0 < 0)
      err = "Marker distance (cm) cannot be negative"

    if(!is.na(rho0) && (rho0 < 0 || rho0 > 0.5))
      err = "The recombination rate must be between 0 and 0.5"

    if(any(is.na(afr1) | afr1 < 0 | afr1 > 1))
      err = "Please check that allele frequencies for M1 are between 0 and 1"

    if(any(is.na(afr2) | afr2 < 0 | afr2 > 1))
      err = "Please check that allele frequencies for M2 are between 0 and 1"

    if(variable != "mutrate" && (is.na(mrate0) || mrate0 < 0 || mrate0 > 1))
      err = "The mutation rate must be between 0 and 1"

    if(variable == "mutrate" && is.na(rho0))
      err = "Please indicate a cM distance or recombination rate"

    if(!is.null(err)) {
      errModal(err)
      req(FALSE)
    }

    # Set up pedigrees
    tryCatch({
      peds = peds()
      if(is.na(mrate0)) mrate0 = 0
      peds[[1]] = peds[[1]] |> setAfreq12(afr1, afr2) |> setMutmod(model = "equal", rate = mrate0)
      peds[[2]] = peds[[2]] |> setAfreq12(afr1, afr2) |> setMutmod(model = "equal", rate = mrate0)
    }, error = errModal, warning = errModal)

    req(peds)

    # Wrapper
    likfun = function(ped, rho) likelihood2(ped, marker1 = 1, marker2 = 2, rho = rho)

    # Compute likelihoods
    dat = switch(variable,
      dist = tibble(
        x = seq(0, xmax, length = N) |> addUser(),
        lik1 = sapply(x, function(cm) likfun(peds[[1]], rho = mapfun(cM = cm))),
        lik2 = sapply(x, function(cm) likfun(peds[[2]], rho = mapfun(cM = cm)))
      ),
      rho = tibble(
        x = seq(0, xmax, length = N) |> addUser(),
        lik1 = sapply(x, function(r) likfun(peds[[1]], rho = r)),
        lik2 = sapply(x, function(r) likfun(peds[[2]], rho = r)),
      ),
      mutrate = tibble(
        x = seq(0, xmax, length = N) |> addUser(),
        lik1 = sapply(x, function(mrate) peds[[1]] |> setMutmod(rate = mrate, update = TRUE) |> likfun(rho = rho0)),
        lik2 = sapply(x, function(mrate) peds[[2]] |> setMutmod(rate = mrate, update = TRUE) |> likfun(rho = rho0)),
      )
    )

    dat$LR = dat$lik1/dat$lik2
    dat$user = FALSE
    if(hasUser)
      dat$user[1] = TRUE

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
    if(debug) message("ibs")
    ibs = getIBS(peds()[[1]])
    HTML(sprintf("# shared alleles: <b>%d</b> and <b>%d</b>", ibs[1], ibs[2]))
  })


  output$graph1 = renderPlotly({
    if(debug) message("graph1")

    datAll = req(LRdat())
    dat0 = datAll[datAll$user, , drop = FALSE]
    dat1 = datAll[!datAll$user, , drop = FALSE]

    xlab = names(VARIABLES)[match(input$variable, VARIABLES)]

    # Hack to avoid max tick at 0.9.
    # (When max LR is around 1.05, default ticks are 0,.3,.6,.9,1.2; 1.2 not shown)
    brks = function(lims) {
      b = scales::extended_breaks()(lims)
      if(max(b) < 1.21)
        b = c(0, 0.25, 0.5, 0.75, 1)
      b
    }

    p1 = ggplot(datAll, aes(x, LR)) +
      theme_classic(base_size = 13) +
      geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5, col = "gray") +
      geom_line(aes(group = 1)) +
      geom_point(data = dat1, size = 0.8) +
      geom_point(data = dat0, size = 2, shape = 21, fill = "yellow", col = "black") +
      labs(x = xlab, y = "LR") +
      scale_y_continuous(limits = c(0, max(1, max(datAll$LR))),
                         breaks = brks, expand = expansion(mult = c(0, 0.08)))

    ggplotly(p1, tooltip = c("x", "y")) |>
      plotly::style(hovertemplate = "x = %{x:.3g}<br>LR = %{y:.4g}<extra></extra>") |>
      plotly::config(displayModeBar = FALSE)
  })

  output$graph2 = renderPlotly({
    if(debug) message("graph2")

    dat = req(LRdat())
    likdat = rbind(data.frame(x = dat$x, lik = dat$lik1, Ped = "Ped 1", user = dat$user),
                   data.frame(x = dat$x, lik = dat$lik2, Ped = "Ped 2", user = dat$user))

    # Extract user points for separate plotting
    likdat0 = likdat[likdat$user, , drop = FALSE]
    likdat1 = likdat[!likdat$user, , drop = FALSE]

    xlab = names(VARIABLES)[match(input$variable, VARIABLES)]

    p2 = ggplot(likdat, aes(x, y = lik, col = as.factor(Ped), group = Ped)) +
      theme_classic(base_size = 13) +
      geom_line() +
      geom_point(data = likdat1, size = 0.8) +
      geom_point(data = likdat0, size = 2, shape = 21, fill = "yellow") +
      labs(x = xlab, y = "Likelihood", col = NULL) +
      scale_y_continuous(limits = c(0, NA),
                         expand = expansion(mult = c(0, 0.08))) +
      scale_colour_manual(values = c("#1F77B4", "#D62780"))

    ggplotly(p2, tooltip = c("x", "y")) |>
      plotly::config(displayModeBar = FALSE) |>
      plotly::style(hovertemplate = "x = %{x:.3g}<br>y = %{y:.4g}<extra></extra>") |>
      plotly::layout(
        legend = list(orientation = "h", x = 1, xanchor = "right", y = 1.09,
                      yanchor = "center", bgcolor = "transparent"),
        showlegend = TRUE)
  })

}


# Run the application
suppressMessages(suppressPackageStartupMessages(
  shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
))
