
IDS = c("A", "B")

stop2 = function (...) {
  a = lapply(list(...), toString)
  a = append(a, list(call. = FALSE))
  do.call(stop, a)
}

errModal = function(..., html = FALSE) {#print("err")
  args = list(...)
  if(length(args) == 1 && inherits(args[[1]], "condition"))
    mess = conditionMessage(args[[1]])
  else
    mess = paste(lapply(args, toString), collapse = "")

  if(html)
    mess = HTML(mess)
  showModal(modalDialog(mess, easyClose = TRUE))
}


addMarkers = function(x) {
  if(hasMarkers(x))
    return(x)
  x |> addMarker(name = "M1", alleles = 1:4) |> addMarker(name = "M2", alleles = 1:4)
}


# Set allele matrix (i.e. genotypes) for both peds
setAmat = function(peds, amat) {
  lapply(peds, function(p)
    setMarkers(p, alleleMatrix = amat, locusAttributes = list(alleles = 1:4)))
}

# Set frequencies for marker 1 and 2
setAfreq12 = function(ped, afr1, afr2) {
  ped |> setAfreq(marker = 1, afreq = afr1) |> setAfreq(marker = 2, afreq = afr2)
}

ibsState = function(gt1, gt2) {
  max(sum(gt1[1] == gt2[1], gt1[2] == gt2[2], na.rm = TRUE),
      sum(gt1[1] == gt2[2], gt1[2] == gt2[1], na.rm = TRUE))
}

getIBS = function(ped) {
  if(!is.ped(ped))
    ped = ped[[1]]

  als = getAlleles(ped, ids = IDS)
  c(ibsState(als["A", 1:2], als["B", 1:2]),
    ibsState(als["A", 3:4], als["B", 3:4]))
}

niceplot = function(ped, title = NULL, fillcol = NULL, cex = 1.4, addbox = TRUE, margins = 2) {
  fill = if(!is.null(fillcol)) setNames(list(IDS), fillcol) else NULL
  plot(ped, title = title, cex = cex, margins = margins, fill = fill, labs = NULL,
       textInside = IDS, marker = 1:2, showEmpty = IDS, keep.par = TRUE)
  if(addbox) box("outer")
}


# UI utils ------------------------------------------------------------------------------------


HR = hr(style = "border-top: 1px solid #BBBBBB; margin: 14px 0px 12px 0px")

alleleInput = function(id, value = "") {
  textInput(id, label = NULL, value = value, width = "100%")
}

alleleRow = function(marker, ids, vals) {
  tagList(
    tags$div(class = "rowname", marker),
    tags$div(class = "allele-cell", alleleInput(ids[1], vals[1])),
    tags$div(class = "allele-cell", alleleInput(ids[2], vals[2])),
    tags$div(class = "allele-cell", alleleInput(ids[3], vals[3])),
    tags$div(class = "allele-cell", alleleInput(ids[4], vals[4]))
  )
}

afInput = function(id, value) {
  numericInput(id, label = NULL, value = value, min = 0, max = 1, step = "any", width = "100%")
}

afRow = function(marker, ids, vals) {
  tagList(
    tags$div(class = "rowname", marker),
    tags$div(class = "afreq-cell", afInput(ids[1], vals[1])),
    tags$div(class = "afreq-cell", afInput(ids[2], vals[2])),
    tags$div(class = "afreq-cell", afInput(ids[3], vals[3])),
    tags$div(class = "afreq-cell", afInput(ids[4], vals[4]))
  )
}
