library(pedtools)
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

fix2 = function(smat) {
  mode(smat) = "character"

  # Entries with exactly two characters
  ch2 = !is.na(smat) & nchar(smat) == 2
  if(!any(ch2))
    return(smat)

  smat[ch2] = sapply(smat[ch2], function(s) {
    ss = strsplit(s, "")[[1]]
    if("/" %in% ss)
      ss[ss == "/"] = "-"
    paste0(ss[1], "/", ss[2])
  })
  smat
}

# Set genotypes for both peds
setGenos = function(peds, genodat) {
  lapply(peds, function(p)
    setMarkers(p, alleleMatrix = genodat, sep = "/",
               locusAttributes = list(alleles = 1:4)))
}

# Set frequencies for marker 1 and 2
setAfreq12 = function(ped, afr)
  ped |> setAfreq(marker = 1, afreq = afr) |> setAfreq(marker = 2, afreq = afr)


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

HR = hr(style = "border-top: 1px solid #BBBBBB; margin: 14px 0px 12px 0px")


builtinPeds = list(
  duo = nuclearPed(father = IDS[1], children = IDS[2]),
  unrelated = list(pedtools::singleton(IDS[1]), pedtools::singleton(IDS[2])),
  sibs = nuclearPed(children = IDS),
  `half-sibs` = halfSibPed() |> relabel(old = 4:5, new = IDS),
  uncle = avuncularPed() |> relabel(old = c(3,6), new = IDS),
  grandparent = linearPed(2) |> relabel(old = c(1,5), new = IDS))
