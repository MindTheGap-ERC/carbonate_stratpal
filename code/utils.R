last_occ = function(t_ext, rate, adm, niche = trivial_niche, gc = trivial_gradient, niche_domain = "time"){
  if (is.function(rate)){
    occ = p3_var_rate(x = rate,
                      from = admtools::min_time(adm),
                      to = admtools::max_time(adm),
                      fmax = 100)
    occ = occ[occ < t_ext]
  } else {
    occ = p3(rate = rate,
             from = admtools::min_time(adm),
             to = admtools::max_time(adm))
    occ = occ[occ < t_ext]
  }
  if (length(occ) <1){return(numeric())}
  if (! niche_domain %in% c("time", "strat")){stop("unknows niche domain, use either \"time\" or \"strat\".")}
  if (niche_domain == "time"){
    occ_strat = occ |>
      apply_niche(niche_def = niche, gc = gc) |>
      time_to_strat(adm, destructive = TRUE)
  }
  if (niche_domain == "strat"){
    occ_strat = occ |>
      time_to_strat(adm, destructive = TRUE) |>
      apply_niche(niche_def = niche, gc = gc)
  }
  # select preserved occurrences
  occ_strat = occ_strat[!is.na(occ_strat)]
  if (length(occ_strat) < 1){return(numeric())}
  highest_occ = max(occ_strat)
  last_occ = highest_occ |>
    strat_to_time(adm)
  return(c("t" = last_occ, "h" = highest_occ))
}

range_offset = function(t_ext, rate, adm, niche, gc, niche_domain = "time"){
  lo = last_occ(t_ext = t_ext, rate = rate, adm = adm, niche = niche, gc = gc, niche_domain = "time")
  t_ro = t_ext - unname(lo["t"])
  h_ext = time_to_strat(t_ext, adm, destructive = FALSE)
  h_ro = h_ext - unname(lo["h"])
  return(c("t" = t_ro, "h" = h_ro))
}


library(admtools)
adm = tp_to_adm(t = scenarioA$t_myr, h = scenarioA$h_m[,1])
plot(adm)
niche1 = function(x) return(rep(1, length(x)))
gc1 = function(x) return(rep(1, length(x)))
t_ext = 1
rate = 500

lo = last_occ(t_ext, rate = rate, adm = adm, niche = niche1, gc = gc1)
lo

range_offset(t_ext, rate = rate, adm = adm, niche = niche1, gc = gc1)

trivial_niche = function(x) rep(1, length(x))
trivial_gradient = function(x) rep(1, length(x))
perfect_preservation = function(x) rep(1, length(x))
indestructible = function(x) rep(1, length(x))
lo = last_occ(t_ext, rate = rate, adm = adm)
lo

gradient_from_data(x) = function(x){
  if (is.data.frame(x)){
    return(approxfun(x = df[,1], y = df[,2], rule = 2))
  }
  if(is.list(x)){
    return(approxfun(x = x[[1]], y = x[[2]], rule = 2))
  }
  stop("Input must be a list or a data frame")
}
