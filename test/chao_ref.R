suppressPackageStartupMessages(library(iNEXT))
cat("iNEXT_version\t", as.character(packageVersion("iNEXT")), "\n", sep = "")

# Three synthetic incidence datasets of increasing coverage.
# Format for datatype="incidence_freq": first element is T (sampling units),
# the rest are per-species incidence counts.
datasets <- list(
  low  = c(40, rep(1, 30), rep(2, 12), rep(3, 6), rep(5, 3), 9, 14, 22),
  mid  = c(120, rep(1, 45), rep(2, 30), rep(3, 20), rep(6, 14), rep(11, 8),
           rep(25, 5), 60, 88),
  high = c(500, rep(1, 20), rep(2, 25), rep(4, 40), rep(9, 50), rep(30, 35),
           rep(120, 18), 300, 410, 455)
)

cat("dataset\tT\tS_obs\tQ1\tQ2\tS_chao2\tSE\tci_lower\tci_upper\n")
for (nm in names(datasets)) {
  x <- datasets[[nm]]
  T <- x[1]
  y <- x[-1]
  res <- ChaoRichness(x, datatype = "incidence_freq")
  cat(sprintf("%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n",
              nm, T, sum(y > 0), sum(y == 1), sum(y == 2),
              res$Estimator, res$Est_s.e.,
              res$`95% Lower`, res$`95% Upper`))
}
