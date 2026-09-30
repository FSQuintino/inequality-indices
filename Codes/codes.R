# NSII - Unbiased estimation under the Gamma distribution
# Reproduces Section 5 (applications)
# Base R only (R >= 4.0). 
#
# Outputs (in output/): fig_gamma_fit.pdf, gof_tests.csv,
#                       table5_gdp.csv, table6_gdp.csv

B_GOF  <- 10000L  # GOF bootstrap replicates
B_BOOT <-  5000L  # CI bootstrap resamples
SEED   <-  2026L
CONF   <-  0.95

dir.create("output", showWarnings = FALSE)

# ---- Data -------------------------------------------------------------------
# GDP per capita, PPP (constant 2021 int. $, thousands), 2023, World Bank
# NY.GDP.PCAP.PP.KD via Our World in Data. Same data in data/gdp_data.csv.

countries <- c(
  "Haiti", "Honduras", "Nicaragua", "Bolivia", "Jamaica",
  "El Salvador", "Guatemala", "Belize", "Ecuador", "Peru",
  "Paraguay", "Grenada", "Dominica", "St. Vincent and the Grenadines",
  "Colombia", "Brazil", "Suriname", "Barbados", "Mexico",
  "Dominican Republic", "Saint Lucia", "Costa Rica", "Argentina",
  "Antigua and Barbuda", "Chile", "St. Kitts and Nevis",
  "Uruguay", "Trinidad and Tobago", "Bahamas", "Panama",
  "Puerto Rico", "Guyana", "Canada", "United States")

gdp <- c(
   2.956,  6.468,  7.487,  9.844, 10.291, 11.404, 12.389, 12.455,
  14.472, 15.294, 15.783, 16.946, 17.420, 18.335, 18.692, 19.018,
  19.044, 19.224, 22.143, 23.088, 23.403, 26.293, 27.105, 28.967,
  29.463, 30.409, 31.019, 31.706, 33.106, 35.864, 42.995, 49.315,
  55.919, 74.578)
names(gdp) <- countries
# gdp <- setNames(read.csv("data/gdp_data.csv")$gdp, countries)  # alternative

n <- length(gdp)
summary(gdp)

# ---- Gamma MLE (profile score equation in alpha) ----------------------------
gamma_mle <- function(x) {
  s  <- log(mean(x)) - mean(log(x))
  f  <- function(a) log(a) - digamma(a) - s
  lo <- 1e-6; hi <- 1
  while (f(hi) > 0) hi <- hi * 2
  while (f(lo) < 0) lo <- lo / 2
  a <- uniroot(f, c(lo, hi), tol = .Machine$double.eps^0.75)$root
  c(alpha = a, lambda = a / mean(x))
}

pars <- gamma_mle(gdp)
al <- pars[["alpha"]]; la <- pars[["lambda"]]
print(pars)

# ---- Figure 1 ---------------------------------------------------------------
plot_fit <- function() {
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
  xs <- seq(0, 80, length.out = 400)
  hist(gdp, probability = TRUE, breaks = 10, main = "", xlab = "gdp",
       cex.axis = 1.3, cex.lab = 1.3)
  lines(xs, dgamma(xs, al, la), col = "red", lwd = 2)
  legend("topright", "Model fitted", lty = 1, col = "red", cex = 1.2)
  plot(ecdf(gdp), main = "ECDF", xlab = "x", ylab = "F(x)", lwd = 2,
       cex.main = 1.6, cex.axis = 1.3, cex.lab = 1.3)
  lines(xs, pgamma(xs, al, la), col = "red", lwd = 2)
  legend("topleft", "Model fitted", lty = 1, col = "red", cex = 1.2)
  res <- qnorm(pgamma(sort(gdp), al, la))
  qqnorm(res, main = "QQ Plot - Gamma Fit", cex.main = 1.6,
         cex.axis = 1.3, cex.lab = 1.3)
  qqline(res, col = "red", lwd = 2)
}
pdf("output/fig_gamma_fit.pdf", width = 12, height = 4); plot_fit(); invisible(dev.off())

# ---- Goodness of fit (parametric bootstrap, both parameters re-estimated) ---
gof_stats <- function(x, alpha, lambda) {
  u <- sort(pgamma(x, shape = alpha, rate = lambda))
  n <- length(u); i <- seq_len(n)
  c(KS  = max(pmax(i/n - u, u - (i - 1)/n)),
    CvM = sum((u - (2*i - 1)/(2*n))^2) + 1/(12*n),
    AD  = -n - mean((2*i - 1) * (log(u) + log1p(-rev(u)))))
}

obs_gof <- gof_stats(gdp, al, la)
set.seed(SEED)
boot_gof <- replicate(B_GOF, {
  xb <- rgamma(n, shape = al, rate = la)
  pb <- gamma_mle(xb)
  gof_stats(xb, pb[["alpha"]], pb[["lambda"]])
})
gof_tab <- data.frame(test = names(obs_gof), statistic = obs_gof,
                      p_value = rowMeans(boot_gof >= obs_gof), row.names = NULL)
print(gof_tab, digits = 4)
write.csv(gof_tab, "output/gof_tests.csv", row.names = FALSE)

# ---- Kernels with m = 2 (Table 5) and sharp constants C_g -------------------
p <- 0.5; th <- 1; bt <- 1
K2 <- list(
  Gini      = list(g = function(a, b) abs(a - b), Cg = 1),
  Power     = list(g = function(a, b) (a^p + b^p) * (a + b)^(1 - p),
                   Cg = max(1, 2^(1 - p))),
  Symmetric = list(g = function(a, b) a^th * b^bt / (a + b)^(th + bt - 1),
                   Cg = th^th * bt^bt / (th + bt)^(th + bt)),
  Maximum   = list(g = function(a, b) pmax(a, b), Cg = 1),
  Minimum   = list(g = function(a, b) pmin(a, b), Cg = 0.5),
  SCV       = list(g = function(a, b) (a - b)^2 / (a + b), Cg = 1),
  Product   = list(g = function(a, b) a * b / (a + b), Cg = 1/4)
)

# ---- Order-statistic kernels (Table 6) --------------------------------------
# sum_{|S|=m} x_{k:S} = sum_i C(i-1,k-1) C(n-i,m-k) x_(i)
os_coef <- function(n, m, a) {
  i <- seq_len(n); cc <- numeric(n)
  for (k in seq_len(m))
    if (a[k] != 0) cc <- cc + a[k] * choose(i - 1, k - 1) * choose(n - i, m - k)
  cc / choose(n, m)
}
hat_I     <- function(x, m, a, Cg) sum(os_coef(length(x), m, a) * sort(x)) /
                                   (Cg * m * mean(x))
hat_I_mat <- function(Xs, m, a, Cg, xbar)
  as.vector(crossprod(os_coef(nrow(Xs), m, a), Xs)) / (Cg * m * xbar)

# Theoretical indices under Gamma(alpha, 1)
aux_A <- function(alpha, r)
  integrate(function(t) 1 - pgamma(t, shape = alpha)^r, 0, Inf,
            rel.tol = 1e-11, subdivisions = 500L)$value
EX <- function(alpha, m, k) {  # E[X_{k:m}]
  r <- k:m
  sum((-1)^(r - k) * choose(r - 1, k - 1) * choose(m, r) *
        vapply(r, function(rr) aux_A(alpha, rr), numeric(1)))
}
I_lin <- function(alpha, a) {
  m <- length(a)
  sum(a * vapply(seq_len(m), function(k) EX(alpha, m, k), numeric(1))) / (m * alpha)
}

cfg <- function(family, a, Cg = 1, conf = "--") {
  list(family = family, m = length(a), a = a, Cg = Cg, conf = conf,
       Ifun = function(alpha) I_lin(alpha, a) / Cg)
}
ext <- function(m, j, k) {
  a <- numeric(m); a[j] <- -1; a[k] <- 1
  cfg("Extended mth Gini", a, 1/(m - k + 1), sprintf("(%d,%d)", j, k))
}
TAB6 <- list(
  cfg("mth Gini", c(-1, 0, 1)), cfg("mth Gini", c(-1, 0, 0, 1)),
  ext(3, 1, 2), ext(3, 2, 3), ext(4, 1, 2), ext(4, 2, 3), ext(4, 3, 4),
  cfg("Linear Order-Statistic", c(-1, 0, 0, 1),    conf = "(-1,0,0,1)"),
  cfg("Linear Order-Statistic", c(-1, 0, 0, 0, 1), conf = "(-1,0,0,0,1)")
)

# ---- Bootstrap resamples (shared by Tables 5 and 6) -------------------------
qs <- c((1 - CONF)/2, 1 - (1 - CONF)/2)
set.seed(SEED)
XB    <- matrix(sample(gdp, n * B_BOOT, replace = TRUE), n, B_BOOT)
XBs   <- apply(XB, 2, sort)
xbarB <- colMeans(XB)
alpha_b <- vapply(seq_len(B_BOOT), function(b)
  gamma_mle(rgamma(n, shape = al, rate = la))[["alpha"]], numeric(1))

# spline of I(alpha) over the bootstrap range (avoids B integrations)
make_spline <- function(Ifun, n_grid = 400L) {
  gr <- seq(max(0.2, min(alpha_b) * 0.98), max(alpha_b) * 1.02, length.out = n_grid)
  splinefun(gr, vapply(gr, Ifun, numeric(1)), method = "natural")
}

# ---- Table 5 ----------------------------------------------------------------
prs <- combn(n, 2L)
Ai  <- XB[prs[1L, ], , drop = FALSE]
Aj  <- XB[prs[2L, ], , drop = FALSE]

tab5 <- do.call(rbind, lapply(names(K2), function(nm) {
  k    <- K2[[nm]]
  est  <- mean(k$g(gdp[prs[1L, ]], gdp[prs[2L, ]])) / (k$Cg * 2 * mean(gdp))
  bval <- colMeans(k$g(Ai, Aj)) / (k$Cg * 2 * xbarB)
  ci   <- quantile(bval, qs, names = FALSE)
  data.frame(index = nm, Cg = k$Cg, est = est, se = sd(bval),
             lo = ci[1], hi = ci[2])
}))
print(tab5, digits = 4)
write.csv(tab5, "output/table5_gdp.csv", row.names = FALSE)

# ---- Table 6 ----------------------------------------------------------------
tab6 <- do.call(rbind, lapply(TAB6, function(cf) {
  b_np <- hat_I_mat(XBs, cf$m, cf$a, cf$Cg, xbarB)
  b_pr <- make_spline(cf$Ifun)(alpha_b)
  ci_np <- quantile(b_np, qs, names = FALSE)
  ci_pr <- quantile(b_pr, qs, names = FALSE)
  data.frame(index = cf$family, m = cf$m, conf = cf$conf,
             np_est = hat_I(gdp, cf$m, cf$a, cf$Cg), np_se = sd(b_np),
             np_lo = ci_np[1], np_hi = ci_np[2],
             par_est = cf$Ifun(al), par_se = sd(b_pr),
             par_lo = ci_pr[1], par_hi = ci_pr[2])
}))
print(tab6, digits = 4)
write.csv(tab6, "output/table6_gdp.csv", row.names = FALSE)
