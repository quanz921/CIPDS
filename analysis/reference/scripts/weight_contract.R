cipds_pooled_mec_weights <- function(dat, n_cycles = 9L) {
  stopifnot(n_cycles == 9L)
  needed <- c("CYCLE", "WTMEC2YR", "WTMEC4YR")
  if (!all(needed %in% names(dat))) stop("Verified original MEC weights missing")
  early <- as.character(dat$CYCLE) %in% c("1999-2000", "2001-2002")
  result <- as.numeric(dat$WTMEC2YR) / n_cycles
  result[early] <- 2 * as.numeric(dat$WTMEC4YR[early]) / n_cycles
  if (any(!is.finite(result) | result <= 0)) stop("Invalid corrected pooled weights")
  result
}
CIPDS_WEIGHT_DESCRIPTION <- "(2/9)*WTMEC4YR for 1999-2002; (1/9)*WTMEC2YR for 2003-2016"
