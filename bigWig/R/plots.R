plot.metaprofile <- function(x, minus.profile = NULL, X0 = plus.profile$X0, draw.error = TRUE, col = c("red", "blue", "lightgrey", "lightgrey"), ylim = NULL, xlim = NULL, xlab = "Distance (bp)", ylab = plus.profile$name, ...) {
  plus.profile = x
  opt.args = list(...)
  
  # check if profiles are all >= 0 (computed with abs.value = TRUE)
  if (!(all(plus.profile$top >= 0) && all(plus.profile$middle >= 0) && all(plus.profile$bottom >= 0)))
    warning("expected plus profile curves to all be >= 0")
  if (!is.null(minus.profile) && !(all(minus.profile$top >= 0) && all(minus.profile$middle >= 0) && all(minus.profile$bottom >= 0)))
    warning("expected minus profile curves to all be >= 0")
  
  # rebuild X coordinates
  N = length(plus.profile$middle)
  
  step = plus.profile$step
  if (is.null(step))
    stop("plus.profile is missing the 'step' element")

  stopifnot(X0 <= length(plus.profile$middle)*step && X0 >= 0)

  if (!is.null(minus.profile)) {
    stopifnot(length(plus.profile$middle) == length(minus.profile$middle))
    stopifnot(plus.profile$step == minus.profile$step)
  }
  x = 1:N * step - X0
  
  # establish xlim if not supplied
  if (is.null(xlim))
    xlim = c(min(x), max(x))
  
  # compute ylim if not supplied
  if (is.null(ylim)) {
    if (draw.error) {
      if (!is.null(minus.profile)) {
        ylim = c(
                 -max(minus.profile$top, minus.profile$middle, minus.profile$bottom),
                 max(plus.profile$top, plus.profile$middle, plus.profile$bottom))
      } else {
        ylim = c(
                 min(plus.profile$top, plus.profile$middle, plus.profile$bottom),
                 max(plus.profile$top, plus.profile$middle, plus.profile$bottom))
      }
    } else {
      if (!is.null(minus.profile)) {
        ylim = c(-max(minus.profile$middle), max(plus.profile$middle))
      } else {
        ylim = c(min(plus.profile$middle), max(plus.profile$middle))
      }
    }
  }
  
  # establish plot area
  plot(x, plus.profile$middle, col = col[1], xlim = xlim, ylim = ylim, type="l", xlab = xlab, ylab = ylab, ...)
  
  # draw error (shared area)
  if (draw.error) {
    polygon(c(x, rev(x)), c(plus.profile$top, rev(plus.profile$bottom)), col=col[3], border=NA, ...)
    
    if (!is.null(minus.profile))
      polygon(c(x, rev(x)), -c(minus.profile$bottom, rev(minus.profile$top)), col=col[4], border=NA, ...)
  }
  
  # draw main plot line
  lines(x, plus.profile$middle, col = col[1], ...)
  if (!is.null(minus.profile)) {
    lines(x, -minus.profile$middle, col = col[2], ...)
    lines(x, rep(0, N), ...)
  }
}
