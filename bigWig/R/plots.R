#
# New bigWig plot code
#

#
# track plotting functions
#

plot.track.bigWig <- function(counts.plus, counts.minus = NULL, newpage = TRUE, ylim = NULL) {
  vec.double <- function(x) {
    as.vector(sapply(x, function(val) c(val, val)))
  }
  
  # get and validate attributes
  attr.plus = attributes(counts.plus)
  stopifnot(length(attr.plus) > 0)
  stopifnot(all(c("chrom", "start", "end", "step") %in% names(attr.plus)))

  if (!is.null(counts.minus)) {
    attr.minus = attributes(counts.minus)

    stopifnot(attr.plus$chrom == attr.minus$chrom)
    stopifnot(attr.plus$start == attr.minus$start)
    stopifnot(attr.plus$end == attr.minus$end)
    stopifnot(attr.plus$step == attr.minus$step)
  }
  
  # define yrange
  if (is.null(ylim)) {
    ylim = NULL
    if (!is.null(counts.minus)) {
      ymax = max(counts.plus, abs(counts.minus))
      ylim = c(-ymax, ymax)
    } else {
      ymax = max(counts.plus)
      ylim = c(0, ymax)
    }
  }

  x.points = attr.plus$start + (1:length(counts.plus)) * attr.plus$step
  y.points.plus = counts.plus
  y.points.minus = counts.minus
  if (attr.plus$step > 1) {
    y.points.plus = c(0, vec.double(counts.plus), 0)
    y.points.minus = c(0, vec.double(-abs(counts.minus)), 0)
    x.points = attr.plus$start + (1:(length(counts.plus) + 1)) * attr.plus$step
    x.points = vec.double(x.points)
  }

  # draw
  #
  if (newpage)
      grid.newpage()

  pushViewport(viewport(xscale=c(attr.plus$start, attr.plus$end), yscale=ylim))

  grid.polygon(x.points, y.points.plus, default.units = "native", gp = gpar(fill = 'red', col=NA))

  if (!is.null(counts.minus)) {
    grid.polygon(x.points, y.points.minus, default.units = "native", gp = gpar(fill = 'blue', col=NA))
    
  }
  popViewport()
}

plot.bed.track.bigWig <- function(bed, chrom, start, end, newpage = TRUE, margin = unit(4, "pt"), height = NULL, draw.direction = TRUE, score.to.color = NULL, mark.color = "white") {
  valid.colors <- function(str) {
    valid = all(grepl("[0-5][0-5]?[0-5]?,[0-5][0-5]?[0-5]?,[0-5][0-5]?[0-5]?", str))
    if (valid) {
      parts = strsplit(str, ",")
      all(sapply(parts, function(color) {
        all(as.integer(color) <= 255)
      }))
    } else FALSE
  }
  
  get.colors <- function(str) {
    parts = strsplit(str, ",")
    sapply(parts, function(color) {
      values = as.integer(color)
      rgb(values[1], values[2], values[3], max = 255)
    })
  }

  stopifnot(!is.null(bed) && is.data.frame(bed) && dim(bed)[1] > 0)
  
  # select content from BED that matches requested region
  overlap.idxs = which(bed[,1] == chrom & bed[,2] < end & bed[,3] > start)
  
  if (length(overlap.idxs) == 0) {
    warning("No BED regions overlap with the requested range.")
    return()
  }
  bed = bed[overlap.idxs,]
  
  has.color = dim(bed)[2] >= 9 && valid.colors(as.character(bed[,9]))
  colors = NULL
  if (has.color)
    colors = get.colors(as.character(bed[,9]))
  
  # if a score.to.color function is given (and we have scores), then
  # override colors
  if (dim(bed)[2] >= 5 && is.function(score.to.color)) {
    colors = score.to.color(bed[,5])
    stopifnot(length(colors) == dim(bed)[1])
    has.color = TRUE
  }
  
  # stack remaining BED regions
  stack.idxs = .Call(stack_bed, bed[,2], bed[,3])
  stack.max = max(stack.idxs)
      
  # draw
  #
  if (newpage)
      grid.newpage()

  if (!is.null(height)) {
    npcTotalHeight = convertHeight(height * stack.max + margin * (stack.max - 1), "npc", valueOnly = TRUE)
    stopifnot(npcTotalHeight <= 1) # can't fit
  }

  # set viewport dimensions
  pushViewport(viewport(xscale=c(start, end)))
  
  if (is.null(height)) {
    # convert everything to 'npc' to do values
    margin.value = convertHeight(margin, "npc", valueOnly = TRUE)
    
    h = (1 - margin.value * (stack.max - 1)) / stack.max
    
    height = unit(h, "npc")
  }
  
  foreach.bed(bed, function(i, chrom, start.i, end.i, strand) {
    stack.level = stack.idxs[i]
    
    # compute coordinates
    x = unit(max(start.i, start), 'native')
    y = convertHeight(unit(1, "npc") - height - (stack.level - 1) * (height + margin), "npc")
    
    #
    w = unit(min(end, end.i) - max(start, start.i), 'native')
    
    # draw rect
    # fill color
    gp = NULL
    if (has.color) {
      gp = gpar(col = colors[i], fill = colors[i])
    } else {
      gp = gpar(fill = "black")
    }
    #
    grid.rect(x, y, width = w, height = height, just = c(0, 0), gp = gp)
    
    # draw visual edge indicators
    #
    soff = 0
    eoff = 0
    if (start.i < start && convertWidth(w, "pt", valueOnly=TRUE) > 36) {
      ybottom = y + unit(4, "pt")
      hi = height - unit(8, "pt")
      
      x0 = convertWidth(x + unit(4, "pt"), "native")
      x1 = convertWidth(x + unit(12, "pt"), "native")
      
      grid.rect(x0, ybottom, width = unit(4, "pt"), height = hi, just = c(0, 0), gp = gpar(col=mark.color, fill=mark.color))
      grid.rect(x1, ybottom, width = unit(4, "pt"), height = hi, just = c(0, 0), gp = gpar(col=mark.color, fill=mark.color))
      
      soff = 16
    }
    if (end.i > end && (convertWidth(w, "pt", valueOnly=TRUE) > 36)) {
      ybottom = y + unit(4, "pt")
      hi = height - unit(8, "pt")
      
      x0 = convertWidth(x + w - unit(8, "pt"), "native")
      x1 = convertWidth(x + w - unit(18, "pt"), "native")
      
      grid.rect(x0, ybottom, width = unit(4, "pt"), height = hi, just = c(0, 0), gp = gpar(col=mark.color, fill=mark.color))
      grid.rect(x1, ybottom, width = unit(4, "pt"), height = hi, just = c(0, 0), gp = gpar(col=mark.color, fill=mark.color))
      
      eoff = 16
    }
    
    # draw direction marks
    #
    if (draw.direction && !is.na(strand)) {
      ybottom = y + unit(4, "pt")
      ytop = y + height - unit(4, "pt")
      
      v1 = convertHeight(ytop, "pt", valueOnly=TRUE)
      v2 = convertHeight(ybottom, "pt", valueOnly=TRUE)
      ymid = unit(mean(c(v1, v2)), "pt")
      
      w0 = convertWidth(w, "pt", valueOnly=TRUE)
      
      if (strand == '+') {
        xs = soff + 8
        xe = w0 - 10 - eoff

        if (xe > xs + 20) {
          for (xi in seq(xs, xe, 20)) {
            grid.lines(x + unit(c(xi, xi + 8, xi), "pt"), unit.c(ybottom, ymid, ytop), gp = gpar(col=mark.color, lwd=2))
          }
        }
      } else if (strand == '-') {
        xs = soff + 16
        xe = w0 - 4 - eoff

        if (xe > xs + 20) {
          for (xi in seq(xs, xe, 20)) {
            grid.lines(x + unit(c(xi, xi - 8, xi), "pt"), unit.c(ybottom, ymid, ytop), gp = gpar(col=mark.color, lwd=2))
          }
        }
      }
    }
  })
    
  popViewport()
}


plot.profile.bigWig <- function(plus.profile, minus.profile = NULL, X0 = 0, draw.error = TRUE, new.page = TRUE, xlim = NULL, ylim = NULL, draw.axis = c(TRUE, TRUE), scalebar = NA, fill = FALSE, smoothed = FALSE, draw = TRUE, col = c("red", "blue", "lightgrey", "lightgrey"), gp = gpar()) {
  # check if profiles are all >= 0 (computed with abs.value = TRUE)
  if (!(all(plus.profile$top >= 0) && all(plus.profile$middle >= 0) && all(plus.profile$bottom >= 0)))
    warning("expected plus profile curves to all be >= 0")
  if (!is.null(minus.profile) && !(all(minus.profile$top >= 0) && all(minus.profile$middle >= 0) && all(minus.profile$bottom >= 0)))
    warning("expected minus profile curves to all be >= 0")

  # rebuild X coordinates
  N = length(plus.profile$middle)
  stopifnot(X0 <= length(plus.profile$middle) && X0 >= 0)
    
  step = plus.profile$step
  if (is.null(step))
    stop("plus.profile is missing the 'step' element")
  
  if (!is.null(minus.profile)) {
    stopifnot(length(plus.profile$middle) == length(minus.profile$middle))
    stopifnot(plus.profile$step == minus.profile$step)
  }

  x = (1:N - X0) * step

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

  # establish xlim if not supplied
  if (is.null(xlim))
    xlim = c(min(x), max(x))
  
  #
  # TODO: create a grob instead of drawing directly!!
  
  # create plot area
  if (new.page)
    grid.newpage()

  pushViewport(plotViewport(margins = par()$mar)) # TODO: This is probably not the 'right' thing to do
  pushViewport(dataViewport(xscale = xlim, yscale = ylim, name = "profilePlot", clip = TRUE))

  # error polygon(s) (see draw.error and smoothed arguments and use colors from col)
  if (draw.error) {
    # TODO: smoothed ...
    
    gperr = gp
    gperr$fill = col[3]
    gperr$col = col[3]
    
    grid.polygon(x = c(x, rev(x)), y = c(plus.profile$top, rev(plus.profile$bottom)), gp = gperr, default.units = "native")
    
    if (!is.null(minus.profile)) {
      gperr = gp
      gperr$fill = col[4]
      gperr$col = col[4]
    
      grid.polygon(x = c(x, rev(x)), y = -c(minus.profile$bottom, rev(minus.profile$top)), gp = gperr, default.units = "native")
    }
  }
  
  # main (see fill, use colors from col)
  gpmain = gp
  gpmain$col = col[1]
  grid.polyline(x = unit(x, "native"), y = unit(plus.profile$middle, "native"), gp = gpmain)
  
  if (!is.null(minus.profile)) {
    gpmain = gp
    gpmain$col = col[2]
    grid.polyline(x = unit(x, "native"), y = unit(-minus.profile$middle, "native"), gp = gpmain)
  }
  
  # stop clipping (not sure if this is the best way to do it ...)
  pushViewport(dataViewport(xscale = xlim, yscale = ylim, name = "profilePlot", clip = "off"))
  
  # axis (see draw.axis and scalebar arguments)
  if (draw.axis[1])
    grid.xaxis(gp = gp)
  if (draw.axis[2])
    grid.yaxis(gp = gp)
  # TODO: scalebar

  # outer rect
  grid.rect(gp = gp)
  
  # TODO: ylab or title ?? (attr of plus.profile; maybe use plus.profile "original argument expression")
  # TODO: "xlab ??
  # TODO: gpar ??
  # TODO: draw ??
  
  upViewport(3)
}

#
# compressed region(s) plots
#
