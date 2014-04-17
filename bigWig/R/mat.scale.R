#
# Functions to scale data matrices
#

# "RPKM stands for Reads Per Kilobase of transcript per Million mapped 
#  reads. FPKM stands for Fragments Per Kilobase of transcript per Million
#  mapped reads. "
#
# from
#  http://cufflinks.cbcb.umd.edu/faq.html
rpkm.scale <- function(mat, step, libSize) {
  # mat is in reads / step
  # m1 = mat * (step / 1000) => reads / kb
  # million mapped := libSize / 10^6
  # m1 / (libSize / 10^6) => reads / kb / million mapped
  scaleFactor = (step / 1000) / (libSize / 10^6)
  
  mat * scaleFactor
}

# rows are scaled to sum to one
densityToOne.scale <- function(mat, na.on.zero = TRUE) {
  t(apply(mat, 1, function(row) {
    if (all(is.na(row)))
      return(row)
    
    total = sum(row, na.rm = TRUE)
    if (total == 0) {
      if (na.on.zero)
        return(rep(NA, length(row)))
      else
        return(rep(0, length(row)))
    }
    row / total
  }))
}

maxToOne.scale <- function(mat) {
  t(apply(mat, 1, function(row) {
    if (all(is.na(row)))
      return(row)
    
    high = max(row, na.rm = TRUE)
    if (high == 0)
      return(rep(0, length(row)))
    
    row / high
  }))
}

zeroToOne.scale <- function(mat) {
  t(apply(mat, 1, function(row) {
    if (all(is.na(row)))
      return(row)
    
    high = max(row, na.rm = TRUE)
    low = min(row, na.rm = TRUE)
    
    if (high == low) {
      if (high == 0)
        return(rep(0, length(row)))
      else
        return(rep(1, length(row)))
    }
    
    (row - low) / (high - low)
  }))
}
