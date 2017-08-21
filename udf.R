# udf - my misc user defined functions in R
# MCC formulas
udf_mcc <- function(TP, FP, FN, TN) {
  num <- (TP*TN) - (FP*FN)
  den <- (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)
  if (den == 0){
    return(0)
  }else {
    return(num / sqrt(den))
  }
}
udf_mcc <- function(y_true, y_pred) {
  y_true <- as.numeric(y_true)
  y_pred <- as.numeric(y_pred)
  
  TP <- sum(as.numeric(y_true == y_pred)[y_true == 1])
  TN <- sum(as.numeric(y_true == y_pred)[y_true == 0])
  FP <- sum(as.numeric(y_true != y_pred)[y_true == 1])
  FN <- sum(as.numeric(y_true != y_pred)[y_true == 0])
  
  num <- (TP*TN) - (FP*FN)
  den <- (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)
  if (den == 0){
    return(0)
  }else {
    return(num / sqrt(den))
  }
}
udf_eval_mcc <- function(y_true, y_prob, show=F) {
  y_true <- as.numeric(y_true)
  y_prob <- as.numeric(y_prob)
  idx <- order(y_prob)
  y_prob_sort <- y_prob[idx]
  y_true_sort <- y_true[idx]
  n <- length(y_true)
  nump <- sum(y_true)
  numn <- n - nump
  
  tn_v <- cumsum(y_true_sort == 0)
  fp_v <- cumsum(y_true_sort == 1)
  fn_v <- numn - tn_v
  tp_v <- nump - fp_v
  sup_v <- tp_v * tn_v - fp_v * fn_v
  inf_v <- sqrt((tp_v + fp_v) * (tp_v + fn_v) * (tn_v + fp_v) * (tn_v + fn_v))
  mcc_v <- sup_v/inf_v
  mcc_v[!is.finite(mcc_v)] <- 0
  best_id <- which.max(mcc_v)
  best_mcc <- mcc_v[best_id]
  best_proba <- y_prob_sort[best_id]
  
  if (show) {
    y_pred <- as.numeric(y_prob > best_proba)
    #score <- udf_mcc(y_true, y_pred)
    #cat("\n", score, best_mcc)
    return(list(p = best_proba, mcc = best_mcc))
  } else {
    return(best_mcc)
  }
}
mcc_eval <- function(y_prob, dtrain) {
  y_true <- getinfo(dtrain, "label")
  best_mcc <- udf_eval_mcc(y_true, y_prob)
  return(list(metric="mcc", value=best_mcc))
}
udf_colMin <- function(dta, mx_cols) {
  dta[, mx_i := sort(rep(1:ceiling(nrow(dta)/10000), length.out = nrow(dta)))]
  dta[, mx_min := 0]
  pb <- txtProgressBar(min = 1, max = max(dta$mx_i), style = 3)
  for (mxi in unique(dta$mx_i)) {
    dta[mx_i == mxi, 
        mx_min := apply(dta[mx_i == mxi, mx_cols, with=F], 1, 
                        FUN = function(x) min(x, na.rm = T))]
    setTxtProgressBar(pb, mxi)
  }
  dta[is.infinite(mx_min), mx_min := NA]
  dta[, mx_i := NULL]
}
udf_colMax <- function(dta, mx_cols) {
  dta[, mx_i := sort(rep(1:ceiling(nrow(dta)/10000), length.out = nrow(dta)))]
  dta[, mx_max := 0]
  pb <- txtProgressBar(min = 1, max = max(dta$mx_i), style = 3)
  for (mxi in unique(dta$mx_i)) {
    dta[mx_i == mxi, 
        mx_max := apply(dta[mx_i == mxi, mx_cols, with=F], 1, 
                        FUN = function(x) max(x, na.rm = T))]
    setTxtProgressBar(pb, mxi)
  }
  dta[is.infinite(mx_max), mx_max := NA]
  dta[, mx_i := NULL]
}
udf_colDigest <- function(dta, mx_cols) {
  dta[, mx_i := sort(rep(1:ceiling(nrow(dta)/10000), length.out = nrow(dta)))]
  dta[, mx_hash := ""]
  pb <- txtProgressBar(min = 1, max = max(dta$mx_i), style = 3)
  for (mxi in unique(dta$mx_i)) {
    dta[mx_i == mxi, mx_hash := apply(dta[mx_i == mxi, mx_cols, with=F], 1, 
                                      FUN = function(x) digest(x,algo="crc32"))]
    setTxtProgressBar(pb, mxi)
  }
  dta[, mx_i := NULL]
}
#
udf_aggByValues <- function(values, Nbins = 100) {
  f_aggBy <- ecdf(quantile(values, prob = (1:1000)/1000))(values)
  f_aggBy <- round(f_aggBy/(1/Nbins))*(1/Nbins)
  f_aggBy <- quantile(values, f_aggBy)
  return(f_aggBy)
}
udf_ecdf <- function(sample, level) {
  return(ecdf(quantile(sample, prob = 1:1000)/1000)(level))
}
udf_quantile <- function(x, q = c(0, 0.01, 0.02, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.98, 0.99, 1)) {
  return(quantile(x, q, na.rm = T))
}
udf_quantile_mtrx <- function(x, q = c(0, 0.01, 0.02, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.98, 0.99, 1)) {
  mx_cols <- colnames(x)
  x <- data.frame(x)
  for (i in 1:length(mxcols))
  {
    mx_tmp <- quantile(x[, i], q, na.rm = T)
    if (i == 1) {
      mx_res <- mx_tmp
    } else {
      mx_res <- rbind(mx_res, mx_tmp)
    }
  }
  rownames(mx_res) <- mx_cols
  return(mx_res)
}
udf_hist_rngcut <- function(x, breaks = 50, main.title = "", cut.min = 0.005, cut.max = 0.995) {
  mx_cut <- quantile(x, c(cut.min, cut.max), na.rm = T)
  x <- x[mx_cut[1] <= x & x <= mx_cut[2]]
  hist(x, breaks, main = paste("hist:", main.title), 
       xlim = mx_cut, col = "grey85", 
       xlab = paste0("range cut: [", cut.min, "-", cut.max, "]"))
  abline(v = mean(x, na.rm = T), col = "red")
  abline(v = median(x, na.rm = T), col = "black")
}
udf_R.squared <- function(y, predy, weights = NULL, intercept = TRUE) {
  # https://stat.ethz.ch/pipermail/r-help/2006-August/111769.html
  if (is.null(weights)) 
  {
    noNA <- !is.na(y) & !is.na(predy) & !is.infinite(y) & !is.infinite(predy)
    m <- 0
    if (intercept) { m <- mean(y[noNA])}
    mss <- sum((y[noNA] - m)^2)
    rss <- sum((y[noNA] - predy[noNA])^2)
  } else {
    noNA <- !is.na(y) & !is.na(predy) & !is.na(weights)
    m <- 0
    if (intercept) { m <- sum(weights[noNA]*y[noNA]/sum(weights[noNA])) }
    mss <- sum(weights[noNA]*(y[noNA] - m)^2)
    rss <- sum(weights[noNA]*(y[noNA] - predy[noNA])^2)
  }
  return(1-rss/mss)
}
#
# F1 score defined for lists (concatenated strings of product_ids)
udf_f1score <- function(list_a, list_b) {
  list_a <- str_split(list_a, ' ')[[1]]
  list_b <- str_split(list_b, ' ')[[1]]
  pr <- length(intersect(list_a,list_b))/length(list_b)
  re <- length(intersect(list_a,list_b))/length(list_a)
  f1 <- 0
  if (pr + re) { f1 <- 2 * pr * re /(pr + re) }
  return(f1)
}
udf_meanF1score <- function(y, predy) {
  colnames(y) <- c("id", "y", "isIn")
  y_zero <- y[, list(y = NA, isIn = sum(isIn)), by = .(id)]
  y <- rbind(y[isIn == 1, ], y_zero[isIn == 0, ])
  y[, isIn := NULL]
  
  colnames(predy) <- c("id", "y", "isIn")
  predy_zero <- predy[, list(y = NA, isIn = sum(isIn)), by = .(id)]
  predy <- rbind(predy[isIn == 1, ], predy_zero[isIn == 0, ])
  predy[, isIn := NULL]
  
  y <- data.table(y)
  predy <- data.table(predy)
  colnames(y) <- c("id", "y")
  colnames(predy) <- c("id", "y")
  
  y <- y[, list(isIn = .N), by = .(id, y)]
  predy <- predy[, list(isIn = .N), by = .(id, y)]
  
  m <- merge(y, predy, by = c("id", "y"), all=T, suffixes=c(".y", ".pred"))
  m[is.na(isIn.y), isIn.y := 0]
  m[is.na(isIn.pred), isIn.pred := 0]
  m[, isIn.both := pmin(isIn.y, isIn.pred)]
  
  agg <- m[, list(p = sum(isIn.both)/sum(isIn.pred), 
                  r = sum(isIn.both)/sum(isIn.y)), keyby = id]
  agg[, f1score := 2*p*r/(p+r)]
  agg[is.na(f1score), f1score := 0]
  
  mean(agg$f1score)
}
udf_meanF1scoreAgg <- function(y, predy) {
  colnames(y) <- c("id", "y", "isIn")
  y_zero <- y[, list(y = NA, isIn = sum(isIn)), by = .(id)]
  y <- rbind(y[isIn == 1, ], y_zero[isIn == 0, ])
  y[, isIn := NULL]
  
  colnames(predy) <- c("id", "y", "isIn")
  predy_zero <- predy[, list(y = NA, isIn = sum(isIn)), by = .(id)]
  predy <- rbind(predy[isIn == 1, ], predy_zero[isIn == 0, ])
  predy[, isIn := NULL]
  
  y <- data.table(y)
  predy <- data.table(predy)
  colnames(y) <- c("id", "y")
  colnames(predy) <- c("id", "y")
  
  y <- y[, list(isIn = .N), by = .(id, y)]
  predy <- predy[, list(isIn = .N), by = .(id, y)]
  
  m <- merge(y, predy, by = c("id", "y"), all=T, suffixes=c(".y", ".pred"))
  m[is.na(isIn.y), isIn.y := 0]
  m[is.na(isIn.pred), isIn.pred := 0]
  m[, isIn.both := pmin(isIn.y, isIn.pred)]
  
  agg <- m[, list(p = sum(isIn.both)/sum(isIn.pred), 
                  r = sum(isIn.both)/sum(isIn.y)), keyby = id]
  agg[, f1score := 2*p*r/(p+r)]
  agg[is.na(f1score), f1score := 0]
  
  return(agg[, .(id, f1score)])
}
#
udf_grepl <- function(patterns, values) {
  mx_filt <- rep(F, length(values))
  for (patt in patterns){ mx_filt <- mx_filt | grepl(patt, values) }
  return(mx_filt)
}
udf_RMSLE <- function(predy, y){
  n <- length(y)
  RMSLE <- (1/n*sum((log(predy + 1) - log(y + 1))^2))^0.5
  return(RMSLE)
}
udf_memorycheck <- function() {
  gc(verbose = F)
  memorycheck. <- sapply(ls(envir = globalenv()), function(x){object.size(get(x))})
  mx_names <- names(memorycheck.)
  memorycheck. <- as.array(as.numeric(memorycheck.))
  mx_order <- order(memorycheck., decreasing = T)
  memorycheck. <- memorycheck.[mx_order]
  memorycheck. <- round(memorycheck./sum(memorycheck., na.rm = T), 2)
  mx_names <- mx_names[mx_order]
  mx_df <- data.frame(obj = mx_names, memory = memorycheck.)
  mx_df <- mx_df[mx_df$memory > 0, ]
  mx_df
}
#udf_memorycheck <- function() {
#  sort(round(sapply(ls(envir = globalenv()),function(x){object.size(get(x))})/sum(sapply(ls(envir = globalenv()),function(x){object.size(get(x))})),2), decreasing=T)
#}
