#' Classify biomes from monthly temperature and precipitation data
#'
#' @param Temp a SpatRaster object with monthly temperature data, with
#'   each month on a different layer.
#' @param Prec a SpatRaster object with montly precipitation data,
#'   with each month on a different layer.
#' @param type the classification type; `type = "broadclass"` returns
#'   a classification into 5 broad categories (polar, cold, temperate,
#'   tropical, and arid), while `type = "class"` will use the
#'   Köppen-Geiger classification (30 categories).
#' @return A SpatRaster with the resulting biome classification.
kg_reclass_raster <- function(Temp, Prec, type = c("broadclass", "class")) {
  type <- match.arg(type)

  if (!identical(dim(Temp), dim(Prec))) {
    stop("Data matrices do not have the same dimensions")
  }

  if (any(Prec[] < 0, na.rm = TRUE)) {
    stop("Precipitation data cannot include negative values")
  }

  Temp[is.na(Temp)] <- -99
  Prec[is.na(Prec)] <- -99

  T_ONDJFM <- terra::mean(Temp[[c(1, 2, 3, 10, 11, 12)]])
  T_AMJJAS <- terra::mean(Temp[[c(4, 5, 6, 7, 8, 9)]])
  tmp <- T_AMJJAS > T_ONDJFM
  SUM_SEL <- c(rep(!tmp, 3), rep(tmp, 6), rep(!tmp, 3))

  Pw <- sum(Prec * !SUM_SEL)
  Ps <- sum(Prec * SUM_SEL)

  Pdry <- min(Prec)

  tmp <- SUM_SEL
  tmp[!tmp] <- NA
  Psdry <- min(Prec * tmp, na.rm = TRUE)
  Pswet <- max(Prec * tmp, na.rm = TRUE)

  tmp <- !SUM_SEL
  tmp[!tmp] <- NA
  Pwdry <- min(Prec * tmp, na.rm = TRUE)
  Pwwet <- max(Prec * tmp, na.rm = TRUE)

  MAT <- terra::mean(Temp)
  MAP <- sum(Prec)
  Tmon10 <- sum(Temp > 10)
  Thot <- max(Temp)
  Tcold <- min(Temp)

  Pthresh <- 2 * MAT + 14 #where temp = -99, this is -184
  Pthresh[Pw > Ps * 2.333] <- 2 * MAT[Pw > Ps * 2.333]
  Pthresh[Ps > Pw * 2.333] <- 2 * MAT[Ps > Pw * 2.333] + 28

  B <- MAP < 10 * Pthresh
  BW <- B & MAP < 5 * Pthresh
  BWh <- BW & MAT >= 18
  BWk <- BW & MAT < 18
  BS <- B & MAP >= 5 * Pthresh
  BSh <- BS & MAT >= 18
  BSk <- BS & MAT < 18

  A <- Tcold >= 18 & !B
  Af <- A & Pdry >= 60
  Am <- A & !Af & Pdry >= 100 - MAP / 25
  Aw <- A & !Af & Pdry < 100 - MAP / 25

  C <- Thot > 10 & Tcold > 0 & Tcold < 18 & !B
  Cs <- C & Psdry < 40 & Psdry < Pwwet / 3
  Cw <- C & Pwdry < Pswet / 10
  overlap <- Cs & Cw
  Cs[overlap & Ps > Pw] <- 0
  Cw[overlap & Ps <= Pw] <- 0
  Csa <- Cs & Thot >= 22
  Csb <- Cs & !Csa & Tmon10 >= 4
  Csc <- Cs & !Csa & !Csb & Tmon10 >= 1 & Tmon10 < 4
  Cwa <- Cw & Thot >= 22
  Cwb <- Cw & !Cwa & Tmon10 >= 4
  Cwc <- Cw & !Cwa & !Cwb & Tmon10 >= 1 & Tmon10 < 4
  Cf <- C & !Cs & !Cw
  Cfa <- Cf & Thot >= 22
  Cfb <- Cf & !Cfa & Tmon10 >= 4
  Cfc <- Cf & !Cfa & !Cfb & Tmon10 >= 1 & Tmon10 < 4

  D <- Thot > 10 & Tcold <= 0 & !B
  Ds <- D & Psdry < 40 & Psdry < Pwwet / 3
  Dw <- D & Pwdry < Pswet / 10
  overlap <- Ds & Dw
  Ds[overlap & Ps > Pw] <- 0
  Dw[overlap & Ps <= Pw] <- 0
  Dsa <- Ds & Thot >= 22
  Dsb <- Ds & !Dsa & Tmon10 >= 4
  Dsd <- Ds & !Dsa & !Dsb & Tcold < (-38)
  Dsc <- Ds & !Dsa & !Dsb & !Dsd

  Dwa <- Dw & Thot >= 22
  Dwb <- Dw & !Dwa & Tmon10 >= 4
  Dwd <- Dw & !Dwa & !Dwb & Tcold < (-38)
  Dwc <- Dw & !Dwa & !Dwb & !Dwd
  Df <- D & !Ds & !Dw
  Dfa <- Df & Thot >= 22
  Dfb <- Df & !Dfa & Tmon10 >= 4
  Dfd <- Df & !Dfa & !Dfb & Tcold < (-38)
  Dfc <- Df & !Dfa & !Dfb & !Dfd

  E <- Thot <= 10 & Thot > (-90) & !B # Added & Thot > (-90)
  ET <- E & Thot > 0
  EF <- E & Thot <= 0 & Thot > (-90) # Added & Thot> (-90)

  if (type == "class") {
    Classes <- list(
      Af, Am, Aw, BWh, BWk, BSh, BSk, Csa, Csb, Csc, Cwa, Cwb,
      Cwc, Cfa, Cfb, Cfc, Dsa, Dsb, Dsc, Dsd, Dwa, Dwb, Dwc, Dwd, Dfa,
      Dfb, Dfc, Dfd, ET, EF
    )
  } else if (type == "broadclass") {
    Classes <- list(A, B, C, D, E)
  }

  Class_cont <- terra::rast(Temp[[1]], vals = NA_integer_)
  for (i in seq_along(Classes)) {
    Class_cont[Classes[[i]]] <- i
  }

  Class_cont
}
