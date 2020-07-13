gausskern.pg <- function (siz, sigma, muadv = 0) {
  if (round(siz) < 1) siz = 1
  x = 1:round(siz)
  mu = c(mean(x), mean(x)) + muadv
  # options(digits = 5)
  fx = exp(-0.5 * ((x - mu[1])/sigma)^2)/sqrt((2 * pi) * sigma**2) 
  fy = exp(-0.5 * ((x - mu[2])/sigma)^2)/sqrt((2 * pi) * sigma**2)
  fx[!is.finite(fx)] = 0
  fy[!is.finite(fy)] = 0
  kern = (fx %*% t(fy))
  kern = kern/(sum(kern, na.rm = T))#;kern
  kern[is.nan(kern)] = 0
  kern
}