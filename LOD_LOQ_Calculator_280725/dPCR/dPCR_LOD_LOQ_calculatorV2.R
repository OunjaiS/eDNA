## -----------------------------
## Load required packages
## -----------------------------
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(drc)) install.packages("drc")

## -----------------------------
## Read input dPCR data from CSV
## -----------------------------
# ปรับ path ให้ตรงกับไฟล์ของคุณ
d <- read.csv('/Users/sarawut/Desktop/LOD_LOQ_Calculator/dPCR/example_data/GC_dPCR_synth_gene_dil.csv', header=TRUE)

# ตัวแปรจากไฟล์
y <- d$npositive
m <- d$ntotal
x.sample <- d$SampleConc   # ความเข้มข้นที่ใส่ใน sample (copies/µL)

# แปลงให้เป็น concentration ใน reaction
x <- x.sample * 2 / 25  # 2 uL sample ใน total 25 uL reaction

# Volume per droplet (µL)
v <- 0.91 / 1000

## -----------------------------
## Define Model for LOD
## -----------------------------
LODmodel.glm = function(xminOfStds, y, m, x, v) {
  std <- x >= xminOfStds
  y.std <- y[std]; m.std <- m[std]; x.std <- x[std]
  y.low <- y[!std]; m.low <- m[!std]
  
  # fit standard curve
  v.offset <- rep(log(v), length(y.std))
  fit.std <- glm(cbind(y.std, m.std - y.std) ~ log(x.std),
                 family = binomial(link = 'cloglog'),
                 offset = v.offset)
  
  # fit background/noise
  v.offset <- rep(log(v), length(y.low))
  fit.low <- glm(cbind(y.low, m.low - y.low) ~ 1,
                 family = binomial(link = 'cloglog'),
                 offset = v.offset)
  
  # fit combined model
  std.ind <- as.integer(std)
  logx <- rep(0, length(x)); logx[std] <- log(x[std])
  X <- cbind(std.ind, logx, 1 - std.ind)
  
  fit.all <- glm(cbind(y, m - y) ~ X - 1,
                 family = binomial(link = 'cloglog'),
                 offset = rep(log(v), length(y)))
  mle <- fit.all$coefficients
  vcv <- vcov(fit.all)
  
  # Calculate LOD
  logLOD <- (mle[3] - mle[1]) / mle[2]
  logLOD.var <- (vcv[3,3] + vcv[1,1] + logLOD^2 * vcv[2,2] -
                   2*vcv[3,1] - 2*logLOD*vcv[3,2] + 2*logLOD*vcv[1,2]) / mle[2]^2
  logLOD.se <- sqrt(logLOD.var)
  
  return(list(beta.mle = mle[1:2],
              beta.se = sqrt(diag(vcv[1:2, 1:2])),
              alpha.mle = mle[3],
              alpha.se = sqrt(vcv[3,3]),
              logLOD = unname(logLOD),
              logLOD.se = unname(logLOD.se),
              aic = fit.all$aic,
              std.ind = std.ind,
              fit.all = fit.all))
}

## -----------------------------
## Print LOD + LOQ summary
## -----------------------------
LODmodel.glm.summary = function(fit, y, m, x, v) {
  zcrit <- qnorm(0.975)
  beta0 <- fit$beta.mle[1]
  beta1 <- fit$beta.mle[2]
  
  beta.lowerCL <- beta0 - zcrit * fit$beta.se[1]
  beta.upperCL <- beta0 + zcrit * fit$beta.se[1]
  alpha.lowerCL <- fit$alpha.mle - zcrit * fit$alpha.se
  alpha.upperCL <- fit$alpha.mle + zcrit * fit$alpha.se
  logLOD.lowerCL <- fit$logLOD - zcrit * fit$logLOD.se
  logLOD.upperCL <- fit$logLOD + zcrit * fit$logLOD.se
  
  LOD <- exp(fit$logLOD)
  
  # คำนวณ LOQ จาก residuals ของ standard curve
  std.x <- x[fit$std.ind == 1]
  lambda.repl <- (-1/v) * log(1 - y/m)
  lambda.std <- lambda.repl[fit$std.ind == 1]
  
  # กรองค่าที่ไม่เหมาะสม
  valid <- !is.infinite(lambda.std) & !is.nan(lambda.std) & lambda.std > 0
  lambda.std.valid <- lambda.std[valid]
  std.x.valid <- std.x[valid]
  
  fitted.log <- beta0 + beta1 * log(std.x.valid)
  residuals <- log(lambda.std.valid) - fitted.log
  sd_resid <- sd(residuals)
  
  LOQ_log <- mean(fitted.log) + 3 * sd_resid
  LOQ <- exp(LOQ_log)
  
  out <- rbind(
    c(beta0, fit$beta.se[1], beta.lowerCL, beta.upperCL),
    c(beta1, fit$beta.se[2], beta.lowerCL, beta.upperCL),
    c(fit$alpha.mle, fit$alpha.se, alpha.lowerCL, alpha.upperCL),
    c(fit$logLOD, fit$logLOD.se, logLOD.lowerCL, logLOD.upperCL),
    c(LOD, NA, exp(logLOD.lowerCL), exp(logLOD.upperCL)),
    c(LOQ_log, NA, NA, NA),
    c(LOQ, NA, NA, NA)
  )
  rownames(out) <- c("beta0", "beta1", "alpha", "logLOD", "LOD", "logLOQ", "LOQ")
  colnames(out) <- c("MLE", "SE", "2.5%", "97.5%")
  print(out)
  
  return(list(LOD=LOD, LOQ=LOQ, residuals=residuals))
}

## -----------------------------
## Select threshold and run model
## -----------------------------
x.unique <- sort(unique(x))
x.test <- x.unique[2:(length(x.unique)-1)]
aics <- sapply(x.test, function(xmin) {
  fit <- LODmodel.glm(xmin, y, m, x, v)
  return(fit$aic)
})

best.xmin <- x.test[which.min(aics)]
fit <- LODmodel.glm(best.xmin, y, m, x, v)

# Show result
res <- LODmodel.glm.summary(fit, y, m, x, v)

# คำนวณ lambda สำหรับ replicate แต่ละจุด
lambda.repl <- (-1/v) * log(1 - y/m)

# สรุป lambda โดยรวม (sum) สำหรับแต่ละค่าความเข้มข้น
ysum <- tapply(y, x, sum)
msum <- tapply(m, x, sum)
lambda.sum <- (-1/v) * log(1 - ysum/msum)

# กำหนดช่วงแกน x, y (log scale)
x.unique <- sort(unique(x))
xylim <- c(min(log(x.unique)) - 1, max(log(x.unique)) + 1)

# สร้าง plot จุดข้อมูลแต่ละ replicate (จุดดำ) และรวม (จุดแดง)
plot(log(x), log(lambda.repl), pch=19, cex=1.2,
     xlab='Log concentration (copies/µL)',
     ylab='Log measured concentration (copies/µL)',
     xlim=xylim, ylim=xylim, las=1)
points(log(x.unique), log(lambda.sum), pch=4, col='red', cex=2, lwd=2)

# เส้น standard curve fit
curve.beta0 <- fit$beta.mle[1]
curve.beta1 <- fit$beta.mle[2]
curve_x <- seq(xylim[1], xylim[2], length.out=100)
curve_y <- curve.beta0 + curve.beta1 * curve_x
lines(curve_x, curve_y, col='blue', lwd=2)

# เส้น LOD
abline(v = fit$logLOD, col='green', lwd=2, lty=2)
abline(h = fit$alpha.mle, col='green', lwd=2, lty=2)

# เส้น LOQ (คำนวณ lambda.LOQ จาก CV=0.25)
CVtarget <- 0.25
lambda.LOQ <- 1/(CVtarget^2 * v)
log.lambda.LOQ <- log(lambda.LOQ)
logLOQ <- (log.lambda.LOQ - curve.beta0) / curve.beta1
abline(v = logLOQ, col='purple', lwd=2, lty=3)

legend("bottomright", legend=c("Replicates", "Summed", "Standard curve", "LOD", "LOQ"),
       pch=c(19,4, NA, NA, NA), lty=c(NA, NA, 1, 2, 3), col=c("black","red","blue","green","purple"))

