# jonashaslbeck@protonmail.com; May 25th, 2023

# -------------------------------------------------------------
# ---------- What is happening here? --------------------------
# -------------------------------------------------------------

# We simulate from a set of multilevel VAR models in which we vary
# 1) the entries of the lagged-effects matrix
# 2) the variance of random intercepts, in a Gaussian distribution
# 3) the length of the time series of each subject

# We evaluate the observed correlation between personwise means


# -------------------------------------------------------------
# ---------- Load Packages ------------------------------------
# -------------------------------------------------------------

# Generate data
library(mlVAR)
library(MASS)

# Plotting
library(RColorBrewer)

# Parallelization
library(foreach)
library(parallel)
library(doParallel)


# -------------------------------------------------------------
# ---------- Function to compute cor btw personwise means ---------
# -------------------------------------------------------------

getBetCor <- function(n_subj, phi_en, Nt, sd_int) {

  # Storage
  l_data <- list()

  # Specify VAR
  phi <- matrix(phi_en, 2, 2)
  Sigma <- diag(2) * sd_int^2 # since we supply cov matrix

  for(i in 1:n_subj) {
    # Draw intercepts
    ints <- mvrnorm(n=1, mu=c(0,0), Sigma=Sigma*.01)

    # Sample data
    l_data[[i]] <- simulateVAR(pars = phi,
                               means = ints,
                               Nt = Nt)
  }

  # Compute Cor
  m_means <- do.call(rbind, lapply(l_data, colMeans))
  s_mean <- cor(m_means)[1,2]

  return(s_mean)

} # eoF

# -------------------------------------------------------------
# ---------- Simulate  ----------------------------------------
# -------------------------------------------------------------

# Size of grid
gsize <- 7

v_phi <- seq(0, 0.45, length=gsize)
v_sd_int <- seq(0, 1, length=gsize)
v_Nt <- c(20, 100, 250, 1000)

# Register clusters
nCores <- gsize
cl <- makeCluster(nCores, outfile = "")
registerDoParallel(cl)

# Loop
out_P <- foreach(i = 1:gsize, .packages = c("mlVAR", "MASS"),
                 .export = c("v_phi", "v_sd_int", "v_Nt", "gsize"),
                 .verbose = TRUE) %dopar% {

                   m_out <- matrix(NA, gsize, 4)

                   for(j in 1:gsize) {
                     for(tt in 1:4) {
                       m_out[j, tt] <- getBetCor(n_subj = 500,
                                                 phi_en = v_phi[i],
                                                 Nt = v_Nt[[tt]],
                                                 sd_int = v_sd_int[j])
                     }
                   }

                   return(m_out)

                 } # end foreach

stopCluster(cl)

# Combine
a_out <- array(NA, dim=c(gsize, gsize, 4))
for(i in 1:gsize) a_out[i, , ] <- out_P[[i]]

# Save results
saveRDS(a_out, "Sim_Out.RDS")


# -------------------------------------------------------------
# ---------- Make Figure --------------------------------------
# -------------------------------------------------------------

# Get color Gradient
color.gradient <- function(x, colors=c("white","blue"), colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}
x <- 1:100
col_grad <- color.gradient(x)
a_out_int <- round(a_out, 2)*100
a_out_int[a_out_int < 0] <- 0

pdf("Fig_SimIllu.pdf", width = 8, height = 8)

par(mfrow=c(2,2))

for(tt in 1:4) {

  par(mar=c(3.5, 3.5, 3,1))
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  axis(1, at=(seq(0, 1, length=gsize+1)+1/(gsize*2))[-(gsize+1)], labels=round(v_phi, 2))
  axis(2, at=(seq(0, 1, length=gsize+1)+1/(gsize*2))[-(gsize+1)], labels=round(v_sd_int, 2), las=2)
  axis_seq <- seq(0, 1, length=gsize+1)
  title(xlab=expression(beta), line=2.5)
  title(ylab=expression(sigma[mu]), line=2.5)
  title(main=bquote(N[t] == .(v_Nt[tt])), font.main=1, cex.main=1.5, line=.5)

  for(i in 1:gsize) {
    for(j in 1:gsize){
      rect(xleft = axis_seq[i],
           ybottom = axis_seq[j],
           xright = axis_seq[i+1],
           ytop = axis_seq[j+1],
           col=col_grad[a_out_int[i,j,tt]])
      text(axis_seq[i]+1/(gsize*2),
           axis_seq[j]+1/(gsize*2),
           label=round(a_out[i,j,tt], 2),
           col=ifelse(a_out[i,j,tt]<0.5, "black", "white"))
    }
  }

}

dev.off()








