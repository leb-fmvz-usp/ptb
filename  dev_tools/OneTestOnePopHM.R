library(rjags)

dataset <- list(tau=0.1,api=c(1.8),bpi=c(26.74),ase=22.5,bse=10.22,asp=88.28,bsp=1.882,Ud=c(92),x=c(1),x1=c(2),n1=c(92),N=c(91),C=1000,k=1,zeros=0)

# initials 1
inits <- list(z=c(0),se=0.7,sp=0.99,dI=c(3))


model{
    for (i in 1:k) {
        n[i] <- n1[i] - 1
        zeros[i] ~ dpois(phi[i])
        phi[i] <- -log(pdfx[i]) + C
        Uy[i] <- min(n[i], d[i])
        Ly[i] <- max(0, n[i] - N[i] + d[i])
        logNn[i] <- logfact(N[i]) - logfact(n[i]) - logfact(N[i] - n[i])
        for (y in 1:Ud[i]) {
            # (y-1) from 0 to (Ud[i]-1) Ud[i]=n[i]+1
            yI[i,y] <- step(Uy[i] - y + 1) * step(y - 1 - Ly[i])
            temy1[i,y] <- yI[i,y] * (logfact(d[i]) - logfact(y - 1) -
                                         logfact(yI[i,y] *(d[i] - y + 1)))
            temy2[i,y] <- yI[i,y] *
                (logfact(N[i] - d[i]) -
                     logfact(n[i] - y + 1) -
                     logfact(yI[i,y] *(N[i] - d[i] - n[i] + y - 1)))
            hypdf[i,y] <-
                yI[i,y] * exp(temy1[i,y] + temy2[i,y] - yI[i,y] * logNn[i])
            Uj[i,y] <- min(x[i],y - 1)
            Lj[i,y] <- max(0,x[i] - n[i] + y - 1)
            for (j in 1:x1[i]) {
                jI[i,y,j] <- step(Uj[i,y] - j + 1) * step(j - 1 - Lj[i,y])
                tem1[i,y,j] <- jI[i,y,j] *
                    (logfact(y - 1) - logfact(j - 1) -
                         logfact(jI[i,y,j] * (y - j)) + (j - 1) *
                         log(se) + (y - j) * log(1 - se))
                tem2[i,y,j] <- jI[i,y,j] *
                    (logfact(n[i] - y + 1) -
                         logfact(x[i] - j + 1) -
                         logfact(jI[i,y,j] *(n[i] - y - x[i] + j)) +
                         (n[i] - y - x[i] + j) * log(sp) + (x[i] - j + 1) *
                         log(1 - sp))
                sj[i,y,j] <- jI[i,y,j] * exp(tem1[i,y,j] + tem2[i,y,j])
            }
            sumj[i,y] <- sum(sj[i,y,])
        }
        pdfx[i] <- inprod(hypdf[i,],sumj[i,])
        #d[i]~dunif(Ld[i],Ud[i])
        d[i] <- z[i] * dI[i]
        dI[i] ~ dcat(pi0[i,])
        pi[i] <- d[i] / N[i]
        for (id in 1:N[i]) {
            p0[i,id] <- pow(id / N[i],api[i] - 1) * pow(1 - id / N[i],bpi[i] - 1)
        }
        p0.sum[i] <- sum(p0[i,])
        for (id in 1:N[i]) {
            pi0[i,id] <- p0[i,id] / p0.sum[i]
        }
        z[i] ~ dbern(tau)
        z0[i] <- 1 - z[i]
    }
    se ~ dbeta(ase,bse)
    sp ~ dbeta(asp,bsp)
}

writeLines(model, 'model.txt')

# Initialize.
jags.mod <- jags.model(file = 'model.txt',
                       data = dataset,
                       inits = inits, 
                       n.chains = 3)

update(jags.mod, n.iter = 1e3)

# Posterior samples.
prev.est <- coda.samples(jags.mod,
                         variable.names = c('pi[1]', 'z0[1]'),
                         n.iter = 10e3)

summary(prev.est)
