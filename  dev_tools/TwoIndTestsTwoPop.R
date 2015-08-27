librarres.pop1(rjags)
model <- 'model{
    res.pop1 ~ dmulti(p1, pop.size1)
    res.pop2 ~ dmulti(p2, pop.size2)
    p1[1,1] <- pi1 * se1 * se2 + (1 - pi1) * (1 - sp1) * (1 - sp2)
    p1[1,2] <- pi1 * se1 * (1 - se2) + (1 - pi1) * (1 - sp1) * sp2
    p1[2,1] <- pi1 * (1 - se1) * se2 + (1 - pi1) * sp1 * (1 - sp2)
    p1[2,2] <- pi1 * (1 - se1) * (1 - se2) + (1 - pi1) * sp1 * sp2
    p2[1,1] <- pi2 * se1 * se2 + (1 - pi2) * (1 - sp1) * (1 - sp2)
    p2[1,2] <- pi2 * se1 * (1 - se2) + (1 - pi2) * (1 - sp1) * sp2
    p2[2,1] <- pi2 * (1 - se1) * se2 + (1 - pi2) * sp1 * (1 - sp2)
    p2[2,2] <- pi2 * (1 - se1) * (1 - se2) + (1 - pi2) * sp1 * sp2
    se1 ~ dbeta(2.82, 2.49)
    sp1 ~ dbeta(15.7, 1.30)
    pi2 ~ dbeta(1.73, 2.71)
    se2 ~ dbeta(8.29, 1.81)
    sp2 ~ dbeta(10.69, 2.71)
    Z ~ dbern(tau1)
    pi1star ~ dbeta(1.27, 9.65)
    pi1 <- Z * pi1star
}'

dataset <- list(pop.size1=132, pop.size2=30, res.pop1=c(3,24,0,3),
                res.pop2=c(0,3,0,129), tau1=0.95)
inits <- list(Z=1, pi1star=0.03, pi2=0.30, se1=0.55,
              sp1=0.98, se2=0.90, sp2=0.85)

writeLines(model, 'model.txt')

# Initialize.
jags.mod <- jags.model(file = 'model.txt',
                       data = dataset,
                       inits = inits, 
                       n.chains = 3)

# Burn in.
update(jags.mod, n.iter = 5e2)

# Posterior samples.
prev.est <- coda.samples(jags.mod,
                         variable.names = c('se1', 'sp1', 'se2', 'sp2',
                                            'pi1', 'pi2'),
                         n.iter = 5e3)

summary(prev.est)
