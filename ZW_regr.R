# zitongs code for mixed spline regr

#`lme4` package
#`glmer` function
mod2r <- glmer(total.spk ~ group*ns(mspeed, knots = c(15, 30), Boundary.knots
                                    = c(3, 50))
               + (1+ns(mspeed,knots = c(15, 30), Boundary.knots = c(3, 50))|rat)
               + (1+ns(mspeed,knots = c(15, 30), Boundary.knots = c(3, 50))|cell.id),
               family = poisson(),
               offset = (log(n)),
               data = dt)
#`lmer` function
mod_er <- lmer(log(rate) ~ group*ns(mdistm15, 2) +ns(mspeedm20, 3)+
                 (1+mdist15|rat),
               data = dt.tot)
#`nlme` package:
#  `lme` function
mod_lme <- lme(log(rate) ~ group*ns(mdistm15, 2),
               random = ~1+mdistm15|rat,
               data = dt.tot,
               control = lmeControl(opt = 'optim'))