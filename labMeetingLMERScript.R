# Lab meeting LMER demo 

library(lme4)
library(splines)
library(ggplot2)
library(lmtest)
library(stringr)

# define data paths and load in data
path <- "E:\\Ratterdam\\labMeetingLMER.csv"
df <- read.csv(path,header=TRUE)


bonf <- 0.05
cipct <- 1-bonf
z <- qnorm(cipct)

#First model 
mod1 <- lm(rate ~ dFactor:cFactor, data=df)
df$fit <- predict(mod1, newdata=df, re.form=NA)

Designmat <- model.matrix(rate ~ dFactor:cFactor, df)
predvar <- diag(Designmat %*% vcov(mod1) %*% t(Designmat))
df$fitCI <- sqrt(predvar)*z

p <- ggplot(df, aes(x=cFactor, y=rate, color=as.factor(dFactor)))+
  geom_line(aes(cFactor, fit))+
  geom_ribbon(aes(y=fit, ymin=fit-fitCI, ymax=fit+fitCI, fill=as.factor(dFactor)), alpha=0.2)+
  geom_point(size=2,alpha=0.5)

print(p)

# model 2
mod2 <- lmer(rate ~ dFactor:ns(cFactor,3)+(1|rat), data=df)
df$fit <- predict(mod2, newdata=df, re.form=NA)

Designmat <- model.matrix(rate ~ dFactor:ns(cFactor,3), df)
predvar <- diag(Designmat %*% vcov(mod2) %*% t(Designmat))
df$fitCI <- sqrt(predvar)*z

p <- ggplot(df, aes(x=cFactor, y=rate, color=as.factor(dFactor)))+
  geom_line(aes(cFactor, fit))+
  geom_ribbon(aes(y=fit, ymin=fit-fitCI, ymax=fit+fitCI, fill=as.factor(dFactor)), alpha=0.2)+
  geom_point(size=2,alpha=0.5)

print(p)

# model 3
mod3 <- lmer(rate ~ dFactor:ns(cFactor,3)+(1+dFactor|rat), data=df)
df$fit <- predict(mod3, newdata=df, re.form=NA)

Designmat <- model.matrix(rate ~ dFactor:ns(cFactor,3), df)
predvar <- diag(Designmat %*% vcov(mod3) %*% t(Designmat))
df$fitCI <- sqrt(predvar)*z

p <- ggplot(df, aes(x=cFactor, y=rate, color=as.factor(dFactor)))+
  geom_line(aes(cFactor, fit))+
  geom_ribbon(aes(y=fit, ymin=fit-fitCI, ymax=fit+fitCI, fill=as.factor(dFactor)), alpha=0.2)+
  geom_point(size=2,alpha=0.5)

print(p)

# model 4
mod4 <- lmer(rate ~ dFactor:ns(cFactor,3)+(1+cFactor|rat), data=df)
df$fit <- predict(mod4, newdata=df, re.form=NA)

Designmat <- model.matrix(rate ~ dFactor:ns(cFactor,3), df)
predvar <- diag(Designmat %*% vcov(mod4) %*% t(Designmat))
df$fitCI <- sqrt(predvar)*z

p <- ggplot(df, aes(x=cFactor, y=rate, color=as.factor(dFactor)))+
  geom_line(aes(cFactor, fit))+
  geom_ribbon(aes(y=fit, ymin=fit-fitCI, ymax=fit+fitCI, fill=as.factor(dFactor)), alpha=0.2)+
  geom_point(size=2,alpha=0.5)

print(p)

# model 4a 
mod4a <- lmer(rate ~ dFactor*cFactor+(1+ns(cFactor,3)|rat), data=df)
df$fit <- predict(mod4a, newdata=df, re.form=NA)

Designmat <- model.matrix(rate ~ dFactor*cFactor, df)
predvar <- diag(Designmat %*% vcov(mod4a) %*% t(Designmat))
df$fitCI <- sqrt(predvar)*z

p <- ggplot(df, aes(x=cFactor, y=rate, color=as.factor(dFactor)))+
  geom_line(aes(cFactor, fit))+
  geom_ribbon(aes(y=fit, ymin=fit-fitCI, ymax=fit+fitCI, fill=as.factor(dFactor)), alpha=0.2)+
  geom_point(size=2,alpha=0.5)

print(p)


# model 5
mod5 <- lmer(rate ~ dFactor*ns(cFactor,3)+(1+cFactor|rat)+ (1+dFactor|rat), data=df)
df$fit <- predict(mod5, newdata=df, re.form=NA)

Designmat <- model.matrix(rate ~ dFactor*ns(cFactor,3), df)
predvar <- diag(Designmat %*% vcov(mod5) %*% t(Designmat))
df$fitCI <- sqrt(predvar)*z

p <- ggplot(df, aes(x=cFactor, y=rate, color=as.factor(dFactor)))+
  geom_line(aes(cFactor, fit))+
  geom_ribbon(aes(y=fit, ymin=fit-fitCI, ymax=fit+fitCI, fill=as.factor(dFactor)), alpha=0.2)+
  geom_point(size=2,alpha=0.5)

print(p)



#Full model 
mod6 <- lmer(rate ~ dFactor*ns(cFactor,3)+(1+dFactor|rat)+(1+cFactor|rat), data=df)
df$fit <- predict(mod6, newdata=df, re.form=NA)

Designmat <- model.matrix(rate ~ dFactor*ns(cFactor,3), df)
predvar <- diag(Designmat %*% vcov(mod6) %*% t(Designmat))
df$fitCI <- sqrt(predvar)*z

p <- ggplot(df, aes(x=cFactor, y=rate, color=as.factor(dFactor)))+
  geom_line(aes(cFactor, fit))+
  geom_ribbon(aes(y=fit, ymin=fit-fitCI, ymax=fit+fitCI, fill=as.factor(dFactor)), alpha=0.2)+
  geom_point(size=2,alpha=0.5)

print(p)

# GLMER model 

modg <- glmer(rate+10 ~ dFactor:ns(cFactor,3)+(1+dFactor|rat)+(1+cFactor|rat),
             family=Gamma, 
             data = df, 
             control = glmerControl(optimizer='Nelder_Mead'),
             nAGQ=0
)

Designmat <- model.matrix(rate ~ dFactor:ns(cFactor,3), df)
predvar <- diag(Designmat %*% vcov(modg) %*% t(Designmat))
df$fitCI <- sqrt(predvar)*z

p <- ggplot(df, aes(x=cFactor, y=rate, color=as.factor(dFactor)))+
  geom_line(aes(cFactor, fit))+
  geom_ribbon(aes(y=fit, ymin=fit-fitCI, ymax=fit+fitCI, fill=as.factor(dFactor)), alpha=0.2)+
  geom_point(size=2,alpha=0.5)

print(p)


# Plots building up the model 

# hist of rates
p <- ggplot(df,aes(x=rate))+
  geom_histogram()
print(p)

# overall C 
(prelim_plot <- ggplot(df, aes(x = cFactor, y = rate)) +
    geom_point(alpha=0.5) +
    geom_smooth(method = "lm"))
print(prelim_plot)

# overall plot of D 
boxplot(rate~dFactor,data=df)


# overall colored by rat / factor
(colour_plot <- ggplot(df, aes(x = cFactor, y = rate, colour = as.factor(rat))) +
    geom_point(size = 2))

(colour_plot <- ggplot(df, aes(x = cFactor, y = rate, colour = as.factor(dFactor))) +
    geom_point(size = 2))

# plot residuals and QQ
basic_c_LM <- lm(rate ~ cFactor, data=df)
basic_d_LM <- lm(rate ~ dFactor, data=df)
plot(basic_c_LM,which=1)
plot(basic_c_LM,which=2)
plot(basic_d_LM,which=1)
plot(basic_d_LM,which=2)

# boxplot by groups
boxplot(rate~rat,data=df)

# plot overall colored by rat / factor
p <- ggplot(df, aes(x=cFactor, y=rate, color=as.factor(rat)))+
  geom_point(aes(cFactor, fit))
print(p)

# plot by rat 
facets <- ggplot(aes(cFactor, rate), data=df)+
  geom_point()+
  facet_wrap(~rat)
print(facets)

# plot by rat/dFactor with fits

(mm_plot <- ggplot(df, aes(x = cFactor, y = rate, colour = as.factor(dFactor))) +
    facet_wrap(~rat, nrow=2) +   # a panel for each mountain range
    geom_point(alpha = 0.5) +
    theme_classic() +
    geom_line(data = cbind(df, pred = predict(mod6)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
    theme(legend.position = "none",
          panel.spacing = unit(2, "lines"))  # adding space between panels
)
print(mm_plot)

# Statistics 
ci <- confint(mod6, method='Wald', level=cipct)
fe <- fixef(mod6)
