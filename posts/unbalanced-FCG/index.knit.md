---
title: "Unbalanced Four Core Genotype Analyses"
description: "You want to run a four core genotype analysis without having all four genotypes"
author: "Jason Lerch and Kamila Szulc-Lerch"
date: "2025-02-08"
draft: true
categories: [R]
format:
  html:
    code-fold: true
    code-summary: "Show the code"
---



This came out of a discussion related to a rat model designed to explore the same question at the Four Core Genotype mouse model: what happens if you don't have all four arms? Are you hopelessly confounded? The short answer is no under some assumptions. We'll explore that with some simple simulations.

To start with we'll create a function to simulate the four core genotype model and then drop one of the arms.



::: {.cell}

```{.r .cell-code}
library(tidyverse)
library(ggplot2)
library(broom)
```
:::

::: {.cell}

```{.r .cell-code}
generateFCGData <- function(
    XXM = 0,
    XXF = 0,
    XYM = 0,
    XYF = 0,
    std = 1,
    ngroup=10) {
  rbind(
    data.frame(group="XXM", gonad="T", chromosomes="XX", 
               volume = rnorm(ngroup, XXM, sd=std)),
    data.frame(group="XXF", gonad="O", chromosomes="XX", 
               volume = rnorm(ngroup, XXF, sd=std)),
    data.frame(group="XYM", gonad="T", chromosomes="XY", 
               volume = rnorm(ngroup, XYM, sd=std)),
    data.frame(group="XYF", gonad="O", chromosomes="XY", 
               volume = rnorm(ngroup, XYF, sd=std))
  )
}
```
:::



Let's assume that having male gonads increases the output variable by 2 standard deviations. That would look as follows:



::: {.cell}

```{.r .cell-code}
ddd <- generateFCGData(XYM=2, XXM=2, ngroup = 500)
  
ggplot(ddd) + 
  aes(x=group, y=volume, fill=chromosomes) + 
  geom_boxplot() + 
  facet_wrap(~gonad, scales="free_x") +
  scale_fill_brewer(palette = "Set1") + 
  theme_classic()
```

::: {.cell-output-display}
![](index_files/figure-html/unnamed-chunk-3-1.png){width=672}
:::
:::



Let's run the stats three ways, first separately fitting gonads and chromosomes and then using an additive model incorporating both chromosomes and gonads:



::: {.cell}

```{.r .cell-code}
summary( ( lfg <- lm(volume ~ gonad, ddd)))
```

::: {.cell-output .cell-output-stdout}

```

Call:
lm(formula = volume ~ gonad, data = ddd)

Residuals:
    Min      1Q  Median      3Q     Max 
-3.7156 -0.6636 -0.0165  0.6493  4.7112 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept) 0.009019   0.031619   0.285    0.775    
gonadT      2.005040   0.044716  44.839   <2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 0.9999 on 1998 degrees of freedom
Multiple R-squared:  0.5016,	Adjusted R-squared:  0.5013 
F-statistic:  2011 on 1 and 1998 DF,  p-value: < 2.2e-16
```


:::

```{.r .cell-code}
summary( ( lfc <- lm(volume ~ chromosomes, ddd)))
```

::: {.cell-output .cell-output-stdout}

```

Call:
lm(formula = volume ~ chromosomes, data = ddd)

Residuals:
    Min      1Q  Median      3Q     Max 
-4.1828 -1.0775 -0.0058  1.0610  3.7782 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept)   1.009075   0.044786  22.531   <2e-16 ***
chromosomesXY 0.004928   0.063337   0.078    0.938    
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 1.416 on 1998 degrees of freedom
Multiple R-squared:  3.03e-06,	Adjusted R-squared:  -0.0004975 
F-statistic: 0.006054 on 1 and 1998 DF,  p-value: 0.938
```


:::

```{.r .cell-code}
summary( ( lfgc <- lm(volume ~ gonad + chromosomes, ddd)))
```

::: {.cell-output .cell-output-stdout}

```

Call:
lm(formula = volume ~ gonad + chromosomes, data = ddd)

Residuals:
    Min      1Q  Median      3Q     Max 
-3.7131 -0.6636 -0.0182  0.6503  4.7088 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept)   0.006555   0.038735   0.169    0.866    
gonadT        2.005040   0.044727  44.828   <2e-16 ***
chromosomesXY 0.004928   0.044727   0.110    0.912    
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 1 on 1997 degrees of freedom
Multiple R-squared:  0.5016,	Adjusted R-squared:  0.5011 
F-statistic:  1005 on 2 and 1997 DF,  p-value: < 2.2e-16
```


:::
:::



The output here is correct - no matter how analyzed, the statistics tell us that having testes increases the output by 2 (everything is in units of standard deviation here), and there is no effect of chromosomes.

Let's now drop one of the arms, and assume that rather than having all four groups we do not have the XY female group. The exact same data would look like this:



::: {.cell}

```{.r .cell-code}
ddd <- ddd %>% filter(group!="XYF")
ggplot(ddd) + 
  aes(x=group, y=volume, fill=chromosomes) + 
  geom_boxplot() + 
  facet_wrap(~gonad, scales="free_x") +
  scale_fill_brewer(palette = "Set1") + 
  theme_classic()
```

::: {.cell-output-display}
![](index_files/figure-html/unnamed-chunk-5-1.png){width=672}
:::
:::



And we run the same statistics:



::: {.cell}

```{.r .cell-code}
summary( ( lrg <- lm(volume ~ gonad, ddd)))
```

::: {.cell-output .cell-output-stdout}

```

Call:
lm(formula = volume ~ gonad, data = ddd)

Residuals:
    Min      1Q  Median      3Q     Max 
-3.7156 -0.6517 -0.0149  0.6648  3.1468 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept) -0.00683    0.04411  -0.155    0.877    
gonadT       2.02089    0.05403  37.405   <2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 0.9864 on 1498 degrees of freedom
Multiple R-squared:  0.4829,	Adjusted R-squared:  0.4826 
F-statistic:  1399 on 1 and 1498 DF,  p-value: < 2.2e-16
```


:::

```{.r .cell-code}
summary( ( lrc <- lm(volume ~ chromosomes, ddd)))
```

::: {.cell-output .cell-output-stdout}

```

Call:
lm(formula = volume ~ chromosomes, data = ddd)

Residuals:
    Min      1Q  Median      3Q     Max 
-3.7430 -0.9230 -0.0218  0.9341  3.6675 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept)    1.00908    0.04077   24.75   <2e-16 ***
chromosomesXY  0.99406    0.07061   14.08   <2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 1.289 on 1498 degrees of freedom
Multiple R-squared:  0.1169,	Adjusted R-squared:  0.1163 
F-statistic: 198.2 on 1 and 1498 DF,  p-value: < 2.2e-16
```


:::

```{.r .cell-code}
summary( ( lrgc <- lm(volume ~ gonad + chromosomes, ddd)))
```

::: {.cell-output .cell-output-stdout}

```

Call:
lm(formula = volume ~ gonad + chromosomes, data = ddd)

Residuals:
    Min      1Q  Median      3Q     Max 
-3.7265 -0.6531 -0.0141  0.6660  3.1468 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept)   -0.00683    0.04413  -0.155    0.877    
gonadT         2.03181    0.06240  32.559   <2e-16 ***
chromosomesXY -0.02184    0.06240  -0.350    0.726    
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 0.9867 on 1497 degrees of freedom
Multiple R-squared:  0.483,	Adjusted R-squared:  0.4823 
F-statistic: 699.2 on 2 and 1497 DF,  p-value: < 2.2e-16
```


:::
:::



Losing one arm of the study means that the design is now unbalanced. Looking at the three different stats models shows that modelling only gonads gives the correct answer of a 2 SD increase, but modelling only chromosomes gives an incorrect answer of an increase of 1SD when having XY chromosomes. What saves us is covarying - modelling gonads + chromosomes again returns the correct answer, since the model simultaneously estimates the effects of gonads and chromosomes while controlling for the other.

Let's show that graphically:



::: {.cell}

```{.r .cell-code}
allStats <- rbind(
  tidy(lfg, conf.int = T) %>% mutate(model = "FCG: G"),
  tidy(lfc, conf.int = T) %>% mutate(model = "FCG: C"),
  tidy(lfgc, conf.int = T) %>% mutate(model = "FCG: G + C"),
  tidy(lrg, conf.int = T) %>% mutate(model = "uFCG: G"),
  tidy(lrc, conf.int = T) %>% mutate(model = "uFCG: C"),
  tidy(lrgc, conf.int = T) %>% mutate(model = "uFCG: G + C")
) %>% filter(term != "(Intercept)")

cols <- RColorBrewer::brewer.pal(2, "Set2")
```

::: {.cell-output .cell-output-stderr}

```
Warning in RColorBrewer::brewer.pal(2, "Set2"): minimal value for n is 3, returning requested palette with 3 different levels
```


:::

```{.r .cell-code}
allStats %>% 
  ggplot() + 
  aes(x=term, y=estimate, colour=term) + 
  geom_point() + 
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high)) + 
  facet_wrap(.~model, scale="free_x", nrow=1) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_color_manual(values = cols[1:2]) + 
  geom_hline(yintercept = 2, color=cols[2]) +
  geom_hline(yintercept = 0, color=cols[1]) +
  labs(title="Model effects", caption="FCG = balanced four core genotypes; uFCG = unbalanced four core genotypes (no XYF); C = chromosomes; G = gonads") 
```

::: {.cell-output-display}
![](index_files/figure-html/unnamed-chunk-7-1.png){width=672}
:::
:::



Now, can't we just do a simple t-test comparing two groups with the desired contrasts which avoids the issue of more complex modelling? Yes, that works:



::: {.cell}

```{.r .cell-code}
summary(lm(volume ~ group, data = ddd %>% filter(group %in% c("XXF", "XXM"))))
```

::: {.cell-output .cell-output-stdout}

```

Call:
lm(formula = volume ~ group, data = ddd %>% filter(group %in% 
    c("XXF", "XXM")))

Residuals:
    Min      1Q  Median      3Q     Max 
-3.7265 -0.6673 -0.0127  0.6692  3.1468 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept) -0.00683    0.04414  -0.155    0.877    
groupXXM     2.03181    0.06242  32.551   <2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 0.9869 on 998 degrees of freedom
Multiple R-squared:  0.515,	Adjusted R-squared:  0.5145 
F-statistic:  1060 on 1 and 998 DF,  p-value: < 2.2e-16
```


:::
:::



But do we pay a price in statistical power? Let's compare the simple two group comparison to the linear model with all three groups (in the unbalanced FCG case) estimating both a gonad and a chromosome effect.



::: {.cell}

```{.r .cell-code}
simNtimes <- function(nsims=100, effectSize=2, ngroup=20) {
  bFCG <- map(1:nsims, ~ generateFCGData(XYM=effectSize, 
                                         XXM=effectSize, 
                                         ngroup=ngroup))
  uFCG <- map(bFCG, ~ .x %>% filter(group!="XYF"))
  tgroup <- map(bFCG, ~ .x %>% filter(group %in% c("XXF", "XXM")))
  
  lbFCG <- map(bFCG, ~ lm(volume ~ gonad + chromosomes, data=.x))
  luFCG <- map(uFCG, ~ lm(volume ~ gonad + chromosomes, data=.x))
  ltgroup <- map(tgroup, ~ lm(volume ~ group, data=.x))
  
  tlbFCG <- map(lbFCG, tidy)
  tluFCG <- map(luFCG, tidy)
  tltgroup <- map(ltgroup, tidy)
  
  data.frame(
    bFCG = map_dbl(tlbFCG, ~.x %>% filter(term == "gonadT") %>% pull(p.value)),
    uFCG = map_dbl(tluFCG, ~.x %>% filter(term == "gonadT") %>% pull(p.value)),
    twogroup = map_dbl(tltgroup, ~.x %>% filter(term == "groupXXM") %>% pull(p.value))
  )
  
}

allSims <- map(seq(0, 2, by=0.1), ~ simNtimes(nsims=200, effectSize=.x, ngroup = 10) %>% mutate(d=.x))
powerDF <- map_dfr(allSims, ~ .x) %>% group_by(d) %>% summarise_at(vars(bFCG:twogroup), function(x) mean(x<0.05))

powerDF %>% 
  pivot_longer(bFCG:twogroup, values_to = "power") %>% 
  ggplot() + aes(x=d, y=power, colour=name) + 
  geom_smooth(method="gam", formula=y~s(x, k=8), se=F) + 
  scale_color_brewer("Model", palette = "Dark2", 
                     label=c("balanced FCG", 
                             "unbalanced FCG", 
                             "two group comparison")) + 
  theme_classic() + 
  xlab("Simulated effect size") + 
  labs(title="Power of gonad effect", caption ="Assuming n=10 per group")
```

::: {.cell-output-display}
![](index_files/figure-html/unnamed-chunk-9-1.png){width=672}
:::
:::



Not really - this shows the power we lose when moving from a balanced FCG dataset to the unbalanced one missing one group. But the power is the same between the two ways of modelling the unbalanced dataset, indicating that in the case of the unbalanced design (i.e. only having three rather than all four groups of the FCG model) the gonad and chromosome terms of the statistical model essentially reduce down to a two-group comparison.

Let's expand this to the FCG-like rat model (male arm only for now). There are four genotypes:

| Genotype | Gonads  | Sex Chromosomes | Sry            |
|----------|---------|-----------------|----------------|
| XX       | Ovaries | XX              | None           |
| XX-Sry   | Testes  | XX              | Transgene      |
| XY       | Testes  | XY              | WT             |
| XY-Sry   | Testes  | XY              | WT + Transgene |

Once again we need a function to simulate from this dataset



::: {.cell}

```{.r .cell-code}
generateFCGRatData <- function(
    XX = 0,
    XXSry = 0,
    XY = 0,
    XYSry = 0,
    std = 1,
    ngroup=10) {
  rbind(
    data.frame(group="XX", gonad="O", chromosomes="XX", sry=0,
               volume = rnorm(ngroup, XX, sd=std)),
    data.frame(group="XXSry", gonad="T", chromosomes="XX", sry=1,
               volume = rnorm(ngroup, XXSry, sd=std)),
    data.frame(group="XY", gonad="T", chromosomes="XY", sry=1,
               volume = rnorm(ngroup, XY, sd=std)),
    data.frame(group="XYSry", gonad="T", chromosomes="XY", sry=2,
               volume = rnorm(ngroup, XYSry, sd=std))
  )
}
```
:::



Let's see what a few different scenarios look like:



::: {.cell}

```{.r .cell-code}
gonads <- generateFCGRatData(XX=0, XXSry=2, XY=2, XYSry=2, ngroup = 500) %>%
  mutate(contrast = "Gonad")
chrom <- generateFCGRatData(XX=0, XXSry=0, XY=2, XYSry=2, ngroup = 500) %>%
  mutate(contrast = "Chromosomes")
sry <- generateFCGRatData(XX=0, XXSry=2, XY=2, XYSry=4, ngroup = 500) %>%
  mutate(contrast = "Sry")
  
rbind(gonads, chrom, sry) %>%
ggplot() + 
  aes(x=group, y=volume, fill=chromosomes) + 
  geom_boxplot() + 
  facet_wrap(~gonad, scales="free_x") +
  scale_fill_brewer(palette = "Set1") + 
  facet_grid(contrast ~ .) + 
  theme_classic()
```

::: {.cell-output-display}
![](index_files/figure-html/unnamed-chunk-11-1.png){width=672}
:::
:::



Can we appropriately detect the right effects in those models, focusing first on analyzing gonads and chromosomes?



::: {.cell}

```{.r .cell-code}
summary(lm(volume ~ gonad + chromosomes, data=gonads))
```

::: {.cell-output .cell-output-stdout}

```

Call:
lm(formula = volume ~ gonad + chromosomes, data = gonads)

Residuals:
    Min      1Q  Median      3Q     Max 
-3.7327 -0.7017  0.0092  0.6782  3.4847 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept)   -0.01110    0.04458  -0.249    0.803    
gonadT         1.99798    0.06305  31.687   <2e-16 ***
chromosomesXY  0.03613    0.05461   0.662    0.508    
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 0.997 on 1997 degrees of freedom
Multiple R-squared:  0.4359,	Adjusted R-squared:  0.4353 
F-statistic: 771.6 on 2 and 1997 DF,  p-value: < 2.2e-16
```


:::

```{.r .cell-code}
summary(lm(volume ~ gonad + chromosomes, data=chrom))
```

::: {.cell-output .cell-output-stdout}

```

Call:
lm(formula = volume ~ gonad + chromosomes, data = chrom)

Residuals:
    Min      1Q  Median      3Q     Max 
-3.2552 -0.6764  0.0075  0.6762  3.9444 

Coefficients:
                Estimate Std. Error t value Pr(>|t|)    
(Intercept)   -0.0002106  0.0452218  -0.005    0.996    
gonadT         0.0081282  0.0639532   0.127    0.899    
chromosomesXY  2.0090115  0.0553851  36.273   <2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 1.011 on 1997 degrees of freedom
Multiple R-squared:  0.4981,	Adjusted R-squared:  0.4976 
F-statistic: 990.8 on 2 and 1997 DF,  p-value: < 2.2e-16
```


:::
:::



Yes. What about treating sry as an ordered factor (i.e. assuming that a dose of 2 is more than a dose of 1, but not assuming that it is purely linear).



::: {.cell}

```{.r .cell-code}
summary(lm(volume ~ sry + chromosomes, data=gonads %>% mutate(sry=as.ordered(sry))))
```

::: {.cell-output .cell-output-stdout}

```

Call:
lm(formula = volume ~ sry + chromosomes, data = gonads %>% mutate(sry = as.ordered(sry)))

Residuals:
    Min      1Q  Median      3Q     Max 
-3.7327 -0.6991  0.0081  0.6783  3.4847 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept)    1.31744    0.03933  33.497   <2e-16 ***
sry.L          1.40547    0.06307  22.285   <2e-16 ***
sry.Q         -0.81990    0.03641 -22.517   <2e-16 ***
chromosomesXY  0.04130    0.06307   0.655    0.513    
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 0.9972 on 1996 degrees of freedom
Multiple R-squared:  0.4359,	Adjusted R-squared:  0.4351 
F-statistic: 514.1 on 3 and 1996 DF,  p-value: < 2.2e-16
```


:::

```{.r .cell-code}
summary(lm(volume ~ sry + chromosomes, data=chrom %>% mutate(sry=as.ordered(sry))))
```

::: {.cell-output .cell-output-stdout}

```

Call:
lm(formula = volume ~ sry + chromosomes, data = chrom %>% mutate(sry = as.ordered(sry)))

Residuals:
    Min      1Q  Median      3Q     Max 
-3.2552 -0.6742  0.0073  0.6762  3.9444 

Coefficients:
               Estimate Std. Error t value Pr(>|t|)    
(Intercept)   -0.004533   0.039890  -0.114    0.910    
sry.L         -0.014917   0.063966  -0.233    0.816    
sry.Q         -0.015249   0.036931  -0.413    0.680    
chromosomesXY  2.023624   0.063966  31.636   <2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 1.011 on 1996 degrees of freedom
Multiple R-squared:  0.4981,	Adjusted R-squared:  0.4974 
F-statistic: 660.4 on 3 and 1996 DF,  p-value: < 2.2e-16
```


:::

```{.r .cell-code}
summary(lm(volume ~ sry + chromosomes, data=sry %>% mutate(sry=as.ordered(sry))))
```

::: {.cell-output .cell-output-stdout}

```

Call:
lm(formula = volume ~ sry + chromosomes, data = sry %>% mutate(sry = as.ordered(sry)))

Residuals:
    Min      1Q  Median      3Q     Max 
-3.5094 -0.6478 -0.0062  0.6706  3.6643 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept)   1.994594   0.038950  51.209   <2e-16 ***
sry.L         2.863991   0.062459  45.854   <2e-16 ***
sry.Q         0.033932   0.036061   0.941    0.347    
chromosomesXY 0.003616   0.062459   0.058    0.954    
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 0.9876 on 1996 degrees of freedom
Multiple R-squared:  0.6786,	Adjusted R-squared:  0.6781 
F-statistic:  1405 on 3 and 1996 DF,  p-value: < 2.2e-16
```


:::
:::



Still works but is a bit more unstable. Some interpretation aids to these models - when treating sry as an ordered factor, there are two possible trends - a linear trends, indicating that with every increase in sry there is an increase in the outcome, and a quadratic trend, indicating that going from 1 to 2 is not the same as going from 0 to 1. So in our third model where we artificially simulated an increase with every dose of sry there is a linear increase (sry.L), but no quadratic increase (sry.Q is non-significant). But in the gonad model there is both a linear and a quadratic effect, since having 2 sry copies (WT + transgene) has the same effect as having one sry copy.

What about power? First just looking at the gonads + chromosomes model:



::: {.cell}

```{.r .cell-code}
simRatNtimes <- function(nsims=100, effectSize=2, ngroup=20) {
  ratG <- map(1:nsims, ~ generateFCGRatData(XX=0,
                                            XXSry=effectSize,
                                            XY=effectSize,
                                            XYSry=effectSize,
                                            ngroup=ngroup))
  ratC <- map(1:nsims, ~ generateFCGRatData(XX=0,
                                            XXSry=0,
                                            XY=effectSize,
                                            XYSry=effectSize,
                                            ngroup=ngroup))
  
  lratG <- map(ratG, ~ lm(volume ~ gonad + chromosomes, data=.x))
  lratC <- map(ratC, ~ lm(volume ~ gonad + chromosomes, data=.x))
  
  tlratG <- map(lratG, tidy)
  tlratC <- map(lratC, tidy)
  
  data.frame(
    gonads = map_dbl(tlratG, ~.x %>% filter(term == "gonadT") %>% pull(p.value)),
    chromosomes = map_dbl(tlratC, ~.x %>% filter(term == "chromosomesXY") %>% pull(p.value))
  )
  
}

allRatSims <- map(seq(0, 2, by=0.1), ~ simRatNtimes(nsims=200, effectSize=.x, ngroup = 10) %>% mutate(d=.x))
powerDF <- map_dfr(allRatSims, ~ .x) %>% group_by(d) %>% summarise_at(vars(gonads:chromosomes), function(x) mean(x<0.05))

powerDF %>% 
  pivot_longer(gonads:chromosomes, values_to = "power") %>% 
  ggplot() + aes(x=d, y=power, colour=name) + 
  geom_smooth(method="gam", formula=y~s(x, k=8), se=F) + 
  scale_color_brewer("Model", palette = "Dark2") + 
  theme_classic() + 
  xlab("Simulated effect size") + 
  labs(title="Power of FCG-like rat study", 
       subtitle="Modelling gonads + chromosomes",
       caption ="Assuming n=10 per group")
```

::: {.cell-output-display}
![](index_files/figure-html/unnamed-chunk-14-1.png){width=672}
:::
:::



Given the unbalanced design for gonads there is some power differential, in that it is easier to detect a chromosome effect than a gonad effect.

And if we model sry as an ordered factor?



::: {.cell}

```{.r .cell-code}
simRatNtimesOF <- function(nsims=100, effectSize=2, ngroup=20) {
  ratG <- map(1:nsims, ~ generateFCGRatData(XX=0,
                                            XXSry=effectSize,
                                            XY=effectSize,
                                            XYSry=effectSize,
                                            ngroup=ngroup) %>%
    mutate(sry = as.ordered(sry))) 
  ratC <- map(1:nsims, ~ generateFCGRatData(XX=0,
                                            XXSry=0,
                                            XY=effectSize,
                                            XYSry=effectSize,
                                            ngroup=ngroup) %>%
    mutate(sry = as.ordered(sry))) 
  ratS <- map(1:nsims, ~ generateFCGRatData(XX=0,
                                            XXSry=effectSize,
                                            XY=effectSize,
                                            XYSry=effectSize*2,
                                            ngroup=ngroup) %>%
    mutate(sry = as.ordered(sry))) 
  
  lratG <- map(ratG, ~ lm(volume ~ sry + chromosomes, data=.x))
  lratC <- map(ratC, ~ lm(volume ~ sry + chromosomes, data=.x))
  lratS <- map(ratS, ~ lm(volume ~ sry + chromosomes, data=.x))
  
  tlratG <- map(lratG, tidy)
  tlratC <- map(lratC, tidy)
  tlratS <- map(lratS, tidy)
  
  data.frame(
    gonads = map_dbl(tlratG, ~.x %>% filter(term == "sry.L") %>% pull(p.value)),
    chromosomes = map_dbl(tlratC, ~.x %>% filter(term == "chromosomesXY") %>% pull(p.value)),
    sryL = map_dbl(tlratS, ~.x %>% filter(term == "sry.L") %>% pull(p.value)),
    sryQ = map_dbl(tlratS, ~.x %>% filter(term == "sry.Q") %>% pull(p.value))
  )
  
}

allRatSims <- map(seq(0, 2, by=0.1), ~ simRatNtimesOF(nsims=200, effectSize=.x, ngroup = 10) %>% mutate(d=.x))
powerDF <- map_dfr(allRatSims, ~ .x) %>% group_by(d) %>% summarise_at(vars(gonads:sryQ), function(x) mean(x<0.05))

powerDF %>% 
  pivot_longer(gonads:sryQ, values_to = "power") %>% 
  ggplot() + aes(x=d, y=power, colour=name) + 
  geom_smooth(method="gam", formula=y~s(x, k=8), se=F) + 
  scale_color_brewer("Model", palette = "Dark2") + 
  theme_classic() + 
  xlab("Simulated effect size") + 
  labs(title="Power of FCG-like rat study", 
       subtitle="Modelling sry + chromosomes",
       caption ="Assuming n=10 per group")
```

::: {.cell-output-display}
![](index_files/figure-html/unnamed-chunk-15-1.png){width=672}
:::
:::



Here the *gonads* model looks at the linear part of the sry term when there is a change in the presence of sry but no dosage effect. *sryL* looks at the linear term when there is a dosage effect, which is now the most powered term (since there are three possible steps), and *sryQ* looks at the quadratic effect when there is a simulated dosage effect which is rightly 0.

