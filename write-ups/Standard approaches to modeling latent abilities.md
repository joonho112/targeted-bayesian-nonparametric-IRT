---
output:
  word_document: default
  html_document: default
---
# Standard approaches to modeling latent abilities (Export)

# Basic setup: the Bayesian Rasch model

In the context of educational assessment, typical data are the responses of a sample of test takers on a set of test items. While a broad set of IRT models is available to analyze such data, this study focuses on the Rasch model (Rasch, 1960) under a Bayesian framework (Fox, 2010; Levy & Mislevy, 2016; Debelak et al., 2022; more citations?). Let $y_{ip}$ denotes person $p$’s response to item $i$ for $p = 1, \dots, N$  and  $i = 1, \dots, I$. $y_{ip}=1$ indicates the correct response for person $p$ at the $i$th item, whereas $y_{ip}=0$ represents the incorrect answer. 

Under a Bayesian hierarchical framework, the first stage of the Rasch model can be written as

$$
\begin{align*} y_{ip} \sim \text{Bernoulli}(\pi_{ip}), \\
\text{logit}(\pi_{ip})=\theta_{p}-\beta_{i}, \\  \tag{1}
\end{align*}
$$

where $\pi_{ip}$ denotes the probability that person $p$ answers item $i$ correctly given the model parameters $\theta_p$ and $\beta_i$, which respectively represent the latent trait of the $p$-th individual and the item difficulty for the $i$-th item. Since the model parameters are both unobservable, researcher only have access to their estimates, $\widehat{\theta}_p$ and $\widehat{\beta}_i$, and the corresponding squared standard errors, $se(\widehat{\theta}_p)^2$ and $se(\widehat{\beta}_i)^2$. In particular, the $\widehat{\theta}_p$ and $se(\widehat{\theta}_p)^2$ can be obtained by maximum likelihood (ML) estimation using only the *within-person* observed response data from person $p$ for the Rasch model. $se(\widehat{\theta}_p)^2$ is estimated as an inverse of the following Fisher information function (Fischer & Molenaar, 1995):  

$$
\mathcal{I}(\theta_p) =1/ se(\widehat{\theta}_p)^2 = \sum^{I}_{i = 1}\{\pi_{ip} \cdot (1-\pi_{ip})\}. \tag{2}
$$

The test information function $\mathcal{I}(\theta_p)$ describes the amount of information provided by the item responses at a specific point on the person parameter continuum (Debelak et al., 2022). Accordingly, larger $\mathcal{I}(\theta_p)$ means smaller sampling error variance of each $\widehat{\theta}_p$ around its true person-specific latent trait $\theta_p$, that is, smaller $se(\widehat{\theta}_p)^2$. 

The second stage of the Bayesian Rasch model assumes that $\theta_p$ are independent and identically distributed with a certain prior distribution $G$, 

$$
\theta_p \overset{iid}{\sim} G \equiv N(\mu_{\theta}, \sigma^2_{\theta}), \quad p = 1, \dots , N.\tag{3}
$$

The prior distribution $G$ is unknown in general, but for the Rasch model, we typically specify $G$ as a Gaussian distribution with two hyperparameters: $\mu_{\theta}$ and $\sigma_{\theta}^2$, the mean and variance of true person-specific latent traits $\theta_{p}$’s, respectively, both defined at the population level. For identification, $\mu_{\theta}$ is constrained to be zero to fix the location of the latent scale. $\sigma_{\theta}^2$ can be assumed to have a prior distribution such as an inverse-gamma distribution (Fox, 2010; Paulon et al., 2018; Paganin et al., 2021) or an inverse-$\chi^2$ distribution (Debelak et al., 2022). 

The second stage model also assumes that the item difficulty parameter $\beta_i$’s are independent draws from a Gaussian distribution with a mean of $\mu_{\beta}$ and a variance of $\sigma_{\beta}^2$:

$$
\beta_i \overset{iid}{\sim}  N(\mu_{\beta}, \sigma^2_{\beta}), \quad i = 1, \dots , I.\tag{4}
$$

A non-informative Gaussian prior with a large variance or an improper uniform prior density can be used for the hyperprior distribution of $\mu_{\beta}$. Again, $\sigma_{\beta}^2$ can have either an inverse-gamma or an inverse-$\chi^2$ prior distribution.  

# Estimation of person-specific latent traits

The estimation of hyperparameters such as $\mu_{\theta}$ and $\sigma_{\theta}$ using Gaussian hierarchical models is found to be robust to misspecification of the prior distribution $G$ (McCulloch & Neuhaus, 2011). The focus of this study is on inferences for $\theta_{p}$, however, which can be sensitive to misspecification of $G$ (citations needed). Under the model’s normality assumptions, the conditional posterior distribution of $\theta_p$ given the person’s observed responses, the item parameters, and the hyperparameters is normal (citation needed): 

$$
\theta_p|\widehat{\theta}_p, \boldsymbol{\beta}, \sigma_{\theta}^2, \mu_{\beta}, \sigma_{\beta}^2 \sim N(\theta_p^{*}, V_p) \qquad p = 1, \dots, N, 
$$

where $\theta_p^{*}$ and $V_p$ are the conditional posterior mean and variance, respectively. When ML estimates for the item parameters are plugged in, $\theta_p^{*}$ becomes an expected a posteriori (EAP) estimate (citation needed). The conditional posterior mean $\theta_p^{*}$ is a weighted average of the prior mean (i.e., $\mu_{\theta}$ fixed at zero) and the ML estimate $\widehat{\theta}_p$ obtained using the person $p$’s observed response data. The weights are given by the so-called *precisions*, the inverse of the $\sigma_{\theta}^2$ and $se(\widehat{\theta}_p)^2$. The posterior variance $V_p$ is proportional to the sum of the two precisions. If $se(\widehat{\theta}_p)^2 = 0$, the ML estimator of $\widehat{\theta}_p$ is perfectly precise and reliable, and thus the posterior mean and the ML estimator are identical ($\theta_p^{*} = \widehat{\theta}_p$). A large $se(\widehat{\theta}_p)^2$ indicates relatively less informative data about the $\theta_p$ than its prior distribution, which results in the posterior mean $\theta_p^{*}$ shrunken more toward the fixed prior mean zero.

The key insight of the Bayesian hierarchical Rasch model is that the observed variation in estiamted person-specific latent traits, $\text{var}(\widehat{\theta}_p)$, reflects two sources of variation: (1) the between-person heterogeneity in true latent traits $\theta_p$ ($\sigma_{\theta}^2$), and (2) the within-person sampling variation of each $\widehat{\theta}_p$ around its $\theta_p$ ($se(\widehat{\theta}_p)^2$). The *shrinkage factor*, $\sigma_{\theta}^2/(\sigma_{\theta}^2 + se(\widehat{\theta}_p)^2)$, is defined as the proportion of the variance of the ML estimator that is due to the underlying heterogeneity across persons. If the shrinkage factor is larger than 0.5 and close to 1.0 for person $p$, it indicates that $se(\widehat{\theta}_p)^2$ is smaller than $\sigma_{\theta}^2$, suggesting more weights on the ML estimate $\widehat{\theta}_p$ obtained using the person $p$’s observed response data. To summarize how informative the $\widehat{\theta}_p$’s are *on average*, we calculate the mean of the squared standard errors $se(\widehat{\theta}_p)^2$ over all $N$ respondents and plug it into the shrinkage factor formula: $I = \sigma_{\theta}^2/(\sigma_{\theta}^2 + \sum_{p = 1}^{P}{ se(\widehat{\theta}_p)^2})$. The average shrinkage factor $I$ is equivalent to the estimate of a test’s reliability in Rasch measurement (Wright & Masters, 1982). A large $I$ value indicates less noisy, more reliable sample-based estimates $\widehat{\theta}_p$’s. In Rasch measurement context, therefore, large $I$ values indicate that a test is more suited to distinguish between respondents in the sample with different $\theta_p$ levels.