

"The likelihood ratio test can be used for fixed or random parameters and is often used for multiparameter tests (Snijders & Bosker, 1999). For tests of random components, restricted maximum likelihood (REML) or full maximum likelihood (ML) estimations methods may be used, but for tests of the fixed components, ML must be used (Snijders & Bosker, 1999)." from Lorah 2014
"### save the result
It is important to note one primary drawback of REML
which lies in model comparisons. Because REML removes
the fixed effects from the estimation, REML uses the
restricted likelihood function that does not contain any
information about the fixed effects (Raudenbush & Bryk,

2002). Therefore, if comparing the fit of different mod-
els using REML, the fixed effects must be identical for the

comparisons of the restricted likelihood to be meaning-
ful. Otherwise, ML must be used for model comparisons

because the full likelihood includes information about
both the fixed effects and the variance components. After
competing models are compared using full ML, the final

model can be estimated with REML (and associated cor-
rections, if necessary) to protect against bias in variance

component estimates and inflated Type-I error rates.

# REML cannot be used for AIC/BIC comparison. :(
#https://stats.stackexchange.com/questions/116770/reml-or-ml-to-compare-two-mixed-effects-models-with-differing-fixed-effects-but
#However, it can be used?
#https://sci-hubtw.hkvisa.net/10.1198/000313006X90396
# refer to the previous email with mason (JAn22)
# but it can use becuase it is same data but different models?
