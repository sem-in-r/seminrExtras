# Congruence Coefficients with SEMinR
library(seminr)
library(seminrExtras)

# Create measurement model ----
mobi_mm <- constructs(
  composite("Image",        multi_items("IMAG", 1:5)),
  composite("Expectation",  multi_items("CUEX", 1:3)),
  composite("Value",        multi_items("PERV", 1:2)),
  composite("Satisfaction", multi_items("CUSA", 1:3)),
  composite("Loyalty",      multi_items("CUSL", 1:3))
)

# Create structural model ----
mobi_sm <- relationships(
  paths(from = "Image",       to = c("Expectation", "Satisfaction", "Loyalty")),
  paths(from = "Expectation", to = c("Value", "Satisfaction")),
  paths(from = "Value",       to = "Satisfaction"),
  paths(from = "Satisfaction", to = "Loyalty")
)

# Estimate the model ----
mobi_pls <- estimate_pls(data = mobi,
                          measurement_model = mobi_mm,
                          structural_model  = mobi_sm)

# ============================================================================
# Congruence coefficients (congruence)
# ============================================================================
# The congruence coefficient rc describes how similarly two constructs relate
# to the other constructs in the model: the cosine similarity of their two
# columns of the construct-correlation matrix, with reliabilities on the
# diagonal (Franke, Sarstedt & Danks, 2021, Eq. 2). A value near 1 says the two
# constructs sit in nearly the same position in the nomological network.
#
# rc is an effect size. It is bounded above by 1, so sampling error can only
# move it downwards, and no significance test for rc has been validated for
# PLS-SEM. Read the magnitude; do not compare it to a cut-off.

congruence(mobi_pls)

# The reliability placed on the diagonal can be changed. Franke et al. (2021)
# specify "the reliabilities" without fixing an estimator, so all four options
# are in specification; they diverge only for Mode B constructs, where rhoA
# returns 1 because internal consistency is undefined for a composite.
congruence(mobi_pls, reliability = "one")

# summary() adds the range and the interpretation note.
summary(congruence(mobi_pls))
