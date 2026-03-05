### Accompanying Code for:
## Partial Least Squares Structural Equation Modeling (PLS-SEM) Using R (Second Edition) - A Workbook (2026)
## Hair, J.F. (Jr), Hult, T.M., Ringle, C.M., Sarstedt, M., Danks, N.P., and Adler, S.

## Chapter 2: Introduction to R and R Studio
# Create a vector of integers
vector <- c(1, 2, 3, 4, 5)

# Install the SEMinR package
install.packages(pkgs = "seminr")

# Load the SEMinR package into the environment
library(seminr)

# Install the textbook packages from CRAN
install.packages(pkgs = c("learnr", "paran", "psych", "seminrExtras"))

# Load the learnr package into the environment
library(learnr)

# Begin learning with learnr
run_tutorial()

# Searching for help using the ? operator
?read.csv

# Check all vignettes available in R
vignette()

# Load the vignette for the reshape() function from the stats package
vignette("reshape", "stats")

# Check all demos available in R
demo()

# Load the SEMinR ECSI demo
demo("seminr-pls-ecsi")
