\name{hybrid}
\alias{hybrid}
\alias{print.hybrid}
\alias{summary.hybrid}
\title{
Hybrid model with structural shift (permits covariate-specific treatment effects)
}
\description{
James Heckman's Hybrid Model with Structural Shift (also known as the Treatment Effects Model). Jointly models outcome regression model and endogenous variable probit model (e.g., outcome associations in the presence of endogenous treatment in observational data). Can handle clustered data. Accommodates treatment effect modification based on observable covariates.
}
\usage{
## S3 method for class "hybrid"
hybrid(outcome, probit, modifiers = NULL, init = NULL, id = NULL, se = "R")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{outcome}{
an object of class "formula" with a numeric vector on the left hand side, and predictors of interest on the right hand side.
}
\item{probit}{
an object of class "formula" with a binary \code{(0/1)} numeric vector on the left hand side (1 indicating medication use), and predictors of medication use on the right hand side (right hand side permitted to contain variables on the right hand side of the outcome equation).
}
\item{modifiers}{
an object of class "formula" with a binary numeric vector indicating medication use on the left hand side, and treatment effect modifiers on the right hand side. If effect modifiers are treatment group specific (e.g., medication dose), set the effect modifier variables to zero for the untreated observations. If any other numeric values are used, they will ultimately be zet to zero. If \code{NULL}, the average treatment effect will be estimated under the assumption of no effect modification.
}
\item{init}{
a vector of initial values. The ordering of subparameters is: \code{alpha} (probit model parameters), \code{beta} (outcome model parameters), \code{eta} (an intercept, with or without effect mofidier paramters), \code{sigmay} (outcome error standard deviation), \code{rho} (error correlation). If \code{NULL}, an initial value will be chosen through OLS linear regression and probit-link GLM without regard to endogeneity.
}
\item{id}{
a numeric vector indicating subject IDs if data are clustered. In the absence of clustered data, this can be left blank (defaults to \code{NULL}).
}
\item{se}{
a string, either \code{"M"} for model-based standard errors (based on inverse observed Fisher information), or \code{"R"} for robust standard errors (based on methods of Huber and White). Defaults to \code{"R"}. If \code{id} is provided for clustered data, the cluster-robust variance estimator (with working independence) will be used even if the user specifies type \code{"M"}.
}
}
\details{
The model is evaluated with numerical minimization of the negative log-likelihood (the \code{BFGS} is used). The probit model and error correlation parameters are weakly identified and hence the error variance is set at unity. The data must be complete (no missing values) and numeric, with the exception of factors, which may be used on the right hand side of equations.
}
\value{
\code{hybrid} prints a summary of the coefficient estimates, standard errors, Wald-based confidence intervals, and p-values for the outcome model, the treatment effects (and potentially effect modifiers), and the medication use probit model.
\item{alpha}{estimate of the medication use probit model parameters.}
\item{beta}{estimate of the outcome model parameters.}
\item{eta}{estimate of the treatment effect, with or without effect modifier parameters.}
\item{sigma}{estimate of the standard deviation of the outcome error.}
\item{rho}{estimate of the correlation between the errors.}
\item{vcov}{entire estimated variance-covariance matrix, provided if the user wishes to perform any more specific hypothesis tests.}
\item{init}{initial value ultimately used, whether specified by the user or generated through the default approach.}
\item{fitted}{vector of fitted outcome values.}
\item{call}{the matched call.}
\item{out.form}{the formula used for the outcome model.}
\item{prob.form}{the formula used for the medication use probit model.}
\item{tx.form}{the formula used for the treatment effects model (potentially with effect modifiers).}
\item{sterr}{the choice of the variance estimate procedure (either model-based or robust).}
\item{labels}{labels for predictors to be passed into output.}
}
\author{
Andrew J. Spieker, Ph.D.
}
\references{
Heckman JJ. Dummy endogenous variables in a simultaneous equation system. \emph{Econometrica} 46(4), 931-959.

Maddala GS. \emph{Limited-dependent and qualitative variables in econometrics.} Cambridgeshire: Cambridge University Press; 1983.

Spieker AJ, Delaney JAC, and McClelland RL. Evaluating the treatment effects model for estimation of cross-sectional associations between risk factors and cardiovascular biomarkers influenced by medication use. \emph{Pharmacoepidemiology and Drug Safety} 24(12), 1286-1296.
}
\examples{
#- Generate Data -#
require(mvtnorm)
set.seed(1)
N <- 2000
X1 <- rnorm(N, 0, 1);
X2 <- rnorm(N, 0, 1);
X3 <- rnorm(N, 0, 1);
errors <- rmvnorm(N, sigma = 50*matrix(c(1, 0.5, 0.5, 1), nrow = 2))
Y <- 50 + X1 + X2 + errors[,1]
Z <- rep(0, N)
Z[(-5 + X1 + X3 + errors[,2]) > 0] <- 1
Y[Z == 1] <- Y[Z == 1] - 0.5*X1[Z == 1]

#- Estimate Model with No Effect Modification -#
hybrid(Y ~ X1 + X2, probit = Z ~ X1 + X3)

#- Estimate Model with Effect Modification -#
hybrid(Y ~ X1 + X2, probit = Z ~ X1 + X3, modifiers = Z ~ X1)

#- Estimate Model with Effect Modification and Model-Based Variance -#
hybrid(Y ~ X1 + X2, probit = Z ~ X1 + X3, modifiers = Z ~ X1, se = "M")
}