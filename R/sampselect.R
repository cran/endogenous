sampselect <- function(outcome, probit, init = NULL, id = NULL, se = "R") {
	negloglik <- function(theta, W, X, Y, Z) {
		N <- dim(X)[1]
		dimw <- dim(W)[2]
		dimx <- dim(X)[2]
		alpha <- theta[1:dimw]
		beta <- theta[(dimw + 1):(dimw + dimx)]
		sigmay <- theta[dimw + dimx + 1]
		rho <- theta[dimw + dimx + 2]
		Ystarhat <- X %*% beta
		Ys0 <- Ystarhat[Z == 0]
		Y0 <- Y[Z == 0]
		Zstarhat <- W %*% alpha
		Z0 <- Zstarhat[Z == 0]
		Z1 <- Zstarhat[Z == 1]
		like01 <- dnorm(Y0, mean = Ys0, sd = sigmay)
		mu.cond0 <- Z0 + (rho/sigmay)*(Y0 - Ys0)
		mu.cond1 <- Z1
		sigma2.cond <- (1 - rho^2)
		utx <- pnorm(0, mean = mu.cond0, sd = sqrt(sigma2.cond))
		tx <- 1 - pnorm(0, mean = mu.cond1, sd = 1)
		like <- rep(0, N)
		like[Z == 0] <- log(like01) + log(utx + 1e-320)
		like[Z == 1] <- log(tx + 1e-320)
		negll = -sum(like)
		negll
	}
	dthetai.part <- function(w, x, y, z, alpha, beta, eta, rho, sigmay) {
		N <- dim(x)[1]
		da.t <- matrix(0, nrow = N, ncol = dim(w)[2])
		db.t <- matrix(0, nrow = N, ncol = dim(x)[2]) 
		ds.t <- matrix(0, nrow = N, 1)
		dr.t <- matrix(0, nrow = N, 1)
		w1 <- w[z == 1,]
		w <- w[z == 0,]
		x <- x[z == 0,]
		y <- y[z == 0]
		a <- (y - x %*% beta)/sigmay
		b1 <- w %*% alpha + (rho/sigmay)*(y - x %*% beta)
		b2 <- sqrt(1 - rho^2)
		b <- -b1/b2
		mr <- dnorm(b)/pnorm(b)
		d <- -mr/b2
		dlda <- matrix(rep(d, dim(w)[2]), ncol = dim(w)[2]) * w
		db.1 <- matrix(rep((a/sigmay), dim(x)[2]), ncol = dim(x)[2]) * x
		e <- mr * rho/(sigmay * b2)
		db.2 <- matrix(rep(e, dim(x)[2]), ncol = dim(x)[2]) * x
		dldb <- db.1 + db.2
		bot <- sqrt(1 - rho^2)
		dbot <- rho/b2
		botsq <- b2^2
		dldr <- -mr * (bot * a + b1 * dbot)/botsq
		ds.1 <- -1/sigmay
		ds.2 <- (a^2)/sigmay
		ds.3 <- -mr * (rho/b2) * (a/sigmay)
		dlds <- ds.1 + ds.2 - ds.3
		dlda1 <- matrix(rep((dnorm(w1 %*% alpha)/pnorm(w1 %*% alpha)), dim(w1)[2]), ncol = dim(w1)[2]) * w1
		da.t[z == 0,] <- dlda
		da.t[z == 1,] <- dlda1
		db.t[z == 0,] <- dldb
		ds.t[z == 0] <- dlds
		dr.t[z == 0] <- dldr
		dldtheta <- matrix(cbind(da.t, db.t, ds.t, dr.t), nrow = N)
		dldtheta
	}
	sampselect.do <- function(outcome, probit, init = NULL, id = NULL, se = "R") {
		if (length(se) != 1) {stop("Argument 'se' must be either 'R' or 'M'")}
		if (se != "R" & se != "M") {stop("Argument 'se' must be either 'R' or 'M'")}
		N <- dim(model.matrix(outcome))[1]
		if (!is.null(id) & length(id) != N) {stop("Argument 'id' must be a numeric vector of length N")}
		N <- dim(model.matrix(outcome))[1]
    	mm.out <- model.frame(outcome)
    	mm.prb <- model.frame(probit)
		Y <- mm.out[,1]
    	X <- matrix(model.matrix(outcome), nrow = N)
    	Z <- mm.prb[,1]
    	W <- matrix(model.matrix(probit), nrow = N)
		dimw <- dim(W)[2]
		dimx <- dim(X)[2]
    	if (is.null(init)) {
    	a.guess <- coef(glm(Z ~ W[,2:dim(W)[2]], family = binomial(link = "probit")))
    	b.model <- lm(Y ~ X[,2:dim(X)[2]], subset = (Z == 0))
    	b.guess <- coef(b.model)
		sigmay.guess <- sd(residuals(b.model))
		rho.guess <- 0.3
		init <- as.numeric(c(a.guess, b.guess, sigmay.guess, rho.guess))
		}
		zz <- suppressWarnings(optim(init, negloglik, method = c("BFGS"), W = W, X = X, Y = Y, Z = Z, hessian = T))
		theta.hat <- zz$par
		a.hat <- theta.hat[1:dimw]
		b.hat <- theta.hat[(dimw + 1):(dimw + dimx)]
		s.hat <- theta.hat[dimw + dimx + 1]
		r.hat <- theta.hat[dimw + dimx + 2]
		zz.hess <- zz$hessian
		if (!is.null(id)) {
			big.meat <- dthetai.part(w = W, x = X, y = Y, z = Z,
    			alpha = a.hat, beta = b.hat, rho = r.hat, sigmay = s.hat)
			uni.id <- unique(id)
			meat.mat <- matrix(0, nrow = (dimx + dimw + 2), ncol = (dimx + dimw + 2))
			for (idx in 1:length(uni.id)) {
				cY <- Y[id == uni.id[idx]]
				idxs <- which(id == uni.id[idx])
				sub.meat <- colSums(matrix(big.meat[idxs,], ncol = (dimx + dimw + 2)))
				meat.mat <- meat.mat + (sub.meat) %*% t(sub.meat)
			}
			varcov <- solve(zz.hess) %*% meat.mat %*% t(solve(zz.hess))
		}
		if (is.null(id) & se == "R") {
    		meat.mat <- dthetai.part(w = W, x = X, y = Y, z = Z,
    			alpha = a.hat, beta = b.hat, rho = r.hat, sigmay = s.hat)
    		meat.mat <- crossprod(meat.mat)
    		varcov <- solve(zz.hess) %*% meat.mat %*% t(solve(zz.hess))
    	}
    	if (is.null(id) & se == "M") {
    		varcov <- solve(zz.hess)
    	}
    	invisible(list(alpha = as.numeric(a.hat), beta = as.numeric(b.hat),
    		sigma = as.numeric(s.hat), rho = as.numeric(r.hat), vcov = varcov, init  = init, fitted = X %*% b.hat))
    	}
    sampselect.obj <- sampselect.do(outcome = outcome, probit = probit, init = init, id = id, se = se)
    sampselect.obj$call <- match.call()
    sampselect.obj$out.form <- outcome
    sampselect.obj$prob.form <- probit
    sampselect.obj$sterr <- se
    sampselect.obj $labels <- c(colnames(model.matrix(outcome)), colnames(model.matrix(probit)))
    class(sampselect.obj) <- "sampselect"
    sampselect.obj
}
print.sampselect <- function(x,...) {
	cat("\nCall:\n")
    print(x$call)
    outcome <- x$out.form
    prb <- x$prob.form
    alpha <- x$alpha
    beta <- x$beta
    sigma <- x$sigma
    rho <- x$rho
    vcov <- x$vcov
    vartype <- x$sterr
    labs <- x$labels
    if (vartype == "R") {varlab <- "Robust St. Err."}
    if (vartype == "M") {varlab <- "St. Err"}
    dimw <- length(alpha)
    dimx <- length(beta)
    st.errs <- sqrt(diag(vcov))
    se.a <- st.errs[1:dimw]
	se.b <- st.errs[(dimw + 1):(dimw + dimx)]
    p.b <- 2*round(1 - pnorm(abs(beta/se.b)), 3)
    p.a <- 2*round(1 - pnorm(abs(alpha/se.a)), 3)
    p.b[p.b < 0.0001] <- "< 0.0001"
    p.a[p.a < 0.0001] <- "< 0.0001"
    main.outcome <- data.frame(beta, se.b, beta - 1.96*se.b, beta + 1.96*se.b, p.b)
    names(main.outcome) <- c("Estimate", varlab, "95% L", "95% U", "P(> |Z|)")
    row.names(main.outcome) <- c(labs[1:dimx])
    cat("\nOutcome Model:\n")
    print(main.outcome)
    main.prob <- data.frame(alpha, se.a, alpha - 1.96*se.a, alpha + 1.96*se.a, p.a)
    names(main.prob) <- c("Estimate", varlab, "95% L", "95% U", "P(> |Z|)")
    row.names(main.prob) <- c(labs[(dimx + 1):(dimx + dimw)])
    cat("\nMedication Use Probit Model:\n")
    print(main.prob)
    rho <- round(rho, 2)
    cat("\nRho = ", rho)
    cat("\nOutcome Error Variance =", sigma^2,"\n")
}
summary.sampselect <- function(object,...) {
	cat("\nCall:\n")
    print(object$call)
    outcome <- object$out.form
    prb <- object$prob.form
    alpha <- object$alpha
    beta <- object$beta
    sigma <- object$sigma
    rho <- object$rho
    vcov <- object$vcov
    vartype <- object$sterr
    labs <- object$labels
    if (vartype == "R") {varlab <- "Robust St. Err."}
    if (vartype == "M") {varlab <- "St. Err"}
    dimw <- length(alpha)
    dimx <- length(beta)
    st.errs <- sqrt(diag(vcov))
    se.a <- st.errs[1:dimw]
	se.b <- st.errs[(dimw + 1):(dimw + dimx)]
    p.b <- 2*round(1 - pnorm(abs(beta/se.b)), 3)
    p.a <- 2*round(1 - pnorm(abs(alpha/se.a)), 3)
    p.b[p.b < 0.0001] <- "< 0.0001"
    p.a[p.a < 0.0001] <- "< 0.0001"
    main.outcome <- data.frame(beta, se.b, beta - 1.96*se.b, beta + 1.96*se.b, p.b)
    names(main.outcome) <- c("Estimate", varlab, "95% L", "95% U", "P(> |Z|)")
    row.names(main.outcome) <- c(labs[1:dimx])
    cat("\nOutcome Model:\n")
    print(main.outcome)
    main.prob <- data.frame(alpha, se.a, alpha - 1.96*se.a, alpha + 1.96*se.a, p.a)
    names(main.prob) <- c("Estimate", varlab, "95% L", "95% U", "P(> |Z|)")
    row.names(main.prob) <- c(labs[(dimx + 1):(dimx + dimw)])
    cat("\nMedication Use Probit Model:\n")
    print(main.prob)
    rho <- round(rho, 2)
    cat("\nRho = ", rho)
    cat("\nOutcome Error Variance =", sigma^2,"\n")
}