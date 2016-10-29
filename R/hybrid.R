hybrid <- function(outcome, probit, modifiers = NULL, init = NULL, id = NULL, se = "R") {
	negloglik <- function(theta, V, W, X, Y, Z) {
		dimv <- dim(V)[2]
		dimw <- dim(W)[2]
		dimx <- dim(X)[2]
		alpha <- theta[1:dimw]
		beta <- theta[(dimw + 1):(dimw + dimx)]
		eta <- theta[(dimw + dimx + 1):(dimw + dimx + dimv)]
		sigmay <- theta[dimw + dimx + dimv + 1]
		rho <- theta[dimw + dimx + dimv + 2]
		Ystarhat <- X %*% beta
		Zstarhat <- W %*% alpha
		if(dimv == 1) {nuisance <- V * eta * Z}
		if (dimv > 1) {nuisance <- V %*% eta * Z}
		all <- dnorm(Y, mean = Ystarhat - nuisance, sd = sigmay)
		mu.cond <- W %*% alpha + (rho/sigmay)*(Y - Ystarhat + nuisance)
		sigma2.cond <- (1 - rho^2)
		tx <- 1 - pnorm(0, mean = mu.cond, sd = sqrt(sigma2.cond))
		utx <- pnorm(0, mean = mu.cond, sd = sqrt(sigma2.cond))
		like <- log(all) + log(utx + 1e-320)
		like[Z==1] <- log(all)[Z==1] + log(tx + 1e-320)[Z==1]
		negll = -sum(like)
		negll
	}
	dthetai.part <- function(v, w, x, y, z, alpha, beta, eta, rho, sigmay) {
		a <- (y - x %*% beta + (v %*% eta)*z)/sigmay
		b1 <- w %*% alpha + (rho/sigmay)*(y - x %*% beta + (v %*% eta)*z)
		b2 <- sqrt(1 - rho^2)
		b <- ((-1)^(1 - z)) * (b1/b2)
		mr <- dnorm(b)/pnorm(b)
		d <- mr * ((-1)^(1 - z))/(b2)
		dlda <- matrix(rep(d, dim(w)[2]), ncol = dim(w)[2]) * w
		db.1 <- matrix(rep((a/sigmay), dim(x)[2]), ncol = dim(x)[2]) * x
		e <- mr * rho * ((-1)^(1 - z))/(sigmay * b2)
		db.2 <- -1 * matrix(rep(e, dim(x)[2]), ncol = dim(x)[2]) * x
		dldb <- db.1 + db.2
		de.1 <- -1 * matrix(rep((a/sigmay), dim(v)[2]), ncol = dim(v)[2]) * v * z
		e <- z * mr * rho * ((-1)^(1 - z))/(sigmay * b2)
		de.2 <- matrix(rep(e, dim(v)[2]), ncol = dim(v)[2]) * v
		dlde <- de.1 + de.2
		bot <- sqrt(1 - rho^2)
		dbot <- rho/b2
		botsq <- b2^2
		dldr <- mr * ((-1)^(1 - z)) * (bot * a + b1 * dbot)/botsq
		ds.1 <- -1/sigmay
		ds.2 <- (a^2)/sigmay
		ds.3 <- mr * ((-1)^(1 - z)) * (rho/b2) * (a/sigmay)
		dlds <- ds.1 + ds.2 - ds.3
		dldtheta <- matrix(cbind(dlda, dldb, dlde, dlds, dldr), nrow = dim(x)[1])
		dldtheta
	}
	hybrid.do <- function(outcome, probit, modifiers = NULL, init = NULL, id = NULL, se = "R") {
		if (length(se) != 1) {stop("Argument 'se' must be either 'R' or 'M'")}
		if (se != "R" & se != "M") {stop("Argument 'se' must be either 'R' or 'M'")}
		N <- dim(model.matrix(outcome))[1]
		if (!is.null(id) & length(id) != N) {stop("Argument 'id' must be a numeric vector of length N")}
    	mm.out <- model.frame(outcome)
    	mm.prb <- model.frame(probit)
    		if (is.null(modifiers)) {V <- matrix(rep(1, N), nrow = N)}
    		if (!is.null(modifiers)) {V <- matrix(model.matrix(modifiers), nrow = N)}
		Y <- mm.out[,1]
    	X <- matrix(model.matrix(outcome), nrow = N)
    	Z <- mm.prb[,1]
    	W <- matrix(model.matrix(probit), nrow = N)
    	dimv <- dim(V)[2]
		dimw <- dim(W)[2]
		dimx <- dim(X)[2]
    	if (is.null(init)) {
    	a.guess <- coef(glm(Z ~ W[,2:dim(W)[2]], family = binomial(link = "probit")))
    	be.model <- lm(Y ~ X[,2:dim(X)[2]] + I(V*Z))
    	be.guess <- coef(be.model)
		sigmay.guess <- sd(residuals(be.model))
		rho.guess <- cor(Z, Y)
		init <- as.numeric(c(a.guess, be.guess, sigmay.guess, rho.guess))
		}
		zz <- suppressWarnings(optim(init, negloglik, method = c("BFGS"), V = V, W = W, X = X, Y = Y, Z = Z, hessian = T))
		theta.hat <- zz$par
		a.hat <- theta.hat[1:dimw]
		b.hat <- theta.hat[(dimw + 1):(dimw + dimx)]
		e.hat <- theta.hat[(dimw + dimx + 1):(dimw + dimx + dimv)]
		s.hat <- theta.hat[dimw + dimx + dimv + 1]
		r.hat <- theta.hat[dimw + dimx + dimv + 2]
		zz.hess <- zz$hessian
		if (!is.null(id)) {
			big.meat <- dthetai.part(v = V, w = W, x = X, y = Y, z = Z,
    			alpha = a.hat, beta = b.hat, eta = e.hat, rho = r.hat, sigmay = s.hat)
			uni.id <- unique(id)
			meat.mat <- matrix(0, nrow = (dimx + dimw + dimv + 2), ncol = (dimx + dimw + dimv + 2))
			for (idx in 1:length(uni.id)) {
				cY <- Y[id == uni.id[idx]]
				idxs <- which(id == uni.id[idx])
				sub.meat <- colSums(matrix(big.meat[idxs,], ncol = (dimx + dimw + dimv + 2)))
				meat.mat <- meat.mat + (sub.meat) %*% t(sub.meat)
			}
			varcov <- solve(zz.hess) %*% meat.mat %*% t(solve(zz.hess))
		}
		if (is.null(id) & se == "R") {
    		meat.mat <- dthetai.part(v = V, w = W, x = X, y = Y, z = Z,
    			alpha = a.hat, beta = b.hat, eta = e.hat, rho = r.hat, sigmay = s.hat)
    		meat.mat <- crossprod(meat.mat)
    		varcov <- solve(zz.hess) %*% meat.mat %*% t(solve(zz.hess))
    	}
    	if (is.null(id) & se == "M") {
    		varcov <- solve(zz.hess)
    	}
    	invisible(list(alpha = as.numeric(a.hat), beta = as.numeric(b.hat), eta = as.numeric(e.hat),
    		sigma = as.numeric(s.hat), rho = as.numeric(r.hat), vcov = varcov, init  = init, fitted = X %*% b.hat))
    	}
    hybrid.obj <- hybrid.do(outcome = outcome, probit = probit, modifiers = modifiers, init = init, id = id, se = se)
    hybrid.obj$call <- match.call()
    hybrid.obj$out.form <- outcome
    hybrid.obj$prob.form <- probit
    hybrid.obj$tx.form <- modifiers
    hybrid.obj$sterr <- se
    if (is.null(modifiers)) {
    	hybrid.obj$labels <- c(colnames(model.matrix(outcome)), "NA", colnames(model.matrix(probit)))
    }
    if (!is.null(modifiers)) {
    	hybrid.obj$labels <- c(colnames(model.matrix(outcome)), colnames(model.matrix(modifiers)), colnames(model.matrix(probit)))
    }
    class(hybrid.obj) <- "hybrid"
    hybrid.obj
}
print.hybrid <- function(x,...) {
	cat("\nCall:\n")
    print(x$call)
    outcome <- x$out.form
    prb <- x$prob.form
    mod <- x$tx.form
    alpha <- x$alpha
    beta <- x$beta
    eta <- x$eta
    sigma <- x$sigma
    rho <- x$rho
    vcov <- x$vcov
    vartype <- x$sterr
    labs <- x$labels
    if (vartype == "R") {varlab <- "Robust St. Err."}
    if (vartype == "M") {varlab <- "St. Err"}
    dimv <- length(eta)
    dimw <- length(alpha)
    dimx <- length(beta)
    st.errs <- sqrt(diag(vcov))
    se.a <- st.errs[1:dimw]
	se.b <- st.errs[(dimw + 1):(dimw + dimx)]
    se.e <- st.errs[(dimw + dimx + 1):(dimx + dimw + dimv)]
    p.b <- 2*round(1 - pnorm(abs(beta/se.b)), 3)
    p.a <- 2*round(1 - pnorm(abs(alpha/se.a)), 3)
    p.e <- 2*round(1 - pnorm(abs(eta/se.e)), 3)
    p.b[p.b < 0.0001] <- "< 0.0001"
    p.a[p.a < 0.0001] <- "< 0.0001"
    p.e[p.e < 0.0001] <- "< 0.0001"
    main.outcome <- data.frame(beta, se.b, beta - 1.96*se.b, beta + 1.96*se.b, p.b)
    names(main.outcome) <- c("Estimate", varlab, "95% L", "95% U", "P(> |Z|)")
    row.names(main.outcome) <- c(labs[1:dimx])
    cat("\nOutcome Model:\n")
    print(main.outcome)
	main.tx <- data.frame(eta, se.e, eta - 1.96*se.e, eta + 1.96*se.e, p.e)
    names(main.tx) <- c("Estimate", varlab, "95% L", "95% U", "P(> |Z|)")
    if (dimv == 1) {
    	lab <- "Average Effect"
    	row.names(main.tx) <- c(lab)
    	cat("\nTreatment Effect:\n")
    	print(main.tx)
    }
    if (dimv > 1) {
    	lab <- "Main Effect"
    	row.names(main.tx) <- c(labs[(dimx + 1):(dimx + dimv)])
    	cat("\nCovariate-Specific Treatment Effects:\n")
    	print(main.tx)
    }
    main.prob <- data.frame(alpha, se.a, alpha - 1.96*se.a, alpha + 1.96*se.a, p.a)
    names(main.prob) <- c("Estimate", varlab, "95% L", "95% U", "P(> |Z|)")
    row.names(main.prob) <- c(labs[(dimx + dimv + 1):(dimx + dimv + dimw)])
    cat("\nMedication Use Probit Model:\n")
    print(main.prob)
    rho <- round(rho, 2)
    cat("\nRho = ", rho)
    cat("\nOutcome Error Variance =", sigma^2,"\n")
}
summary.hybrid <- function(object,...) {
	cat("\nCall:\n")
    print(object$call)
    outcome <- object$out.form
    prb <- object$prob.form
    mod <- object$tx.form
    alpha <- object$alpha
    beta <- object$beta
    eta <- object$eta
    sigma <- object$sigma
    rho <- object$rho
    vcov <- object$vcov
    vartype <- object$sterr
    labs <- object$labels
    if (vartype == "R") {varlab <- "Robust St. Err."}
    if (vartype == "M") {varlab <- "St. Err"}
    dimv <- length(eta)
    dimw <- length(alpha)
    dimx <- length(beta)
    st.errs <- sqrt(diag(vcov))
    se.a <- st.errs[1:dimw]
	se.b <- st.errs[(dimw + 1):(dimw + dimx)]
    se.e <- st.errs[(dimw + dimx + 1):(dimx + dimw + dimv)]
    p.b <- 2*round(1 - pnorm(abs(beta/se.b)), 3)
    p.a <- 2*round(1 - pnorm(abs(alpha/se.a)), 3)
    p.e <- 2*round(1 - pnorm(abs(eta/se.e)), 3)
    p.b[p.b < 0.0001] <- "< 0.0001"
    p.a[p.a < 0.0001] <- "< 0.0001"
    p.e[p.e < 0.0001] <- "< 0.0001"
    main.outcome <- data.frame(beta, se.b, beta - 1.96*se.b, beta + 1.96*se.b, p.b)
    names(main.outcome) <- c("Estimate", varlab, "95% L", "95% U", "P(> |Z|)")
    row.names(main.outcome) <- c(labs[1:dimx])
    cat("\nOutcome Model:\n")
    print(main.outcome)
	main.tx <- data.frame(eta, se.e, eta - 1.96*se.e, eta + 1.96*se.e, p.e)
    names(main.tx) <- c("Estimate", varlab, "95% L", "95% U", "P(> |Z|)")
    if (dimv == 1) {
    	lab <- "Average Effect"
    	row.names(main.tx) <- c(lab)
    	cat("\nTreatment Effect:\n")
    	print(main.tx)
    }
    if (dimv > 1) {
    	lab <- "Main Effect"
    	row.names(main.tx) <- c(labs[(dimx + 1):(dimx + dimv)])
    	cat("\nCovariate-Specific Treatment Effects:\n")
    	print(main.tx)
    }
    main.prob <- data.frame(alpha, se.a, alpha - 1.96*se.a, alpha + 1.96*se.a, p.a)
    names(main.prob) <- c("Estimate", varlab, "95% L", "95% U", "P(> |Z|)")
    row.names(main.prob) <- c(labs[(dimx + dimv + 1):(dimx + dimv + dimw)])
    cat("\nMedication Use Probit Model:\n")
    print(main.prob)
    rho <- round(rho, 2)
    cat("\nRho = ", rho)
    cat("\nOutcome Error Variance =", sigma^2,"\n")
}