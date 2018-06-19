
# From mice 2.25
function (object, method = "smallsample") 
{
    call <- match.call()
    if (!is.mira(object)) 
    stop("The object must have class 'mira'")
    m <- length(object$analyses)
    fa <- getfit(object, 1)
    if (m == 1) {
        warning("Number of multiple imputations m=1. No pooling done.")
        return(fa)
    }
    analyses <- getfit(object)
    if (class(fa)[1] == "lme" & !requireNamespace("nlme", quietly = TRUE)) 
    stop("Package 'nlme' needed fo this function to work. Please install it.", 
        call. = FALSE)
    if ((class(fa)[1] == "mer" | class(fa)[1] == "lmerMod") & 
        !requireNamespace("lme4", quietly = TRUE)) 
    stop("Package 'lme4' needed fo this function to work. Please install it.", 
        call. = FALSE)
    mess <- try(coef(fa), silent = TRUE)
    if (inherits(mess, "try-error")) 
    stop("Object has no coef() method.")
    mess <- try(vcov(fa), silent = TRUE)
    if (inherits(mess, "try-error")) 
    stop("Object has no vcov() method.")
    if (class(fa)[1] == "mer" | class(fa)[1] == "lmerMod") {
        k <- length(lme4::fixef(fa))
        names <- names(lme4::fixef(fa))
    }
    else if (class(fa)[1] == "polr") {
        k <- length(coef(fa)) + length(fa$zeta)
        names <- c(names(coef(fa)), names(fa$zeta))
    }
    else {
        k <- length(coef(fa))
        names <- names(coef(fa))
    }
    qhat <- matrix(NA, nrow = m, ncol = k, dimnames = list(1:m, 
                                                        names))
    u <- array(NA, dim = c(m, k, k), dimnames = list(1:m, names, 
                                                    names))
    for (i in 1:m) {
        fit <- analyses[[i]]
        if (class(fit)[1] == "mer") {
            qhat[i, ] <- lme4::fixef(fit)
            ui <- as.matrix(vcov(fit))
            if (ncol(ui) != ncol(qhat)) 
            stop("Different number of parameters: class mer, fixef(fit): ", 
                ncol(qhat), ", as.matrix(vcov(fit)): ", ncol(ui))
            u[i, , ] <- array(ui, dim = c(1, dim(ui)))
        }
        else if (class(fit)[1] == "lmerMod") {
            qhat[i, ] <- lme4::fixef(fit)
            ui <- vcov(fit)
            if (ncol(ui) != ncol(qhat)) 
            stop("Different number of parameters: class lmerMod, fixed(fit): ", 
                ncol(qhat), ", vcov(fit): ", ncol(ui))
            u[i, , ] <- array(ui, dim = c(1, dim(ui)))
        }
        else if (class(fit)[1] == "lme") {
            qhat[i, ] <- fit$coefficients$fixed
            ui <- vcov(fit)
            if (ncol(ui) != ncol(qhat)) 
            stop("Different number of parameters: class lme, fit$coefficients$fixef: ", 
                ncol(qhat), ", vcov(fit): ", ncol(ui))
            u[i, , ] <- array(ui, dim = c(1, dim(ui)))
        }
        else if (class(fit)[1] == "polr") {
            qhat[i, ] <- c(coef(fit), fit$zeta)
            ui <- vcov(fit)
            if (ncol(ui) != ncol(qhat)) 
            stop("Different number of parameters: class polr, c(coef(fit, fit$zeta): ", 
                ncol(qhat), ", vcov(fit): ", ncol(ui))
            u[i, , ] <- array(ui, dim = c(1, dim(ui)))
        }
        else if (class(fit)[1] == "survreg") {
            qhat[i, ] <- coef(fit)
            ui <- vcov(fit)
            parnames <- dimnames(ui)[[1]]
            select <- !(parnames %in% "Log(scale)")
            ui <- ui[select, select]
            if (ncol(ui) != ncol(qhat)) 
            stop("Different number of parameters: class survreg, coef(fit): ", 
                ncol(qhat), ", vcov(fit): ", ncol(ui))
            u[i, , ] <- array(ui, dim = c(1, dim(ui)))
        }
        else {
            qhat[i, ] <- coef(fit)
            ui <- vcov(fit)
            ui <- expandvcov(qhat[i, ], ui)
            if (ncol(ui) != ncol(qhat)) 
            stop("Different number of parameters: coef(fit): ", 
                ncol(qhat), ", vcov(fit): ", ncol(ui))
            u[i, , ] <- array(ui, dim = c(1, dim(ui)))
        }
    }
    qbar <- apply(qhat, 2, mean)
    ubar <- apply(u, c(2, 3), mean)
    e <- qhat - matrix(qbar, nrow = m, ncol = k, byrow = TRUE)
    b <- (t(e) %*% e)/(m - 1)
    t <- ubar + (1 + 1/m) * b
    r <- (1 + 1/m) * diag(b/ubar)
    lambda <- (1 + 1/m) * diag(b/t)
    dfcom <- df.residual(object)
    df <- mice.df(m, lambda, dfcom, method)
    fmi <- (r + 2/(df + 3))/(r + 1)
    names(r) <- names(df) <- names(fmi) <- names(lambda) <- names
    fit <- list(call = call, call1 = object$call, call2 = object$call1, 
                nmis = object$nmis, m = m, qhat = qhat, u = u, qbar = qbar, 
                ubar = ubar, b = b, t = t, r = r, dfcom = dfcom, df = df, 
                fmi = fmi, lambda = lambda)
    oldClass(fit) <- c("mipo", oldClass(object))
    return(fit)
}
 
