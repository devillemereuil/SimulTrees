require(dplyr)

# Tidier and glance for gls() models
tidy.gls <- function(x, 
                     conf.int = FALSE,
                     conf.level = 0.95,
                     ...) 
{
    summary(x)[["tTable"]] %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var="term") %>%
        dplyr::rename(estimate=Value,
                      std.error=Std.Error,
                      statistic=`t-value`,
                      p.value=`p-value`)
}

glance.gls <- function(x) 
{
    summary(x) %>%
        # Should it be tibble?
        {data.frame(
            # R-squared, deviance, statistics and p.value are not available for gls() AFAIK
            sigma   = .[["sigma"]],
            df      = .[["dims"]][["p"]],
            logLik  = .[["logLik"]],
            AIC     = .[["AIC"]],
            BIC     = .[["BIC"]],
            df.residual = .[["dims"]][["N"]] - .[["dims"]][["p"]]
        )}
}
