limma.one.sided <- function (fit, lower = FALSE) 
{
    se.coef <- sqrt(fit$s2.post) * fit$stdev.unscaled
    df.total <- fit$df.prior + fit$df.residual
    pt(fit$t, df = df.total, lower.tail = lower)[, 1]
}

limma.one.sided.use.coef <- function (fit, coef, lower = FALSE) {
    se.coef <- sqrt(fit$s2.post) * fit$stdev.unscaled
    df.total <- fit$df.prior + fit$df.residual
    pt(fit$t, df = df.total, lower.tail = lower)[, coef]
}
permute <- function(l.t, g.t , t.t, time = Sys.time(), maxit=50000,  nPerm,threads, quietly = FALSE){
    
    set.seed(as.numeric(time))
    
    shuffTrans = sapply(1:nPerm, function(x){
        res = rep(Inf, length(t.t))
        for(j in unique(l.t)){
            res[l.t==j] = sample(t.t[l.t==j])
        }
        return(res)
    })
    
    l_index = 1
    c_index = 1
    t_index = 1:nPerm
    
    trios = cbind(l_index,c_index,t_index)
    cit = citpp::cit(L = l.t,G = g.t,T = shuffTrans,trios,threads=threads, quietly = quietly)
    
    return(cit)
    
}
permute_reactive <- function(l.t, g.t , t.t, time = Sys.time(), maxit=50000,  nPerm,threads, quietly = FALSE){
    
    set.seed(as.numeric(time))
    
    shuffTrans = sapply(1:nPerm, function(x){
        res = rep(Inf, length(t.t))
        for(j in unique(l.t)){
            res[l.t==j] = sample(t.t[l.t==j])
        }
        return(res)
    })
    
    l_index = 1
    c_index = 1
    t_index = 1:nPerm
    
    trios = cbind(l_index,c_index,t_index)
    
    cit = citpp::cit(L = l.t,T = g.t,G = shuffTrans,trios[,c(1,3,2)],threads=threads, quietly = quietly)
    
    return(cit)
    
}