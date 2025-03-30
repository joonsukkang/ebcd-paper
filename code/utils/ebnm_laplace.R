ebnm_laplace <- function(x,
                         s = 1,
                         mode = 0,
                         scale = "estimate",
                         g_init = NULL,
                         fix_g = FALSE,
                         output = ebnm::ebnm_output_default(),
                         optmethod = NULL,
                         control = NULL) {
  
  # Create a custom initialization function that forces alpha to be very large
  custom_initpar <- function(g_init, mode, scale, pointmass, x, s) {
    par <- ebnm:::pl_initpar(g_init, mode, scale, pointmass, x, s)
    # Force alpha to be very large to make point mass weight effectively zero
    par$alpha <- 30
    
    return(par)
  }
  
  # Create a custom parameter checking function that keeps alpha fixed
  custom_checkg <- function(g_init, fix_g, mode, scale, pointmass, call) {
    result <- ebnm:::pl_checkg(g_init, fix_g, mode, scale, pointmass, call)
    
    # Always fix alpha (first parameter) to keep point mass at zero
    # But allow beta (scale parameter) to be optimized
    result$fix_par[1] <- TRUE
    
    return(result)
  }
  
  # Set up call
  call <- match.call()
  
  # Ensure control is a list
  control <- if (is.null(control)) ebnm:::nlm_control_defaults() else control
  
  # Use the point_laplace machinery, but with our custom functions
  retlist <- ebnm:::parametric_workhorse(
    x = x,
    s = s,
    mode = mode,
    scale = scale,
    pointmass = TRUE,
    g_init = g_init,
    fix_g = fix_g,
    output = output,
    optmethod = optmethod,
    control = control,
    checkg_fn = custom_checkg,
    initpar_fn = custom_initpar,
    scalepar_fn = ebnm:::pl_scalepar,
    precomp_fn = ebnm:::pl_precomp,
    nllik_fn = ebnm:::pl_nllik,
    postcomp_fn = ebnm:::pl_postcomp,
    summres_fn = ebnm:::pl_summres,
    partog_fn = ebnm:::pl_partog,
    postsamp_fn = ebnm:::pl_postsamp,
    call = call
  )
  
  return(ebnm:::as_ebnm(retlist, call))
}