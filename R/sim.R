
#' Generate Poisson Point Pattern
#'
#' \code{sim_pois} generates a random point pattern using the Poisson process.
#'
#' @param lambda_int An integer specifying the intensity of the Poisson process; the expected
#' number of points *per unit area*.  The total number of points in
#' the simulated pattern will be random with expectation value ‘mu =
#' lambda * a’ where ‘a’ is the area of the window in which the pattern is simulated.
#' @param win_len A numeric specifying the window side-length of a quadratic window in which the pattern is simulated.
#' 
#' @return A point-pattern 
#'
#' @examples
#' pp = sim_pois(100)
#' 
#' @export
sim_pois <- function(lambda_int, win_len = 1){
    
    pp = spatstat::rpoispp(lambda_int, win = spatstat::owin(c(0, win_len), c(0, win_len)))

    return(pp)
}

#' Add step-gradient mark distributions to a point pattern
#'
#' \code{add_markdist_step} adds step-gradient mark distributions to a point pattern.
#'
#' @param pp A point-pattern.
#' @param low_marks A numeric or numeric vector specifying the lower value of the step mark distribution. If a vector is given 'n' mark distributions are added where 'n' is the length of the vector.
#' @param high_marks A numeric or numeric vector specifying the upper value of the step mark distribution. If a vector is given 'n' mark distributions are added where 'n' is the length of the vector.
#' @param step_border A numeric specifying the relative x-position within the point-pattern window where the values of the marks will change from low to high.
#' 
#' @return A point-pattern with added mark distributions.
#'
#' @examples
#' low_expr = c(10, 10)
#' high_expr = c(15, 20)
#' pp = sim_pois(100)
#' pp = add_markdist_step(pp, low_expr, high_expr)
#' 
#' @export
add_markdist_step <- function(pp, low_marks, high_marks, step_border = 0.5){
###add step mark distribution
    
    ##mark dist init
    npoints = pp[['n']]
    nmarks = length(low_marks)
    marx = as.data.frame(matrix(NA, nrow = npoints, ncol = nmarks))
    rownames(marx) = paste('p', 1:npoints, sep = '')
    colnames(marx) = paste('step', 1:nmarks, sep = '_')

    ##step mark dist
    x = pp[['x']]
    x_max = pp[['window']][['xrange']][2]
    x_border = x_max * step_border
    low_ind = which(x < x_border)
    high_ind = setdiff(1:npoints, low_ind)
    
    for(j_gene in 1:nmarks){
        marx[low_ind, j_gene] = low_marks[j_gene]
        marx[high_ind, j_gene] = high_marks[j_gene]
    }

    ##add to previous mark dist if exists
    marx_prev = pp[['marks']]
    if(!is.null(marx_prev)){
        marx = cbind(marx_prev, marx)
    }

    ##set mark dist
    pp[['marks']] = marx

    return(pp)
}

#' Add hot-spot mark distributions to a point pattern
#'
#' \code{add_markdist_hotspot} adds squared hot-spot mark distributions to a point pattern.
#'
#' @param pp A point-pattern.
#' @param low_marks A numeric or numeric vector specifying the lower value of the mark distribution. If a vector is given 'n' mark distributions are added where 'n' is the length of the vector.
#' @param high_marks A numeric or numeric vector specifying the upper value of the mark distribution. If a vector is given 'n' mark distributions are added where 'n' is the length of the vector.
#' @param hotspot_size A numeric specifying the relative size of the point-pattern window where the values of the marks will change from low to high.
#' 
#' @return A point-pattern with added mark distributions.
#'
#' @examples
#' low_expr = c(10, 10)
#' high_expr = c(15, 20)
#' pp = sim_pois(100)
#' pp = add_markdist_hotspot(pp, low_expr, high_expr)
#' 
#' @export
add_markdist_hotspot <- function(pp, low_marks, high_marks, hotspot_size = 0.2){
###add step mark distribution

    ##mark dist init
    npoints = pp[['n']]
    nmarks = length(low_marks)
    marx = as.data.frame(matrix(NA, nrow = npoints, ncol = nmarks))
    rownames(marx) = paste('p', 1:npoints, sep = '')
    colnames(marx) = paste('hotspot', 1:nmarks, sep = '_')

    ##hotspot mark dist
    x = pp[['x']]
    y = pp[['y']]
    x_max = pp[['window']][['xrange']][2]
    half_hot_len = x_max * hotspot_size / 2
    x_half = x_max * 0.5
    hot_min = x_half - half_hot_len
    hot_max = x_half + half_hot_len
    
    high_ind = which(x >= hot_min & x <= hot_max & y >= hot_min & y <= hot_max)
    low_ind = setdiff(1:npoints, high_ind)
    for(j_gene in 1:nmarks){
        marx[low_ind, j_gene] = low_marks[j_gene]
        marx[high_ind, j_gene] = high_marks[j_gene]
    }
    
    ##add to previous mark dist if exists
    marx_prev = pp[['marks']]
    if(!is.null(marx_prev)){
        marx = cbind(marx_prev, marx)
    }

    ##set mark dist
    pp[['marks']] = marx

    return(pp)
}

#' Add streak mark distributions to a point pattern
#'
#' \code{add_markdist_streak} adds rectangular streak mark distributions to a point pattern.
#'
#' @param pp A point-pattern.
#' @param low_marks A numeric or numeric vector specifying the lower value of the mark distribution. If a vector is given 'n' mark distributions are added where 'n' is the length of the vector.
#' @param high_marks A numeric or numeric vector specifying the upper value of the mark distribution. If a vector is given 'n' mark distributions are added where 'n' is the length of the vector.
#' @param streak_len_frac A numeric specifying the relative length of the point-pattern window where the values of the marks will change from low to high.
#' @param streak_height_frac A numeric specifying the relative height of the point-pattern window where the values of the marks will change from low to high.
#' 
#' @return A point-pattern with added mark distributions.
#'
#' @examples
#' low_expr = c(10, 10)
#' high_expr = c(15, 20)
#' pp = sim_pois(100)
#' pp = add_markdist_streak(pp, low_expr, high_expr)
#' 
#' @export
add_markdist_streak <- function(pp, low_marks, high_marks, streak_len_frac = 0.9, streak_height_frac = 0.05){
###add step mark distribution

    ##mark dist init
    npoints = pp[['n']]
    nmarks = length(low_marks)
    marx = as.data.frame(matrix(NA, nrow = npoints, ncol = nmarks))
    rownames(marx) = paste('p', 1:npoints, sep = '')
    colnames(marx) = paste('streak', 1:nmarks, sep = '_')

    ##streak mark dist
    x = pp[['x']]    
    x_max = pp[['window']][['xrange']][2]
    half_streak_len = x_max * streak_len_frac / 2
    x_half = x_max * 0.5
    streak_x_min = x_half - half_streak_len
    streak_x_max = x_half + half_streak_len

    y = pp[['y']]
    y_max = pp[['window']][['yrange']][2]
    y_half = y_max * 0.5
    half_streak_height = y_max * streak_height_frac / 2
    streak_y_min = y_half - half_streak_height
    streak_y_max = y_half + half_streak_height
        
    high_ind = which(x >= streak_x_min & x <= streak_x_max & y >= streak_y_min & y <= streak_y_max)
    low_ind = setdiff(1:npoints, high_ind)
    for(j_gene in 1:nmarks){
        marx[low_ind, j_gene] = low_marks[j_gene]
        marx[high_ind, j_gene] = high_marks[j_gene]
    }
    
    ##add to previous mark dist if exists
    marx_prev = pp[['marks']]
    if(!is.null(marx_prev)){
        marx = cbind(marx_prev, marx)
    }

    ##set mark dist
    pp[['marks']] = marx

    return(pp)
}

##set all marks as bg
##identify cells in selection window
##note that this will reflect sampling of cells as the number of positive cells in the window will be stochastic
add_markdist_windowed <- function(pp, bg_marks, spike_marks, x_frac = 0.1, y_frac = x_frac, x_start_frac = 0.5, y_start_frac = 0.5){
    
    ##x
    x = pp[['x']]    
    x_max = pp[['window']][['xrange']][2]
    w_x_middle = x_max * x_start_frac

    ##window borders
    w_half_width = x_max * x_frac / 2
    w_x_min = w_x_middle - w_half_width
    w_x_max = w_x_middle + w_half_width

    ##y
    y = pp[['y']]
    y_max = pp[['window']][['yrange']][2]
    w_y_middle = y_max * y_start_frac

    ##window borders
    w_half_height = y_max * y_frac / 2
    w_y_min = w_y_middle - w_half_height
    w_y_max = w_y_middle + w_half_height

    ##get indices of points in window
    high_ind = which(x >= w_x_min & x <= w_x_max & y >= w_y_min & y <= w_y_max)
    n_spiked = length(high_ind)
    
    print(sprintf('Number of cells in spike-window: %i', n_spiked))
    
    ##set spike marks to points in window
    marx = bg_marks
    marx[high_ind, ] = spike_marks[high_ind, ]
    rownames(marx) = paste('p', 1:nrow(marx), sep = '')
    colnames(marx) = paste('g', 1:ncol(marx), sep = '')
    pp[['marks']] = marx

    return(list(pp = pp, n_spiked = n_spiked, spiked_ind = high_ind))
}
