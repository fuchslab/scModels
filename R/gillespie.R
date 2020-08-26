#' Gillespie algorithm for mRNA generating processes
#'
#' Gillespie algorithms allow synthetic data simulation via three different
#' underlying mRNA generating processes: the basic process consists of a
#' simple death-birth model of mRNA transcription and degradation; the
#' switching process considers additionally gene activation and deactivation,
#' with mRNA transcription only happening in active gene states; the
#' bursting process, transcribes mRNA in bursts with geometrically distributed burst sizes.
#' The basic_burst model combines both the basic and the burst model.
#'

#' @param n Number of observations
#' @param r.degr mRNA degradation rate (all models)
#' @param r.act DNA activation rate (Switching Model)
#' @param r.deact DNA deactivation rate (Switching Model)
#' @param r.on Transcription rate during gene activation (Switching model)
#' @param r.burst Bursty transcription rate (Bursting model)
#' @param s.burst Mean burst size (Bursting Model)
#' @name gmRNA
#' @rdname gmRNA
#' @export
#' @examples
#' x <- gmRNA_basic(100, 0.75, 0.001)
#' plot(density(x))
gmRNA_basic <- function(n, r.on, r.degr) {
  cpp_gmRNA_basic(n, r.on, r.degr)
}


#' @rdname gmRNA
#' @export
#' @examples
#' x <- gmRNA_switch(100, 0.23, 0.15, 0.75, 0.001)
#' plot(density(x))
gmRNA_switch <- function(n, r.act, r.deact, r.on, r.degr) {
  cpp_gmRNA_switch(n, r.act, r.deact, r.on, r.degr)
}


#' @rdname gmRNA
#' @export
#' @examples
#' x <- gmRNA_burst(10, 0.15, 0.75, 0.001)
#' plot(density(x))
gmRNA_burst <- function(n, r.burst, s.burst, r.degr) {
  cpp_gmRNA_burst(n, r.burst, s.burst, r.degr)
}


#' @rdname gmRNA
#' @export
#' @examples
#' x <- gmRNA_basic_burst(10, 0.75, 0.15, 0.5, 0.001)
#' plot(density(x))
gmRNA_basic_burst <- function(n, r.on, r.burst, s.burst, r.degr) {
    cpp_gmRNA_basic_burst(n, r.on, r.burst, s.burst, r.degr)
}
