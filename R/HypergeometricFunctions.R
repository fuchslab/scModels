
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR Description.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


#' The kummer function
#'
#' A function implemented by Diethelm Wuertz
#' Calculate the Confluent Hypergeometric Function of the First Kind
#' The kummer_gsl uses the GNU Scientific Library implementation.
#'
#' @param x complex function argument
#' @param a complex index a
#' @param b complex index b
#' @param lnchf an integer indicating if the output is logarithmic or not
#' @param ip an integer indicating the number of terms to compute in
#'  the kummer series
#'
#' @name KummerM
#' @rdname KummerM
#' @useDynLib mpb2
#' @export


kummerM =
function(x, a, b, lnchf = 0, ip = 0)
{

    # You can also input real arguments:
    if (!is.complex(x)) x = complex(real = x, imaginary = 0*x)
    if (!is.complex(a)) a = complex(real = a, imaginary = 0)
    if (!is.complex(b)) b = complex(real = b, imaginary = 0)

    # Calculate KummerM:
    chm = rep(complex(real = 0, imaginary = 0), length = length(x))
    for(i in 1:length(x)) {
      if(Re(x[i]) < 0 && Im(x[i]) == 0) {
        value = .Fortran("chfm",
                       as.double(-Re(x[i])),
                       as.double(Im(x[i])),
                       as.double(Re(b)-Re(a)),
                       as.double(Im(a)),
                       as.double(Re(b)),
                       as.double(Im(b)),
                       as.double(Re(chm[i])),
                       as.double(Im(chm[i])),
                       as.integer(length(x[i])),
                       as.integer(lnchf),
                       as.integer(ip),
                       PACKAGE = "mpb2")
        value[[7]] = value[[7]] + Re(x[i])
      } else {
        value = .Fortran("chfm",
                         as.double(Re(x[i])),
                         as.double(Im(x[i])),
                         as.double(Re(a)),
                         as.double(Im(a)),
                         as.double(Re(b)),
                         as.double(Im(b)),
                         as.double(Re(chm[i])),
                         as.double(Im(chm[i])),
                         as.integer(length(x[i])),
                         as.integer(lnchf),
                         as.integer(ip),
                         PACKAGE = "mpb2")
      }
    }
    result = complex(real = value[[7]], imaginary = value[[8]])

    # Return Value:
    result
}
