#Color Scheme from Marcos Longo that mimics the color bar in GrADS
#Paleta de cores parecida com a do GrADS
#paleta <- function (n)
#
#  $Id: col.grads.r,v 1.2 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

col.grads <- function (n)
{
    if ((n <- as.integer(n[1])) > 0) {
        k <- as.integer(n%/%2)
        h <- c(240/360,   42/360,  30/360)
        s <- c(   1.00,     0.46,    1.00)
        v <- c(   0.51,     0.99,    1.00)
        c(hsv(h = seq(h[1], h[2], length = k), s = seq(s[1],
            s[2], length = k), v = seq(v[1], v[2], length = k)),
            hsv(h = seq(h[2], h[3], length = n - k + 1)[-1],
                s = seq(s[2], s[3], length = n - k + 1)[-1],
                v = seq(v[2], v[3], length = n - k + 1)[-1]))
    }
    else character(0)
}

#A mesma paleta, so que na ordem inversa...
#ipaleta <- function (n)
icol.grads <- function (n)
{
    if ((n <- as.integer(n[1])) > 0) {
        k <- as.integer(n%/%2)
        h <- c( 30/360,   42/360, 240/360)
        s <- c(   1.00,     0.46,    1.00)
        v <- c(   1.00,     0.99,    0.51)
        c(hsv(h = seq(h[1], h[2], length = k), s = seq(s[1],
            s[2], length = k), v = seq(v[1], v[2], length = k)),
            hsv(h = seq(h[2], h[3], length = n - k + 1)[-1],
                s = seq(s[2], s[3], length = n - k + 1)[-1],
                v = seq(v[2], v[3], length = n - k + 1)[-1]))
    }
    else character(0)
}

