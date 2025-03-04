#!/usr/bin/env Rscript

# Copyright Genome Research Limited (c) 2025
#
# Author Martin Pollard <mp15@sanger.ac.uk>
#
# This file is part of arc.chaos.
#
# arc.chaos is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# arc.chaos is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with arc.chaos If not, see <https://www.gnu.org/licenses/>.



#  Borrowed from yilong li plot.rearrangements
arc = function(x0, x1, y, xr, yr, col, lwd) {
    x = (x0 + x1)/2  # Center of arc
    xr = x - x0 	 # x-radius of arc

    apply(
        cbind(x, y, xr, yr, col),
        1,
        function(z) {
            x   = as.numeric(z[1])
            y   = as.numeric(z[2])
            xr  = as.numeric(z[3])
            yr  = as.numeric(z[4])
            col = z[5]
            x_points = seq(x - xr, x + xr, length.out = 200)
            y_points = y + yr * sqrt( 1  -  ( (x_points-x) / xr )^2 )
            
            lines(
                x_points,
                y_points,
                col = col,
                lwd = lwd
            )
        }
    )

    return()
}

#' @param blocks is a data frame containing the CN blocks and their depth
#' @param arcs is a data from containing links, start, end, strand
#' @export
plot_arc_chaos <- function(blocks, links, chr, start, end){
    ylen <- end - start
    graphics::plot(c(), xlim=c(start,end), ylim=c(0,50))
#    plotrix::gap.plot(c(0),c(0), gap=c(81000000,90500000,91000000,108000000), gap.axis = "x", xlim = c(start,end), ylim = c(0,50))

    apply(blocks, 1,
        function(z)
        {
            start <- as.numeric(z[3])
            end <- as.numeric(z[4])
            graphics::rect(xleft = start, ybottom = 0, xright = end, ytop = 30)
        }
    )
    apply(links, 1,
        function(line)
        {
            start1 <- as.numeric(line[2])
            start2 <- as.numeric(line[5])
            strand1 = line[9]
            strand2 = line[10]

            strand1dir <- if (strand1 == "+") ylen/50 else -(ylen/50)
            strand2dir <- if (strand2 == "+") ylen/50 else -(ylen/50)
            arc(x0 = start1, x1 = start2, y = 40, xr = 1, yr = 5, col = "pink", lwd = 1)
            graphics::lines(rbind(c(start1, 30),c(start1, 40)), col ="yellow")
            graphics::lines(rbind(c(start1, 37),c(start1 + strand1dir, 37)))
            graphics::lines(rbind(c(start2, 30),c(start2, 40)), col ="yellow")
            graphics::lines(rbind(c(start2, 35),c(start2 + strand2dir, 35)))
        }
    )
}

chr <- "chr5"
start <- 80200000
end <- 126000000
blocks <- read.csv('lineage_blocks.csv')
links <- read.csv('lineage_arcs.csv')
plot_arc_chaos(blocks, links, chr, start, end)