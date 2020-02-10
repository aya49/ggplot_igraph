require(purrr)
require(igraph)
require(RColorBrewer)
require(ggplot2)
require(ggrepel)
require(plotly)

#' @title Determines cell hierarchy layout.
#' @description Determines cell hierarchy layout and returns the X, Y coordinate
#'  of each cell population.
#' @param gr A list containing data frames \code{e} and \code{v}.
#' @param layout_fun A function from the \code{igraph} package that
#'  indicates what layout should be used if a cell hierarchy is to be ploted;
#'  all such functions have prefix \code{layout_}.
#' @return A list containing data frames \code{e} and \code{v}; each data frame
#'  contains an X, Y column or coordinate for each node and edge.
#' @examples
#'
#'
#' @seealso
#'  \code{\link[flowGraph]{set_layout}}
#'  \code{\link[igraph]{as_data_frame}}
#'  \code{\link[igraph]{layout.reingold.tilford}}
#'  \code{\link[igraph]{V}}
#' @rdname set_layout_graph
#' @export
#' @importFrom igraph graph_from_data_frame layout.reingold.tilford V
#' @importFrom purrr map_int
set_layout_graph <- function(gr, layout_fun=NULL) {
    # layout.circle NULL <-
    # layout.reingold.tilford FUN is a layout
    # function from the igraph package assume
    # graph is connected, used internally

    gr_e <- gr$e
    gr_v <- gr$v

    # edit layout
    gr0 <- igraph::graph_from_data_frame(gr_e)
    if (base::is.null(layout_fun)) {
        gr_vxy_ <- igraph::layout.reingold.tilford(gr0)
    } else {
        gr_vxy_ <- layout_fun(gr0)
    }
    gr_vxy <- base::as.data.frame(gr_vxy_)
    if (base::is.null(layout_fun)) {
        ## edit layout manually
        gys <- base::sort(base::unique(gr_vxy[, 2]))
        gxns <- purrr::map_int(gys, function(y) sum(gr_vxy[, 2] == y))
        maxlayern <- max(gxns)
        maxlayer <- gys[base::which.max(gxns)]
        maxlayertf <- gr_vxy[, 2] == maxlayer
        gr_vxy[maxlayertf, 1] <- base::rank(gr_vxy[maxlayertf, 1]) - 1
        for (gy in gys) {
            if (gy == maxlayer)
                next()
            layertf <- gr_vxy[, 2] == gy
            gr_vxy[layertf, 1]=base::rank(gr_vxy[layertf, 1]) - 1 +
                floor((maxlayern - sum(layertf))/2)
        }
        # turn plot sideways
        gr_vxy <- gr_vxy[, 2:1]
        xmax <- max(gr_vxy[, 1])
        if (gr_vxy[1, 1] == xmax)
            gr_vxy[, 1] <- xmax - gr_vxy[, 1]
    }

    # get node
    base::colnames(gr_vxy) <- c("x", "y")
    gr_v_xy <- gr_vxy[base::match(gr_v$phenotype, names(igraph::V(gr0)[[]])), ]
    gr_v$x <- gr_v_xy$x
    gr_v$y <- gr_v_xy$y

    # get edge
    e_match <- base::match(gr_e$from, gr_v$phenotype)
    gr_e$from.x <- gr_v$x[e_match]
    gr_e$from.y <- gr_v$y[e_match]
    e_match2 <- base::match(gr_e$to, gr_v$phenotype)
    gr_e$to.x <- gr_v$x[e_match2]
    gr_e$to.y <- gr_v$y[e_match2]

    return(base::list(e=gr_e, v=gr_v))
}

#' @title Creates a cell hierarchy plot.
#' @description Creates a cell hierarchy plot given a list of nodes and edges.
#' @param gr A list containing data frames \code{e} and \code{v}.
#' @param main A string containing the plot title. If this is set to NULL, the
#'  function will look for a plot title in the \code{main} slot of \code{gr};
#'  otherwise, this defaults to "".
#' @param bgedges A logical variable indicating whether or not
#'  edges not specified for plotting should be plotted as light grey
#'  in the background. If this is \code{NULL}, the function will look for a
#'  \code{bgedges} in the \code{bgedges} slot of \code{gr};
#'  otherwise, this defaults to \code{TRUE}.
#' @param interactive A logical variable indicating whether the plot should be
#'  an interacttive plot; see package \code{plotly}.
#' @param colour_palette A colour palette e.g. the default palette if the user
#'  sets this to \code{NULL} is \code{rev(RColorBrewer::brewer.pal(10,"RdBu"))}.
#' @param ... Other parameters for \code{ggplot} if \code{interactive}
#'  is set to \code{FALSE}; other parameters for \code{plot_ly}
#'  if \code{interactive} is set to \code{TRUE}.
#' @return  A \code{ggplot} object if \code{interactive} is set to \code{FALSE};
#'  a \code{plotly} object if \code{interactive} is set to \code{TRUE}.
#' @examples
#'
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_plot}}
#'  \code{\link[flowGraph]{get_phen_meta}}
#'  \code{\link[flowGraph]{ggdf}}
#'  \code{\link[flowGraph]{fg_get_feature}}
#'  \code{\link[flowGraph]{fg_get_summary}}
#'  \code{\link[ggplot2]{labs}},
#'  \code{\link[ggplot2]{scale_colour_brewer}},
#'  \code{\link[ggplot2]{geom_segment}},
#'  \code{\link[ggplot2]{aes}},
#'  \code{\link[ggplot2]{geom_point}}
#'  \code{\link[ggrepel]{geom_label_repel}}
#'  \code{\link[plotly]{plot_ly}}
#' @rdname gggraph
#' @export
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous
#'  theme element_blank ggtitle scale_fill_brewer geom_segment geom_point
#'  scale_colour_gradientn
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggrepel geom_label_repel
#' @importFrom plotly plot_ly layout
#' @importFrom purrr map
gggraph <- function(gr, main=NULL, bgedges=NULL, interactive=FALSE,
                    colour_palette=NULL, ...) {
    # gr_v: name x y label size colour sizeb colourb
    # gr_e: from to from.x from.y to.x to.y colour

    # place holder; all variables specified in plotting functions
    # refers to columns names in the given data frame
    size <- 0

    gr_v <- gr$v
    gr_e <- gr$e
    gr_e$e_ind <- gr$e$from%in%gr$v$phenotype[gr$v$v_ind] &
        gr$e$to%in%gr$v$phenotype[gr$v$v_ind]

    if (base::is.null(main))
        if (base::is.null(gr$main)) {
            main <- ""
        } else {
            main <- gr$main
        }
    if (base::is.null(bgedges))
        if (base::is.null(gr$bgedges)) {
            bgedges <- TRUE
        } else {
            bgedges <- gr$bgedges
        }

    if (!interactive) {
        if(base::is.null(colour_palette))
            colour_palette <- rev(RColorBrewer::brewer.pal(10,"RdBu"))

        gp <- ggplot2::ggplot(gr_v, ggplot2::aes(x=x, y=y, ...)) +
            ggplot2::scale_x_continuous(expand=c(0,1)) +  # expand x limits
            ggplot2::scale_y_continuous(expand=c(0,1)) + # expand y limits
            # theme_bw()+  # use the ggplot black and white theme
            ggplot2::theme(
                axis.text.x=ggplot2::element_blank(),  # rm x-axis text
                axis.text.y=ggplot2::element_blank(), # rm y-axis text
                axis.ticks=ggplot2::element_blank(),  # rm axis ticks
                axis.title.x=ggplot2::element_blank(), # rm x-axis labels
                axis.title.y=ggplot2::element_blank(), # rm y-axis labels
                panel.background=ggplot2::element_blank(),
                panel.border=ggplot2::element_blank(),
                panel.grid.major=ggplot2::element_blank(),  #rm grid labels
                panel.grid.minor=ggplot2::element_blank(),  #rm grid labels
                plot.background=ggplot2::element_blank())

        # base graph
        gp <- gp + ggplot2::ggtitle(main)

        if (bgedges)  # keep greyed out edges on
            gp <- gp +
            ggplot2::geom_segment(
                data=gr_e[!gr_e$e_ind,], color="grey",
                ggplot2::aes(x=from.x,xend=to.x, y=from.y,yend=to.y))

        gp <- gp +
            ggplot2::geom_segment(
                data=gr_e[gr_e$e_ind,],
                color="grey50",
                ggplot2::aes(x=from.x,xend=to.x, y=from.y,yend=to.y)) +
            ggplot2::geom_segment(
                data=gr_e[gr_e$e_ind & gr_e$colour!=0,],
                ggplot2::aes(x=from.x,xend=to.x, y=from.y,yend=to.y,
                             color=colour)) +
            ggplot2::geom_point(
                data=gr_v[gr_v$v_ind,],
                ggplot2::aes(x=x,y=y, color=colour, size=size)) +
            ggrepel::geom_label_repel(
                data=gr_v[gr_v$label_ind,],
                ggplot2::aes(x=x, y=y,label=label, color=colour),
                nudge_x=-.1, direction="y", hjust=1,
                segment.size=0.2,force=1.5) +
            ggplot2::scale_colour_gradientn(colours=colour_palette)

    } else {
        vx <- gr_v$x[gr_v$v_ind]
        vy <- gr_v$y[gr_v$v_ind]
        gp <- plotly::plot_ly(x=~vx, y=~vy, mode="markers", opacity=1,
                              color=gr_v$colour[gr_v$v_ind],
                              size=gr_v$size[gr_v$v_ind],
                              text=gr_v$label_long[gr_v$v_ind],
                              hoverinfo="text", ...)

        gr_e$colour_new <- gr_e$colour
        gr_e$colour_new[gr_e$colour<0] <- "turquoise"
        gr_e$colour_new[gr_e$colour>0] <- "cherry"
        gr_e$colour_new[gr_e$colour==0] <- "grey50"
        if (!bgedges) {
            el <- purrr::map(which(gr_e$e_ind), function(j) {
                list(type="line", line=list(color=gr_e$colour[j], width=.5),
                     x0=gr_e$from.x[j], x1=gr_e$to.x[j],
                     y0=gr_e$from.y[j], y1=gr_e$to.y[j])
            })
        } else {
            el <- purrr::map(seq_len(nrow(gr_e)), function(j) {
                list(type="line", line=list(color=gr_e$colour[j], width=.5),
                     x0=gr_e$from.x[j], x1=gr_e$to.x[j],
                     y0=gr_e$from.y[j], y1=gr_e$to.y[j])
            })
        }

        axis <- list(title="", showgrid=FALSE, showticklabels=FALSE,
                     zeroline=FALSE)

        gp <- plotly::layout(gp, title=main, shapes=el, xaxis=axis, yaxis=axis)
    }

    return(gp)
}


################### START HEREs ###################

load(gr)

gr <- set_layout_graph(gr) # set x y coordinates (layout)
gr <- ggdf(gr) # add plotting parameters to data frames, editted below
gr$v$colour <- rnorm(nrow(gr$v),0,2) # random node colours
gr$v$label <- gr$v$phenotype # give some node labels
gr$v$v_ind <- gr$v$label_ind <- TRUE # which nodes to colour in and label
gr$e$e_ind <- !grepl("[-]",gr$e$from) & !grepl("[-]",gr$e$to) # which edges to show

# create and save ggplot
gch <- gggraph(gr, main="Example cell hierarchy with markers A, B, C, D")
ggplot2::ggsave("gr.png",
                plot=gch, scale=1, width=9, height=9,
                units="in", dpi=500, limitsize=TRUE)
