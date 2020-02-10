# prepare packages
libr = function(pkgs) {
    if (length(setdiff(pkgs, rownames(installed.packages()))) > 0)
        install.packages(setdiff(pkgs, rownames(installed.packages())), verbose=F)
    sapply(pkgs, require, character.only=T)
}
libr(c("purrr", "igraph", "RColorBrewer", "ggplot2", "ggrepel", "plotly"))

# load example graph; v contains nodes named by "phenotype" & e contains edges.
load("gr.Rdata")

# source the functions
source("gr.R")

# set plotting parameters
gr <- set_layout_graph(gr, layout_fun=NULL) # set x y coordinates (layout)
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
