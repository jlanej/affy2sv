

# plotMad
# -----------------------------------------------------------------------------
# 24/04/2014

#' @export
plotMad <- function(object, colChr, colPos, colLRR, colGType, colBAF, chr, 
    save.as, bp.from, bp.to, main="", legend=FALSE, rt=FALSE, pch=20, 
    pch.cex=2, wixh=par("din")[1], height=par("din")[2]) {
    
    if(length(class(object)) == 1) {
        if(class(object) == "character") {
            object <- fread(object)
        }
    } else {
        if(length(class(object)) == 2) {
            if(class(object)[1] != "data.table") {
                stop("Invalid object type. Allowed: 'character' or 'data.table'.")
            }
        } else {
            stop("Invalid object type. Allowed: 'character' or 'data.table'.")
        }
    }

    oldnames <- colnames(object)
    setnames(object, colnames(object), c("Name", "Chr", "Position", "Log.R.Ratio", "GType", "B.Allele.Freq"))

    color.baf <- c("red", "salmon", "red", "skyblue")
    names(color.baf) <- c("AA", "AB", "BB", "NA")

    # filter data
    object$Chr <- as.character(object$Chr)
    object <- object[object$Chr == chr, ]
    object$Position <- as.numeric(as.character(object$Position))
    object <- object[order(object$Position), ]
    object$Log.R.Ratio[object$Log.R.Ratio < -1] <- -1
    object$Log.R.Ratio[object$Log.R.Ratio > 1] <- 1

    if (!missing(bp.from)) {
        object <- object[object$Position >= bp.from, ]
    }
    if (!missing(bp.to)) {
        object <- object[object$Position <= bp.to, ]
    }
    # /

    if (!rt) {
        grid.newpage()
    }

    # create LRR and BAF plots
    plot.lrr <- ggplot(data = object, aes(x = Position, y = Log.R.Ratio, color="LRR"))
    plot.lrr <- plot.lrr + geom_point(shape = pch, size=pch.cex) + ggtitle(main)
    plot.lrr <- plot.lrr + theme(plot.title = element_text(lineheight = 0.8, face = "bold"))
    plot.lrr <- plot.lrr + ylim(-1, 1)
    if (legend) {
        plot.lrr <- plot.lrr + scale_color_manual(name = "", labels = c("LRR"), values = c("black"))
        plot.lrr <- plot.lrr + theme(legend.position = "bottom")
    } else {
        plot.lrr <- plot.lrr + scale_color_manual(name = "", labels = c("LRR"), values = c("black"), guide = FALSE)
    }

    plot.baf <- ggplot(data = object, aes(x = Position, y = B.Allele.Freq, color=GType))
    plot.baf <- plot.baf + geom_point(shape = pch, size=pch.cex)
    plot.baf <- plot.baf + theme(panel.background = element_rect(fill = NA), panel.grid = element_blank())
    plot.baf <- plot.baf + scale_color_manual(name = "GType", labels = c("AA", "AB", "BB", "NA"), values = color.baf)
    plot.baf <- plot.baf + ylim(0, 1)
    if (legend) {
        plot.baf <- plot.baf + theme(legend.position = "bottom")
    }
    # /

    # extract gtable from each plot
    gtab.lrr <- ggplot_gtable(ggplot_build(plot.lrr))
    gtab.baf <- ggplot_gtable(ggplot_build(plot.baf))
    # /

    # overlap the panel of BAF plot on the one of LRR plot
    pp <- c(subset(gtab.lrr$layout, name == "panel", se = t:r))
    plot <- gtable_add_grob(gtab.lrr, gtab.baf$grobs[[which(gtab.baf$layout$name == "panel")]], pp$t,
        pp$l, pp$b, pp$l)
    # /

    # mixing axis
    ia <- which(gtab.baf$layout$name == "axis-l")
    ga <- gtab.baf$grobs[[ia]]
    ax <- ga$children[[2]]
    ax$wixhs <- rev(ax$wixhs)
    ax$grobs <- rev(ax$grobs)
    ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
    plot <- gtable_add_cols(plot, gtab.baf$wixhs[gtab.baf$layout[ia, ]$l], length(plot$wixhs) - 1)
    plot <- gtable_add_grob(plot, ax, pp$t, length(plot$wixhs) - 1, pp$b)
    # /

    # second label for y axis (BAF)
    ia <- which(gtab.baf$layout$name == "ylab")
    ylab <- gtab.baf$grobs[[ia]]
    plot <- gtable_add_cols(plot, plot$wixhs[plot$layout[ia, ]$l], length(plot$wixhs) - 1)
    plot <- gtable_add_grob(plot, ylab, pp$t, length(plot$wixhs) - 1, pp$b)
    # /

    # seconds legend for y axis (BAF)
    if (legend) {
        ia <- which(gtab.baf$layout$name == "guide-box")
        #yleg <- gtab.baf$grobs[[ia]]
        plot <- gtable_add_rows(plot, gtab.baf$heights[gtab.baf$layout[ia, ]$b] - unit(0.15, "cm"), length(gtab.baf$heights) - 1)
        #plot <- gtable_add_grob(plot, yleg, length(plot$heights) - 1, pp$l, pp$b)
    }
    # /

    if (rt) {
        return(plot)
    }

    if (!missing(save.as)) {
        ggsave(filename = save.as, plot = plot, wixh = wixh, height = height)
    } else {
        grid.draw(plot)
    }
}
