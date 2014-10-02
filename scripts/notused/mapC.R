function (x, y = NULL, view = 1, giblocs = NULL, minrange = NA, 
    maxrange = NA, trim.range = 0.98, names = FALSE, value = FALSE, 
    show.na = FALSE, log.data = FALSE, col.pos = c("white", NA, 
        "red"), col.neg = c("white", NA, "blue"), col.na = "gray80", 
    mask.data = NULL, grid = FALSE, title = NULL) 
{
    if (missing(x)) {
        stop("Missing x input")
    }
    if (inherits(x, "HTCexp")) 
        xdata <- intdata(x)
    else if (is.matrix(x)) 
        xdata <- x
    else stop("Input has to belong to 'HTCexp' or 'matrix' classes")
    if (!is.null(y)) {
        view = 2
        if (inherits(y, "HTCexp")) {
            ydata <- intdata(y)
        }
        else if (is.matrix(y)) 
            ydata <- y
        else stop("Input2 has to belong to 'HTCexp' or 'matrix' classes")
    }
    if (log.data) {
        xdata[which(xdata < 0)] <- NA
        xdata[which(xdata > 0)] <- log2(xdata[which(xdata > 0)])
        if (!is.null(y)) {
            ydata[which(ydata < 0)] <- NA
            ydata[which(ydata > 0)] <- log2(ydata[which(ydata > 
                0)])
        }
    }
    if (!is.null(y)) {
        if (!isBinned(x) || !isBinned(y)) 
            stop("x and y have to be binned to plot them on the same scale")
    }
    if (trim.range < 1 && is.na(maxrange) && is.na(minrange)) {
        xmaxrange <- quantile(abs(xdata[xdata != 0]), probs = trim.range, 
            na.rm = TRUE)
        xminrange <- quantile(abs(xdata[xdata != 0]), probs = 1 - 
            trim.range, na.rm = TRUE)
        if (!is.null(y)) {
            ymaxrange <- quantile(abs(ydata[ydata != 0]), probs = trim.range, 
                na.rm = TRUE)
            yminrange <- quantile(abs(ydata[ydata != 0]), probs = 1 - 
                trim.range, na.rm = TRUE)
        }
    }
    else {
        if (is.na(maxrange)) {
            xmaxrange <- max(abs(xdata[xdata != 0]), na.rm = TRUE)
            if (!is.null(y)) 
                ymaxrange <- max(abs(ydata[ydata != 0]), na.rm = TRUE)
        }
        else {
            xmaxrange = maxrange
            if (!is.null(y)) 
                ymaxrange = maxrange
        }
        if (is.na(minrange)) {
            xminrange <- min(abs(xdata[xdata != 0]), na.rm = TRUE)
            if (!is.null(y)) 
                yminrange <- min(abs(ydata[ydata != 0]), na.rm = TRUE)
        }
        else {
            xminrange = minrange
            if (!is.null(y)) 
                yminrange = minrange
        }
    }
    print(paste("minrange=", round(xminrange, 6), " - maxrange=", 
        round(xmaxrange, 6)))
    xdata[which(xdata <= xminrange & xdata > 0)] <- xminrange
    xdata[which(xdata >= -xminrange & xdata < 0)] <- -xminrange
    xdata[which(xdata >= xmaxrange & xdata > 0)] <- xmaxrange
    xdata[which(xdata <= -xmaxrange & xdata < 0)] <- -xmaxrange
    if (!is.null(y)) {
        print(paste("minrange=", round(yminrange, 6), " - maxrange=", 
            round(ymaxrange, 6)))
        ydata[which(ydata <= yminrange & ydata > 0)] <- yminrange
        ydata[which(ydata >= -yminrange & ydata < 0)] <- -yminrange
        ydata[which(ydata >= ymaxrange & ydata > 0)] <- ymaxrange
        ydata[which(ydata <= -ymaxrange & ydata < 0)] <- -ymaxrange
    }
    if (!is.null(giblocs)) {
        if (!inherits(x, "HTCexp")) {
            warning("Cannot diplay genomeIntervals blocs. 'x' has to be a HTCexp object.")
        }
        else {
            if (names) {
                names <- FALSE
                warning("Cannot diplay names and genomeIntervals blocs.")
            }
            if (!is.list(giblocs)) 
                giblocs <- list(giblocs)
            stopifnot(unlist(lapply(giblocs, inherits, "Genome_intervals")))
            ntrack <- length(giblocs)
            sizeblocs <- 0.1
            ygi <- y_intervals(x)
            xgi <- x_intervals(x)
        }
    }
    if (view == 1) {
        if (!is.null(giblocs)) {
            design <- matrix(1:4, 2, 2, byrow = TRUE)
            layout(design, widths = c(sizeblocs * ntrack, 1 - 
                sizeblocs * ntrack), heights = c(sizeblocs * 
                ntrack, 1 - sizeblocs * ntrack))
            par(mar = c(0, 0, 0, 0))
            plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
            addImageTracks(x, giblocs, orientation = "h")
            addImageTracks(x, giblocs, orientation = "v")
        }
        else {
            layout(matrix(1, 1, 1, byrow = TRUE), heights = c(1))
        }
    }
    else if (view == 2) {
        if (!isIntraChrom(x)) 
            stop("The triangle view is available for intrachromosomal data only")
        rx <- range(x)
        if (!is.null(y)) {
            if (seq_name(x) != seq_name(y)) 
                stop("Objects x and y are not on the same chromosome")
            if (range(x)[1] < range(y)[1] && range(x)[2] < range(y)[2]) 
                y <- extractRegion(y, chr = seq_name(x), from = range(y)[1], 
                  to = range(x)[2], exact = TRUE)
            else if (range(x)[1] > range(y)[1] && range(x)[2] > 
                range(y)[2]) 
                y <- extractRegion(y, chr = seq_name(x), from = range(x)[1], 
                  to = range(y)[2], exact = TRUE)
            else if (range(y)[1] < range(x)[1] && range(y)[2] < 
                range(x)[2]) 
                y <- extractRegion(y, chr = seq_name(x), from = range(x)[1], 
                  to = range(y)[2], exact = TRUE)
            else if (range(y)[1] > range(x)[1] && range(y)[2] > 
                range(x)[2]) 
                y <- extractRegion(y, chr = seq_name(x), from = range(y)[1], 
                  to = range(x)[2], exact = TRUE)
            ry <- range(y)
        }
        if (!is.null(giblocs)) {
            if (!is.null(y)) {
                if (diff(rx) >= diff(ry)) {
                  design <- rbind(rep(2, 3), rep(1, 3), c(4, 
                    3, 5))
                  lmar <- (ry[1] - rx[1])/diff(rx)
                  rmar <- (rx[2] - ry[2])/diff(rx)
                  layout(design, heights = c((1 - sizeblocs * 
                    ntrack)/2, sizeblocs * ntrack, (1 - sizeblocs * 
                    ntrack)/2), widths = c(lmar, (1 - rmar - 
                    lmar), rmar))
                  addImageTracks(x, giblocs, orientation = "h")
                }
                else {
                  design <- rbind(c(4, 2, 5), rep(1, 3), rep(3, 
                    3))
                  lmar <- (rx[1] - ry[1])/diff(ry)
                  rmar <- (ry[2] - rx[2])/diff(ry)
                  layout(design, heights = c((1 - sizeblocs * 
                    ntrack)/2, sizeblocs * ntrack, (1 - sizeblocs * 
                    ntrack)/2), widths = c(lmar, (1 - rmar - 
                    lmar), rmar))
                  addImageTracks(x, giblocs, orientation = "h")
                }
            }
            else {
                design <- matrix(2:1, 2, 1, byrow = TRUE)
                layout(design, heights = c(1 - sizeblocs * ntrack, 
                  sizeblocs * ntrack))
                addImageTracks(x, giblocs, orientation = "h")
            }
        }
        else {
            if (!is.null(y)) {
                if (diff(rx) >= diff(ry)) {
                  design <- matrix(c(1, 1, 1, 3, 2, 4), ncol = 3, 
                    byrow = TRUE)
                  lmar <- (ry[1] - rx[1])/diff(rx)
                  rmar <- (rx[2] - ry[2])/diff(rx)
                  layout(design, widths = c(lmar, (1 - rmar - 
                    lmar), rmar), heights = c(1, 1))
                }
                else {
                  design <- matrix(c(3, 1, 4, 2, 2, 2), ncol = 3, 
                    byrow = TRUE)
                  lmar <- (rx[1] - ry[1])/diff(ry)
                  rmar <- (ry[2] - rx[2])/diff(ry)
                  layout(design, widths = c(lmar, (1 - rmar - 
                    lmar), rmar), heights = c(1, 1))
                }
            }
            else {
                layout(matrix(1, 1, 1, byrow = TRUE), heights = c(1))
            }
        }
    }
    if (!names) 
        par(mar = c(0, 0, 0, 0))
    else par(mar = c(mean(sapply(rownames(xdata), nchar))/2, 
        mean(sapply(colnames(xdata), nchar))/2, 0, 0))
    if (view == 1) 
        heatmapC(xdata, names = names, value = value, show.na = show.na, 
            col.pos = col.pos, col.neg = col.neg, col.na = col.na, 
            mask.data = mask.data, grid = grid, title = title)
    else {
        triViewC(xdata, show.na = show.na, col.pos = col.pos, 
            col.neg = col.neg, col.na = col.na)
        if (!is.null(y)) 
            triViewC(ydata, flip = TRUE, show.na = show.na, col.pos = col.pos, 
                col.neg = col.neg, col.na = col.na)
    }
}

