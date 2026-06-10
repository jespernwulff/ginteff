## ginteffplot.R
## R port of Stata's ginteffplot command -- visualises the AIE point
## estimate with capped-spike CI on a horizontal axis, and (optionally)
## the per-observation interaction effects as a thin strip below.
##
## Layout matches Stata's default: x-axis is the interaction-effect
## value, AIE row sits on top with a capped spike, individual effects
## as a tick strip below.
##
## Uses ggplot2 if available, otherwise falls back to base graphics.


#' Plot a ginteff result
#'
#' Draws the average interaction effect (AIE) on a horizontal axis with
#' its capped-spike confidence interval, optionally overlaid with the
#' per-observation interaction effects as a thin strip below.
#'
#' @param x         A `ginteff` object.
#' @param output    Integer index of the row of `aie` to plot when there
#'                  is more than one result (e.g. multi-level factor).
#'                  Default 1. Pass `NULL` (or `"all"`) to draw every row
#'                  in one figure -- stacked panels sharing a common
#'                  x-axis, the analogue of Stata's `xcommon(*)`.
#' @param obseff    Logical. If TRUE, the per-observation interaction
#'                  effects are plotted as a tick strip below the AIE.
#'                  Requires that `ginteff()` was called with `obseff = TRUE`.
#' @param against   Optional name of a covariate. If supplied alongside
#'                  `obseff = TRUE`, the per-obs effects are plotted
#'                  against that covariate (a scatter on the y-axis is
#'                  produced instead of the strip layout). This is an R
#'                  extension and is not part of Stata's `ginteffplot`.
#' @param data      Data frame to draw `against` from. Required when
#'                  `against` is set.
#' @param zeroline  Logical. Draw a vertical reference line at x = 0.
#'                  Default TRUE.
#' @param xlab      Optional x-axis label override.
#' @param point_color,range_color,obs_color  Colours.
#' @param ...       Ignored (reserved for future options).
#'
#' @return A ggplot object (if ggplot2 is available) or NULL (base graphics).
#' @export
ginteffplot <- function(x,
                        output     = 1L,
                        obseff     = FALSE,
                        against    = NULL,
                        data       = NULL,
                        zeroline   = TRUE,
                        xlab       = NULL,
                        point_color = "black",
                        range_color = "black",
                        obs_color   = "grey40",
                        ...) {
    stopifnot(inherits(x, "ginteff"))

    has_obs <- isTRUE(obseff) && !is.null(x$obseff)
    if (isTRUE(obseff) && is.null(x$obseff))
        warning("`obseff = TRUE` ignored: ginteff() was called without ",
                "`obseff = TRUE`.")

    ## output = NULL / "all": every row as a stacked panel on a common
    ## x-axis (Stata's ginteffplot with xcommon(*)).
    if (is.null(output) || identical(output, "all")) {
        if (!is.null(against))
            stop("`against` is not supported with output = NULL; ",
                 "plot one output at a time.")
        if (length(x$aie) == 1L) {
            output <- 1L
        } else {
            return(.ginteffplot_all(x, has_obs, zeroline, xlab,
                                    point_color, range_color, obs_color))
        }
    }

    if (output < 1L || output > length(x$aie))
        stop("`output` out of range; ginteff produced ",
             length(x$aie), " result(s).")

    aie <- as.numeric(x$aie[output])
    se  <- as.numeric(x$se[output])
    lo  <- as.numeric(x$ci[output, 1])
    hi  <- as.numeric(x$ci[output, 2])
    lab <- x$labels[output]

    obs_v <- if (has_obs) as.numeric(x$obseff[, output]) else NULL
    obs_x <- NULL
    if (has_obs && !is.null(against)) {
        if (is.null(data))
            stop("Pass `data = ` so the `against` covariate can be retrieved.")
        if (!against %in% names(data))
            stop("Variable `", against, "` not found in data.")
        obs_x <- data[[against]]
        if (length(obs_x) != length(obs_v))
            stop("`data` has ", length(obs_x), " rows but ginteff stored ",
                 length(obs_v), " obs effects.")
    }

    if (is.null(xlab)) xlab <- "Interaction effect"

    if (requireNamespace("ggplot2", quietly = TRUE)) {
        return(.ginteffplot_gg(aie, lo, hi, lab, obs_v, obs_x, against,
                               zeroline, xlab,
                               point_color, range_color, obs_color,
                               level = x$level))
    }
    .ginteffplot_base(aie, lo, hi, lab, obs_v, obs_x, against,
                      zeroline, xlab,
                      point_color, range_color, obs_color,
                      level = x$level)
}


# ---- ggplot2 backend -------------------------------------------------------

.ginteffplot_gg <- function(aie, lo, hi, lab, obs_v, obs_x, against,
                            zeroline, xlab,
                            point_color, range_color, obs_color,
                            level) {

    ## Mode 2 (R extension): per-obs effects vs a covariate.
    if (!is.null(obs_x)) {
        return(.ginteffplot_gg_against(aie, lo, hi, lab, obs_v, obs_x,
                                       against, zeroline,
                                       point_color, range_color, obs_color,
                                       level))
    }

    ## Mode 1 (Stata default): horizontal AIE row + obs strip below.
    ## Y positions: 2 = AIE row, 1 = obs strip.
    df_aie <- data.frame(x = aie, y = 2, lo = lo, hi = hi)

    p <- ggplot2::ggplot()

    if (zeroline)
        p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dotted",
                                     colour = "grey50")

    if (!is.null(obs_v)) {
        df_obs <- data.frame(x = obs_v, y = 1)
        p <- p + ggplot2::geom_point(
            data = df_obs,
            mapping = ggplot2::aes(x = x, y = y),
            shape = "|", size = 4, alpha = 0.5,
            colour = obs_color)
    }

    p <- p +
        ## capped-spike CI for the AIE
        ggplot2::geom_errorbarh(
            data = df_aie,
            mapping = ggplot2::aes(xmin = lo, xmax = hi, y = y),
            height = 0.15, colour = range_color, linewidth = 0.6) +
        ## AIE marker (filled square, like Stata's default)
        ggplot2::geom_point(
            data = df_aie,
            mapping = ggplot2::aes(x = x, y = y),
            shape = 15, size = 3.5, colour = point_color) +
        ## "AIE" text label above the marker
        ggplot2::geom_text(
            data = df_aie,
            mapping = ggplot2::aes(x = x, y = y, label = "AIE"),
            vjust = -1.1, size = 3.5, colour = point_color)

    ## Y axis: hide ticks, reserve space top/bottom
    y_breaks <- if (!is.null(obs_v)) c(1, 2) else 2
    y_labels <- if (!is.null(obs_v)) c("Individual", "AIE") else "AIE"

    p <- p +
        ggplot2::scale_y_continuous(
            breaks = y_breaks, labels = y_labels,
            limits = c(0.5, 2.6),
            expand = ggplot2::expansion(add = c(0.1, 0.3))) +
        ggplot2::labs(
            x = xlab, y = NULL,
            title = sprintf("Average interaction effect: %s", lab),
            subtitle = sprintf("Capped spike: %.0f%% CI (delta method)",
                               100 * level)) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            axis.line.y        = ggplot2::element_blank(),
            axis.ticks.y       = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_blank())

    p
}


.ginteffplot_gg_against <- function(aie, lo, hi, lab, obs_v, obs_x,
                                    against, zeroline,
                                    point_color, range_color, obs_color,
                                    level) {
    df <- data.frame(x = obs_x, y = obs_v)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_point(colour = obs_color, alpha = 0.4, size = 1) +
        ggplot2::geom_smooth(method = "loess", se = FALSE,
                             formula = y ~ x, colour = "steelblue") +
        ggplot2::geom_hline(yintercept = aie, colour = point_color,
                            linetype = "dashed") +
        ggplot2::annotate(
            "text",
            x = max(df$x, na.rm = TRUE),
            y = aie,
            label = "AIE", hjust = 1.1, vjust = -0.5,
            colour = point_color) +
        ggplot2::labs(x = against, y = "Interaction effect",
                      title = sprintf("Per-observation interaction effects: %s",
                                      lab),
                      subtitle = sprintf("Dashed line: AIE; %.0f%% CI not shown",
                                         100 * level)) +
        ggplot2::theme_bw()

    if (zeroline)
        p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dotted",
                                     colour = "grey50")
    p
}


# ---- base graphics fallback ------------------------------------------------

.ginteffplot_base <- function(aie, lo, hi, lab, obs_v, obs_x, against,
                              zeroline, xlab,
                              point_color, range_color, obs_color,
                              level) {

    if (!is.null(obs_x)) {
        ## "against" view: scatter
        graphics::plot.default(obs_x, obs_v,
             pch = 20, col = obs_color,
             xlab = against, ylab = "Interaction effect",
             main = sprintf("Per-observation interaction effects: %s", lab))
        graphics::abline(h = aie, col = point_color, lty = 2)
        if (zeroline) graphics::abline(h = 0, col = "grey50", lty = 3)
        return(invisible(NULL))
    }

    ## Stata-style horizontal layout
    xs   <- c(lo, hi)
    if (!is.null(obs_v)) xs <- c(xs, obs_v)
    pad  <- 0.04 * diff(range(xs, na.rm = TRUE))
    xlim <- range(xs, na.rm = TRUE) + c(-pad, pad)
    ylim <- c(0.5, 2.6)
    has_obs <- !is.null(obs_v)

    graphics::plot.new()
    graphics::plot.window(xlim = xlim, ylim = ylim)
    graphics::axis(1)
    graphics::title(main = sprintf("Average interaction effect: %s", lab),
                    sub  = sprintf("Capped spike: %.0f%% CI", 100 * level),
                    xlab = xlab)
    if (has_obs) {
        graphics::axis(2, at = c(1, 2),
                       labels = c("Individual", "AIE"),
                       las = 1, lwd = 0)
    } else {
        graphics::axis(2, at = 2, labels = "AIE", las = 1, lwd = 0)
    }
    if (zeroline) graphics::abline(v = 0, col = "grey50", lty = 3)

    if (has_obs) {
        ## thin tick strip at y = 1
        graphics::segments(obs_v, 1 - 0.12, obs_v, 1 + 0.12,
                           col = grDevices::adjustcolor(obs_color, 0.5))
    }

    ## AIE capped spike at y = 2
    graphics::segments(lo, 2, hi, 2, col = range_color, lwd = 1.5)
    cap <- 0.10
    graphics::segments(lo, 2 - cap, lo, 2 + cap,
                       col = range_color, lwd = 1.5)
    graphics::segments(hi, 2 - cap, hi, 2 + cap,
                       col = range_color, lwd = 1.5)
    graphics::points(aie, 2, pch = 15, col = point_color, cex = 1.4)
    graphics::text(aie, 2 + 0.18, "AIE", col = point_color, cex = 0.9)

    invisible(NULL)
}


# ---- all-outputs (faceted / stacked) view ----------------------------------

.ginteffplot_all <- function(x, has_obs, zeroline, xlab,
                             point_color, range_color, obs_color) {
    k    <- length(x$aie)
    labs <- x$labels
    if (is.null(xlab)) xlab <- "Interaction effect"

    if (requireNamespace("ggplot2", quietly = TRUE)) {
        df_aie <- data.frame(
            label = factor(labs, levels = labs),
            x  = as.numeric(x$aie),
            lo = as.numeric(x$ci[, 1]),
            hi = as.numeric(x$ci[, 2]),
            y  = 2
        )
        p <- ggplot2::ggplot()
        if (zeroline)
            p <- p + ggplot2::geom_vline(xintercept = 0,
                                         linetype = "dotted",
                                         colour = "grey50")
        if (has_obs) {
            df_obs <- data.frame(
                label = factor(rep(labs, each = nrow(x$obseff)),
                               levels = labs),
                x = as.numeric(x$obseff),
                y = 1
            )
            p <- p + ggplot2::geom_point(
                data = df_obs,
                mapping = ggplot2::aes(x = x, y = y),
                shape = "|", size = 4, alpha = 0.5, colour = obs_color)
        }
        p <- p +
            ggplot2::geom_errorbarh(
                data = df_aie,
                mapping = ggplot2::aes(xmin = lo, xmax = hi, y = y),
                height = 0.15, colour = range_color, linewidth = 0.6) +
            ggplot2::geom_point(
                data = df_aie,
                mapping = ggplot2::aes(x = x, y = y),
                shape = 15, size = 3, colour = point_color) +
            ggplot2::facet_wrap(~label, ncol = 1) +
            ggplot2::scale_y_continuous(
                breaks = if (has_obs) c(1, 2) else 2,
                labels = if (has_obs) c("Individual", "AIE") else "AIE",
                limits = c(0.5, 2.6)) +
            ggplot2::labs(
                x = xlab, y = NULL,
                title = "Average interaction effects",
                subtitle = sprintf(
                    "Capped spikes: %.0f%% CIs (delta method); common x-axis",
                    100 * x$level)) +
            ggplot2::theme_classic() +
            ggplot2::theme(
                axis.line.y        = ggplot2::element_blank(),
                axis.ticks.y       = ggplot2::element_blank(),
                panel.grid.major.y = ggplot2::element_blank(),
                strip.background   = ggplot2::element_blank(),
                strip.text         = ggplot2::element_text(hjust = 0))
        return(p)
    }

    ## base-graphics fallback: stacked panels via mfrow, common xlim
    op <- graphics::par(mfrow = c(k, 1), mar = c(4, 6, 2, 1))
    on.exit(graphics::par(op), add = TRUE)
    xs <- c(x$ci)
    if (has_obs) xs <- c(xs, as.numeric(x$obseff))
    pad  <- 0.04 * diff(range(xs, na.rm = TRUE))
    xlim <- range(xs, na.rm = TRUE) + c(-pad, pad)
    for (i in seq_len(k)) {
        obs_i <- if (has_obs) as.numeric(x$obseff[, i]) else NULL
        graphics::plot.new()
        graphics::plot.window(xlim = xlim, ylim = c(0.5, 2.6))
        graphics::axis(1)
        graphics::title(main = x$labels[i], xlab = xlab, cex.main = 0.9)
        if (zeroline) graphics::abline(v = 0, col = "grey50", lty = 3)
        if (!is.null(obs_i))
            graphics::segments(obs_i, 1 - 0.12, obs_i, 1 + 0.12,
                               col = grDevices::adjustcolor(obs_color, 0.5))
        lo <- x$ci[i, 1]; hi <- x$ci[i, 2]
        graphics::segments(lo, 2, hi, 2, col = range_color, lwd = 1.5)
        graphics::segments(c(lo, hi), 2 - 0.1, c(lo, hi), 2 + 0.1,
                           col = range_color, lwd = 1.5)
        graphics::points(x$aie[i], 2, pch = 15, col = point_color,
                         cex = 1.4)
    }
    invisible(NULL)
}

#' @export
plot.ginteff <- function(x, ...) ginteffplot(x, ...)
