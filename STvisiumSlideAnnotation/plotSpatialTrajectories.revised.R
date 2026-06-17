#' @title Plot spatial trajectories
#'
#' @description Displays the spatial course of spatial trajectory that was
#' drawn with \code{createSpatialTrajectories()} on a surface plot.
#' Increase the transparency via argument \code{pt_alpha} to highlight
#' the trajectory's course.
#'
#' @param pt_alpha2 Numeric value. Specifies the transparency of the spots
#' that fall into the trajectory's reach.
#' @inherit argument_dummy params
#' @inherit check_color_to params
#' @inherit check_display params
#' @inherit check_method params
#' @inherit check_pt params
#' @inherit check_sample params
#' @inherit check_smooth params
#' @inherit check_trajectory params
#' @inherit check_uniform_genes params
#'
#' @param sgmt_size The size of the segment arrrow specified as a numeric value.
#'
#' @inherit ggplot_family return
#'
#' @export
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF269T_diet
#'
#' plotSpatialTrajectories(object, color_by = "histology", pt_clrp = "npg", ids = "horizontal_mid")
#'
plotSpatialTrajectories <- function(object,
                                    ids = NULL,
                                    color_by = NULL,
                                    alpha_by = NULL,
                                    method_gs = NULL,
                                    smooth = NULL,
                                    smooth_span = NULL,
                                    width = NULL,
                                    pt_size = NULL,
                                    pt_alpha = 0.5,
                                    pt_alpha2 = 0.9,
                                    pt_clr = NULL,
                                    pt_clrp = NULL,
                                    clrp_adjust = NULL,
                                    pt_clrsp = NULL,
                                    sgmt_clr = NULL,
                                    sgmt_size = NULL,
                                    display_facets = FALSE,
                                    display_image = NULL,
                                    display_title = NULL,
                                    uniform_genes = NULL,
                                    arrow = ggplot2::arrow(length = ggplot2::unit(x = 0.125, "inches")),
                                    nrow = NULL,
                                    ncol = NULL,
                                    xrange = getCoordsRange(object)[["x"]],
                                    yrange = getCoordsRange(object)[["y"]],
                                    verbose = NULL,
                                    ...){

  check_object(object)
  hlpr_assign_arguments(object)

  if(base::is.null(ids)){ ids <- getSpatialTrajectoryIds(object) }

  df <-
    purrr::map_df(
      .x = ids,
      .f = function(id){

        background_df <-
          getCoordsDfST(object, id = id, width = width) %>%
          dplyr::mutate(ids = {{id}})

        if(base::is.character(color_by)){

          background_df <-
            hlpr_join_with_color_by(
              object = object,
              df = background_df,
              color_by = color_by,
              smooth = smooth,
              smooth_span = smooth_span,
              method_gs = method_gs,
              verbose = FALSE
            )

        }

        return(background_df)

      }
    )

  if(base::isTRUE(display_image)){

    base_plot <-
      ggplot2::ggplot() +
      ggpLayerImage(object = object)

  } else {

    base_plot <-
      ggplot2::ggplot(
        data = df,
        mapping = ggplot2::aes_string(x = "x", y = "y")
      )

  }


  if(base::length(ids) == 1){ display_facets <- FALSE }

  if(base::isTRUE(display_facets)){

    facet_add_on <- ggplot2::facet_wrap(facets = . ~ ids, nrow = nrow, ncol = ncol)

  } else {

    facet_add_on <- NULL

  }

  # scale spot size to plot frame
  mx_range <- base::max(c(base::diff(xrange), base::diff(yrange)))

  if(containsImage(object)){

    mx_dims <- base::max(getImageDims(object))

  } else {

    mx_dims <-
      purrr::map_dbl(df[,c("x", "y")], .f = base::max) %>%
      base::max()

  }

  pt_size <- (mx_dims/mx_range)*pt_size

  params <- SPATA2:::adjust_ggplot_params(params = list(color = pt_clr, size = pt_size))

  base_plot +
    geom_point_fixed(
      params,
      data = df,
      mapping = ggplot2::aes(x = x, y = y, alpha = rel_loc, color = .data[[color_by]])
    ) +
    ggpLayerSpatialTrajectories(
      object = object,
      ids = ids,
      arrow = arrow,
      size = sgmt_size,
      color = sgmt_clr
    ) +
    facet_add_on +
    ggplot2::theme_void() +
    ggplot2::scale_alpha_manual(values = c("inside" = pt_alpha2, "outside" = pt_alpha), guide = "none") +
    scale_color_add_on(
      aes = "color",
      variable = SPATA2:::pull_var(df, color_by),
      clrp = pt_clrp,
      clrsp = pt_clrsp,
      clrp.adjust = clrp_adjust,
      ...
    ) +
    ggplot2::coord_equal(xlim = xrange, ylim = yrange)

}