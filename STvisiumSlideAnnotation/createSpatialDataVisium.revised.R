createSpatialDataVisium <- function(dir,
                                    sample,
                                    img_ref = "lowres",
                                    img_active = "lowres",
                                    resize_images = NULL,
                                    unload = TRUE,
                                    meta = list(),
                                    misc = list(),
                                    verbose = TRUE){

  # get all files in folder and subfolders
  files <- base::list.files(dir, full.names = TRUE, recursive = TRUE)

  # check required image availability
  req_images <- base::unique(c(img_ref, img_active))

  confuns::check_one_of(
    input = req_images,
    against = c("lowres", "hires"),
    ref.input = "required images"
  )

  lowres_path <- base::file.path(dir, "spatial", "tissue_lowres_image.png")
  hires_path <- base::file.path(dir, "spatial", "tissue_hires_image.png")

  if("lowres" %in% req_images){

    if(!lowres_path %in% files){

      stop(glue::glue("'{lowres_path}' is missing."))

    }

  }

  if("hires" %in% req_images){

    if(!hires_path %in% files){

      stop(glue::glue("'{hires_path}' is missing."))

    }

  }

  # load in data

  # check and load tissue positions for different space ranger versions
  v1_coords_path <- base::file.path(dir, "spatial", "tissue_positions_list.csv")
  v2_coords_path <- base::file.path(dir, "spatial", "tissue_positions.csv")

  if(v2_coords_path %in% files){

    space_ranger_version <- 2
    coords_df <- read_coords_visium(dir_coords = v2_coords_path)
    misc[["dirs"]][["coords"]] <- v2_coords_path

  } else if(v1_coords_path %in% files){

    space_ranger_version <- 1
    coords_df <- read_coords_visium(dir_coords = v1_coords_path)
    misc[["dirs"]][["coords"]] <- v1_coords_path

  }

  xmean <- base::mean(coords_df$x_orig, na.rm = TRUE)
  ymean <- base::mean(coords_df$y_orig, na.rm = TRUE)

  method <- spatial_methods[["VisiumSmall"]]

  # load scalefactors
  scale_factors <-
    jsonlite::read_json(path = base::file.path(dir, "spatial", "scalefactors_json.json"))

  # load images
  # reference image
  img_list <- list()

  if("hires" %in% req_images){

    img_list[["hires"]] <-
      createHistoImage(
        dir = hires_path,
        sample = sample,
        img_name ="hires",
        scale_factors =
          list(
            image = scale_factors$tissue_hires_scalef
          ),
        reference = img_ref == "hires",
        verbose = verbose
      )

  }

  if("lowres" %in% req_images){

    img_list[["lowres"]] <-
      createHistoImage(
        dir = lowres_path,
        sample = sample,
        img_name ="lowres",
        scale_factors =
          list(
            image = scale_factors$tissue_lowres_scalef
          ),
        reference = img_ref == "lowres",
        verbose = verbose
      )
  }

  # compute spot size
  spot_size <-
    scale_factors$fiducial_diameter_fullres *
    scale_factors[[stringr::str_c("tissue", img_ref, "scalef", sep = "_")]] /
    base::max(getImageDims(img_list[[img_ref]]))*100 # Visium * 100

  method@method_specifics[["spot_size"]] <- spot_size * 1.1

  # create output
  sp_data <-
    createSpatialData(
      sample = sample,
      hist_img_ref = img_list[[img_ref]],
      hist_imgs = img_list[req_images[req_images != img_ref]],
      active = img_active,
      unload = unload,
      resize_images = resize_images,
      coordinates = coords_df,
      method = method,
      meta = meta,
      misc = misc,
      verbose = verbose
    )

  # compute pixel scale factor to

  sp_data <- computePixelScaleFactor(sp_data, verbose = verbose)

  sp_data <- computeCaptureArea(sp_data)

  return(sp_data)


}