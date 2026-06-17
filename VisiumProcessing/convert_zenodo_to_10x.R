#!/usr/bin/env Rscript
# ==============================================================================
# zenodo14204217  ST-TNBC  →  10x Visium 标准格式转换（修订版）
# 数据源：
#   - 计数矩阵：rawCountsMatrices/*.tsv.bz2 (spots × genes, 需转置)
#   - 坐标信息：byArray/{slide}/{subslide}/spot_data-all-*.tsv 或 allSpots.RDS
#   - 图像文件：byArray/.../small.jpg 或任意 .jpg/.png
#   - 样本映射：Clinical/ids.RDS  (id 列对应 ST_TNBC_ID)
#   - 临床信息：Clinical/Clinical.xlsx
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(DropletUtils)
  library(jsonlite)
  library(readxl)
  library(data.table)
})

# ──────────────────────────────────────────────────────────────────────────────
# 命令行参数
# ──────────────────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--raw_dir",   type="character", default=NULL,
              help="原始数据根目录（包含 byArray/, rawCountsMatrices/, Clinical/）"),
  make_option("--out_dir",   type="character", default=NULL,
              help="输出根目录"),
  make_option("--sample_id", type="integer",   default=NULL,
              help="仅处理单个 ST_TNBC_ID（示例：--sample_id 1）"),
  make_option("--ncores",    type="integer",   default=1,
              help="并行核心数（仅当 --sample_id 未指定时生效）"),
  make_option("--overwrite", action="store_true", default=FALSE,
              help="强制覆盖已有输出"),
  make_option("--dry_run",   action="store_true", default=FALSE,
              help="只打印计划，不写文件")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$raw_dir)) stop("--raw_dir 是必填参数")
if (is.null(opt$out_dir)) stop("--out_dir 是必填参数")

# ──────────────────────────────────────────────────────────────────────────────
# 辅助函数
# ──────────────────────────────────────────────────────────────────────────────
say <- function(level, ...) {
  icons <- c(INFO="·", WARN="⚠", ERROR="✗", OK="✓", STEP="▶")
  msg <- paste0("[", level, "] ", icons[level], " ", paste0(...))
  cat(msg, "\n")
  flush.console()
}
info  <- function(...) say("INFO",  ...)
warn  <- function(...) say("WARN",  ...)
err   <- function(...) say("ERROR", ...)
ok    <- function(...) say("OK",    ...)
step  <- function(...) say("STEP",  ...)

# 1) 读取计数矩阵（tsv.bz2, spots × genes → genes × spots）
read_count_matrix <- function(count_file) {
  info("  读取计数矩阵: ", basename(count_file))
  dt <- fread(cmd = paste("bzcat", shQuote(count_file)), header = TRUE, sep = "\t", data.table = FALSE)
  spot_names <- dt[[1]]
  dt[[1]] <- NULL
  # 基因名去掉 Ensembl 版本号
  gene_names <- sub("\\.\\d+$", "", colnames(dt))
  mat_dense <- as.matrix(dt)
  rownames(mat_dense) <- spot_names
  colnames(mat_dense) <- gene_names
  mat <- t(as(mat_dense, "dgCMatrix"))   # genes × spots
  info("    基因数: ", nrow(mat), "  spots: ", ncol(mat))
  return(mat)
}

# 2) 读取空间坐标（优先 spot_data-all-*.tsv，后备 allSpots.RDS）
read_spot_coords <- function(array_dir) {
  tsv_files <- list.files(array_dir, pattern = "^spot_data-all-.*\\.tsv$", full.names = TRUE)
  if (length(tsv_files) > 0) {
    info("  读取坐标: ", basename(tsv_files[1]))
    df <- read.table(tsv_files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # 假设至少包含 x, y, pixel_x, pixel_y 列
    coords <- data.frame(
      barcode            = paste0(df$x, "x", df$y),
      in_tissue          = 1L,
      array_row          = as.integer(round(df$y)),
      array_col          = as.integer(round(df$x)),
      pxl_row_in_fullres = as.integer(round(df$pixel_y)),
      pxl_col_in_fullres = as.integer(round(df$pixel_x)),
      stringsAsFactors   = FALSE
    )
    return(coords)
  }
  # 后备 allSpots.RDS
  rds_path <- file.path(array_dir, "allSpots.RDS")
  if (file.exists(rds_path)) {
    warn("  未找到 spot_data TSV，使用 allSpots.RDS")
    spots <- readRDS(rds_path)
    df <- as.data.frame(spots)
    # 自动匹配列名（不区分大小写）
    find_col <- function(patterns, df) {
      for (p in patterns) {
        idx <- grep(p, names(df), ignore.case = TRUE)
        if (length(idx) == 1) return(names(df)[idx])
      }
      return(NULL)
    }
    x_col <- find_col(c("^x$", "^array_col$", "^col$"), df)
    y_col <- find_col(c("^y$", "^array_row$", "^row$"), df)
    px_col <- find_col(c("pixel_x", "pxl_col"), df)
    py_col <- find_col(c("pixel_y", "pxl_row"), df)
    if (is.null(x_col) || is.null(y_col))
      stop("allSpots.RDS 中无坐标列，实际列名: ", paste(names(df), collapse=", "))
    coords <- data.frame(
      barcode            = paste0(df[[x_col]], "x", df[[y_col]]),
      in_tissue          = 1L,
      array_row          = as.integer(round(as.numeric(df[[y_col]]))),
      array_col          = as.integer(round(as.numeric(df[[x_col]]))),
      pxl_row_in_fullres = if (!is.null(py_col)) as.integer(round(as.numeric(df[[py_col]]))) else as.integer(round(as.numeric(df[[y_col]]) * 10)),
      pxl_col_in_fullres = if (!is.null(px_col)) as.integer(round(as.numeric(df[[px_col]]))) else as.integer(round(as.numeric(df[[x_col]]) * 10)),
      stringsAsFactors   = FALSE
    )
    return(coords)
  }
  stop("坐标文件缺失：既无 spot_data-all-*.tsv 也无 allSpots.RDS")
}

# 3) 对齐计数矩阵与坐标（barcode 匹配）
align_matrix_coords <- function(mat, coords) {
  mat_barcodes <- colnames(mat)
  coord_barcodes <- coords$barcode
  common <- intersect(mat_barcodes, coord_barcodes)
  if (length(common) == 0) {
    # 尝试去除小数（"1.0x1.0" → "1x1"）
    clean <- function(v) sub("\\.0$", "", sub("\\.0x", "x", v))
    mat_clean <- clean(mat_barcodes)
    coord_clean <- clean(coord_barcodes)
    common2 <- intersect(mat_clean, coord_clean)
    if (length(common2) > 0) {
      info("  去小数后匹配成功，交集: ", length(common2))
      mat_idx <- which(mat_clean %in% common2)
      coord_idx <- which(coord_clean %in% common2)
      mat <- mat[, mat_idx, drop=FALSE]
      coords <- coords[coord_idx, ]
      ord <- match(mat_clean[mat_idx], coord_clean[coord_idx])
      coords <- coords[ord, ]
      colnames(mat) <- coords$barcode
    } else {
      warn("  barcode 完全不匹配，按行数顺序对齐")
      n <- min(ncol(mat), nrow(coords))
      mat <- mat[, 1:n, drop=FALSE]
      coords <- coords[1:n, ]
      colnames(mat) <- coords$barcode
    }
  } else {
    mat <- mat[, common, drop=FALSE]
    coords <- coords[match(common, coord_barcodes), ]
  }
  list(mat = mat, coords = coords)
}

# 4) 图像缩放（依赖 magick，若缺失则复制原图）
process_images <- function(array_dir, spatial_out, coords) {
  hires_out <- file.path(spatial_out, "tissue_hires_image.png")
  lowres_out <- file.path(spatial_out, "tissue_lowres_image.png")
  # 查找图像文件（优先 small.jpg）
  img_files <- list.files(array_dir, pattern = "\\.(jpg|jpeg|png)$", ignore.case = TRUE, full.names = TRUE)
  small_jpg <- img_files[grepl("small\\.jpg$", img_files, ignore.case = TRUE)]
  src_img <- if (length(small_jpg) > 0) small_jpg[1] else if (length(img_files) > 0) img_files[1] else NULL
  
  if (!is.null(src_img) && requireNamespace("magick", quietly = TRUE)) {
    img <- magick::image_read(src_img)
    info <- magick::image_info(img)
    hires_scale <- min(1, 2000 / info$width)
    lowres_scale <- min(1, 600 / info$width)
    magick::image_write(magick::image_scale(img, "2000x>"), hires_out, format = "png")
    magick::image_write(magick::image_scale(img, "600x>"),  lowres_out, format = "png")
    hires_scale <- magick::image_info(magick::image_read(hires_out))$width / info$width
    lowres_scale <- magick::image_info(magick::image_read(lowres_out))$width / info$width
  } else if (!is.null(src_img)) {
    warn("  magick 未安装，直接复制图像（不缩放）")
    file.copy(src_img, hires_out, overwrite = TRUE)
    file.copy(src_img, lowres_out, overwrite = TRUE)
    # 估算缩放因子（全分辨率尺寸从坐标推算）
    full_w <- max(coords$pxl_col_in_fullres) + 200
    hires_scale <- 2000 / full_w
    lowres_scale <- 600 / full_w
  } else {
    warn("  未找到图像，生成空白图")
    full_w <- max(coords$pxl_col_in_fullres, na.rm = TRUE) + 200
    full_h <- max(coords$pxl_row_in_fullres, na.rm = TRUE) + 200
    blank <- function(path, w, h) {
      png(path, width = w, height = h, bg = "white")
      par(mar = c(0,0,0,0)); plot.new(); dev.off()
    }
    blank(hires_out, 2000, as.integer(2000 * full_h / full_w))
    blank(lowres_out, 600, as.integer(600 * full_h / full_w))
    hires_scale <- 2000 / full_w
    lowres_scale <- 600 / full_w
  }
  list(hires_scale = hires_scale, lowres_scale = lowres_scale)
}

# 5) 写 scalefactors_json.json
write_scalefactors <- function(path, hires_scale, lowres_scale, coords) {
  # 估算 spot 直径（取相邻 spot 中位距离 * 0.65）
  spot_diam <- tryCatch({
    n <- min(200, nrow(coords))
    sub <- coords[1:n, ]
    d <- as.matrix(dist(cbind(sub$pxl_col_in_fullres, sub$pxl_row_in_fullres)))
    diag(d) <- NA
    median(apply(d, 1, min, na.rm = TRUE), na.rm = TRUE) * 0.65
  }, error = function(e) 89.0)
  sf <- list(
    spot_diameter_fullres = round(spot_diam, 4),
    tissue_hires_scalef   = round(hires_scale, 6),
    fiducial_diameter_fullres = round(spot_diam * 2.5, 4),
    tissue_lowres_scalef  = round(lowres_scale, 6)
  )
  writeLines(jsonlite::toJSON(sf, auto_unbox = TRUE, pretty = TRUE), path)
}

# 6) 写 H5 文件
write_h5 <- function(path, mat) {
  mat <- as(mat, "dgCMatrix")
  DropletUtils::write10xCounts(
    path = path, x = mat, barcodes = colnames(mat),
    gene.id = rownames(mat), gene.symbol = rownames(mat),
    gene.type = rep("Gene Expression", nrow(mat)),
    overwrite = TRUE, type = "HDF5", genome = "GRCh38",
    version = "3", chemistry = "Spatial 3' v1"
  )
}

# ──────────────────────────────────────────────────────────────────────────────
# 主处理函数（单个样本）
# ──────────────────────────────────────────────────────────────────────────────
process_one_sample <- function(row, raw_dir, out_dir, overwrite, dry_run, clin_df) {
  st_id <- row$id
  slide <- row$slide
  subslide <- row$subslide
  count_file <- row$count
  has_annot <- row$hasAnnot
  
  array_dir <- file.path(raw_dir, "byArray", slide, subslide)
  count_path <- file.path(raw_dir, "rawCountsMatrices", count_file)
  sample_out <- file.path(out_dir, paste0("zenodo14204217_", st_id))
  spatial_out <- file.path(sample_out, "spatial")
  h5_out <- file.path(sample_out, "filtered_feature_bc_matrix.h5")
  pos_out <- file.path(spatial_out, "tissue_positions_list.csv")
  sf_out <- file.path(spatial_out, "scalefactors_json.json")
  clin_out <- file.path(sample_out, "clinical_metadata.csv")
  
  step(sprintf("── ST_TNBC_ID=%d [%s/%s] hasAnnot=%s ──", st_id, slide, subslide, has_annot))
  if (dry_run) {
    info("[DRY RUN] 将处理 ", count_path)
    return(list(status = "dry_run", id = st_id))
  }
  if (!overwrite && file.exists(h5_out) && file.exists(pos_out)) {
    info("输出已存在，跳过（使用 --overwrite 强制重跑）")
    return(list(status = "skipped", id = st_id))
  }
  if (!file.exists(count_path)) {
    err("计数文件不存在: ", count_path)
    return(list(status = "error", id = st_id, msg = "missing count file"))
  }
  if (!dir.exists(array_dir)) {
    err("阵列目录不存在: ", array_dir)
    return(list(status = "error", id = st_id, msg = "missing array dir"))
  }
  dir.create(spatial_out, recursive = TRUE, showWarnings = FALSE)
  
  # 1) 计数矩阵
  mat <- tryCatch(read_count_matrix(count_path), error = function(e) NULL)
  if (is.null(mat)) return(list(status = "error", id = st_id, msg = "read count failed"))
  
  # 2) 坐标
  coords <- tryCatch(read_spot_coords(array_dir), error = function(e) NULL)
  if (is.null(coords)) return(list(status = "error", id = st_id, msg = "read coords failed"))
  
  # 3) 对齐
  aligned <- tryCatch(align_matrix_coords(mat, coords), error = function(e) NULL)
  if (is.null(aligned)) return(list(status = "error", id = st_id, msg = "alignment failed"))
  mat <- aligned$mat
  coords <- aligned$coords
  if (ncol(mat) == 0 || nrow(coords) == 0) {
    err("对齐后无有效 spot")
    return(list(status = "error", id = st_id, msg = "no common spots"))
  }
  
  # 4) 图像处理
  img_info <- tryCatch(process_images(array_dir, spatial_out, coords), error = function(e) {
    warn("图像处理失败，使用默认缩放: ", e$message)
    list(hires_scale = 0.15, lowres_scale = 0.05)
  })
  
  # 5) 写出文件
  write.csv(coords, pos_out, row.names = FALSE, quote = FALSE)
  write_scalefactors(sf_out, img_info$hires_scale, img_info$lowres_scale, coords)
  write_h5(h5_out, mat)
  ok("  H5 文件写入成功 (", nrow(mat), " genes × ", ncol(mat), " spots)")
  
  # 临床元数据（可选）
  if (!is.null(clin_df) && st_id %in% clin_df$ST_TNBC_ID) {
    clin_sub <- clin_df[clin_df$ST_TNBC_ID == st_id, , drop = FALSE]
    clin_sub$selected_slide <- slide
    clin_sub$selected_subslide <- subslide
    clin_sub$selected_count_file <- count_file
    clin_sub$hasAnnot_selected <- has_annot
    write.csv(clin_sub, clin_out, row.names = FALSE)
    ok("  临床元数据已保存")
  }
  ok(sprintf("ST_TNBC_ID=%d 完成 → %s", st_id, sample_out))
  return(list(status = "success", id = st_id))
}

# ──────────────────────────────────────────────────────────────────────────────
# 主流程
# ──────────────────────────────────────────────────────────────────────────────
step("====== 读取 ids.RDS ======")
ids_path <- file.path(opt$raw_dir, "Clinical", "ids.RDS")
if (!file.exists(ids_path)) stop("找不到 ids.RDS: ", ids_path)
ids <- readRDS(ids_path)
ids$id <- as.integer(ids$id)          # 确保 ST_TNBC_ID 为整数
ids$hasAnnot[is.na(ids$hasAnnot)] <- FALSE
info("ids.RDS 行数: ", nrow(ids), ", ST_TNBC_ID 范围: ", min(ids$id, na.rm=TRUE), "~", max(ids$id, na.rm=TRUE))

# 为每个 ID 选择最佳阵列（优先 hasAnnot==TRUE）
best_list <- split(ids, ids$id)
best_rows <- lapply(best_list, function(sub) {
  ann <- sub[sub$hasAnnot == TRUE, ]
  if (nrow(ann) > 0) ann[1, ] else sub[1, ]
})
best <- do.call(rbind, best_rows)
rownames(best) <- NULL
info("共 ", nrow(best), " 个 ST_TNBC_ID 有对应阵列")

# 如果指定了单个样本，则只处理该样本
if (!is.null(opt$sample_id)) {
  best <- best[best$id == opt$sample_id, ]
  if (nrow(best) == 0) stop("ST_TNBC_ID ", opt$sample_id, " 不在 ids.RDS 中")
  info("仅处理样本: ", opt$sample_id)
}

# 读取临床信息
clin_df <- NULL
clin_path <- file.path(opt$raw_dir, "Clinical", "Clinical.xlsx")
if (file.exists(clin_path)) {
  clin_df <- as.data.frame(read_excel(clin_path))
  clin_df$ST_TNBC_ID <- as.integer(clin_df$ST_TNBC_ID)
  info("临床数据: ", nrow(clin_df), " 行")
} else {
  warn("Clinical.xlsx 不存在，跳过临床元数据输出")
}

# 执行转换
step("====== 开始转换 ======")
if (!opt$dry_run) dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

results <- if (opt$ncores > 1 && is.null(opt$sample_id)) {
  parallel::mclapply(seq_len(nrow(best)), function(i) {
    process_one_sample(best[i, ], opt$raw_dir, opt$out_dir, opt$overwrite, opt$dry_run, clin_df)
  }, mc.cores = opt$ncores)
} else {
  lapply(seq_len(nrow(best)), function(i) {
    process_one_sample(best[i, ], opt$raw_dir, opt$out_dir, opt$overwrite, opt$dry_run, clin_df)
  })
}

# 汇总
step("====== 汇总 ======")
statuses <- sapply(results, `[[`, "status")
ok("成功: ", sum(statuses == "success"))
info("跳过: ", sum(statuses == "skipped"))
info("空跑: ", sum(statuses == "dry_run"))
if (any(statuses == "error")) {
  err("失败: ", sum(statuses == "error"))
  for (r in results[statuses == "error"]) {
    err("  ST_TNBC_ID=", r$id, "  ", r$msg)
  }
  quit(status = 1)
}
step("====== 全部完成 ======")