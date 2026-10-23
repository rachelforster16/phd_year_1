# Build one long-format CSV from every Sharkipedia trend dataset.
#
# Setup:   install.packages(c("rvest", "httr", "dplyr", "purrr", "readr", "stringr", "tibble"))
# Output:  sharkipedia_trends_combined.csv (one row per dataset per year)
#
# Raw HTML is cached in ./sharkipedia_cache, so you can stop and re-run
# without re-downloading anything.

library(rvest)
library(httr)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)

BASE       <- "https://www.sharkipedia.org"
LAST_PAGE  <- 44      # from the pagination on /trends
DELAY      <- 0.7     # seconds between uncached requests (be polite)
CACHE_DIR  <- "sharkipedia_cache"
OUT_FILE   <- "sharkipedia_trends_combined.csv"
UA         <- "sharkipedia-research-scraper (personal research use)"

dir.create(CACHE_DIR, showWarnings = FALSE)

# output column -> key as shown in the Details page tables
META_FIELDS <- c(
  species                 = "Species",
  family                  = "Family",
  order                   = "Order",
  subclass                 = "Subclass",
  unit                    = "Unit",
  source_observation      = "Source observation",
  general_data_type       = "General data type",
  unit_freeform           = "Unit freeform",
  unit_transformation     = "Unit transformation",
  latitude                = "Latitude",
  longitude               = "Longitude",
  sampling_method_general = "Sampling Method general",
  comments                = "Comments",
  location_name           = "Name",
  ocean                   = "Ocean",
  author_year             = "AuthorYear"
)

get_page <- function(url, cache_name) {
  path <- file.path(CACHE_DIR, cache_name)
  if (file.exists(path)) return(read_html(path))
  for (attempt in 1:4) {
    r <- tryCatch(GET(url, user_agent(UA), timeout(30)), error = function(e) NULL)
    if (!is.null(r) && status_code(r) == 200) {
      writeLines(content(r, "text", encoding = "UTF-8"), path, useBytes = TRUE)
      Sys.sleep(DELAY)
      return(read_html(path))
    }
    message("  ", url, " failed, retrying")
    Sys.sleep(2 * attempt)
  }
  stop("Failed to fetch ", url)
}

collect_ids <- function() {
  ids <- character()
  for (page in seq_len(LAST_PAGE)) {
    doc   <- get_page(sprintf("%s/trends?page=%d", BASE, page), sprintf("list_%d.html", page))
    hrefs <- html_attr(html_elements(doc, "a"), "href")
    found <- str_match(hrefs, "^(?:https://www\\.sharkipedia\\.org)?/trends/(\\d+)$")[, 2]
    found <- setdiff(unique(found[!is.na(found)]), ids)
    ids   <- c(ids, found)
    message("list page ", page, ": ", length(found), " new datasets")
  }
  ids
}

cells_of <- function(tr) html_text2(html_elements(tr, "th,td"))

parse_detail <- function(doc) {
  meta <- list()
  rows <- tibble(year = character(), measurement = character())
  
  for (tbl in html_elements(doc, "table")) {
    trs <- html_elements(tbl, "tr")
    if (length(trs) == 0) next
    header <- cells_of(trs[[1]])
    
    if (length(header) >= 2 && tolower(header[1]) == "year" && tolower(header[2]) == "value") {
      for (tr in trs[-1]) {
        cells <- cells_of(tr)
        if (length(cells) >= 2) {
          rows <- add_row(rows, year = cells[1], measurement = cells[2])
        }
      }
    } else {
      # key/value metadata table (first row may itself be a key/value pair)
      for (tr in trs) {
        cells <- cells_of(tr)
        if (length(cells) == 2 && nzchar(cells[1]) && is.null(meta[[cells[1]]])) {
          meta[[cells[1]]] <- cells[2]
        }
      }
    }
  }
  list(meta = meta, rows = rows)
}

scrape_dataset <- function(id, i, n) {
  doc <- get_page(sprintf("%s/trends/%s", BASE, id), sprintf("detail_%s.html", id))
  p   <- parse_detail(doc)
  
  if (nrow(p$rows) == 0) {
    message("[", i, "/", n, "] ", id, ": no Year/Value rows found")
    return(NULL)
  }
  
  meta_vals <- map(META_FIELDS, ~ p$meta[[.x]] %||% "")
  out <- bind_cols(
    tibble(dataset_id = id),
    p$rows,
    as_tibble(meta_vals)
  )
  if (i %% 25 == 0) message("[", i, "/", n, "] done")
  out
}

ids <- collect_ids()
message("\n", length(ids), " datasets found. Fetching details...\n")

result <- imap_dfr(ids, ~ scrape_dataset(.x, .y, length(ids)))

# column order: year + measurement up front, then metadata
result <- result |>
  select(dataset_id, year, measurement, everything())

write_excel_csv(result, OUT_FILE)  # UTF-8 with BOM so Excel reads accents correctly
message("\nWrote ", nrow(result), " rows to ", OUT_FILE)