read_pesticides <- function(
    frp_ss, frp_tab, frp_data, frp_lims, frp_pestName, frp_noPivot = 1:4, frp_limUnits) {

  library(tidyverse)  # Easier data manipulation
  library(readxl)     # Import Excel files directly
  library(janitor)    # Provides renaming functions
  library(data.table) # For transpose() function

  # Build Excel ranges
  frp_data     <- paste(frp_tab, frp_data,     sep = "!")
  frp_lims     <- paste(frp_tab, frp_lims,     sep = "!")
  frp_pestName <- paste(frp_tab, frp_pestName, sep = "!")

  # Import the raw pesticide data.
  data.Raw <- read_excel(frp_ss, range = frp_data) %>%
    rename_with(~ make_clean_names(.x, case="none", sep_out = "."), all_of(frp_noPivot)) %>%
    rename(Sample.grams = Mass.g)

  # Import Limits (LOD, LOQ, etc.)
  meta.Limits <- read_excel(
      frp_ss, range = frp_lims,
      col_names = names(read_excel(frp_ss, range = frp_pestName))) %>%
    rename(Limit = 1) %>%
    transpose(keep.names = "Cmpd", make.names = "Limit")

  if(frp_limUnits == "ng") {
    # If limits in raw data are reported in ng, we want to convert to sample-specific limits.
    # We'll retrieve average sample mass, then divide limits by that number.
    av.mass <- data.Raw %>% pull(Sample.grams) %>% mean() %>% round(digits = 4)

    meta.Limits <- meta.Limits %>%
      # Convert limits to ppb
      mutate(across(!Cmpd, ~ .x / av.mass)) %>%
      rename_with(~ sub(" .*", ".ppb", .x))

  } else if (frp_limUnits == "ppb") {
    # In case limits are already reported as sample-specific ppb limits.
    meta.Limits <- meta.Limits %>%
      rename_with(~ sub(" .*", ".ppb", .x))
  }

  meta.Limits <- meta.Limits %>% mutate(across(!Cmpd, ~ round(.x, digits = 2)))

  # if(!any(grepl("ULOQ", names(meta.Limits)))) {
  #   # If dataset reports ULOL vs ULOQ.
  #   # Mutate() functions below won't work unless we have an ULOQ value.
  #   # So I set this to NA just in case it's different from ULOQ
  #   ULOQ.col <- paste0("ULOQ.", frp_limUnits)
  #   meta.Limits[[ULOQ.col]] <- NA
  # }

  # Pivot data to longer format -----------------------------------------------
  data_Long <- data.Raw %>%
    # Convert all pesticide columns to character.(Any column)
    # We'll get an error when we pivot_longer if we don't.
    mutate(across(-all_of(frp_noPivot), ~as.character(.x))) %>%

    # Convert columns to rows
    pivot_longer(
      cols      = -all_of(frp_noPivot),
      names_to  = "Cmpd",
      values_to = "ppb") %>%

    # Add on LOQ values for each pesticide
    left_join(meta.Limits, by = "Cmpd") %>%

    # Convert every detection to a number. Pesticides not found should be 0ppb
    mutate(ppb = ifelse(ppb %in% c("n.d.", "N/F"), 0, ppb)) ## %>%

    # <LOQs or >ULOQs get replaced with LOQ or ULOQ value, respectively
    # mutate(ppb = case_match(ppb, "<LOQ"   ~ as.character(LOD.ppb),
    #                              ">ULOQ"  ~ as.character(ULOQ.ppb),
    #                              .default = as.character(ppb)))

  # Now we can drop the placeholder ULOQ column
  # if(exists("ULOQ.col")) data_Long <- data_Long %>% select(-starts_with("ULOQ"))

  # Return dataframes (in a list) ---------------------------------------------
  return(data_Long)
}
