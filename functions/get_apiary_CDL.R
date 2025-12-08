get_apiary_CDL <- function(gc_apiaries, gc_rad, gc_yr, gc_dTol,
                         gc_shp = F, gc_tiff = F, gc_CDL = T, gc_log = T) {

  library(tidyverse)
  library(sf)         # Work with spatial data formats
  library(units)      # Support units in R vectors, for set_units()
  library(CropScapeR) # Retrieve CDL Data

  # Loop through each apiary. ----------------------------------------------------
  # Create a blank dataframe to hold CDL observations for each apiary
  CDL.apiaries <- data.frame(
    Apiary = factor(),
    Year   = factor(),
    Rad.km = factor(),
    Crop   = factor(),
    Pct    = double()
  )

  for (i in 1:nrow(gc_apiaries)) {

    # Set the origin of the circle
    lat    <- gc_apiaries$Lat[i]
    lon    <- gc_apiaries$Lon[i]

    # Define points of the circle. See https://gis.stackexchange.com/a/322432
    circle <- st_point(c(lon, lat), dim = "XY") %>%

      # GMaps is WGS84, or ESPG 4326.
      st_sfc(crs = 4326) %>%

      # Create a circle. This will be pixellated since it traces a circle
      # following "s2 cell boundaries".
      # See https://r-spatial.github.io/sf/articles/sf7.html#buffers
      st_buffer(set_units(gc_rad, km), ) %>%

      # Simplify to speed up computation.
      st_simplify(dTolerance = gc_dTol)  %>%

      # Convert buffer to coordinates to feed into CDLStat. The CDL uses the Albers
      # projection. EPSG value from https://gis.stackexchange.com/a/59750
      st_transform(crs = 5070)

    # Optionally: Save the buffer as an ESRI shapefile. You can check the
    # result by uploading the shapefile to https://mapshaper.org/
    if (gc_shp) {

      setwd("./functions") # Move to functions folder

      shp <- paste0("buffer-", substr(gc_apiaries$Apiary[i], 1,4),
        "-", gc_yr, "-", gc_rad, "km-", gc_dTol, "dTol")

      st_write(circle, paste0(shp, ".shp"), append = F)

      # Get files we just created, zip them, then remove now-zipped files.
      shpfiles <- paste0(grep(paste0("^", shp, ".*"), dir(getwd()), value=T))
      zip(zipfile = paste0(shp, ".zip"), files=shpfiles)
      file.remove(shpfiles)

      setwd("../") # Return to main directory
    }

    # To get data for a precise buffer and not a bounding box, GetCDL() aoi
    # needs to be a vector of coordinates or a URL to a zipped shapefile.
    # Uploading a .zip to the internet is out of the scope here so we'll convert
    # coordinates to a vector. IMPORTANT: the coordinates have a max length
    polygon <- circle[[1]][[1]]
    nCoords <- length(polygon)/2

    # The CropScape package takes a vector of buffer coordinates.
    # Here I convert the dataframe of coordinates to a vector.
    coords <- NULL
    for (n in 1:nCoords) {
      coords <- c(coords, as.vector(polygon[n,]))
    }

    # Check if coordinates are too long. Otherwise, throw error.
    if(nchar(paste0(coords, collapse = ',')) > 1955) {
      stop("The generated buffer contains too many coordinates to submit to the CDL as a string of points. Try increasing dTolerance.")}

    # Get the CDL data for current apiary.
    # If we wanted to save a TIFF:
    if(gc_tiff) {
      CDLData <- GetCDLData(aoi = coords, year = gc_yr, type = 'ps',
                            save_path = paste0(gc_apiaries$Apiary[i], "-CDL.tif"))
    }

    if(gc_CDL) {
      # Get acerage:
      CDLStat <- GetCDLStat(aoi = coords, year = gc_yr, type = 'ps')
      CDLStat <- CDLStat %>%
        rename("Crop" = "Category") %>%  mutate(Crop = as.factor(Crop))

      # Convert acreage to Percentages (proportions were small numbers)
      # The pull() function separately computes acreage of all landscapes in buffer
      CDLStat <- CDLStat %>% mutate(
        # We need the farm name AND apiary # to join onto Apiary df
        Apiary = as.factor(gc_apiaries$Apiary[i]),
        Year   = as.factor(gc_yr),
        Rad.km = as.factor(gc_rad),
        Pct = Acreage /
          pull(summarize(CDLStat, totalBuffAc = sum(Acreage))) * 100) %>%
        select(Apiary, Year, Rad.km, Crop, Pct)

      # Add the current observation onto the data for all apiaries.
      CDL.apiaries <- bind_rows(CDL.apiaries, CDLStat)

      # Log to the console current status of the function
      if(gc_log) {
        print(paste("Finished apiary", i, "of", nrow(gc_apiaries)))
      } # End: if(gc_log)
    }   # End: if(gc_cdl)
  }     # End: for()

  # Return the dataframe
  return(CDL.apiaries)
}