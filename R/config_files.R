##  --initialization_method arg           initialization method
##                                        (domainAverageInitialization,
##                                        pointInitialization,
##                                        wxModelInitialization)
##  --time_zone arg                       time zone (common choices are:
##                                        America/New_York, America/Chicago,
##                                        America/Denver, America/Phoenix,
##                                        America/Los_Angeles, America/Anchorage;
##                                        use 'auto-detect' to try and find the
##                                        time zone for the dem.  All choices are
##                                        listed in date_time_zonespec.csv)
##  --wx_model_type arg                   type of wx model to download
##                                        (UCAR-NAM-CONUS-12-KM,
##                                        UCAR-NAM-ALASKA-11-KM,
##                                        UCAR-NDFD-CONUS-2.5-KM,
##                                        UCAR-RAP-CONUS-13-KM,
##                                        UCAR-GFS-GLOBAL-0.5-DEG,
##                                        NOMADS-GFS-GLOBAL-0.25-DEG,
##                                        NOMADS-HIRES-ARW-ALASKA-5-KM,
##                                        NOMADS-HIRES-NMM-ALASKA-5-KM,
##                                        NOMADS-HIRES-ARW-CONUS-5-KM,
##                                        NOMADS-HIRES-NMM-CONUS-5-KM,
##                                        NOMADS-NAM-ALASKA-11.25-KM,
##                                        NOMADS-NAM-CONUS-12-KM,
##                                        NOMADS-NAM-NORTH-AMERICA-32-KM,
##                                        NOMADS-NAM-NEST-ALASKA-3-KM,
##                                        NOMADS-NAM-NEST-CONUS-3-KM,
##                                        NOMADS-HRRR-ALASKA-3-KM,
##                                        NOMADS-HRRR-CONUS-3-KM,
##                                        NOMADS-HRRR-CONUS-SUBHOURLY-3-KM,
##                                        NOMADS-HRRR-ALASKA-SUBHOURLY-3-KM,
##                                        NOMADS-RAP-CONUS-13-KM,
##                                        NOMADS-RAP-NORTH-AMERICA-32-KM)
##  --forecast_duration arg               forecast duration to download (in
##                                        hours)
##  --forecast_filename arg               path/filename of an already downloaded
##                                        wx forecast file
##  --forecast_time arg                   specific time to run in wx model
##  --match_points arg (=1)               match simulation to points(true, false)
##  --input_speed arg                     input wind speed
##  --input_speed_units arg               units of input wind speed (mps, mph,
##                                        kph, kts)
##  --output_speed_units arg (=mph)       units of output wind speed (mps, mph,
##                                        kph, kts)
##  --input_direction arg                 input wind direction
##  --input_speed_grid arg                path/filename of input raster speed
##                                        file (*.asc)
##  --input_dir_grid arg                  path/filename of input raster dir file
##                                        (*.asc)
##  --uni_air_temp arg                    surface air temperature
##  --air_temp_units arg                  surface air temperature units (K, C, R,
##                                        F)
##  --uni_cloud_cover arg                 cloud cover
##  --cloud_cover_units arg               cloud cover units (fraction, percent,
##                                        canopy_category)
##  --fetch_station arg (=0)              download a station file from an
##                                        internet server (Mesonet API)
##                                        (true/false)
##  --start_year arg                      point initialization: start year for
##                                        simulation
##  --start_month arg                     point initialization: start month for
##                                        simulation
##  --start_day arg                       point initialization: start day for
##                                        simulation
##  --start_hour arg                      point initialization: start hour for
##                                        simulation
##  --start_minute arg                    point initialization: start minute for
##                                        simulation
##  --end_year arg                        point initialization: end year for
##                                        simulation
##  --end_month arg                       point initialization: end month for
##                                        simulation
##  --end_day arg                         point initialization: end day for
##                                        simulation
##  --end_hour arg                        point initialization: end hour for
##                                        simulation
##  --end_minute arg                      point initialization: end minute for
##                                        simulation
##  --number_time_steps arg               point initialization: number of
##                                        timesteps for simulation
##  --fetch_metadata arg (=0)             get weather station metadata for a
##                                        domain
##  --metadata_filename arg               filename for weather station metadata
##  --fetch_type arg                      fetch weather station from bounding box
##                                        (bbox) or by station ID (stid)
##  --fetch_current_station_data arg (=0) fetch the latest weather station data
##                                        (true) or fetch a timeseries (false)
##                                        (true/false)
##  --station_buffer arg (=0)             distance around dem to fetch station
##                                        data
##  --station_buffer_units arg (=km)      Units of distance around DEM
##  --fetch_station_name arg              list of stations IDs to fetch
##  --wx_station_filename arg             path/filename of input wx station file
##  --write_wx_station_kml arg (=0)        point initialization: write a Google
##                                        Earth kml file for the input wx
##                                        stations (true, false)
##  --write_wx_station_csv arg (=0)       point initialization: write a csv of
##                                        the interpolated weather data
##                                        (true,false)
##  --units_input_wind_height arg         units of input wind height (ft, m)
##  --output_wind_height arg              height of output wind speed above the
##                                        vegetation
##  --units_output_wind_height arg        units of output wind height (ft, m)
##  --diurnal_winds arg (=0)              include diurnal winds in simulation
##                                        (true, false)
##  --year arg                            year of simulation
##  --month arg                           month of simulation
##  --day arg                             day of simulation
##  --hour arg                            hour of simulation
##  --minute arg                          minute of simulation
##  --mesh_choice arg                     mesh resolution choice (coarse, medium,
##                                        fine)
##  --mesh_resolution arg                 mesh resolution
##  --units_mesh_resolution arg           mesh resolution units (ft, m)
##  --output_buffer_clipping arg (=0)     percent to clip buffer on output files
##  --write_wx_model_goog_output arg (=0) write a Google Earth kmz output file
##                                        for the raw wx model forecast (true,
##                                        false)
##  --write_goog_output arg (=0)          write a Google Earth kmz output file
##                                        (true, false)
##  --goog_out_resolution arg (=-1)       resolution of Google Earth output file
##                                        (-1 to use mesh resolution)
##  --units_goog_out_resolution arg (=m)  units of Google Earth resolution (ft,
##                                        m)
##  --goog_out_color_scheme arg (=default)
##                                        Sets the color scheme for kml outputs,
##                                        available options:
##                                         default (ROYGB), oranges, blues,
##                                        greens,pinks, magic_beans,
##                                        pink_to_green,ROPGW
##  --goog_out_vector_scaling arg (=0)    Enable Vector Scaling based on Wind
##                                        speed
##  --write_wx_model_shapefile_output arg (=0)
##                                        write a shapefile output file for the
##                                        raw wx model forecast (true, false)
##  --write_shapefile_output arg (=0)     write a shapefile output file (true,
##                                        false)
##  --shape_out_resolution arg (=-1)      resolution of shapefile output file (-1
##                                        to use mesh resolution)
##  --units_shape_out_resolution arg (=m) units of shapefile resolution (ft, m)
##  --write_wx_model_ascii_output arg (=0)
##                                        write ascii fire behavior output files
##                                        for the raw wx model forecast (true,
##                                        false)
##  --write_ascii_output arg (=0)         write ascii fire behavior output files
##                                        (true, false)
##  --ascii_out_resolution arg (=-1)      resolution of ascii fire behavior
##                                        output files (-1 to use mesh
##                                        resolution)
##  --units_ascii_out_resolution arg (=m) units of ascii fire behavior output
##                                        file resolution (ft, m)
##  --write_vtk_output arg (=0)           write VTK output file (true, false)
##  --write_farsite_atm arg (=0)          write a FARSITE atm file (true, false)
##  --write_pdf_output arg (=0)           write PDF output file (true, false)
##  --pdf_out_resolution arg (=-1)        resolution of pdf output file (-1 to
##                                        use mesh resolution)
##  --units_pdf_out_resolution arg (=m)   units of PDF resolution (ft, m)
##  --pdf_linewidth arg (=1)              width of PDF vectors (in pixels)
##  --pdf_basemap arg (=topofire)         background image of the geospatial pdf,
##                                        default is topo map
##  --pdf_height arg                      height of geospatial pdf
##  --pdf_width arg                       width of geospatial pdf
##  --pdf_size arg (=letter)              pre-defined pdf sizes (letter, legal,
##                                        tabloid)
##  --output_path arg                     path to where output files will be
##                                        written
##  --non_neutral_stability arg (=0)      use non-neutral stability (true, false)
##  --alpha_stability arg                 alpha value for atmospheric stability
##  --input_points_file arg               input file containing lat,long,z for
##                                        requested output points (z in m above
##                                        ground)
##  --output_points_file arg              file to write containing output for
##                                        requested points
##  --existing_case_directory arg         path to an existing OpenFOAM case
##                                        directory
##  --momentum_flag arg (=0)              use momentum solver (true, false)
##  --number_of_iterations arg (=300)     number of iterations for momentum
##                                        solver
##  --mesh_count arg                      number of cells in the mesh


#' Configuration options that are common across different WindNinja configuration types
#'
#' @param num_threads numeric: number of threads to use
#' @param elevation string or Raster: a raster layer, or the filename of an elevation file (*.asc, *.lcp, *.tif, *.img)
#' @param vegetation string: dominant type of vegetation ("grass", "brush", "trees")
#' @param mesh_choice string: mesh resolution choice ("coarse", "medium", "fine")
#' @param output_formats character: one or more of "asc", "kmz", "shp"
#' @param ascii_out_resolution numeric: resolution of ascii output files (-1 to use mesh resolution)
#' @param units_ascii_out_resolution string: units of ascii output file resolution ("ft", "m")
#' @param shape_out_resolution numeric: resolution of shapefile output files (-1 to use mesh resolution)
#' @param units_shape_out_resolution string: units of shapefile output file resolution ("ft", "m")
#' @param goog_out_resolution numeric: resolution of Google Earth output file (-1 to use mesh resolution)
#' @param units_goog_out_resolution string: units of Google Earth resolution ("ft", "m")
#'
#' @return A named list of options
#'
#' @seealso \code{\link{wn_config_domain_average}}
#'
#' @examples
#' opts <- wn_config_common(elevation = wn_demo_file("missoula_valley_elevation"),
#'                          ascii_out_resolution = 100)
#'
#' @export
wn_config_common <- function(num_threads = max(1, availableCores()-1),
                             elevation, vegetation = "grass", mesh_choice = "coarse",
                             output_formats = "asc",
                             ascii_out_resolution = -1L, units_ascii_out_resolution = "m",
                             shape_out_resolution = -1L, units_shape_out_resolution = "m",
                             goog_out_resolution = -1L, units_goog_out_resolution = "m") {
    ## elevation
    if (is.character(elevation) && file.exists(elevation)) {
        ## is file
        if (!tolower(fs::path_ext(elevation)) %in% c("asc", "tif", "lcp", "img")) stop("elevation file should be of type asc, lcp, tif, or img")
        elevation_file <- elevation
    } else if (inherits(elevation, "RasterLayer")) {
        ## raster layer, write to file
        elevation_file <- tempfile(fileext = ".tif")
        raster::writeRaster(elevation, elevation_file)
    }
    output_formats <- match_arg(output_formats, c("asc", "kmz", "shp"), several_ok = TRUE)
    write_goog_output <- "kmz" %in% output_formats
    write_shapefile_output <- "shp" %in% output_formats
    write_ascii_output <- "asc" %in% output_formats
    ##write_farsite_atm        = false

    list(num_threads = num_threads,
         elevation_file = elevation_file,
         vegetation = vegetation, mesh_choice = mesh_choice,
         write_ascii_output = write_ascii_output,
         ascii_out_resolution = ascii_out_resolution, units_ascii_out_resolution = units_ascii_out_resolution,
         write_shapefile_output = write_shapefile_output,
         shape_out_resolution = shape_out_resolution, units_shape_out_resolution = units_shape_out_resolution,
         write_goog_output = write_goog_output,
         goog_out_resolution = goog_out_resolution, units_goog_out_resolution = units_goog_out_resolution
         )
}
##mesh_resolution          = 250.0
##units_mesh_resolution    = m

##  --goog_out_color_scheme arg (=default)
##                                        Sets the color scheme for kml outputs,
##                                        available options:
##                                         default (ROYGB), oranges, blues,
##                                        greens,pinks, magic_beans,
##                                        pink_to_green,ROPGW
##  --goog_out_vector_scaling arg (=0)    Enable Vector Scaling based on Wind
##                                        speed
##  --write_wx_model_shapefile_output arg (=0)
##                                        write a shapefile output file for the
##                                        raw wx model forecast (true, false)
##  --write_shapefile_output arg (=0)     write a shapefile output file (true,
##                                        false)



## check entries in opts
## internal
wn_check_config <- function(opts) {
    checked <- as.list(rep(FALSE, length(opts)))
    names(checked) <- names(opts)
    if ("num_threads" %in% names(opts)) {
        num_threads <- opts$num_threads
        assert_that(is.numeric(num_threads), is.scalar(num_threads), num_threads > 0)
        checked$num_threads <- TRUE
    }
    if ("initialization_method" %in% names(opts)) {
        opts$initialization_method <- match.arg(opts$initialization_method, c("domainAverageInitialization"))
        checked$initialization_method <- TRUE
    }
    if ("momentum_flag" %in% names(opts)) {
        momentum_flag <- opts$momentum_flag
        assert_that(is.flag(momentum_flag), !is.na(momentum_flag))
        checked$momentum_flag <- TRUE
    }
    if ("number_of_iterations" %in% names(opts)) {
        number_of_iterations <- opts$number_of_iterations
        assert_that(is.numeric(number_of_iterations), is.scalar(number_of_iterations), number_of_iterations > 0)
        checked$number_of_iterations <- TRUE
    }
    if ("input_speed" %in% names(opts)) {
        input_speed <- opts$input_speed
        assert_that(is.numeric(input_speed), all(input_speed > 0))
        if (length(input_speed) != 1) stop("only a single input_speed is currently supported")
        checked$input_speed <- TRUE
    }
    if ("input_direction" %in% names(opts)) {
        input_direction <- opts$input_direction
        assert_that(is.numeric(input_direction), all(input_direction >= 0 & input_direction <= 360)) ## assuming they must be in this range, TODO check
        if (length(input_direction) != 1) stop("only a single input_direction is currently supported")
        checked$input_direction <- TRUE
    }
    if (all(c("input_speed", "input_direction") %in% names(opts))) {
        assert_that(length(input_direction) == length(input_speed))
    }
    if ("elevation_file" %in% names(opts)) {
        assert_that(file.exists(opts$elevation_file))
        checked$elevation_file <- TRUE
    } else {
        stop("elevation_file entry is missing from configuration")
    }
    if ("vegetation" %in% names(opts)) {
        vegetation <- opts$vegetation
        opts$vegetation <- match_arg(opts$vegetation, c("grass", "trees", "brush"))
        checked$vegetation <- TRUE
    }

    speed_units <- c("mps", "mph", "kph", "kts")
    if ("input_speed_units" %in% names(opts)) {
        opts$input_speed_units <- match_arg(opts$input_speed_units, speed_units)
        checked$input_speed_units <- TRUE
    }
    if ("output_speed_units" %in% names(opts)) {
        opts$output_speed_units <- match_arg(opts$output_speed_units, speed_units)
        checked$output_speed_units <- TRUE
    }

    height_units <- c("m", "ft")
    if ("units_input_wind_height" %in% names(opts)) {
        opts$units_input_wind_height <- match_arg(opts$units_input_wind_height, height_units)
        checked$units_input_wind_height <- TRUE
    }
    if ("units_output_wind_height" %in% names(opts)) {
        opts$units_output_wind_height <- match.arg(opts$units_output_wind_height, height_units)
        checked$units_output_wind_height <- TRUE
    }
    if ("input_wind_height" %in% names(opts)) {
        input_wind_height <- opts$input_wind_height
        assert_that(is.numeric(input_wind_height), is.scalar(input_wind_height), input_wind_height > 0)
        checked$input_wind_height <- TRUE
    }
    if ("output_wind_height" %in% names(opts)) {
        output_wind_height <- opts$output_wind_height
        assert_that(is.numeric(output_wind_height), is.scalar(output_wind_height), output_wind_height > 0)
        checked$output_wind_height <- TRUE
    }
    if ("mesh_choice" %in% names(opts)) {
        opts$mesh_choice <- match_arg(opts$mesh_choice, c("coarse", "medium", "fine"))
        checked$mesh_choice <- TRUE
    }

    if ("write_ascii_output" %in% names(opts)) {
        write_ascii_output <- opts$write_ascii_output
        assert_that(is.flag(write_ascii_output), !is.na(write_ascii_output))
        checked$write_ascii_output <- TRUE
    }

    if ("output_path" %in% names(opts)) {
        if (!fs::dir_exists(opts$output_path)) stop("output_path directory ", opts$output_path, " does not exist")
        checked$output_path <- TRUE
    }
    ## TODO “ascii_out_resolution”, “units_ascii_out_resolution”, “write_shapefile_output”, “shape_out_resolution”, “units_shape_out_resolution”, “write_goog_output”, “goog_out_resolution”, “units_goog_out_resolution”

    if (!all(as.logical(checked))) {
        notchecked <- names(checked)[!as.logical(checked)]
        warning(sprintf("some options not checked: %s", paste(dQuote(notchecked), collapse = ", ")))
    }
    opts
}


#' Generate a configuration file with domain-average wind initialization
#'
#' @param input_speed numeric: vector of wind speeds
#' @param input_speed_units string: units of input wind speed ("mps", "mph", "kph", "kts")
#' @param input_direction numeric: vector of input directions in degrees (direction that wind is blowing towards - e.g. 135 means that the wind is blowing towards the northwest)
#' @param input_wind_height numeric: height of input wind speed above the vegetation
#' @param units_input_wind_height string: units of input wind height ("ft", "m")
#' @param momentum_flag logical: use momentum solver?
#' @param number_of_iterations numeric: number of iterations for momentum solver
#' @param output_wind_height numeric: height of input wind speed above the vegetation
#' @param units_output_wind_height string: units of output wind height ("ft", "m")
#' @param output_speed_units string: units of input wind speed ("mps", "mph", "kph", "kts")
#' @param ... : additional parameters, see \code{\link{wn_config_common}}
#'
#'
#' @return The path to the configuration file
#'
#' @seealso \code{\link{wn_config_common}}
#'
#' @examples
#' \dontrun{
#'   my_config <- wn_config_domain_average(elevation = wn_demo_file("missoula_valley_elevation"),
#'                                         input_speed = 10, input_direction = 270)
#' }
#'
#' @export
wn_config_domain_average <- function(input_speed, input_speed_units = "mps", input_direction, input_wind_height = 10, units_input_wind_height = "m",
                                     momentum_flag = FALSE, number_of_iterations = 300L,
                                     output_wind_height = 10, units_output_wind_height = "m", output_speed_units = "mps", ...) {

    assert_that(is.numeric(input_speed))
    assert_that(is.numeric(input_direction))
    if (is.scalar(input_direction) && length(input_speed) > 1) {
        input_direction <- rep(input_direction, length(input_speed))
    } else if (is.scalar(input_speed) && length(input_direction) > 1) {
        input_speed <- rep(input_speed, length(input_direction))
    }
    assert_that(length(input_direction) == length(input_speed))


##output_buffer_clipping   = 5.0

    wn_check_config(c(list(initialization_method = "domainAverageInitialization",
                           momentum_flag = momentum_flag, number_of_iterations = number_of_iterations,
                           input_speed = input_speed, input_speed_units = input_speed_units,
                           input_direction = input_direction,
                           input_wind_height = input_wind_height,
                           units_input_wind_height = units_input_wind_height,
                           output_wind_height = output_wind_height,
                           units_output_wind_height = units_output_wind_height,
                           output_speed_units = output_speed_units),
                      wn_config_common(...)))
}


## write the config file
## internal, not exported
wn_write_config_file <- function(filename, config_args) {
    assert_that(is.string(filename))
    assert_that(is.list(config_args))
    wn_check_config(config_args)
    totext <- function(z) {
        if (is.logical(z)) {
            tolower(as.character(z)) ## "true" or "false"
        } else {
            as.character(z)
        }
    }
    config_args <- lapply(config_args, totext)
    writeLines(sprintf("%s = %s", names(config_args), config_args), con = filename)
    filename
}
