// =============================================================================
// Project 3: Cropping Cycles from Space
// Sentinel-2 NDVI/NDWI time series and CHIRPS rainfall, Bojonegoro, East Java
// -----------------------------------------------------------------------------
// Self-initiated portfolio project using public data.
// Author: Utik Tri Wulan Cahya, PhD
//
// Question: How many crop cycles per year are grown on Bojonegoro's cropland,
//           and how does the cropping calendar follow rainfall (2019-2024)?
//
// Outputs (Google Drive folder P3_cropping_cycles):
//   1. P3_timeseries_10day.csv  - 10-day NDVI, NDWI, rainfall, cloud-free share
//   2. P3_ndvi_monthly_2023.tif - 12 monthly NDVI composites (one band per month)
//   3. P3_boundary.geojson      - study area boundary
//
// How to use: paste into the Earth Engine Code Editor, click Run, then open
// the Tasks tab and click RUN on each of the three exports.
// =============================================================================


// ----- 1. SETTINGS -----------------------------------------------------------

var POINT      = ee.Geometry.Point([111.88, -7.15]);  // inside Bojonegoro Regency
var START      = '2019-01-01';
var END        = '2025-01-01';     // end date is exclusive
var STEP_DAYS  = 10;               // length of each composite period
var MAP_YEAR   = 2023;             // year for the monthly NDVI maps
var CLEAR_THR  = 0.60;             // Cloud Score+ threshold (higher = stricter)
var STAT_SCALE = 60;               // m; resolution for regional averages
var FOLDER     = 'P3_cropping_cycles';
var CRS        = 'EPSG:32749';     // WGS 84 / UTM 49S


// ----- 2. STUDY AREA AND CROPLAND MASK ---------------------------------------

// Regency boundary from geoBoundaries v6 (CC BY 4.0)
var region = ee.FeatureCollection('WM/geoLab/geoBoundaries/600/ADM2')
               .filter(ee.Filter.eq('shapeGroup', 'IDN'))
               .filterBounds(POINT);
var regionGeom = region.geometry();
print('Study area:', region.aggregate_array('shapeName'));
print('Area (km2):', regionGeom.area(100).divide(1e6));

// Cropland from ESA WorldCover 2021 (class 40 = cropland, CC BY 4.0)
var cropland = ee.ImageCollection('ESA/WorldCover/v200').first()
                 .select('Map').eq(40).selfMask().clip(regionGeom);


// ----- 3. SENTINEL-2 WITH CLOUD SCORE+ MASKING -------------------------------

var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
           .filterBounds(regionGeom)
           .filterDate(START, END)
           .linkCollection(csPlus, ['cs_cdf']);

// NDVI = (NIR - Red) / (NIR + Red)            -> vegetation greenness
// NDWI = (NIR - SWIR1) / (NIR + SWIR1)        -> canopy / surface wetness
//        (Gao 1996; also called LSWI; useful for flooded rice fields)
function addIndices(img) {
  var clear = img.select('cs_cdf').gte(CLEAR_THR);
  var ndvi = img.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var ndwi = img.normalizedDifference(['B8', 'B11']).rename('NDWI');
  // toFloat() gives every image the same band type, so composites can merge
  return ndvi.addBands(ndwi).toFloat().updateMask(clear)
             .copyProperties(img, ['system:time_start']);
}
var indices = s2.map(addIndices);

// A fully masked placeholder keeps the median valid when a period has no images
var empty = ee.Image.constant([0, 0]).rename(['NDVI', 'NDWI']).toFloat().updateMask(0);

var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').select('precipitation');


// ----- 4. 10-DAY TIME SERIES -------------------------------------------------

var nPeriods = ee.Date(END).difference(ee.Date(START), 'day').divide(STEP_DAYS).floor();
var periodIdx = ee.List.sequence(0, nPeriods.subtract(1));

var series = ee.FeatureCollection(periodIdx.map(function (i) {
  var t0 = ee.Date(START).advance(ee.Number(i).multiply(STEP_DAYS), 'day');
  var t1 = t0.advance(STEP_DAYS, 'day');

  var col = indices.filterDate(t0, t1);
  var comp = col.merge(ee.ImageCollection([empty])).median().updateMask(cropland);

  var stats = comp.reduceRegion({
    reducer: ee.Reducer.mean(), geometry: regionGeom,
    scale: STAT_SCALE, maxPixels: 1e10, tileScale: 8
  });

  // Share of cropland pixels with at least one cloud-free observation
  var validShare = comp.select('NDVI').mask().rename('valid')
                       .updateMask(cropland)
                       .reduceRegion({
                         reducer: ee.Reducer.mean(), geometry: regionGeom,
                         scale: STAT_SCALE, maxPixels: 1e10, tileScale: 8
                       }).get('valid');

  // Rainfall: sum over the period, averaged over the regency
  var rain = chirps.filterDate(t0, t1).sum().reduceRegion({
    reducer: ee.Reducer.mean(), geometry: regionGeom, scale: 5566, maxPixels: 1e9
  }).get('precipitation');

  return ee.Feature(null, {
    date_start: t0.format('YYYY-MM-dd'),
    date_mid:   t0.advance(STEP_DAYS / 2, 'day').format('YYYY-MM-dd'),
    NDVI:       stats.get('NDVI'),
    NDWI:       stats.get('NDWI'),
    valid_share: validShare,
    n_images:   col.size(),
    rain_mm:    rain
  });
}));


// ----- 5. MONTHLY NDVI MAPS FOR ONE YEAR -------------------------------------

var months = ee.List.sequence(1, 12);
var monthly = ee.ImageCollection.fromImages(months.map(function (m) {
  var t0 = ee.Date.fromYMD(MAP_YEAR, m, 1);
  return indices.filterDate(t0, t0.advance(1, 'month'))
                .select('NDVI')
                .merge(ee.ImageCollection([empty.select('NDVI')]))
                .median()
                .updateMask(cropland)
                .rename(ee.String('m').cat(ee.Number(m).format('%02d')));
}));
var monthlyStack = monthly.toBands()
  .rename(['m01','m02','m03','m04','m05','m06','m07','m08','m09','m10','m11','m12']);


// ----- 6. PREVIEW ------------------------------------------------------------

Map.centerObject(region, 10);
var ndviVis = {min: 0.1, max: 0.85, palette: ['#f3efe0', '#c9d9a0', '#7fb069', '#2f7a3a', '#0f4a22']};
Map.addLayer(cropland, {palette: ['#d9b44a']}, 'Cropland (WorldCover 2021)', false);
Map.addLayer(monthlyStack.select('m03'), ndviVis, 'NDVI March ' + MAP_YEAR);
Map.addLayer(monthlyStack.select('m08'), ndviVis, 'NDVI August ' + MAP_YEAR, false);
Map.addLayer(region.style({color: '#000000', fillColor: '00000000', width: 2}), {}, 'Bojonegoro');
// The time series is not previewed here (it is heavy); it is computed in the export.


// ----- 7. EXPORTS ------------------------------------------------------------

Export.table.toDrive({
  collection: series, description: 'P3_timeseries_10day', folder: FOLDER,
  fileFormat: 'CSV',
  selectors: ['date_start', 'date_mid', 'NDVI', 'NDWI', 'valid_share', 'n_images', 'rain_mm']
});

Export.image.toDrive({
  image: monthlyStack.toFloat(), description: 'P3_ndvi_monthly_' + MAP_YEAR, folder: FOLDER,
  region: regionGeom.bounds(), scale: 60, crs: CRS, maxPixels: 1e10, fileFormat: 'GeoTIFF'
});

Export.table.toDrive({
  collection: region.select(['shapeName']), description: 'P3_boundary', folder: FOLDER,
  fileFormat: 'GeoJSON'
});

// Data credits:
// Sentinel-2: contains modified Copernicus Sentinel data (2019-2024).
// Cloud Score+: Pasquarella et al. (2023).
// CHIRPS: Funk et al. (2015), Climate Hazards Center, UC Santa Barbara.
// ESA WorldCover 2021 (CC BY 4.0). geoBoundaries v6 (CC BY 4.0).
