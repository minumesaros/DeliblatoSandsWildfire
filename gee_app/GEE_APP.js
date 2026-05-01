var roi = ee.FeatureCollection("projects/ee-minucher/assets/DeliblatoROI_SQFIN"),
    roisrp = ee.FeatureCollection("projects/ee-minucher/assets/SRPDeliblatoWGS84"),
    stat4 = ee.Image("projects/ee-minucher/assets/StaticV4"),
    fc = ee.ImageCollection("ECMWF/NRT_FORECAST/IFS/OPER");

// CiROCCO FireRiskApp V007

/**** USER INPUTS ****/
// imports: roi, roisrp, stat4

/**** DATASETS ****/
var era5LandHourly = ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY');
var era5Hourly = ee.ImageCollection('ECMWF/ERA5/HOURLY');
var ifsForecast = ee.ImageCollection('ECMWF/NRT_FORECAST/IFS/OPER');
var cmip6 = ee.ImageCollection('NASA/GDDP-CMIP6');

/**** HELPERS ****/
function norm(img, minVal, maxVal) {
  return img.subtract(minVal)
            .divide(ee.Number(maxVal).subtract(minVal))
            .clamp(0, 1);
}

function pad2(n) {
  return (n < 10 ? '0' : '') + n;
}

function daysInMonth(year, month) {
  return new Date(Date.UTC(year, month, 0)).getUTCDate();
}

function meanOverArea(img, geom, scale) {
  return ee.Number(img.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: geom,
    scale: scale,
    bestEffort: true,
    maxPixels: 1e9
  }).values().get(0));
}

function paintBoundary() {
  return ee.Image().paint(roisrp, 1, 2);
}

function round1(x) {
  return (x === null || x === undefined) ? 'n/a' : (Math.round(x * 10) / 10).toFixed(1);
}

/**** METADATA FOR UI ****/
// Previous measurements availability from intersection of ERA5-Land and ERA5.
var pastMinMs = Math.max(
  ee.Number(era5LandHourly.aggregate_min('system:time_start')).getInfo(),
  ee.Number(era5Hourly.aggregate_min('system:time_start')).getInfo()
);

var pastMaxMs = Math.min(
  ee.Number(era5LandHourly.aggregate_max('system:time_start')).getInfo(),
  ee.Number(era5Hourly.aggregate_max('system:time_start')).getInfo()
);

var PAST_MIN_DATE = new Date(pastMinMs);
var PAST_MAX_DATE = new Date(pastMaxMs);

var PAST_MIN_YEAR = PAST_MIN_DATE.getUTCFullYear();
var PAST_MIN_MONTH = PAST_MIN_DATE.getUTCMonth() + 1;
var PAST_MIN_DAY = PAST_MIN_DATE.getUTCDate();

var PAST_MAX_YEAR = PAST_MAX_DATE.getUTCFullYear();
var PAST_MAX_MONTH = PAST_MAX_DATE.getUTCMonth() + 1;
var PAST_MAX_DAY = PAST_MAX_DATE.getUTCDate();

// Forecast availability from the latest run.
var latestForecastCreation = ee.Number(ifsForecast.aggregate_max('creation_time'));
var latestForecastRun = ifsForecast.filter(
  ee.Filter.eq('creation_time', latestForecastCreation)
);

var latestForecastRunLabel = ee.Date(latestForecastCreation)
  .format('YYYY-MM-dd HH:mm')
  .getInfo();

var forecastDates = ee.List(
  latestForecastRun.aggregate_array('forecast_time')
).map(function(t) {
  return ee.Date(ee.Number(t)).format('YYYY-MM-dd');
}).distinct().sort().getInfo();

var forecastDefault = forecastDates.length > 1 ? forecastDates[1] : forecastDates[0];

/**** PREVIOUS MEASUREMENTS BRANCH ****/
function buildPastRiskImage(dateString) {
  var target = ee.Date(dateString + 'T13:00:00');
  var targetEnd = target.advance(1, 'hour');

  var eraLand = era5LandHourly.filterDate(target, targetEnd).first();
  var era = era5Hourly.filterDate(target, targetEnd).first();

  var tC = eraLand.select('temperature_2m').subtract(273.15).rename('tC');
  var tdC = eraLand.select('dewpoint_temperature_2m').subtract(273.15).rename('tdC');

  var solarW = eraLand.select('surface_solar_radiation_downwards_hourly')
    .divide(3600)
    .rename('solarW');

  var rainMm = eraLand.select('total_precipitation_hourly')
    .multiply(1000)
    .rename('rainMm');

  var u10 = eraLand.select('u_component_of_wind_10m');
  var v10 = eraLand.select('v_component_of_wind_10m');
  var wind = u10.pow(2).add(v10.pow(2)).sqrt().rename('wind');

  var soilMoist = eraLand.select('volumetric_soil_water_layer_1').rename('soilMoist');

  var soilTempC = eraLand.select('soil_temperature_level_1')
    .subtract(273.15)
    .rename('soilTempC');

  var cloud = era.select('total_cloud_cover').rename('cloud');

  var es = tC.expression(
    '0.6108 * exp(17.27 * T / (T + 237.3))', {T: tC}
  ).rename('es');

  var ea = tdC.expression(
    '0.6108 * exp(17.27 * Td / (Td + 237.3))', {Td: tdC}
  ).rename('ea');

  var rh = ea.divide(es).multiply(100).clamp(0, 100).rename('RH');
  var vpd = es.subtract(ea).max(0).rename('VPD');

  var tRisk          = norm(tC, 15, 38);
  var vpdRisk        = norm(vpd, 0.5, 3.5);
  var windRisk       = norm(wind, 1, 10);
  var solarRisk      = norm(solarW, 150, 850);
  var rainDryRisk    = ee.Image(1).subtract(norm(rainMm, 0.2, 2.0));
  var soilDryRisk    = ee.Image(1).subtract(norm(soilMoist, 0.12, 0.38));
  var soilTempRisk   = norm(soilTempC, 10, 32);
  var cloudClearRisk = ee.Image(1).subtract(norm(cloud, 0.3, 1.0));

  var weatherIndex = tRisk.multiply(0.20)
    .add(vpdRisk.multiply(0.20))
    .add(windRisk.multiply(0.20))
    .add(solarRisk.multiply(0.10))
    .add(rainDryRisk.multiply(0.10))
    .add(soilDryRisk.multiply(0.10))
    .add(soilTempRisk.multiply(0.05))
    .add(cloudClearRisk.multiply(0.05))
    .rename('weatherIndex');

  var modifier = ee.Image(0.2).add(weatherIndex.multiply(1.6)).rename('modifier');
  modifier = modifier.where(rainMm.gt(2.0), 0.15);
  modifier = modifier.where(rh.gt(95), 0.25);
  modifier = modifier.where(rainMm.gt(0.5).and(rh.gt(85)), 0.35);

  var fireRisk = stat4.multiply(modifier)
    .clamp(0, 100)
    .rename('fireRisk');

  return fireRisk.updateMask(fireRisk.gt(0));
}

/**** SHORT TERM FORECAST BRANCH ****/
function buildForecastRiskResult(dateString) {
  var target = ee.Date(dateString + 'T13:00:00');

  var img = ee.Image(
    latestForecastRun.map(function(im) {
      var diff = ee.Number(im.get('forecast_time')).subtract(target.millis()).abs();
      return im.set('timeDiff', diff);
    }).sort('timeDiff').first()
  );

  var tC = img.select('temperature_2m_sfc').rename('tC');
  var tdC = img.select('dewpoint_temperature_2m_sfc').rename('tdC');

  var u10 = img.select('u_component_of_wind_10m_sfc');
  var v10 = img.select('v_component_of_wind_10m_sfc');
  var wind = u10.pow(2).add(v10.pow(2)).sqrt().rename('wind');

  var rainMmH = img.select('total_precipitation_rate_sfc')
    .multiply(3600 * 1000)
    .rename('rainMmH');

  var solarJ = img.select('surface_solar_radiation_downwards_sfc').rename('solarJ');

  var soilTempC = img.select('soil_temperature_sol1')
    .subtract(273.15)
    .rename('soilTempC');

  var soilMoist = img.select('volumetric_soil_moisture_sol1').rename('soilMoist');

  var es = tC.expression(
    '0.6108 * exp(17.27 * T / (T + 237.3))', {T: tC}
  ).rename('es');

  var ea = tdC.expression(
    '0.6108 * exp(17.27 * Td / (Td + 237.3))', {Td: tdC}
  ).rename('ea');

  var rh = ea.divide(es).multiply(100).clamp(0, 100).rename('RH');
  var vpd = es.subtract(ea).max(0).rename('VPD');

  var scale = 28000;

  var tMean = meanOverArea(tC.rename('x'), roisrp, scale);
  var tdMean = meanOverArea(tdC.rename('x'), roisrp, scale);
  var rhMean = meanOverArea(rh.rename('x'), roisrp, scale);
  var vpdMean = meanOverArea(vpd.rename('x'), roisrp, scale);
  var windMean = meanOverArea(wind.rename('x'), roisrp, scale);
  var rainMean = meanOverArea(rainMmH.rename('x'), roisrp, scale);
  var solarMean = meanOverArea(solarJ.rename('x'), roisrp, scale);
  var soilMoistMean = meanOverArea(soilMoist.rename('x'), roisrp, scale);
  var soilTempMean = meanOverArea(soilTempC.rename('x'), roisrp, scale);

  var tConst = ee.Image.constant(tMean).clip(roi);
  var rhConst = ee.Image.constant(rhMean).clip(roi);
  var vpdConst = ee.Image.constant(vpdMean).clip(roi);
  var windConst = ee.Image.constant(windMean).clip(roi);
  var rainConst = ee.Image.constant(rainMean).clip(roi);
  var solarConst = ee.Image.constant(solarMean).clip(roi);
  var soilMoistConst = ee.Image.constant(soilMoistMean).clip(roi);
  var soilTempConst = ee.Image.constant(soilTempMean).clip(roi);

  var tRisk        = norm(tConst, 15, 38);
  var vpdRisk      = norm(vpdConst, 0.5, 3.5);
  var windRisk     = norm(windConst, 1, 12);
  var solarRisk    = norm(solarConst, 0, 3.0e7);
  var rainDryRisk  = ee.Image(1).subtract(norm(rainConst, 0.1, 2.0));
  var soilDryRisk  = ee.Image(1).subtract(norm(soilMoistConst, 0.12, 0.38));
  var soilTempRisk = norm(soilTempConst, 10, 32);

  var weatherIndex = tRisk.multiply(0.22)
    .add(vpdRisk.multiply(0.22))
    .add(windRisk.multiply(0.22))
    .add(solarRisk.multiply(0.12))
    .add(rainDryRisk.multiply(0.10))
    .add(soilDryRisk.multiply(0.07))
    .add(soilTempRisk.multiply(0.05))
    .rename('weatherIndex');

  var modifier = ee.Image(0.2).add(weatherIndex.multiply(1.6)).rename('modifier');
  modifier = modifier.where(rainConst.gt(2.0), 0.15);
  modifier = modifier.where(rhConst.gt(95), 0.25);
  modifier = modifier.where(rainConst.gt(0.5).and(rhConst.gt(85)), 0.35);

  var fireRisk = stat4.multiply(modifier)
    .clamp(0, 100)
    .rename('fireRisk')
    .updateMask(stat4.gt(0));

  var meteo = ee.Dictionary({
    tC: tMean,
    tdC: tdMean,
    rh: rhMean,
    vpd: vpdMean,
    wind: windMean,
    rain: rainMean,
    solarJ: solarMean,
    soilTempC: soilTempMean,
    soilMoist: soilMoistMean
  });

  return {
    fireRisk: fireRisk,
    meteo: meteo
  };
}

/**** LONG TERM CLIMATE BRANCH ****/
function bandMonthlyEnsembleMedian(scenarioName, startDate, endDate, bandName, monthNum) {
  var base = cmip6
    .filterDate(startDate, endDate)
    .filter(ee.Filter.eq('scenario', scenarioName))
    .filter(ee.Filter.calendarRange(monthNum, monthNum, 'month'))
    .filter(ee.Filter.listContains('system:band_names', bandName))
    .select(bandName);

  var models = ee.List(base.aggregate_array('model')).distinct();

  var perModelMeans = ee.ImageCollection.fromImages(
    models.map(function(model) {
      model = ee.String(model);
      return base
        .filter(ee.Filter.eq('model', model))
        .mean()
        .rename(bandName)
        .set('model', model);
    })
  );

  return perModelMeans.median().rename(bandName);
}

function buildClimateRiskImage(monthNum, scenario) {
  var histStart = '1995-01-01';
  var histEnd   = '2014-12-31';
  var futStart  = '2031-01-01';
  var futEnd    = '2041-12-31';
  var scale = 27830;

  var histTas  = bandMonthlyEnsembleMedian('historical', histStart, histEnd, 'tas', monthNum).subtract(273.15);
  var histPr   = bandMonthlyEnsembleMedian('historical', histStart, histEnd, 'pr', monthNum).multiply(86400);
  var histHurs = bandMonthlyEnsembleMedian('historical', histStart, histEnd, 'hurs', monthNum);
  var histRsds = bandMonthlyEnsembleMedian('historical', histStart, histEnd, 'rsds', monthNum);
  var histWind = bandMonthlyEnsembleMedian('historical', histStart, histEnd, 'sfcWind', monthNum);

  var futTas  = bandMonthlyEnsembleMedian(scenario, futStart, futEnd, 'tas', monthNum).subtract(273.15);
  var futPr   = bandMonthlyEnsembleMedian(scenario, futStart, futEnd, 'pr', monthNum).multiply(86400);
  var futHurs = bandMonthlyEnsembleMedian(scenario, futStart, futEnd, 'hurs', monthNum);
  var futRsds = bandMonthlyEnsembleMedian(scenario, futStart, futEnd, 'rsds', monthNum);
  var futWind = bandMonthlyEnsembleMedian(scenario, futStart, futEnd, 'sfcWind', monthNum);

  var histTasMean  = meanOverArea(histTas.rename('x'), roisrp, scale);
  var histPrMean   = meanOverArea(histPr.rename('x'), roisrp, scale);
  var histHursMean = meanOverArea(histHurs.rename('x'), roisrp, scale);
  var histRsdsMean = meanOverArea(histRsds.rename('x'), roisrp, scale);
  var histWindMean = meanOverArea(histWind.rename('x'), roisrp, scale);

  var futTasMean  = meanOverArea(futTas.rename('x'), roisrp, scale);
  var futPrMean   = meanOverArea(futPr.rename('x'), roisrp, scale);
  var futHursMean = meanOverArea(futHurs.rename('x'), roisrp, scale);
  var futRsdsMean = meanOverArea(futRsds.rename('x'), roisrp, scale);
  var futWindMean = meanOverArea(futWind.rename('x'), roisrp, scale);

  var dTas  = futTasMean.subtract(histTasMean);
  var dHurs = histHursMean.subtract(futHursMean);
  var dPr   = histPrMean.subtract(futPrMean);
  var dRsds = futRsdsMean.subtract(histRsdsMean);
  var dWind = futWindMean.subtract(histWindMean);

  var tTerm    = ee.Image.constant(dTas.divide(4)).clamp(-1, 1);
  var hursTerm = ee.Image.constant(dHurs.divide(20)).clamp(-1, 1);
  var prTerm   = ee.Image.constant(dPr.divide(3)).clamp(-1, 1);
  var rsdsTerm = ee.Image.constant(dRsds.divide(40)).clamp(-1, 1);
  var windTerm = ee.Image.constant(dWind.divide(2)).clamp(-1, 1);

  var climateSignal = tTerm.multiply(0.35)
    .add(hursTerm.multiply(0.20))
    .add(prTerm.multiply(0.20))
    .add(windTerm.multiply(0.15))
    .add(rsdsTerm.multiply(0.10))
    .rename('climateSignal');

  var climateModifier = ee.Image(1)
    .add(climateSignal.multiply(0.5))
    .clamp(0.7, 1.5)
    .rename('climateModifier');

  var fireRisk = stat4.multiply(climateModifier)
    .clamp(0, 100)
    .rename('fireRisk');

  return fireRisk.updateMask(stat4.gt(0));
}

/**** UI HELPERS ****/
function getPastDateString() {
  return pastYearSelect.getValue() + '-' +
         pastMonthSelect.getValue() + '-' +
         pastDaySelect.getValue();
}

function refreshPastMonthOptions() {
  var year = parseInt(pastYearSelect.getValue(), 10);
  var minMonth = (year === PAST_MIN_YEAR) ? PAST_MIN_MONTH : 1;
  var maxMonth = (year === PAST_MAX_YEAR) ? PAST_MAX_MONTH : 12;

  var currentMonth = parseInt(pastMonthSelect.getValue() || minMonth, 10);
  var monthItems = [];
  for (var m = minMonth; m <= maxMonth; m++) {
    monthItems.push(pad2(m));
  }

  pastMonthSelect.items().reset(monthItems);

  if (currentMonth < minMonth) currentMonth = minMonth;
  if (currentMonth > maxMonth) currentMonth = maxMonth;

  pastMonthSelect.setValue(pad2(currentMonth), false);
  refreshPastDayOptions();
}

function refreshPastDayOptions() {
  var year = parseInt(pastYearSelect.getValue(), 10);
  var month = parseInt(pastMonthSelect.getValue(), 10);

  var minDay = 1;
  var maxDay = daysInMonth(year, month);

  if (year === PAST_MIN_YEAR && month === PAST_MIN_MONTH) minDay = PAST_MIN_DAY;
  if (year === PAST_MAX_YEAR && month === PAST_MAX_MONTH) maxDay = PAST_MAX_DAY;

  var currentDay = parseInt(pastDaySelect.getValue() || minDay, 10);
  var dayItems = [];
  for (var d = minDay; d <= maxDay; d++) {
    dayItems.push(pad2(d));
  }

  pastDaySelect.items().reset(dayItems);

  if (currentDay < minDay) currentDay = minDay;
  if (currentDay > maxDay) currentDay = maxDay;

  pastDaySelect.setValue(pad2(currentDay), false);
}

function setPastDate(jsDate) {
  if (jsDate < PAST_MIN_DATE) jsDate = new Date(pastMinMs);
  if (jsDate > PAST_MAX_DATE) jsDate = new Date(pastMaxMs);

  pastYearSelect.setValue(String(jsDate.getUTCFullYear()), false);
  refreshPastMonthOptions();
  pastMonthSelect.setValue(pad2(jsDate.getUTCMonth() + 1), false);
  refreshPastDayOptions();
  pastDaySelect.setValue(pad2(jsDate.getUTCDate()), false);
}

function shiftPastDay(delta) {
  var current = new Date(getPastDateString() + 'T00:00:00Z');
  current.setUTCDate(current.getUTCDate() + delta);
  setPastDate(current);
}

function shiftForecastDay(delta) {
  var current = forecastDateSelect.getValue();
  var idx = forecastDates.indexOf(current);
  if (idx < 0) idx = 0;
  idx = Math.max(0, Math.min(forecastDates.length - 1, idx + delta));
  forecastDateSelect.setValue(forecastDates[idx], false);
}

function shiftClimateMonth(delta) {
  var current = parseInt(climateMonthSelect.getValue(), 10);
  current = current + delta;
  if (current < 1) current = 12;
  if (current > 12) current = 1;
  climateMonthSelect.setValue(pad2(current), false);
}

function updateModeUI() {
  var mode = timeframeSelect.getValue();

  pastPanel.style().set('shown', mode === 'Previous measurements');
  forecastPanel.style().set('shown', mode === 'Short term Forecast (15 days)');
  climatePanel.style().set('shown', mode === 'Long term climate (2031-2041)');

  if (mode === 'Long term climate (2031-2041)') {
    prevButton.setLabel('Previous month');
    nextButton.setLabel('Next month');
  } else {
    prevButton.setLabel('Previous day');
    nextButton.setLabel('Next day');
  }
}

function shiftSelection(delta) {
  var mode = timeframeSelect.getValue();
  if (mode === 'Previous measurements') shiftPastDay(delta);
  if (mode === 'Short term Forecast (15 days)') shiftForecastDay(delta);
  if (mode === 'Long term climate (2031-2041)') shiftClimateMonth(delta);
  renderCurrentSelection();
}

function renderCurrentSelection() {
  var mode = timeframeSelect.getValue();
  statusLabel.setValue('Loading...');
  meteoLabel.setValue('');
  map.layers().reset();

  var fireRisk;
  var layerName;
  var statusText;

  if (mode === 'Previous measurements') {
    var pastDate = getPastDateString();
    fireRisk = buildPastRiskImage(pastDate);
    layerName = 'Previous measurements fire risk 13:00';
    statusText = 'Previous measurements: ' + pastDate + ' at 13:00 UTC';
  }

  if (mode === 'Short term Forecast (15 days)') {
    var fcDate = forecastDateSelect.getValue();
    var forecastResult = buildForecastRiskResult(fcDate);
    fireRisk = forecastResult.fireRisk;
    layerName = 'Short term forecast fire risk 13:00';
    statusText = 'Short term Forecast (15 days): ' + fcDate + ' at 13:00 UTC | latest run: ' + latestForecastRunLabel + ' UTC';

    forecastResult.meteo.evaluate(function(vals) {
      if (timeframeSelect.getValue() !== 'Short term Forecast (15 days)') return;
      meteoLabel.setValue(
        'Air temperature: ' + round1(vals.tC) + ' °C\n' +
        'Dew point temperature: ' + round1(vals.tdC) + ' °C\n' +
        'Relative humidity: ' + round1(vals.rh) + ' %\n' +
        'VPD: ' + round1(vals.vpd) + ' kPa\n' +
        'Wind speed: ' + round1(vals.wind) + ' m/s\n' +
        'Solar radiation: ' + round1(vals.solarJ) + ' J/m²\n' +
        'Precipitation rate: ' + round1(vals.rain) + ' mm/h\n' +
        'Soil temperature: ' + round1(vals.soilTempC) + ' °C\n' +
        'Soil moisture: ' + round1(vals.soilMoist)
      );
    });
  }

  if (mode === 'Long term climate (2031-2041)') {
    var climateMonth = parseInt(climateMonthSelect.getValue(), 10);
    var climateScenario = climateScenarioSelect.getValue();
    fireRisk = buildClimateRiskImage(climateMonth, climateScenario);
    layerName = 'Long term climate fire risk';
    statusText = 'Long term climate (2031-2041): ' + climateScenario + ' | month ' + climateMonthSelect.getValue() + ' | 2031–2041 ensemble median';
  }

  map.addLayer(
    fireRisk.clip(roi),
    {min: 0, max: 100, palette: riskPalette},
    layerName,
    true
  );

  map.addLayer(
    paintBoundary(),
    {palette: ['00FFFF']},
    'Protected area boundary',
    true
  );

  statusLabel.setValue(statusText);
}

/**** UI ****/
ui.root.clear();

var map = ui.Map();
map.centerObject(roisrp, 10);
map.setOptions('HYBRID');

var title = ui.Label('CiRROCO wildfire risk', {
  fontWeight: 'bold',
  fontSize: '20px',
  margin: '0 0 8px 0'
});

var subtitle = ui.Label('Select timeframe and date', {
  fontSize: '14px',
  margin: '0 0 10px 0'
});

var timeframeLabel = ui.Label('Timeframe:', {margin: '0 0 4px 0'});
var timeframeSelect = ui.Select({
  items: [
    'Previous measurements',
    'Short term Forecast (15 days)',
    'Long term climate (2031-2041)'
  ],
  value: 'Previous measurements',
  style: {stretch: 'horizontal'}
});

var dateLabel = ui.Label('Selection:', {margin: '8px 0 4px 0'});

/*** Previous measurements controls ***/
var pastYearItems = [];
for (var y = PAST_MIN_YEAR; y <= PAST_MAX_YEAR; y++) {
  pastYearItems.push(String(y));
}

var pastYearSelect = ui.Select({
  items: pastYearItems,
  value: String(PAST_MAX_YEAR),
  style: {width: '80px'}
});

var pastMonthSelect = ui.Select({
  items: [],
  style: {width: '70px'}
});

var pastDaySelect = ui.Select({
  items: [],
  style: {width: '70px'}
});

pastYearSelect.onChange(function() { refreshPastMonthOptions(); });
pastMonthSelect.onChange(function() { refreshPastDayOptions(); });

var pastPanel = ui.Panel(
  [pastYearSelect, pastMonthSelect, pastDaySelect],
  ui.Panel.Layout.Flow('horizontal')
);

/*** Short term forecast controls ***/
var forecastDateSelect = ui.Select({
  items: forecastDates,
  value: forecastDefault,
  style: {stretch: 'horizontal'}
});

var forecastInfo = ui.Label('Latest run: ' + latestForecastRunLabel + ' UTC', {
  color: 'gray',
  margin: '4px 0 0 0',
  fontSize: '11px'
});

var forecastPanel = ui.Panel({
  widgets: [forecastDateSelect, forecastInfo],
  style: {shown: false}
});

/*** Long term climate controls ***/
var climateMonthItems = [];
for (var m = 1; m <= 12; m++) {
  climateMonthItems.push(pad2(m));
}

var climateMonthSelect = ui.Select({
  items: climateMonthItems,
  value: '08',
  style: {width: '70px'}
});

var climateScenarioSelect = ui.Select({
  items: ['ssp245', 'ssp585'],
  value: 'ssp245',
  style: {width: '90px'}
});

var climateInfo = ui.Label('Projection window: 2031–2041 ensemble median', {
  color: 'gray',
  margin: '4px 0 0 0',
  fontSize: '11px'
});

var climatePanel = ui.Panel({
  widgets: [
    ui.Panel([climateMonthSelect, climateScenarioSelect], ui.Panel.Layout.Flow('horizontal')),
    climateInfo
  ],
  style: {shown: false}
});

var prevButton = ui.Button({
  label: 'Previous day',
  style: {width: '120px', margin: '6px 4px 0 0'},
  onClick: function() { shiftSelection(-1); }
});

var nextButton = ui.Button({
  label: 'Next day',
  style: {width: '120px', margin: '6px 0 0 4px'},
  onClick: function() { shiftSelection(1); }
});

var showButton = ui.Button({
  label: 'Show map',
  style: {stretch: 'horizontal', margin: '8px 0 0 0'},
  onClick: renderCurrentSelection
});

var statusLabel = ui.Label('', {
  color: 'gray',
  margin: '8px 0 0 0'
});

var riskPalette = [
  '#006400', '#228B22', '#7FBF3F',
  '#ADFF2F', '#FFFF00', '#FFD700',
  '#FFA500', '#FF7F00', '#FF4500',
  '#8B0000'
];

var legendTitle = ui.Label('Fire risk', {
  fontWeight: 'bold',
  margin: '12px 0 6px 0'
});

var legend = ui.Panel();

function addLegendRow(colors, text) {
  var row = ui.Panel({layout: ui.Panel.Layout.Flow('horizontal')});
  colors.forEach(function(c) {
    row.add(ui.Label('', {
      backgroundColor: c,
      padding: '8px',
      margin: '0 2px 4px 0'
    }));
  });
  row.add(ui.Label(text, {margin: '0 0 4px 6px'}));
  legend.add(row);
}

addLegendRow(['#006400', '#228B22', '#7FBF3F'], '0–30   Low');
addLegendRow(['#ADFF2F', '#FFFF00', '#FFD700'], '30–60   Medium');
addLegendRow(['#FFA500', '#FF7F00', '#FF4500'], '60–90   High');
addLegendRow(['#8B0000'], '90–100  Very high');

var meteoLabel = ui.Label('', {
  color: 'gray',
  margin: '8px 0 0 0',
  fontSize: '11px',
  whiteSpace: 'pre-wrap'
});

var buttonRow = ui.Panel(
  [prevButton, nextButton],
  ui.Panel.Layout.Flow('horizontal')
);

var panel = ui.Panel({
  widgets: [
    title,
    subtitle,
    timeframeLabel,
    timeframeSelect,
    dateLabel,
    pastPanel,
    forecastPanel,
    climatePanel,
    buttonRow,
    showButton,
    statusLabel,
    legendTitle,
    legend,
    meteoLabel
  ],
  style: {
    width: '360px',
    padding: '12px'
  }
});

timeframeSelect.onChange(function() {
  updateModeUI();
});

ui.root.add(panel);
ui.root.add(map);

// Initial state
setPastDate(PAST_MAX_DATE);
updateModeUI();
renderCurrentSelection();