/*  Global class for simulating the movement of particle through a 1km wind grid

 credit: All the credit for this work goes to: https://github.com/cambecc for creating the repo:
 https://github.com/cambecc/earth. The majority of this code is directly take nfrom there, since its awesome.

 This class takes a canvas element and an array of data (1km GFS from http://www.emc.ncep.noaa.gov/index.php?branch=GFS)
 and then uses a mercator (forward/reverse) projection to correctly map wind vectors in "map space".

 The "start" method takes the bounds of the map at its current extent and starts the whole gridding,
 interpolation and animation process.
 */

var GradientCanvas = function(params) {
  var MIN_INTENSITY = params.minIntensity || 0;
  var MAX_INTENSITY = params.maxIntensity || 10;
  var OPACITY = params.opacity || 0.97;
  var reverseX = params.reverseX || false;
  var reverseY = params.reverseY || false;
  var DPX = params.dpx || 2;

  var defaulColorScale = [
    "rgba(  0,  0,  0, 0.0)", // transparent for no data
    "rgba( 36,104,180, 0.8)",
    "rgba( 60,157,194, 0.8)",
    "rgba(128,205,193, 0.8)",
    "rgba(151,218,168, 0.8)",
    "rgba(198,231,181, 0.8)",
    "rgba(238,247,217, 0.8)",
    "rgba(255,238,159, 0.8)",
    "rgba(252,217,125, 0.8)",
    "rgba(255,182,100, 0.8)",
    "rgba(252,150, 75, 0.8)",
    "rgba(250,112, 52, 0.8)",
    "rgba(245, 64, 32, 0.8)",
    "rgba(237, 45, 28, 0.8)",
    "rgba(220, 24, 32, 0.8)",
    "rgba(180,  0, 35, 0.8)"
  ];

  const colorScale = params.colorScale || defaulColorScale;

  var interpolateFn = interpolateBI;
  if (typeof params.interpolateType === 'string') {
    interpolateFn = (params.interpolateType.toLowerCase() == 'nearestneighbor')? interpolateNN : interpolateBI;
  }
  if (typeof params.interpolateType === 'function') {
    interpolateFn = params.interpolateType;
  }

  var builder;
  var grid;
  var gridData = params.data;
  var date;
  var λ0, φ0, Δλ, Δφ, ni, nj;

  var setData = function(data) {
    gridData = data;
  };

  var setOptions = function(options) {
    if (options.hasOwnProperty("minIntensity"))
      MIN_INTENSITY = options.minIntensity;

    if (options.hasOwnProperty("maxIntensity"))
      MAX_INTENSITY = options.maxIntensity;

    if (options.hasOwnProperty("opacity")) OPACITY = +options.opacity;

    if (options.hasOwnProperty("dpx")) DPX = +options.dpx;
  };

  // interpolation for value
  var bilinearInterpolateVector = function(x, y, g00, g10, g01, g11) {
    var rx = 1 - x;
    var ry = 1 - y;
    var a = rx * ry,
      b = x * ry,
      c = rx * y,
      d = x * y;
    var v = g00 * a + g10 * b + g01 * c + g11 * d;
    return v;
  };

  var createBuilder = function(data) {
    if(!data.dx) data.dx = (data.lo2 - data.lo1) / (data.nx-1);
    if(!data.dy) data.dy = (data.la1 - data.la2) / (data.ny-1);
    var nx = data.nx;
    var ny = data.ny;
    return {
      header: data,
      data: function(i) {
        return data.d['v'][i];
      }
    };
  };

  var buildGrid = function(data, callback) {

    builder = createBuilder(data);
    var header = builder.header;

    λ0 = header.lo1;
    φ0 = header.la1; // the grid's origin (e.g., 0.0E, 90.0N)

    Δλ = header.dx;
    Δφ = header.dy; // distance between grid points (e.g., 2.5 deg lon, 2.5 deg lat)

    ni = header.nx;
    nj = header.ny; // number of grid points W-E and N-S (e.g., 144 x 73)

    grid = [];
    var p = 0;
    var isContinuous = Math.floor(ni * Δλ) >= 360;

    for (var j = 0; j < nj; j++) {
      var row = [];
      for (var i = 0; i < ni; i++, p++) {
//        row[i] = builder.data(p);
        if (reverseX) {
          row.unshift(builder.data(p));
        } else {
          row.push(builder.data(p));
        }
      }
      if (isContinuous) {
        // For wrapped grids, duplicate first column as last column to simplify interpolation logic
//        row.push(row[0]);
        if (reverseX) {
          row.unshift(row[row.length-1]);
        } else {
          row.push(row[0]);
        }
      }
//      grid[j] = row;
      if (reverseY) {
        grid.unshift(row);
      } else {
        grid.push(row);
      }
    }

    callback(interpolateFn);
  };

  /**
   * Get interpolated grid value from Lon/Lat position by bilinear
   * @param λ {Float} Longitude
   * @param φ {Float} Latitude
   * @param grid {2D-Array} data grid
   * @returns {Object}
   */
  function interpolateBI(λ, φ, grid) {
    if (!grid) return null;

    var i = floorMod(λ - λ0, 360) / Δλ; // calculate longitude index in wrapped range [0, 360)
    var j = (φ0 - φ) / Δφ; // calculate latitude index in direction +90 to -90

    var fi = Math.floor(i),
      ci = fi + 1;
    var fj = Math.floor(j),
      cj = fj + 1;

    var row;
    if ((row = grid[fj])) {
      var g00 = row[fi];
      var g10 = row[ci];
      if (isValue(g00) && isValue(g10) && (row = grid[cj])) {
        var g01 = row[fi];
        var g11 = row[ci];
        if (isValue(g01) && isValue(g11)) {
          // All four points found, so interpolate the value.
          return bilinearInterpolateVector(i - fi, j - fj, g00, g10, g01, g11);
        }
      }
    }
    return null;
  };

  /**
   * Get interpolated grid value from Lon/Lat position by nearest neighbor
   * @param λ {Float} Longitude
   * @param φ {Float} Latitude
   * @param grid {2D-Array} data grid
   * @returns {Object}
   */
  function interpolateNN(λ, φ, grid) {
    if (!grid) return null;

    var i = floorMod(λ - λ0, 360) / Δλ; // calculate longitude index in wrapped range [0, 360)
    var j = (φ0 - φ) / Δφ; // calculate latitude index in direction +90 to -90

    var fi = Math.round(i);
    var fj = Math.round(j);

    var row = grid[fj]
    if (row) {
      var g00 = row[fi]
      if (isValue(g00)) return g00;
    }
    return null;
  };

  /**
   * @returns {Boolean} true if the specified value is not null and not undefined.
   */
  var isValue = function(x) {
    return x !== null && x !== undefined && (typeof x === 'number');
  };

  /**
   * @returns {Number} returns remainder of floored division, i.e., floor(a / n). Useful for consistent modulo
   *          of negative numbers. See http://en.wikipedia.org/wiki/Modulo_operation.
   */
  var floorMod = function(a, n) {
    return a - n * Math.floor(a / n);
  };

  /**
   * @returns {Number} the value x clamped to the range [low, high].
   */
  var clamp = function(x, range) {
    return Math.max(range[0], Math.min(x, range[1]));
  };

  /**
   * @returns {Boolean} true if agent is probably a mobile device. Don't really care if this is accurate.
   */
  var isMobile = function() {
    return /android|blackberry|iemobile|ipad|iphone|ipod|opera mini|webos/i.test(
      navigator.userAgent
    );
  };

  var buildBounds = function(bounds, width, height) {
    var upperLeft = bounds[0];
    var lowerRight = bounds[1];
    var x = Math.round(upperLeft[0]); //Math.max(Math.floor(upperLeft[0], 0), 0);
    var y = Math.max(Math.floor(upperLeft[1], 0), 0);
    var xMax = Math.min(Math.ceil(lowerRight[0], width), width - 1);
    var yMax = Math.min(Math.ceil(lowerRight[1], height), height - 1);
    return {
      x: x,
      y: y,
      xMax: width,
      yMax: yMax,
      width: width,
      height: height
    };
  };

  var deg2rad = function(deg) {
    return (deg / 180) * Math.PI;
  };

  var invert = function(x, y) {
    var latlon = params.map.containerPointToLatLng(L.point(x, y));
    return [latlon.lng, latlon.lat];
  };

  var project = function(lat, lon) {
    var xy = params.map.latLngToContainerPoint(L.latLng(lat, lon));
    return [xy.x, xy.y];
  };

  var gradient = function (grad) {
        // create a 256x1 gradient color
        var canvas = document.createElement('canvas');
        var ctx = canvas.getContext('2d');
        var gradient = ctx.createLinearGradient(0, 0, 0, 256);

        canvas.width = 1;
        canvas.height = 256;

	var max = grad.length
        for (var i in grad) {
            gradient.addColorStop(i/max, grad[i]);
        }

        ctx.fillStyle = gradient;
        ctx.fillRect(0, 0, 1, 256);

        return ctx.getImageData(0, 0, 1, 256).data;
    };

  var createField = function(columns, bounds, callback) {
    /**
     * @returns value at the point (x, y), or null if value is undefined at that point.
     */
    function field(x, y) {
      var column = columns[Math.round(x)];
      return (column && column[Math.round(y)]) || null;
    }

    // Frees the massive "columns" array for GC. Without this, the array is leaked (in Chrome) each time a new
    // field is interpolated because the field closure's context is leaked, for reasons that defy explanation.
    field.release = function() {
      columns = [];
    };

    callback(bounds, field, columns);
  };

  var interpolateField = function(gridInterpolateFn, bounds, extent, callback) {
    var columns = [];
    var x = bounds.x;

    function interpolateColumn(x) {
      var w = bounds.width;
      var h = bounds.height;
      var column = [];
      for (var y = bounds.y; y <= bounds.yMax; y += DPX) {
        var coord = invert(x, y);
        if (coord) {
          var λ = coord[0],
            φ = coord[1];
          if (isFinite(λ)) {
            var value = gridInterpolateFn(λ, φ, grid);
            for(var k=0; k<DPX; k++) column[y + k] = value;
          }
        }
      }
      for(var k=0; k<DPX; k++) columns[x + k] = column;
    }

    (function batchInterpolate() {
      var start = Date.now();
      while (x < bounds.width) {
        interpolateColumn(x);
        x += DPX;
        if (Date.now() - start > 500) { // not to block too long
          setTimeout(batchInterpolate, 25);
          return;
        }
      }
      createField(columns, bounds, callback);
    })();
  };

  var draw = function(bounds, field, columns) {
    var colorStyles = gradient(colorScale);
    var intensityColorScale = function (val, out) { // map value to a style
      var idx = Math.max(
        0,
        Math.min(
          255,
          Math.round(((val - MIN_INTENSITY) / (MAX_INTENSITY - MIN_INTENSITY)) * 255)
        )
      );
      idx = idx * 4;
      //return [colorStyles[idx], colorStyles[idx+1], colorStyles[idx+2], colorStyles[idx+3]];
      out[0] = colorStyles[idx + 0];
      out[1] = colorStyles[idx + 1];
      out[2] = colorStyles[idx + 2];
      out[3] = colorStyles[idx + 3];
    };

    var w = bounds.width;
    var h = bounds.height;
    var g = params.canvas.getContext("2d");
    var imageData = g.getImageData(0, 0, w, h);
    var px = imageData.data;

    var color = [0,0,0,0];
//console.log("[gradient]draw", bounds, w, h, px.length / 4, grid, columns);
    for(var i=0; i<w; i++) {
      var col = columns[i];
      for(var j=0; j<h; j++) {
        var k = (j*w + i) * 4;
        intensityColorScale(col[j], color);
        px[k + 0] = color[0]; // red
        px[k + 1] = color[1]; // green
        px[k + 2] = color[2]; // blue
        px[k + 3] = color[3]; // alpha
      }
    }
    g.putImageData(imageData, 0, 0);
  };



  var start = function(bounds, width, height, extent) {
    var mapBounds = {
      south: deg2rad(extent[0][1]),
      north: deg2rad(extent[1][1]),
      east: deg2rad(extent[1][0]),
      west: deg2rad(extent[0][0]),
      width: width,
      height: height
    };

    stop();


    // build grid
    buildGrid(gridData, function(gridInterpolateFn) {
      interpolateField(
        gridInterpolateFn,
        buildBounds(bounds, width, height),
        mapBounds,
        function(bounds, field, columns) {
          obj.field = field;
          draw(bounds, field, columns);
        }
      );
    });
  };

  var stop = function() {
    if (obj.field) obj.field.release();
  };

  var obj = {
    params: params,
    start: start,
    stop: stop,
    createField: createField,
    interpolatePoint: interpolateFn,
    setData: setData,
    setOptions: setOptions
  };

  return obj;
};

