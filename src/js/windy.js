/*  Global class for simulating the movement of particle through a 1km wind grid

 credit: All the credit for this work goes to: https://github.com/cambecc for creating the repo:
 https://github.com/cambecc/earth. The majority of this code is directly take nfrom there, since its awesome.

 This class takes a canvas element and an array of data (1km GFS from http://www.emc.ncep.noaa.gov/index.php?branch=GFS)
 and then uses a mercator (forward/reverse) projection to correctly map wind vectors in "map space".

 The "start" method takes the bounds of the map at its current extent and starts the whole gridding,
 interpolation and animation process.
 */

var Windy = function(params) {
	var self = this;

  var MIN_VELOCITY_INTENSITY = params.minVelocity || 0; // velocity at which particle intensity is minimum (m/s)
  var MAX_VELOCITY_INTENSITY = params.maxVelocity || 10; // velocity at which particle intensity is maximum (m/s)
  var VELOCITY_SCALE =
    (params.velocityScale || 0.005) *
    (Math.pow(window.devicePixelRatio, 1 / 3) || 1); // scale for wind velocity (completely arbitrary--this value looks nice)
  var MAX_PARTICLE_AGE = params.particleAge || 90; // max number of frames a particle is drawn before regeneration
  var MIN_PARTICLE_AGE = params.particleMinAge || 10; // min number of frames a particle is drawn before regeneration
  var PARTICLE_LINE_WIDTH = params.lineWidth || 1; // line width of a drawn particle
  var PARTICLE_MULTIPLIER = params.particleMultiplier || 1 / 300; // particle count scalar (completely arbitrary--this values looks nice)
  var PARTICLE_REDUCTION = Math.pow(window.devicePixelRatio, 1 / 3) || 1.6; // multiply particle count for mobiles by this amount
  var FRAME_RATE = params.frameRate || 15;
  var FRAME_TIME = 1000 / FRAME_RATE; // desired frames per second
  var OPACITY = params.opacity || 0.97;
  var reverseX = params.reverseX || false;
  var reverseY = params.reverseY || false;
  var waveStyle = params.waveStyle || false; // particle color set by particle age for true, by intensity for false
  var dataFn = params.dataFn || null; // function for data convert

  var defaulColorScale = [
    "rgb( 36,104,180)",
    "rgb( 60,157,194)",
    "rgb(128,205,193)",
    "rgb(151,218,168)",
    "rgb(198,231,181)",
    "rgb(238,247,217)",
    "rgb(255,238,159)",
    "rgb(252,217,125)",
    "rgb(255,182,100)",
    "rgb(252,150, 75)",
    "rgb(250,112, 52)",
    "rgb(245, 64, 32)",
    "rgb(237, 45, 28)",
    "rgb(220, 24, 32)",
    "rgb(180,  0, 35)"
  ];

  const colorScale = params.colorScale || defaulColorScale;

  var NULL_WIND_VECTOR = [NaN, NaN, null]; // singleton for no wind in the form: [u, v, magnitude]

  self.grid = null;
  var gridData = params.data;
  var date;

  var setData = function(data) {
    gridData = data;
  };

  var setOptions = function(options) {
    if (options.hasOwnProperty("minVelocity"))
      MIN_VELOCITY_INTENSITY = options.minVelocity;

    if (options.hasOwnProperty("maxVelocity"))
      MAX_VELOCITY_INTENSITY = options.maxVelocity;

    if (options.hasOwnProperty("velocityScale"))
      VELOCITY_SCALE =
        (options.velocityScale || 0.005) *
        (Math.pow(window.devicePixelRatio, 1 / 3) || 1);

    if (options.hasOwnProperty("particleAge"))
      MAX_PARTICLE_AGE = options.particleAge;

    if (options.hasOwnProperty("particleMinAge"))
      MIN_PARTICLE_AGE = options.particleMinAge;

    if (options.hasOwnProperty("lineWidth"))
      PARTICLE_LINE_WIDTH = options.lineWidth;

    if (options.hasOwnProperty("particleMultiplier"))
      PARTICLE_MULTIPLIER = options.particleMultiplier;

    if (options.hasOwnProperty("opacity")) OPACITY = +options.opacity;

    if (options.hasOwnProperty("frameRate")) FRAME_RATE = options.frameRate;
    FRAME_TIME = 1000 / FRAME_RATE;

    if (options.hasOwnProperty("waveStyle")) waveStyle = options.waveStyle;

    if (options.hasOwnProperty("dataFn")) dataFn = options.dataFn;
  };


	/**
	 * @returns {Boolean} true if the specified value is not null and not undefined.
	 */
	var isValue = function(x) {
		return x !== null && x !== undefined;// && !isNaN(x);
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

	/**
	 * Calculate distortion of the wind vector caused by the shape of the projection at point (x, y). The wind
	 * vector is modified in place and returned by this function.
	 */
	var distort = function(projection, λ, φ, x, y, scale, wind) {
		var u = wind[0] * scale;
		var v = wind[1] * scale;
		var d = distortion(projection, λ, φ, x, y);

		// Scale distortion vectors by u and v, then add.
		wind[0] = d[0] * u + d[2] * v;
		wind[1] = d[1] * u + d[3] * v;
		return wind;
	};

	var distortion = function(projection, λ, φ, x, y) {
		var τ = 2 * Math.PI;
		var H = Math.pow(10, -5.2); // 0.00000630957344480193
		//var H = 0.0000360;          // 0.0000360°φ ~= 4m  (from https://github.com/cambecc/earth/blob/master/public/libs/earth/1.0.0/micro.js#L13)
		//var H = 5; // ToDo:   Why does this work?
		var hλ = λ < 0 ? H : -H;
		var hφ = φ < 0 ? H : -H;

		//var pλ = project(φ, λ + hλ);
		//var pφ = project(φ + hφ, λ);
		const pλ = mapToCanvas(λ + hλ, φ);
		const pφ = mapToCanvas(λ, φ + hφ);

		// Meridian scale factor (see Snyder, equation 4-3), where R = 1. This handles issue where length of 1º λ
		// changes depending on φ. Without this, there is a pinching effect at the poles.
		var k = Math.cos((φ / 360) * τ);
		return [
			(pλ[0] - x) / hλ / k,
			(pλ[1] - y) / hλ / k,
			(pφ[0] - x) / hφ,
			(pφ[1] - y) / hφ
		];
	};

	/**
	* Project a point on the map
	* @param λ Longitude
	* @param φ Latitude
	* @return [x, y]
	*/
	var mapToCanvas = function(λ, φ) {
		const ymin = mercY(mapBounds.south);
		const ymax = mercY(mapBounds.north);
		const xFactor = canvasBound.width / ( mapBounds.east - mapBounds.west );
		const yFactor = canvasBound.height / ( ymax - ymin );

		let y = mercY(deg2rad(φ) );
		const x = (deg2rad(λ) - mapBounds.west) * xFactor;
		y = (ymax - y) * yFactor;
		return [x, y];
	};

	var mercY = function(φ) {
		return Math.log( Math.tan( φ / 2 + Math.PI / 4 ) );
	}

	function Header(data) {
		var self = this;
		self.uData;
		self.vData;
		self.scalar;

		var tmp;
		data.forEach(function(record) {
			switch (
				record.header.parameterCategory +
				"," +
				record.header.parameterNumber
			) {
			case "1,2":
			case "2,2":
				tmp = record;
				self.uData = record.data;
				break;
			case "1,3":
			case "2,3":
				self.vData = record.data;
				break;
			default:
				self.scalar = record.data;
			}
		});

		var supported = true;

		if (data.length < 2 ) supported = false;
		if (!supported) console.log("Windy Error: data must have at least two components (u,v)");

		var header = tmp.header;
		if (header.hasOwnProperty("gridDefinitionTemplate") && header.gridDefinitionTemplate != 0 ) supported = false;
		if (!supported) {
			console.log("Windy Error: Only data with Latitude_Longitude coordinates is supported");
		}
		supported = true;  // reset for futher checks

		this.λ0 = header.lo1;
		this.φ0 = header.la1; // the grid's origin (e.g., 0.0E, 90.0N)

		this.Δλ = header.dx;
		this.Δφ = header.dy; // distance between grid points (e.g., 2.5 deg lon, 2.5 deg lat)

		this.ni = header.nx;
		this.nj = header.ny; // number of grid points W-E and N-S (e.g., 144 x 73)

		if (header.hasOwnProperty("scanMode")) {
			var scanModeMask = header.scanMode.toString(2)
			scanModeMask = ('0'+scanModeMask).slice(-8);
			var scanModeMaskArray = scanModeMask.split('').map(Number).map(Boolean);

			if (scanModeMaskArray[0]) this.Δλ = -this.Δλ;
			if (scanModeMaskArray[1]) this.Δφ = -this.Δφ;
			if (scanModeMaskArray[2]) supported = false;
			if (scanModeMaskArray[3]) supported = false;
			if (scanModeMaskArray[4]) supported = false;
			if (scanModeMaskArray[5]) supported = false;
			if (scanModeMaskArray[6]) supported = false;
			if (scanModeMaskArray[7]) supported = false;
			if (!supported) console.log("Windy Error: Data with scanMode: "+header.scanMode+ " is not supported.");
		}
		this.date = new Date(header.refTime);
		this.date.setHours(this.date.getHours() + header.forecastTime);
	}

	Header.prototype.data = function(i, uData, vData) {
		return [this.uData[i], this.vData[i]];
	}


	function Grid(header, canvasBound, mapBounds, callback) {
		this.header = header;
		this.bounds = canvasBound;
		this.grid = [];

		this.λ0 = header.λ0;
		this.φ0 = header.φ0; // the grid's origin (e.g., 0.0E, 90.0N)

		this.Δλ = header.Δλ;
		this.Δφ = header.Δφ; // distance between grid points (e.g., 2.5 deg lon, 2.5 deg lat)

		this.ni = header.ni;
		this.nj = header.nj; // number of grid points W-E and N-S (e.g., 144 x 73)

		var ni = header.ni;
		var nj = header.nj;

		var grid = this.grid;
		var p = 0;
		var isContinuous = Math.floor(ni * this.Δλ) >= 360;

		const uData = header.uData;
		const vData = header.vData;

		for (var j = 0; j < nj; j++) {
			var row = [];
			for (var i = 0; i < ni; i++, p++) {
				//if (reverseX) {
				//	row.unshift(builder.data(p, uData, vData));
				//} else {
					row.push(header.data(p, uData, vData));
					//var v = [uData[p], vData[p]];
					//row.push(v);
				//}
			}
			if (isContinuous) {
				// For wrapped grids, duplicate first column as last column to simplify interpolation logic
				//if (reverseX) {
				//	row.unshift(row[row.length-1]);
				//} else {
					row.push(row[0]);
				//}
			}
			//if (reverseY) {
				grid.unshift(row);
			//} else {
			//	grid.push(row);
			//}
		}

		callback(this);
	}
	Grid.prototype.release = function() {
		this.grid = null;
		this.header = null;
	}

	Grid.prototype.get = function(x, y) {
		var row = this.grid[Math.round(x)];
		return (row && row[Math.round(y)]) || NULL_WIND_VECTOR;
	}

	Grid.prototype.interpolate = function(λ, φ) {
		if (!this.grid) return null;

		var i = this.floorMod(λ - this.λ0, 360) / this.Δλ; // calculate longitude index in wrapped range [0, 360)
		var j = (this.φ0 - φ) / this.Δφ; // calculate latitude index in direction +90 to -90

		var fi = Math.floor(i),
			ci = fi + 1;
		var fj = Math.floor(j),
			cj = fj + 1;

		var grid = this.grid;
		var row;
		if ((row = grid[fj])) {
			var g00 = row[fi];
			var g10 = row[ci];
			if (isValue(g00) && isValue(g10) && (row = grid[cj])) {
				var g01 = row[fi];
				var g11 = row[ci];
				if (isValue(g01) && isValue(g11)) {
					// All four points found, so interpolate the value.
					return this.bilinearInterpolateVector(i - fi, j - fj, g00, g10, g01, g11);
				}
			}
		}
		return null;
	};

	Grid.prototype.bilinearInterpolateVector = function(x, y, g00, g10, g01, g11) {
		var rx = 1 - x;
		var ry = 1 - y;
		var a = rx * ry,
		b = x * ry,
		c = rx * y,
		d = x * y;
		var u = g00[0] * a + g10[0] * b + g01[0] * c + g11[0] * d;
		var v = g00[1] * a + g10[1] * b + g01[1] * c + g11[1] * d;
		return [u, v, Math.sqrt(u * u + v * v)];
	}

	/**
	 * @returns {Number} returns remainder of floored division, i.e., floor(a / n). Useful for consistent modulo
	 *          of negative numbers. See http://en.wikipedia.org/wiki/Modulo_operation.
	 */
	Grid.prototype.floorMod = function(a, n) {
		return a - n * Math.floor(a / n);
	};

	Grid.prototype.randomize = function(p) {
		// UNDONE: this method is terrible
		var x, y;
		var bounds = this.bounds;
		var v;
		var safetyNet = 0;
		do {
			x = Math.round(Math.floor(Math.random() * bounds.width) + bounds.x);
			y = Math.round(Math.floor(Math.random() * bounds.height) + bounds.y);
			v = this.interpolate(x, y);
		} while ((v === null || v[2] === null) && safetyNet++ < 30);
		p.x = x;
		p.y = y;
		return p;
	};

	function Particule() {
		this.x;
		this.y;
		this.age;
		//this.maxAge;

		this.xt;
		this.yt;
		this.intensity;
	}
	Particule.prototype.maxAge = MAX_PARTICLE_AGE;



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

  var invert = function(x, y, windy) {
    var latlon = params.map.containerPointToLatLng(L.point(x, y));
    return [latlon.lng, latlon.lat];
  };

  var project = function(lat, lon, windy) {
    var xy = params.map.latLngToContainerPoint(L.latLng(lat, lon));
    return [xy.x, xy.y];
  };

	var context2D; // cache params.canvas.getContext("2d");
	var mapBounds;
	var canvasBound;
	var colorStyles; // for runtime
	var buckets; // for runtime
	var particles; // for runtime


	var animate = function(grid2) {
		var grid = self.grid || grid2;
		self.grid = grid;
		var particleCount = Math.round(
			canvasBound.width * canvasBound.height * PARTICLE_MULTIPLIER
		);
		if (isMobile()) {
			particleCount *= PARTICLE_REDUCTION;
		}

		particles = [];
		for (var i = 0; i < particleCount; i++) {
			particles.push(
				grid.randomize({
					age: Math.floor((Math.random() * (MAX_PARTICLE_AGE-MIN_PARTICLE_AGE)) + MIN_PARTICLE_AGE) + 0
				})
			);
		}

		frame();
	};

	function evolve() {
		var projection = {}; // map.crs used instead
		var mapArea = (mapBounds.south - mapBounds.north) * (mapBounds.west - mapBounds.east);
		var velocityScale = VELOCITY_SCALE * Math.pow(mapArea, 0.4);

		buckets.forEach(function(bucket) {
			bucket.length = 0;
		});

		particles.forEach(function(particle) {
			if (particle.age > MAX_PARTICLE_AGE) {
				self.grid.randomize(particle).age = 0;
			}

			var x = particle.x;
			var y = particle.y;
			var coord = invert(x, y);
			var λ = coord[0];
			var φ = coord[1];
			var v = self.grid.interpolate(λ, φ); // vector at current position
			if (v) {
				v = distort(projection, λ, φ, x, y, velocityScale, v);
				var m = v[2];
				if (m === null || m == 0) {
					particle.age = MAX_PARTICLE_AGE; // particle has escaped the grid, never to return...
				} else {
//console.log(self, x, y, v)
					var xt = x + v[0];
					var yt = y + v[1];

					// Path from (x,y) to (xt,yt) is visible, so add this particle to the appropriate draw bucket.
					particle.xt = xt;
					particle.yt = yt;
					if(waveStyle) m = particle.age;
					buckets[colorStyles.indexFor(m)].push(particle);
				}
			}
			particle.age += 1;
		});
	}

	var draw = function (g) {
		// Fade existing particle trails.
		var prev = "lighter";
		g.globalCompositeOperation = "destination-in";
		g.fillRect(canvasBound.x, canvasBound.y, canvasBound.width, canvasBound.height);
		g.globalCompositeOperation = prev;
		g.globalAlpha = OPACITY === 0 ? 0 : OPACITY * 0.9;

		// Draw new particle trails.
		buckets.forEach(function(bucket, i) {
			if (bucket.length > 0) {
				g.beginPath();
				g.strokeStyle = colorStyles[i];
				bucket.forEach(function(particle) {
					g.moveTo(particle.x, particle.y);
					g.lineTo(particle.xt, particle.yt);
					particle.x = particle.xt;
					particle.y = particle.yt;
				});
				g.stroke();
			}
		});
	}

	var animationLoop;
	var then = Date.now();
	var frame = function () {
		animationLoop = requestAnimationFrame(frame);
		var now = Date.now();
		var delta = now - then;
		if (delta > FRAME_TIME) {
			then = now - (delta % FRAME_TIME);
			evolve();
			draw(context2D);
		}
	};

	var start = function(bounds, width, height, extent) {
		mapBounds = {
			south: deg2rad(extent[0][1]),
			north: deg2rad(extent[1][1]),
			east: deg2rad(extent[1][0]),
			west: deg2rad(extent[0][0]),
			width: width,
			height: height
		};

		stop();

		function windIntensityColorScale(min, max) {
			colorScale.indexFor = function(m) {
				// map velocity speed to a style
				return Math.max(
					0,
					Math.min(
						colorScale.length - 1,
						Math.round(((m - min) / (max - min)) * (colorScale.length - 1))
					)
				);
			};

			return colorScale;
		}

		colorStyles = windIntensityColorScale(
			MIN_VELOCITY_INTENSITY,
			MAX_VELOCITY_INTENSITY
		);
		if (waveStyle) {
			colorStyles = windIntensityColorScale(
				MIN_PARTICLE_AGE,
				MAX_PARTICLE_AGE
			);
		}

		buckets = colorStyles.map(function() {
			return [];
		});

		var fadeFillStyle = `rgba(0, 0, 0, ${OPACITY})`;
		context2D = params.canvas.getContext("2d");
		context2D.lineWidth = PARTICLE_LINE_WIDTH;
		context2D.fillStyle = fadeFillStyle;
		context2D.globalAlpha = 0.6;

		canvasBound = buildBounds(bounds, width, height);

		// build grid
		var header = new Header(gridData);
		if(dataFn && typeof dataFn === 'function') header.data = dataFn;

		new Grid(header, canvasBound, mapBounds, function(grid) {
			self.grid = grid;
			animate(grid);
		})
		//frame();
	};

	var stop = function() {
		//if (self.grid) self.grid.release();
		if (animationLoop) L.Util.cancelAnimFrame(animationLoop);
	};

	self.params = params;
	self.start = start;
	self.stop = stop;
	self.draw = draw;
	self.setData = setData;
	self.setOptions = setOptions;

	return self;
};





