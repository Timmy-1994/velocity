L.GradientLayer = (L.Layer ? L.Layer : L.Class).extend({
  options: {
    colorScale: null,
    data: null,
    reverseX: false,
    reverseY: false,
    dpx: 2, // 2x2 px
    getValue: null, // function(lat, lon) return value
  },

  _map: null,
  _canvasLayer: null,
  _context: null,
  _mouseControl: null,
  _init: null,
  _timer: 0, // delay draw
  _width: 0,
  _height: 0,

  initialize: function(options) {
    L.setOptions(this, options);
  },

  onAdd: function(map) {
    // determine where to add the layer
    this._paneName = this.options.paneName || "overlayPane";

    // fall back to overlayPane for leaflet < 1
    let pane = map._panes.overlayPane;
    if (map.getPane) {
      // attempt to get pane first to preserve parent (createPane voids this)
      pane = map.getPane(this._paneName);
      if (!pane) {
        pane = map.createPane(this._paneName);
      }
    }
    // create canvas, add to map pane
    this._canvasLayer = L.canvasLayer({ pane: pane }).delegate(this);
    this._canvasLayer.addTo(map);

    this._map = map;
  },

  onRemove: function(map) {
    this._destroy();
  },

  setData: function(data) {
    this.options.data = data;
    if (this._init) {
      this._init.setData(data);
      this._clearAndRestart();
    }
    this.fire("load");
  },

  setOpacity: function(opacity) {
    console.log("this._canvasLayer", this._canvasLayer);
    this._canvasLayer.setOpacity(opacity);
  },

  setOptions: function(options) {
    this.options = Object.assign(this.options, options);
    if (options.hasOwnProperty("data")) this.options.data = options.data;
    if (this._init) {
      this._init.setOptions(options);
      if (options.hasOwnProperty("data")) this._init.setData(options.data);
      this._clearAndRestart();
    }
    //this._initMouseHandler(true);

    this.fire("load");
  },

  /*------------------------------------ PRIVATE ------------------------------------------*/

  onDrawLayer: function(overlay, params, doneCb) {
    var self = this;

    if (!this._init) {
      this._initLayer(this);
      return;
    }

    if (!this.options.data) {
      return;
    }

    if (this._timer) clearTimeout(self._timer);

    this._timer = setTimeout(function() {
      self._draw();
      doneCb();
    }, 250); // draw data is delayed
  },

  _initLayer: function(self) {
    self._resize();

    // copy options
    const options = Object.assign(
      { canvas: self._canvasLayer._canvas, map: this._map },
      self.options
    );
    this._init = new GradientCanvas(options);

    // prepare context global var, start drawing
    this._context = this._canvasLayer._canvas.getContext("2d");
    this._canvasLayer._canvas.classList.add("gradient-overlay");
    this.onDrawLayer();

    this._map.on("dragstart", self._init.stop);
    this._map.on("dragend", self._clearAndRestart);
    this._map.on("zoomstart", self._init.stop);
    this._map.on("zoomend", self._clearAndRestart);
    this._map.on("resize", self._resize);

    //this._initMouseHandler(false);
  },

  _initMouseHandler: function(voidPrevious) {
    if (voidPrevious) {
      this._map.removeControl(this._mouseControl);
      this._mouseControl = false;
    }
    if (!this._mouseControl && this.options.displayValues) {
      var options = this.options.displayOptions || {};
      options["leafletGradient"] = this;
      this._mouseControl = L.control.velocity(options).addTo(this._map);
    }
  },

  _clearAndRestart: function() {
    if (this._context) this._context.clearRect(0, 0, this._width, this._height);
    if (this._init) this._draw();
  },

  _draw: function() {
    var bounds = this._map.getBounds();
    var size = this._map.getSize();

    // bounds, width, height, extent
    this._init.start(
      [
        [0, 0],
        [size.x, size.y]
      ],
      size.x,
      size.y,
      [
        [bounds._southWest.lng, bounds._southWest.lat],
        [bounds._northEast.lng, bounds._northEast.lat]
      ]
    );
  },

  _resize: function(e) {
    var size = (e) ? e.newSize : this._map.getSize();
    this._width = size.x;
    this._height = size.y;
    if (this._init) this._init.stop();
    if (this._context) this._context.clearRect(0, 0, this._width, this._height);
  },

  _destroy: function() {
    if (this._timer) clearTimeout(this._timer);
    if (this._init) this._init.stop();
    if (this._context) this._context.clearRect(0, 0, this._width, this._height);
    if (this._mouseControl) this._map.removeControl(this._mouseControl);
    this._mouseControl = null;
    this._init = null;
    this._map.removeLayer(this._canvasLayer);
  }
});

L.gradientLayer = function(options) {
  return new L.GradientLayer(options);
};
