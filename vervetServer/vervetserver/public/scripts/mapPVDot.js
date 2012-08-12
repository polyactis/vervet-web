

function mapPVDot(lat, lon, value, title, objIndex, minValue, maxValue) {
	this.dotColorScale = pv.Scale.linear(minValue, maxValue).range("blue", "red");
	this.dotSizeScale = pv.Scale.linear(minValue, maxValue).range(50, 500);
	
	this.lat = lat;
	this.lon = lon;
	this.value = value;
	this.title = title;
	this.objIndex = objIndex;
	this.minValue = minValue;
	this.maxValue = maxValue;
	this.radius = this.dotSizeScale(this.value);
}

mapPVDot.prototype = pv.extend(GOverlay);

mapPVDot.prototype.initialize = function(map) {
  this.map = map;
  this.canvas = document.createElement("div");
  this.canvas.setAttribute("class", "canvas");
  //Call instance method instanceFoo() on x
  var pane = map.@com.google.gwt.maps.client.MapWidget::getPane(G_MAP_MAP_PANE);
  pane.parentNode.appendChild(this.canvas);
  //map.getPane(G_MAP_MAP_PANE).parentNode.appendChild(this.canvas);
};

mapPVDot.prototype.redraw = function(force) {
  if (!force) return;
  var c = this.canvas, m = this.map, r = this.radius;

  /* Get the pixel locations of the crimes. */
  var pixels = [m.@com.google.gwt.maps.client.MapWidget::convertLatLngToDivPixel(new GLatLng(this.lat, this.lon))];
  //var pixels = [m.fromLatLngToDivPixel(new GLatLng(this.lat, this.lon))];
  /* Update the canvas bounds. Note: may be large. */
  function x(p) {return p.x};
  function y(p) {return p.y};
  var x = { min: pv.min(pixels, x) - r, max: pv.max(pixels, x) + r };
  var y = { min: pv.min(pixels, y) - r, max: pv.max(pixels, y) + r };
  c.style.width =  this.radius*2 + "px";
  c.style.height = this.radius*2 + "px";
  c.style.left = x.min + "px";
  c.style.top = y.min + "px";

  /* Render the visualization. */
  dotColorScale = this.dotColorScale;
  dotSizeScale = this.dotSizeScale;
  title = this.title;
  
	new pv.Panel()
	      .canvas(c)
	      .left(-x.min)
	      .top(-y.min)
	    .add(pv.Panel)
	      .data({value:this.value})
	    .add(pv.Dot)
	      .left(function() {return pixels[this.parent.index].x})
	      .top(function() {return pixels[this.parent.index].y})
	      .strokeStyle(function(x, d) {return dotColorScale[d.value]})
	      .fillStyle(function(x, d) {return dotColorScale[d.value]})
	      .size(r)
	    .anchor("center").add(pv.Label)
	  .textStyle("black")
	  .text('')
	  .title(title)
	.root.render();
};

mapPVDot.prototype.remove = function() {
	this.canvas.parentNode.removeChild( this.canvas );
};