var map = null;
var overlayGroups = new Object();
var defaultColor = "#FF0000";
var algorithmParameters = [];
var clickGroup = -1;
var clickPos = null;
var sourcePos = null;
var targetPos = null;
var groupCount = 0;
var blockFlag = false;
var openStreetMapTypeOptions = {
    getTileUrl: function(coord, zoom) {
        var normalizedCoord = getNormalizedCoord(coord, zoom);
        if (!normalizedCoord) {
            return null;
        }
        var bound = Math.pow(2, zoom);
        return "http://tile.openstreetmap.org/" + zoom + "/" + normalizedCoord.x + "/" + (bound - normalizedCoord.y - 1) + ".png";
    },
    tileSize: new google.maps.Size(256, 256),
    maxZoom: 18,
    minZoom: 0,
    name: "OpenStreetMap",
    alt: "OpenStreetMap"
};
var openStreetMapType = new google.maps.ImageMapType(openStreetMapTypeOptions);

window.onpopstate = function(event) {    
    // Check for PermaLink and execute if neccessary
    var tokens = document.URL.split("index.html?");
    if (tokens.length > 1) {
    	ReRun(tokens[1]);
    }
};

window.setTimeout(manipulateMapElements , 100);
function manipulateMapElements() {
	// This is a bad hack!
	// But I realy want the Map-Button to disply 'GoogleMap'
    var osm = getElementByXpath("//div[@title='OpenStreetMap']");
    var gm = getElementByXpath("//div[@title='Show street map']") || getElementByXpath("//div[@title='Stadtplan anzeigen']");
    if (osm) {
        osm.style.minWidth = "20px";
        gm.style.minWidth = "20px";
        gm.innerHTML = "GoogleMap";
    } else {
    	window.setTimeout(manipulateMapElements, 100);
    }
}

function getElementByXpath(path) {
    return document.evaluate(path, document, null, 9, null).singleNodeValue;
};

// Initializes the map.
function InitializeMap() {
	
	// Set up right-click-menus
    HideContextMenu();
    HideEntryMenu();
    HideFunctionCanvas();

    // Set some default map options.
    var kitLatLng = new google.maps.LatLng(49.013871, 8.419561);
    var initialMapOptions = {
        zoom: 8,
        center: kitLatLng,
        streetViewControl: false,
        disableDefaultUI: true,
        mapTypeControl: true,
        mapTypeControlOptions: {
            //mapTypeIds: ["OpenStreetMap", google.maps.MapTypeId.ROADMAP]
    		mapTypeIds: ["OpenStreetMap", google.maps.MapTypeId.ROADMAP, google.maps.MapTypeId.TERRAIN]
        }
    };
    
    // Load the map into the canvas.
    map = new google.maps.Map(document.getElementById("mapCanvas"), initialMapOptions);
    map.mapTypes.set('OpenStreetMap', openStreetMapType);
    map.setMapTypeId('OpenStreetMap');
    
    // Shift to bounding box.
    new Ajax.Request('/bounding_box', {
        onSuccess: function(response) {
            var data = response.responseText.evalJSON();
            var bounds = new google.maps.LatLngBounds(new google.maps.LatLng(data.Min.Latitude, data.Min.Longitude), new google.maps.LatLng(data.Max.Latitude, data.Max.Longitude));
            map.fitBounds(bounds);
        }
    });
       
    // Add listener for clicking onto the map (place source and target markers)
    google.maps.event.addListener(map, 'click', function(e) {
        clickPos = e.latLng;
        HideContextMenu();    
        SetSource()
        Run();
    });
    google.maps.event.addListener(map, 'rightclick', function(e) {
        ShowContextMenu(e.latLng);
        e.stop();
    });
    google.maps.event.addListener(map, 'mousemove', function(e) {
        HideContextMenu();
        HideEntryMenu();
    });
    
    // Retrieve algorithms.
    new Ajax.Request('/algorithms_information', {
        onSuccess: function(response) {
            var data = response.responseText.evalJSON();
            for (var i = 0; i < data.NumberOfAlgorithms; ++i) {
                $('select_algorithm').options.add(new Option(data.Algorithms[i].Name, data.Algorithms[i].Key));
                if (data.Algorithms[i].Parameter) {
                    algorithmParameters[data.Algorithms[i].Key] = data.Algorithms[i].Parameter;
                    algorithmParameters[data.Algorithms[i].Key].sort(function(a, b){return a.Index < b.Index});
                }    
            }
            LoadParameter();
                        
            // Check for PermaLink and execute if neccessary
            var tokens = document.URL.split("index.html?");
            if (tokens.length > 1) {
            	ReRun(tokens[1]);
            }
        }
    });

}

function SetSource() {
    sourcePos = clickPos;
}

function SetTarget() {
    targetPos = clickPos;
}

function ShowContextMenu(latlng) {
    clickPos = latlng;
    var contextMenu = $('contextMenu');
    var pos = LatLngToPixel(latlng);
    contextMenu.style.left = (pos.x - 10) + 'px';
    contextMenu.style.top = (pos.y - 10) + 'px';    
    contextMenu.style.display = '';
}

function HideContextMenu() {
    var contextMenu = $('contextMenu');
    contextMenu.style.display = 'none';    
}

function ShowEntryMenu(e) {
    var entryMenu = $('entryMenu');
    entryMenu.style.left = (e.clientX - 10) + 'px';
    entryMenu.style.top = (e.clientY - 10) + 'px';    
    entryMenu.style.display = '';
    return false;
}

function HideEntryMenu() {
    var entryMenu = $('entryMenu');
    entryMenu.style.display = 'none';    
}

function ShowFunctionCanvas() {
    var entryMenu = $('functionCanvas');  
    entryMenu.style.display = '';
}

function HideFunctionCanvas() {
    var entryMenu = $('functionCanvas');
    entryMenu.style.display = 'none';    
}

function ResetQuery() {
    sourcePos = null;
    targetPos = null;
    ClearOverlays();
    window.history.pushState({}, "Map Server", "/");
    HideContextMenu();
}

function Block() {
	blockFlag = true;
    $('all').setAttribute("class", "block");
}

function Unblock() {
	blockFlag = false;
    $('all').setAttribute("class", "unblock");
}

function Run() {
    if (sourcePos != null || targetPos != null) {
        ClearOverlays();
        var key = $('select_algorithm').value
        var ajaxRequest = "/" + key;
        if (sourcePos != null) {
            ajaxRequest = ajaxRequest + "?SourceLat=" + sourcePos.lat();
            ajaxRequest = ajaxRequest + "&SourceLon=" + sourcePos.lng();
        }
        if (targetPos != null) {
            var delimiter = ((sourcePos == null) ? "?" : "&");
            ajaxRequest = ajaxRequest + delimiter + "TargetLat=" + targetPos.lat();
            ajaxRequest = ajaxRequest + "&TargetLon=" + targetPos.lng();
        }
        if (algorithmParameters[key]) {
            for (var i = 0; i < algorithmParameters[key].length; i++) {
                if (algorithmParameters[key][i].Type == "bool") {
                    ajaxRequest = ajaxRequest + "&" + algorithmParameters[key][i].Name + "=" + $(algorithmParameters[key][i].Name).checked;
                } else if (algorithmParameters[key][i].Type == "time") {
                    ajaxRequest = ajaxRequest + "&" + algorithmParameters[key][i].Name + "=" + TimeValueToSeconds($(algorithmParameters[key][i].Name).value);
                } else if (algorithmParameters[key][i].Type == "string") {
                    ajaxRequest = ajaxRequest + "&" + algorithmParameters[key][i].Name + "=" + Encode($(algorithmParameters[key][i].Name).value);
                } else {
                    ajaxRequest = ajaxRequest + "&" + algorithmParameters[key][i].Name + "=" + $(algorithmParameters[key][i].Name).value;
                }
            }
        }
        window.history.pushState({request: ajaxRequest}, "Map Server", "index.html?" + ajaxRequest);
        if (!blockFlag) {
            Block();
            new Ajax.Request(ajaxRequest, {
                onSuccess: function(response) {
                    DisplayResult(response.responseText.evalJSON());
                    Unblock();
                },
                onComplete: function(response) {
        			if (response.status == 0) {
        				DisplayResult({"ResultSet":[{"Info":"Timeout", "Color":"FF0000", "Pathes":[], "Nodes":[]}]});
        				Unblock();
        			}
        		}
            });
        } 
    }
}

function ReRun(ajaxRequest) {
	var tokens = ajaxRequest.split("?", 2);
	var key = tokens[0].substr(1);

	$('select_algorithm').value = key;
	LoadParameter();
	
	var params = tokens[1].split("&");
	var SourceLat = 1000;
	var SourceLon = 1000;
	var TargetLat = 1000;
	var TargetLon = 1000;
	for (var i = 0; i < params.length; i++) {
		tokens = params[i].split("=", 2);
		if (tokens[0] == "SourceLat") {
			SourceLat = tokens[1];
		} else if (tokens[0] == "SourceLon") {
			SourceLon = tokens[1];
		} else if (tokens[0] == "TargetLat") {
			TargetLat = tokens[1];
		} else if (tokens[0] == "TargetLon") {
			TargetLon = tokens[1];
		} else {
		    if (algorithmParameters[key]) {
		        for (var j = 0; j < algorithmParameters[key].length; j++) {
		        	if (tokens[0] == algorithmParameters[key][j].Name) {
		        		if (algorithmParameters[key][j].Type == "bool") {
			        		if (tokens[1] == "true") {
			        			$(algorithmParameters[key][j].Name).checked = true;
			        		} else {
			        			$(algorithmParameters[key][j].Name).checked = false;
			        		}
		        		} else if (algorithmParameters[key][j].Type == "time") {
			        		$(algorithmParameters[key][j].Name).value = SecondsToTimeValue(tokens[1]);
		        		} else {
			        		$(algorithmParameters[key][j].Name).value = tokens[1];
		        		}
		        	}
		        }
		    }			
		}
	}
	if (SourceLat < 1000 && SourceLon < 1000) {
	    sourcePos = new google.maps.LatLng(SourceLat, SourceLon);
	}
	if (TargetLat < 1000 && TargetLon < 1000) {
	    targetPos = new google.maps.LatLng(TargetLat, TargetLon);
	}

    ClearOverlays();
    if (!blockFlag) {
        Block();
        new Ajax.Request(ajaxRequest, {
            onSuccess: function(response) {
                DisplayResult(response.responseText.evalJSON());
                Unblock();
            },
            onComplete: function(response) {
    			if (response.status == 0) {
    				DisplayResult({"ResultSet":[{"Info":"Timeout", "Color":"FF0000", "Pathes":[], "Nodes":[]}]});
    				Unblock();
    			}
    		}
        }); 
    }
}

function Encode(str) {
    return str.replace(/#/g, "%23");
}

function ToggleAll() {
    var state;
    for (var group in overlayGroups) {
    	state = !overlayGroups[group][0].getVisible();
    	break;
    }
    for (var group in overlayGroups) {
        for (var i = 0; i < overlayGroups[group].length; i++) {
            overlayGroups[group][i].setVisible(state);
        }
    }
}

function ToggleGroup() {
    if (overlayGroups[clickGroup]) {
        for (var i = 0; i < overlayGroups[clickGroup].length; i++) {
            overlayGroups[clickGroup][i].setVisible(!overlayGroups[clickGroup][i].getVisible());
        }
    }
}

function PinGroup() {
	$('tr_' + clickGroup).setAttribute('pinned', true);
}

function DeleteGroup() {
	for (var i = 0; i < overlayGroups[clickGroup].length; i++) {
		overlayGroups[clickGroup][i].setMap(null);
	}
	$('infoArea').removeChild($('tr_' + clickGroup));
	delete overlayGroups[clickGroup];
}

function SetCurrentGroup(i) {
	clickGroup = i;
}

function LoadParameter() {
    viaNode = null;
    viaRegion = null;
    var htmlString = "";
    var key = $('select_algorithm').value
    if (algorithmParameters[key]) {
        for (var i = 0; i < algorithmParameters[key].length; i++) {
            var onChangeString = "";
            if (algorithmParameters[key][i].RecalculateOnChange) {
            	if ((algorithmParameters[key][i].Type == "time")) {
                    onChangeString = " onblur='Run();'";
            	} else {
                    onChangeString = " onchange='Run();'";
            	}
            }
            var defaultValue = algorithmParameters[key][i].DefaultValue
            var oldInput = $(algorithmParameters[key][i].Name);
            if (oldInput) {
            	if (algorithmParameters[key][i].Type == "bool") {
                	defaultValue = "" + oldInput.checked;
                } else if (algorithmParameters[key][i].Type == "time") {
                	defaultValue = TimeValueToSeconds(oldInput.value);
                } else {
                	defaultValue = oldInput.value;
                }
            }
            if (algorithmParameters[key][i].Type == "bool") {
                if (defaultValue == "true") {
                    htmlString = htmlString + "<p><input id='" + algorithmParameters[key][i].Name + "'" + onChangeString + " type='checkbox' name='boxes' checked='checked' style='position: relative; top: 3px; left: -1px;'> " + algorithmParameters[key][i].Display + "</p>";
                } else {
                    htmlString = htmlString + "<p><input id='" + algorithmParameters[key][i].Name + "'" + onChangeString + " type='checkbox' name='boxes' style='position: relative; top: 3px; left: -1px;'> " + algorithmParameters[key][i].Display + "</p>";
                }
            } else if (algorithmParameters[key][i].Type == "string") {
                htmlString = "<p>" + algorithmParameters[key][i].Display + ":<br/><textarea wrap='off' cols='25' rows='10' id='" + algorithmParameters[key][i].Name + "'" + onChangeString + ">" + defaultValue + "</textarea></p>" + htmlString;
            } else if (algorithmParameters[key][i].Type == "time") {
            	var timeInput = BuildTimeInput(algorithmParameters[key][i], defaultValue, onChangeString)
                htmlString = "<p>" + algorithmParameters[key][i].Display + ":<br/>" + timeInput + "</p>" + htmlString;
            } else {
                var limits = "";
                if (algorithmParameters[key][i].MinValue) {
                	limits = limits + " min='" + algorithmParameters[key][i].MinValue + "'";
                }
                if (algorithmParameters[key][i].MaxValue) {
                	limits = limits + " max='" + algorithmParameters[key][i].MaxValue + "'";
                }
                if (algorithmParameters[key][i].Step) {
                	limits = limits + " step='" + algorithmParameters[key][i].Step + "'";
                }
                htmlString = "<p>" + algorithmParameters[key][i].Display + ":<br/><input id='" + algorithmParameters[key][i].Name + "'" + onChangeString + " value='" + defaultValue + "' type='number'" + limits + " size='30' maxlength='30'></input></p>" + htmlString;
            }
        }
    }
    $('parameterArea').innerHTML = htmlString;
    ClearOverlays();
}

function DisplayResult(data) {
	var oldInfos = $('infoArea').innerHTML;
    infoHTML = "";
    if (!data.ResultSet) return;
    if (data.ResultSet.length > 2) {
        infoHTML = infoHTML + "<tr id='tr_top' onclick='ToggleAll();'><td class='ColorMarker' style='background-color: #222222;'></td>";
        infoHTML = infoHTML + "<td class='TextMarker'>Toggle All</td></tr>";
    }
    for (var Group = 0; Group < data.ResultSet.length; Group++) {
    	var GroupID = Group + groupCount;
        infoHTML = infoHTML + "<tr id='tr_" + GroupID + "' onclick='ToggleGroup();' onmouseover='SetCurrentGroup(" + GroupID + ")'><td class='ColorMarker' style='background-color: #" + data.ResultSet[Group].Color + ";'></td>";
        infoHTML = infoHTML + "<td class='TextMarker'>" + data.ResultSet[Group].Info + "</td></tr>";
        overlayGroups[GroupID] = [];
        if (data.ResultSet[Group].Pathes.length > 0) {
            DrawPathes(data.ResultSet[Group].Pathes, defaultColor, GroupID);
        }
        if (data.ResultSet[Group].Nodes.length > 0) {
            DrawNodes(data.ResultSet[Group].Nodes, GroupID);
        }
        if (data.ResultSet[Group].Parameter) {
        	SetParameter(data.ResultSet[Group].Parameter);
        }
    }
    if (data.ResultSet.length == 0) {
        infoHTML = infoHTML + "<tr id='tr_top'><td class='ColorMarker' style='background-color: #222222;'></td>";
        infoHTML = infoHTML + "<td class='TextMarker'>No Result</td></tr>";
    }
    $('infoArea').innerHTML = infoHTML + oldInfos;
    for (var group in overlayGroups) {
        $('tr_' + group).oncontextmenu = ShowEntryMenu;
    }
    groupCount += data.ResultSet.length;
}

function ClearOverlays() {
	if ($('tr_top')) {
    	$('infoArea').removeChild($('tr_top'));
	}
    for (var Group in overlayGroups) {
    	var info = $('tr_' + Group);
    	if (!info.getAttribute('pinned')) {
        	for (var i = 0; i < overlayGroups[Group].length; i++) {
        		overlayGroups[Group][i].setMap(null);
        	}
        	$('infoArea').removeChild(info);
        	delete overlayGroups[Group];
    	}
    }
}

function ParameterType(name) {
    var key = $('select_algorithm').value
    if (algorithmParameters[key]) {
        for (var i = 0; i < algorithmParameters[key].length; i++) {
        	if (algorithmParameters[key][i].Name == name) {
        		return algorithmParameters[key][i].Type;
        	}
        }
    }
    return "";
}

function SetParameter(parameter) {
    var key = $('select_algorithm').value
    if (algorithmParameters[key]) {
        for (var i = 0; i < algorithmParameters[key].length; i++) {
        	if (parameter[algorithmParameters[key][i].Name]) {
            	if (algorithmParameters[key][i].Type == "bool") {
                    if (parameter[algorithmParameters[key][i].Name] == "true") {
                    	$(algorithmParameters[key][i].Name).checked = "checked";
                    } else {
                    	$(algorithmParameters[key][i].Name).removeAttribute("checked");
                    }
                } else if (algorithmParameters[key][i].Type == "time") {
                	$(algorithmParameters[key][i].Name).value = SecondsToTimeValue(parameter[algorithmParameters[key][i].Name]);
                } else {
                	$(algorithmParameters[key][i].Name).value = parameter[algorithmParameters[key][i].Name];
                }
        	}
        }
    }  
}

function SecondsToTimeValue(seconds) {
	if (HasTimeInput()) {
		seconds = Number(seconds);
		var h = Math.floor(seconds / 3600);
		var m = Math.floor((seconds % 3600) / 60);
		var s = Math.floor(seconds % 60);
		return (h < 10 ? "0" : "") + h + ":" + (m < 10 ? "0" : "") + m + ":" + (s < 10 ? "0" : "") + s;
	} else {
		return seconds;
	}
}

function TimeValueToSeconds(timeValue) {
	if (HasTimeInput()) {
		var tokens = timeValue.split(':');
		var seconds = 0;
		seconds += tokens[0] ? ((+tokens[0]) * 3600) : 0;
		seconds += tokens[1] ? ((+tokens[1]) * 60) : 0;
		seconds += tokens[2] ? (+tokens[2]) : 0;
		return seconds;	
	} else {
		return timeValue;
	}
}

function BuildTimeInput(parameter, defaultValue, onChangeString) {
	var step = parameter.Step ? parameter.Step : "1";
	var timeValue = SecondsToTimeValue(defaultValue);
	if (HasTimeInput()) {
		return "<input id='" + parameter.Name + "'" + onChangeString + " value='" + timeValue + "' type='time' step='" + step + "' size='30' maxlength='30'></input>"
	} else {
		return "<input id='" + parameter.Name + "'" + onChangeString + " value='" + timeValue + "' type='number' step='" + step + "' min='0' size='30' maxlength='30'></input>"
	}
}

function HasTimeInput() {
    var input = document.createElement('input');
    input.setAttribute('type','time');
    var notADateValue = 'not-a-date';
    input.setAttribute('value', notADateValue); 
    return (input.value !== notADateValue);
}

function DrawNodes(nodes, Group) {
    for (var i = 0; i < nodes.length; i++) {
        var options = {
            position: new google.maps.LatLng(nodes[i].Pos.Lat, nodes[i].Pos.Lon),
            title: (nodes[i].Pos.Lat + ", " + nodes[i].Pos.Lon),
            visible: true,
            map: map
        };
        if (nodes[i].Info.length > 0) {
        	if (nodes[i].Info.indexOf('.png', nodes[i].Info.length - 4) !== -1) {
            	var image = {url: nodes[i].Info, size: new google.maps.Size(20, 32), anchor: new google.maps.Point(9, 9)};
            	options.icon = image; 
        	} else {
                options.icon = ("http://chart.apis.google.com/chart?chst=d_map_pin_letter&chld=" + nodes[i].Info.charAt(0) + "|" + nodes[i].Color + "|000000");  
        	}
        }
        var marker = new google.maps.Marker(options)
        overlayGroups[Group].push(marker);
    }    
}

function DrawPathes(pathes, color, Group) {
    for (var i = 0; i < pathes.length; i++) {
        if (pathes[i].Color) {
            DrawPath(pathes[i].Nodes, "#" + pathes[i].Color, pathes[i], Group);
        } else {
            DrawPath(pathes[i].Nodes, color, pathes[i], Group);
        }
    }
}

function DrawPath(nodes, color, path, Group) {
    var line = new google.maps.Polyline({
        path: PolylineToPath(nodes),
        visible: true,
        strokeColor: color,
        strokeOpacity: 0.6,
        strokeWeight: 6,
        map: map
    });
    overlayGroups[Group].push(line);
}

function PolylineToPath(nodes) {
    var path = new Array(nodes.length);    
    for (var j = 0; j < nodes.length; j++) {
        path[j] = new google.maps.LatLng(nodes[j].Lat, nodes[j].Lon);
    }    
    return path;
}

function getNormalizedCoord(coord, zoom) {
    var y = coord.y;
    var x = coord.x;

    var tileRange = 1 << zoom;
    y = tileRange - y - 1;
    
    if (y < 0 || y >= tileRange) {
        return null;
    }
    if (x < 0 || x >= tileRange) {
        x = (x % tileRange + tileRange) % tileRange;
    }

    return {x: x, y: y};
}

function LatLngToPixel(latlng) {
    var scale = Math.pow(2, map.getZoom());
    var nw = new google.maps.LatLng(
        map.getBounds().getNorthEast().lat(),
        map.getBounds().getSouthWest().lng()
    );
    var worldCoordinateNW = map.getProjection().fromLatLngToPoint(nw);
    var worldCoordinate = map.getProjection().fromLatLngToPoint(latlng);
    var caurrentLatLngOffset = new google.maps.Point(
        Math.floor((worldCoordinate.x - worldCoordinateNW.x) * scale),
        Math.floor((worldCoordinate.y - worldCoordinateNW.y) * scale)
    );
    return caurrentLatLngOffset;
}
