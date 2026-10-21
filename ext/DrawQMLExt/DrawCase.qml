// The drawing window of draw_case (DrawQMLExt): a view of the Julia `Editor`. Everything it shows
// comes from Julia.background()/overlay()/panel() (JSON); every action calls one Julia function, then refresh().
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Dialogs
import jlqml

ApplicationWindow {
    id: win
    width: 1400; height: 900; visible: true
    title: "ElastoPlasm — draw a case"
    property var panel: JSON.parse(Julia.panel())
    function refresh() { panel = JSON.parse(Julia.panel()); canvas.requestPaint(); overlay.requestPaint() }
    function act(f) { f(); refresh() }                       // a Julia call that changes the editor, then the view
    onClosing: Julia.cancel()

    Rectangle {
        id: root; anchors.fill: parent; color: palette.window
        RowLayout {
            anchors.fill: parent; spacing: 0

            // ---- the drawing area, in two layers: the background (grid, domain, the points to be created) repaints
            // after edits, zoom and resizing; the overlay (zones, the polygon being drawn, the cursor) at every move
            Canvas {
                id: canvas
                Layout.fillWidth: true; Layout.fillHeight: true
                focus: true
                function seg(c, l, color, w) { c.strokeStyle = color; c.lineWidth = w; c.beginPath()
                    for (var i = 0; i < l.length; i++) { c.moveTo(l[i][0], l[i][1]); c.lineTo(l[i][2], l[i][3]) } c.stroke() }
                function path(c, p, close) { c.beginPath(); for (var i = 0; i < p.length; i++) { i == 0 ? c.moveTo(p[i][0], p[i][1]) : c.lineTo(p[i][0], p[i][1]) }
                    if (close) c.closePath() }
                onPaint: {
                    var s = JSON.parse(Julia.background(width, height))
                    var c = getContext("2d")
                    c.reset(); c.fillStyle = "white"; c.fillRect(0, 0, width, height)
                    // the picture to trace over, under everything (loaded once; painted again when it arrives)
                    if (s.backdrop) { if (!isImageLoaded(s.backdrop.url)) loadImage(s.backdrop.url)
                        else { var b = s.backdrop.rect; c.globalAlpha = s.backdrop.opacity; c.drawImage(s.backdrop.url, b[0], b[1], b[2], b[3]); c.globalAlpha = 1.0 } }
                    // axes box and ticks
                    c.strokeStyle = "#333"; c.lineWidth = 1; c.strokeRect(s.box[0], s.box[1], s.box[2], s.box[3])
                    c.fillStyle = "#333"; c.font = "11px sans-serif"
                    for (var i = 0; i < s.ticks.x.length; i++) { var t = s.ticks.x[i]; c.fillText(t.label, t.px - 8, s.box[1] + s.box[3] + 14) }
                    for (var j = 0; j < s.ticks.z.length; j++) { var u = s.ticks.z[j]; c.fillText(u.label, s.box[0] - 40, u.py + 4) }
                    c.fillText("x [m]", s.box[0] + s.box[2]/2 - 10, s.box[1] + s.box[3] + 30)
                    c.fillText("z [m]", 4, s.box[1] - 8)
                    c.font = "bold 12px sans-serif"; c.fillText(s.title, s.box[0] + s.box[2]/2 - 150, 18)
                    // grid, domain, sides
                    seg(c, s.sub, "rgba(128,128,128,0.15)", 0.5)
                    seg(c, s.el, "rgba(128,128,128,0.45)", 0.7)
                    c.strokeStyle = "black"; c.lineWidth = 1.5
                    c.strokeRect(s.domain[0], s.domain[3], s.domain[2] - s.domain[0], s.domain[1] - s.domain[3])
                    for (var k = 0; k < s.sides.length; k++) seg(c, [s.sides[k].seg], s.sides[k].color, 4)
                    // the zones' fills, then the points to be created over them
                    c.globalAlpha = 0.15
                    for (var f = 0; f < s.fills.length; f++) { path(c, s.fills[f].pts, true); c.fillStyle = s.fills[f].color; c.fill() }
                    c.globalAlpha = 1.0
                    for (var q = 0; q < s.points.length; q++) { c.fillStyle = s.points[q].color; var xy = s.points[q].xy; c.beginPath()
                        for (var n = 0; n < xy.length; n += 2) { c.moveTo(xy[n] + 2.2, xy[n + 1]); c.arc(xy[n], xy[n + 1], 2.2, 0, 2*Math.PI) }
                        c.fill() }                               // one fill per zone: fast on fine meshes
                }
                onImageLoaded: requestPaint()
                onPainted: shots.painted()
                Canvas {
                    id: overlay
                    anchors.fill: parent
                    onPaint: {
                        var s = JSON.parse(Julia.overlay(width, height))
                        var c = getContext("2d")
                        c.reset()
                        // the zones' outlines (their fills are on the background)
                        for (var z = 0; z < s.zones.length; z++) { var zn = s.zones[z]
                            canvas.path(c, zn.pts, true); c.strokeStyle = zn.color; c.lineWidth = zn.selected ? 3 : 2; c.stroke() }
                        // vertices
                        for (var y = 0; y < s.zones.length; y++) { c.fillStyle = s.zones[y].color; var hs = s.zones[y].handles
                            for (var v = 0; v < hs.length; v++) c.fillRect(hs[v][0] - 4, hs[v][1] - 4, 8, 8) }
                        if (s.picked) { c.strokeStyle = "#111"; c.lineWidth = 2; c.strokeRect(s.picked[0] - 7, s.picked[1] - 7, 14, 14) }   // its row in the panel
                        // the polygon being drawn, the cursor, the domain handle, the cursor's coordinates
                        if (s.sketch.length > 0) { c.setLineDash([6, 4]); canvas.path(c, s.sketch, false); c.strokeStyle = "black"; c.lineWidth = 1.5; c.stroke(); c.setLineDash([]) }
                        if (s.cursor) { var cx = s.cursor.xy[0], cy = s.cursor.xy[1], r = s.cursor.r; c.beginPath()
                            if (s.cursor.shape === "square") c.rect(cx - r, cy - r, 2*r, 2*r)
                            else if (s.cursor.shape === "diamond") { c.moveTo(cx, cy - r - 2); c.lineTo(cx + r + 2, cy); c.lineTo(cx, cy + r + 2); c.lineTo(cx - r - 2, cy); c.closePath() }
                            else c.arc(cx, cy, r, 0, 2*Math.PI)
                            c.strokeStyle = s.cursor.color; c.lineWidth = r > 6 ? 2.5 : 2; c.stroke() }
                        if (s.handle) { c.fillStyle = "white"; c.fillRect(s.handle[0] - 7, s.handle[1] - 7, 14, 14)
                            c.strokeStyle = "black"; c.lineWidth = 2; c.strokeRect(s.handle[0] - 7, s.handle[1] - 7, 14, 14) }
                        c.fillStyle = "#333"; c.font = "11px sans-serif"; c.fillText(s.readout, s.box[0] + s.box[2] - 260, s.box[1] + s.box[3] + 30)
                        // calibration points clicked on the picture: a cross and its number
                        c.strokeStyle = "#7a5cc7"; c.fillStyle = "#7a5cc7"; c.lineWidth = 2; c.font = "bold 12px sans-serif"
                        for (var m = 0; m < s.calib.length; m++) { var k = s.calib[m].xy; c.beginPath()
                            c.moveTo(k[0] - 8, k[1]); c.lineTo(k[0] + 8, k[1]); c.moveTo(k[0], k[1] - 8); c.lineTo(k[0], k[1] + 8); c.stroke()
                            c.fillText(s.calib[m].n, k[0] + 6, k[1] - 6) }
                        if (s.calib.length === 2) { c.setLineDash([4, 3]); c.beginPath(); c.moveTo(s.calib[0].xy[0], s.calib[0].xy[1])
                            c.lineTo(s.calib[1].xy[0], s.calib[1].xy[1]); c.stroke(); c.setLineDash([]) }      // the scale bar measured
                    }
                    onPainted: shots.painted()
                }
                MouseArea {
                    anchors.fill: parent; hoverEnabled: true; acceptedButtons: Qt.LeftButton | Qt.RightButton
                    cursorShape: win.panel.zones.drawing ? Qt.CrossCursor : Qt.ArrowCursor      // a crosshair while a zone is drawn
                    onPressed: (m) => { canvas.forceActiveFocus(); win.act(() => Julia.press(m.x, m.y, m.button === Qt.RightButton, width, height)) }
                    // moving and leaving change only the overlay; zooming both layers; none of them the panel (its checks)
                    onPositionChanged: (m) => { Julia.move(m.x, m.y, width, height); overlay.requestPaint() }
                    onReleased: win.act(() => Julia.release())
                    onExited: { Julia.leave(); overlay.requestPaint() }
                    onWheel: (w) => { Julia.wheel(w.x, w.y, w.angleDelta.y, width, height); canvas.requestPaint(); overlay.requestPaint() }
                }
                Keys.onPressed: (k) => {
                    var names = {}; names[Qt.Key_Return] = "enter"; names[Qt.Key_Enter] = "enter"; names[Qt.Key_Backspace] = "backspace"
                    names[Qt.Key_Escape] = "escape"; names[Qt.Key_Z] = "z"; names[Qt.Key_Y] = "y"
                    if (names[k.key] === undefined) return
                    Julia.key(names[k.key], (k.modifiers & (Qt.ControlModifier | Qt.MetaModifier)) !== 0, (k.modifiers & Qt.ShiftModifier) !== 0)
                    win.refresh(); k.accepted = true
                }
            }

            // ---- the side panel
            ColumnLayout {
                Layout.preferredWidth: 380; Layout.maximumWidth: 380; Layout.fillHeight: true; Layout.margins: 8
                TabBar {
                    id: tabs; Layout.fillWidth: true
                    currentIndex: win.panel.tab
                    onCurrentIndexChanged: if (currentIndex !== win.panel.tab) win.act(() => Julia.set_step(currentIndex))
                    TabButton { text: "1 Domain" } TabButton { text: "2 Zones" } TabButton { text: "3 Options" } TabButton { text: "4 Run" }
                }
                Label { text: win.panel.status; wrapMode: Text.WordWrap; Layout.fillWidth: true; color: "#444" }
                StackLayout {
                    currentIndex: win.panel.tab; Layout.fillWidth: true; Layout.fillHeight: true

                    // 1 Domain
                    ScrollView { contentWidth: availableWidth; clip: true
                        GridLayout { width: parent.width; columns: 2; columnSpacing: 8
                            Button { text: "Start from a case file…"; Layout.columnSpan: 2; onClicked: loadDialog.open() }
                            Label { text: "Domain"; font.bold: true; Layout.columnSpan: 2; Layout.topMargin: 6 }
                            Label { text: "Lx [m]" } TextField { text: win.panel.domain.Lx; Layout.fillWidth: true; onEditingFinished: win.act(() => Julia.set_field("Lx", text)) }
                            Label { text: "Lz [m]" } TextField { text: win.panel.domain.Lz; Layout.fillWidth: true; onEditingFinished: win.act(() => Julia.set_field("Lz", text)) }
                            Label { text: "element size h [m]" } TextField { text: win.panel.domain.h; Layout.fillWidth: true; onEditingFinished: win.act(() => Julia.set_field("h", text)) }
                            Label { text: "width y [m], 3D" } TextField { text: win.panel.domain.width; placeholderText: "none: 2D"; Layout.fillWidth: true
                                onEditingFinished: win.act(() => Julia.set_field("width", text)) }
                            Label { text: win.panel.domain.nel; Layout.columnSpan: 2 }
                            Label { text: "points per element" } ComboBox { model: win.panel.domain.nis; currentIndex: win.panel.domain.ni; Layout.fillWidth: true
                                onActivated: (i) => win.act(() => Julia.set_ni(i)) }
                            // a picture to draw on: load it, give its scale (its width, or a scale bar), the domain follows it
                            Label { text: "Picture to draw on"; font.bold: true; Layout.columnSpan: 2; Layout.topMargin: 8 }
                            RowLayout { Layout.columnSpan: 2
                                Button { text: win.panel.backdrop ? "Replace…" : "Picture…"; onClicked: imageDialog.open() }
                                Button { text: "Scale bar"; visible: win.panel.backdrop !== null; highlighted: win.panel.calib !== null
                                    onClicked: { win.act(() => Julia.scale_bar()); canvas.forceActiveFocus() } }
                                Button { text: "Remove"; visible: win.panel.backdrop !== null; onClicked: win.act(() => Julia.remove_backdrop()) } }
                            Label { visible: win.panel.backdrop !== null; text: "picture width [m]" }
                            TextField { visible: win.panel.backdrop !== null; text: win.panel.backdrop ? win.panel.backdrop.width : ""; Layout.fillWidth: true
                                onEditingFinished: win.act(() => Julia.set_backdrop_field("width", text)) }
                            Label { visible: win.panel.calib !== null && win.panel.calib.clicked === 2
                                text: "scale bar length [m]\n(" + (win.panel.calib ? win.panel.calib.shown : "") + " m as shown)" }
                            RowLayout { visible: win.panel.calib !== null && win.panel.calib.clicked === 2; Layout.fillWidth: true
                                TextField { text: win.panel.calib ? win.panel.calib.length : ""; Layout.fillWidth: true
                                    onEditingFinished: win.act(() => Julia.set_scale_length(text)) }
                                Button { text: "Apply"; onClicked: win.act(() => Julia.apply_scale()) } }
                            Repeater { model: win.panel.backdrop ? [["x0", "picture left x [m]"], ["z0", "picture bottom z [m]"], ["opacity", "opacity"]] : []
                                delegate: RowLayout { Layout.columnSpan: 2
                                    Label { text: modelData[1]; Layout.preferredWidth: 130 }
                                    TextField { text: win.panel.backdrop ? win.panel.backdrop[modelData[0]] : ""; Layout.fillWidth: true
                                        onEditingFinished: win.act(() => Julia.set_backdrop_field(modelData[0], text)) } } }
                            Button { visible: win.panel.backdrop !== null; text: "Fit domain to picture"; Layout.columnSpan: 2
                                onClicked: win.act(() => Julia.fit_domain()) }
                            Label { text: "Boundary conditions"; font.bold: true; Layout.columnSpan: 2; Layout.topMargin: 8 }
                            Repeater { model: win.panel.domain.bcs
                                delegate: RowLayout { Layout.columnSpan: 2
                                    Label { text: modelData.side; Layout.preferredWidth: 130 }
                                    ComboBox { model: win.panel.domain.kinds; currentIndex: modelData.index; Layout.fillWidth: true
                                        onActivated: (i) => win.act(() => Julia.set_bc(modelData.side, i)) } } }
                            Label { Layout.columnSpan: 2; Layout.fillWidth: true; wrapMode: Text.WordWrap; color: "gray"; Layout.topMargin: 8
                                text: "Drag the white handle (upper-right corner) to size the domain, or type the values; a width extrudes the drawing along y (blank: 2D). Boundaries: fixed (red), roller (blue), free (grey). Scroll to zoom." }
                        } }

                    // 2 Zones
                    ScrollView { contentWidth: availableWidth; clip: true
                        GridLayout { width: parent.width; columns: 3; columnSpacing: 8
                            Label { text: "Zones"; font.bold: true; Layout.columnSpan: 3 }
                            ComboBox { model: win.panel.zones.names; currentIndex: win.panel.zones.selected; Layout.columnSpan: 3; Layout.fillWidth: true
                                onActivated: (i) => win.act(() => Julia.select_zone(i)) }
                            Label { text: "name" } TextField { text: win.panel.zones.name; Layout.columnSpan: 2; Layout.fillWidth: true
                                onEditingFinished: win.act(() => Julia.rename_zone(text)) }
                            Repeater { model: win.panel.zones.props
                                delegate: RowLayout { Layout.columnSpan: 3
                                    Label { text: modelData.label; Layout.preferredWidth: 90 }
                                    TextField { text: modelData.value; placeholderText: "default"; Layout.preferredWidth: 110
                                        onEditingFinished: win.act(() => Julia.set_prop(modelData.key, text)) }
                                    Label { text: modelData.default; color: "gray" } } }
                            Label { text: "circle [m]"; font.bold: true; Layout.columnSpan: 3; Layout.topMargin: 6
                                visible: win.panel.zones.circle !== null }
                            RowLayout { Layout.columnSpan: 3; visible: win.panel.zones.circle !== null   // a Ball zone: centre and radius
                                Repeater { model: ["x", "z", "r"]
                                    delegate: RowLayout { Label { text: modelData }
                                        TextField { text: win.panel.zones.circle ? win.panel.zones.circle[modelData] : ""; Layout.preferredWidth: 80
                                            onEditingFinished: win.act(() => Julia.set_circle(modelData, text)) } } } }
                            Label { text: "vertex [m]"; font.bold: true; Layout.columnSpan: 3; Layout.topMargin: 6
                                visible: win.panel.zones.vertex !== null }
                            RowLayout { Layout.columnSpan: 3; visible: win.panel.zones.vertex !== null   // the vertex last clicked or dragged
                                Repeater { model: ["x", "z"]
                                    delegate: RowLayout { Label { text: modelData }
                                        TextField { text: win.panel.zones.vertex ? win.panel.zones.vertex[modelData] : ""; Layout.preferredWidth: 95
                                            onEditingFinished: win.act(() => Julia.set_picked(modelData, text)) } } }
                                Button { text: "×"; implicitWidth: 28; onClicked: win.act(() => Julia.remove_picked())
                                    ToolTip.visible: hovered; ToolTip.text: "delete this vertex" } }
                            RowLayout { Layout.columnSpan: 3
                                Button { id: newZone; text: "+ New zone ▾"; highlighted: win.panel.zones.drawing
                                    onClicked: shapes.open()
                                    Menu { id: shapes; y: newZone.height                     // the shape of the new zone
                                        Repeater { model: [["polygon", "Draw (polygon)"], ["rectangle", "Rectangle"], ["triangle", "Triangle"], ["circle", "Circle"]]
                                            MenuItem { text: modelData[1]
                                                onTriggered: { win.act(() => Julia.new_zone(modelData[0])); canvas.forceActiveFocus() } } } } } }
                            RowLayout { Layout.columnSpan: 3
                                Button { text: "↑"; onClicked: win.act(() => Julia.swap_zone(-1)) }
                                Button { text: "↓"; onClicked: win.act(() => Julia.swap_zone(1)) }
                                Button { text: "delete"; onClicked: win.act(() => Julia.delete_zone()) } }
                            RowLayout { Layout.columnSpan: 3
                                Button { text: "undo"; onClicked: win.act(() => Julia.undo()) }
                                Button { text: "redo"; onClicked: win.act(() => Julia.redo()) } }
                            RowLayout { Layout.columnSpan: 3
                                CheckBox { text: "snap grid"; checked: win.panel.zones.grid; onToggled: win.act(() => Julia.toggle("grid", checked)) }
                                CheckBox { text: "snap shapes"; checked: win.panel.zones.shapes; onToggled: win.act(() => Julia.toggle("shapes", checked)) }
                                CheckBox { text: "points"; checked: win.panel.zones.points; onToggled: win.act(() => Julia.toggle("points", checked)) } }
                            Label { Layout.columnSpan: 3; Layout.fillWidth: true; wrapMode: Text.WordWrap; color: "gray"
                                text: "+ New zone: draw a polygon (click vertices; first vertex or Enter: close), a rectangle (two corners), a triangle (three corners) or a circle (centre, then a point on it). Drag a vertex (its coordinates are then typed above), a circle's centre or its handle; click an edge to insert a vertex, right-click one to delete it; Backspace/Esc: undo the last vertex / cancel; Cmd/Ctrl+Z: undo. Later zones win; blank = default." }
                            Label { Layout.columnSpan: 3; Layout.fillWidth: true; wrapMode: Text.WordWrap; text: win.panel.zones.summary }
                        } }

                    // 3 Options
                    ScrollView { contentWidth: availableWidth; clip: true
                        ColumnLayout { width: parent.width
                            Label { text: "Solver options"; font.bold: true }
                            ComboBox { model: win.panel.options.sections; currentIndex: win.panel.options.section; Layout.fillWidth: true
                                onActivated: (i) => win.act(() => Julia.set_section(i)) }
                            Label { text: "Ticked fields override the default; picking or typing a value ticks the box."; color: "gray"; wrapMode: Text.WordWrap; Layout.fillWidth: true }
                            Repeater { model: win.panel.options.rows
                                delegate: ColumnLayout { Layout.fillWidth: true
                                    property var row: modelData
                                    RowLayout { Layout.fillWidth: true
                                        CheckBox { checked: row.checked; onToggled: win.act(() => Julia.tick_option(row.key, checked)) }
                                        Label { text: row.label; font.bold: row.kind === "list"; Layout.preferredWidth: 110 }
                                        ComboBox { visible: row.kind === "menu"; model: row.values; currentIndex: row.index; Layout.fillWidth: true
                                            onActivated: (i) => win.act(() => Julia.set_option(row.key, i)) }
                                        TextField { visible: row.kind === "numeric"; text: row.text; Layout.fillWidth: true
                                            onEditingFinished: win.act(() => Julia.type_option(row.key, text)) }
                                        ComboBox { visible: row.kind === "numeric"; model: row.values; displayText: "presets"; Layout.preferredWidth: 90
                                            onActivated: (i) => win.act(() => Julia.set_option(row.key, i)) } }
                                    GridLayout { visible: row.kind === "list"; columns: 2; Layout.leftMargin: 30
                                        Repeater { model: row.kind === "list" ? row.values : []
                                            delegate: CheckBox { text: modelData; checked: row.picks[index]
                                                onToggled: win.act(() => Julia.tick_value(row.key, index, checked)) } } } } }
                        } }

                    // 4 Run
                    ScrollView { contentWidth: availableWidth; clip: true
                        GridLayout { width: parent.width; columns: 2; columnSpacing: 8
                            Label { text: "Run"; font.bold: true; Layout.columnSpan: 2 }
                            ComboBox { model: win.panel.run.workflows; currentIndex: win.panel.run.workflow; Layout.columnSpan: 2; Layout.fillWidth: true
                                onActivated: (i) => win.act(() => Julia.set_workflow(i)) }
                            Repeater { model: win.panel.run.fields
                                delegate: RowLayout { Layout.columnSpan: 2
                                    Label { text: modelData.label; Layout.preferredWidth: 160 }
                                    TextField { text: modelData.value; Layout.fillWidth: true; onEditingFinished: win.act(() => Julia.set_run(modelData.key, text)) } } }
                            Label { text: "Checks"; font.bold: true; Layout.columnSpan: 2; Layout.topMargin: 8 }
                            Label { text: win.panel.run.errors; color: "firebrick"; visible: text !== ""; wrapMode: Text.WordWrap; Layout.columnSpan: 2; Layout.fillWidth: true }
                            Label { text: win.panel.run.warnings; color: "darkorange"; visible: text !== ""; wrapMode: Text.WordWrap; Layout.columnSpan: 2; Layout.fillWidth: true }
                            Label { text: win.panel.run.info; wrapMode: Text.WordWrap; Layout.columnSpan: 2; Layout.fillWidth: true }
                            Label { text: "folder: " + win.panel.run.folder; visible: win.panel.run.folder !== ""; wrapMode: Text.WrapAnywhere
                                color: "gray"; Layout.columnSpan: 2; Layout.fillWidth: true; Layout.topMargin: 8 }
                            RowLayout { Layout.columnSpan: 2; Layout.topMargin: 8
                                Button { text: "Create…"; highlighted: true; onClicked: { if (Julia.ready()) folderDialog.open(); win.refresh() } }
                                Button { text: "Cancel"; onClicked: { Julia.cancel(); Qt.quit() } } }
                            Label { text: "Create…: choose or make the case's folder. The case is written there, and its runs go there: execute!(setup_simulation(case))."
                                color: "gray"; wrapMode: Text.WordWrap; Layout.columnSpan: 2; Layout.fillWidth: true }
                        } }
                }
            }
        }
    }

    Image { visible: false; source: win.panel.backdrop ? win.panel.backdrop.url : ""            // reports the picture's size
        onStatusChanged: if (status === Image.Ready && !win.panel.backdrop.loaded) win.act(() => Julia.backdrop_size(implicitWidth, implicitHeight)) }
    FileDialog { id: imageDialog; fileMode: FileDialog.OpenFile; currentFolder: startdir
        nameFilters: ["Images (*.png *.jpg *.jpeg *.bmp *.gif *.tif *.tiff)"]
        onAccepted: win.act(() => Julia.set_backdrop(selectedFile.toString())) }
    FolderDialog { id: folderDialog; currentFolder: startdir                               // Create…: the case's folder
        onAccepted: { if (Julia.create(selectedFolder.toString())) Qt.quit(); else win.refresh() } }
    FileDialog { id: loadDialog; fileMode: FileDialog.OpenFile; nameFilters: ["Case files (*.toml)"]; currentFolder: startdir
        onAccepted: win.act(() => Julia.load(selectedFile.toString())) }

    // screenshots for tests and docs: save each listed step as a PNG, then close
    QtObject {
        id: shots
        property int i: 0
        property int waiting: 0                                  // layers still to paint before the grab
        function next() { if (i >= shotsteps.length) { Qt.quit(); return }
            Julia.set_step(shottabs[i]); waiting = 2; win.refresh() }
        function painted() { if (waiting > 0 && --waiting === 0) grabNext() }
        function grabNext() { root.grabToImage(function(r) { r.saveToFile(shotdir + "/step" + shotsteps[shots.i] + ".png"); shots.i++; shots.next() }) }
    }
    Timer { interval: 400; running: shotsteps.length > 0; onTriggered: shots.next() }
}
