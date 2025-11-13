import QtQuick
import QtQuick.Controls

ApplicationWindow {
    visible: true
    width: 600
    height: 600
    title: "system palette test"

    SystemPalette { id: sys; colorGroup: SystemPalette.Active }

    Component.onCompleted: {
        console.log("Window:", sys.window)
        console.log("Text:", sys.text)
        console.log("Base:", sys.base)
        console.log("Button:", sys.button)
        console.log("Highlight:", sys.highlight)
        console.log("Highlighted Text:", sys.highlightedText)
    }

    Rectangle {
        anchors.fill: parent
        color: "#515151"
        Column {
            anchors.centerIn: parent
            spacing: 8
            Repeater {
                model: [
                    { name: "Window", color: sys.window },
                    { name: "WindowText", color: sys.windowText },
                    { name: "Base", color: sys.base },
                    { name: "AlternateBase", color: sys.alternateBase },
                    { name: "Text", color: sys.text },
                    { name: "Button", color: sys.button },
                    { name: "ButtonText", color: sys.buttonText },
                    { name: "BrightText", color: sys.brightText },
                    { name: "Light", color: sys.light },
                    { name: "Midlight", color: sys.midlight },
                    { name: "Mid", color: sys.mid },
                    { name: "Dark", color: sys.dark },
                    { name: "Shadow", color: sys.shadow },
                    { name: "Highlight", color: sys.highlight },
                    { name: "HighlightedText", color: sys.highlightedText }
                ]
                delegate: Row {
                    spacing: 10
                    Rectangle {
                        width: 40; height: 20
                        color: modelData.color
                        border.width: 1
                        border.color: "black"
                    }
                    Text { text: modelData.name; color: "#FFFFFF" }
                }
            }
        }
    }
}
