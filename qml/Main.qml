import QtQuick 2.15
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.VectorImage
import QtQuick.Dialogs

ApplicationWindow {
    width: 1000
    height: 700
    visible: true
    color: "#1A1111"
    title: "QSynthesizer"

    FontLoader { id: boxicons; source: "assets/fonts/boxicons.ttf" }
    SystemPalette { id: sys; colorGroup: SystemPalette.Active }

    GridLayout {
        anchors.fill: parent
        columns: 2
        rows: 2
        columnSpacing: 0
        rowSpacing: 0

        anchors.bottom: parent.bottom
        anchors.bottomMargin: 15
        anchors.right: parent.right
        anchors.rightMargin: 15

        // top bar for top projucers
        Rectangle {
            Layout.row: 0
            Layout.column: 0
            Layout.columnSpan: 2
            Layout.fillWidth: true
            Layout.preferredHeight: 48
            color: sys.base
            radius: 0

            Text {
                text: "QSynthesizer"
                anchors.centerIn: parent
                color: sys.windowText
                font.pixelSize: 20
                font.bold: true
                font.family: "Inter"
            }

            Rectangle {
                width: 32; height: 32
                anchors.right: parent.right
                anchors.rightMargin: 12
                anchors.verticalCenter: parent.verticalCenter
                radius: 12
                color: sys.base
                Text {
                    text: "×"
                    anchors.centerIn: parent
                    color: sys.windowText
                    font.pixelSize: 32
                }
                MouseArea {
                    anchors.fill: parent
                    cursorShape: Qt.PointingHandCursor
                    onClicked: Qt.quit()
                }
            }
        }

        // sidebar
        Rectangle {
            Layout.row: 1
            Layout.column: 0
            Layout.fillHeight: true
            Layout.preferredWidth: 64
            color: sys.base
            radius: 0

            ColumnLayout {
                anchors.horizontalCenter: parent.horizontalCenter
                anchors.bottom: parent.bottom
                anchors.top: parent.top
                spacing: 12

                Rectangle {
                    Layout.alignment: Qt.AlignTop
                    width: 48; height: 48; radius: 12
                    color: sys.highlightedText
                    Text {
                        anchors.centerIn: parent
                        text: "+"
                        color: sys.windowText
                        font.pixelSize: 32
                    }
                    MouseArea {
                        anchors.fill: parent
                        cursorShape: Qt.PointingHandCursor
                        onClicked: fileDialog.open()
                    }
                    FileDialog {
                            id: fileDialog
                            currentFolder: StandardPaths.standardLocations(StandardPaths.PicturesLocation)[0]
                    }
                }
                Rectangle {
                    width: 48; height: 48; radius: 12
                    color: sys.base
                    Layout.alignment: Qt.AlignTop
                    Text {
                        anchors.centerIn: parent
                        font.pixelSize: 32
                        font.family: boxicons.font.family
                        text: "\uEDCB"
                        color: sys.windowText
                    }
                    MouseArea {
                        anchors.fill: parent
                        cursorShape: Qt.PointingHandCursor
                    }
                }

                Item {
                    Layout.fillHeight: true
                }

                Rectangle {
                    width: 48; height: 48; radius: 12
                    color: sys.base
                    Layout.alignment: Qt.AlignBottom
                    Text {
                        anchors.centerIn: parent
                        font.family: boxicons.font.family
                        text: "\uED52"
                        color: sys.windowText
                        font.pixelSize: 32
                    }
                }
                MouseArea {
                    anchors.fill: parent
                    cursorShape: Qt.PointingHandCursor
                }
            }
        }

        // main area
        Rectangle {
            Layout.row: 1
            Layout.column: 1
            Layout.fillWidth: true
            Layout.fillHeight: true
            color: sys.mid
            radius: 16
            antialiasing: true
            clip: true

            ColumnLayout {
                anchors.fill: parent
                anchors.margins: 16
                spacing: 10

                // waveform box
                Rectangle {
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    radius: 12
                    color: sys.mid
                    gradient: Gradient {
                        orientation: Gradient.Horizontal
                        GradientStop { position: 0.0; color: sys.mid }
                        GradientStop { position: 1.0; color: sys.light }
                    }

                    Text {
                        anchors.centerIn: parent
                        text: "waveform or something"
                        color: sys.text
                    }
                }

                // bottom buttons
                RowLayout {
                    Text {
                        font.family: boxicons.font.family
                        text: "\uED70"
                        font.pixelSize: 32
                        font.bold: true
                        color: sys.text
                    }

                    Text {
                        text: " Synthesize"
                        font.family: "Inter"
                        font.bold: true
                        font.pixelSize: 24
                        color: sys.text
                    }
                }

                RowLayout {
                    spacing: 8

                    ButtonGroup { id: waveGroup }

                    Rectangle {
                        id: sineButton
                        radius: 16
                        height: 32
                        color: sys.light
                        Layout.preferredWidth: 80

                        RowLayout {
                            VectorImage {
                                preferredRendererType: VectorImage.CurveRenderer
                                source: "assets/svg/sine.svg"
                                scale: 0.7
                            }

                            Text {
                                anchors.centerIn: parent
                                text: "    Sine"
                                color: sys.text
                                font.pixelSize: 16
                                font.family: "Inter"
                            }
                        }

                        MouseArea {
                            anchors.fill: parent
                            onClicked: waveGroup.checkedButton = sineButton
                            cursorShape: Qt.PointingHandCursor
                        }
                    }

                    Rectangle {
                        id: squareButton
                        radius: 16
                        height: 32
                        color: sys.light
                        Layout.preferredWidth: 100

                        RowLayout {
                            VectorImage {
                                preferredRendererType: VectorImage.CurveRenderer
                                source: "assets/svg/square.svg"
                                scale: 0.7
                            }

                            Text {
                                anchors.centerIn: parent
                                text: "    Square"
                                color: sys.text
                                font.pixelSize: 16
                                font.family: "Inter"
                            }
                        }

                        MouseArea {
                            anchors.fill: parent
                            onClicked: waveGroup.checkedButton = squareButton
                            cursorShape: Qt.PointingHandCursor
                        }
                    }

                    Rectangle {
                        id: sawButton
                        radius: 16
                        height: 32
                        color: sys.light
                        Layout.preferredWidth: 80

                        RowLayout {
                            VectorImage {
                                preferredRendererType: VectorImage.CurveRenderer
                                source: "assets/svg/saw.svg"
                                scale: 0.7
                            }

                            Text {
                                anchors.centerIn: parent
                                text: "    Saw"
                                color: sys.text
                                font.pixelSize: 16
                                font.family: "Inter"
                            }
                        }

                        MouseArea {
                            anchors.fill: parent
                            onClicked: waveGroup.checkedButton = sawButton
                            cursorShape: Qt.PointingHandCursor
                        }
                    }
                }

                // synth controls
                Rectangle {
                    Layout.fillWidth: true
                    Layout.preferredHeight: 240
                    radius: 12
                    color: sys.base
                    border.color: sys.midlight
                    border.width: 1

                    Text {
                        anchors.centerIn: parent
                        text: "insert sine waves here"
                        color: sys.text
                    }
                }
            }
        }
    }
}
