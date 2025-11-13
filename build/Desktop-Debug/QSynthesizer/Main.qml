import QtQuick 2.15
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.VectorImage

ApplicationWindow {
    width: 1000
    height: 700
    visible: true
    color: "#1A1111"
    title: "QSynthesizer"

    FontLoader { id: boxicons; source: "assets/fonts/boxicons.ttf" }

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
            color: "#1A1111"
            radius: 0

            Text {
                text: "QSynthesizer"
                anchors.centerIn: parent
                color: "#f8d0d0"
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
                color: "#1A1111"
                Text {
                    text: "×"
                    anchors.centerIn: parent
                    color: "#F0DDDC"
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
            color: "#1A1111"
            radius: 0

            ColumnLayout {
                anchors.horizontalCenter: parent.horizontalCenter
                anchors.bottom: parent.bottom
                anchors.top: parent.top
                spacing: 12

                Rectangle {
                    Layout.alignment: Qt.AlignTop
                    width: 48; height: 48; radius: 12
                    color: "#723331"
                    Text {
                        anchors.centerIn: parent
                        text: "+"
                        color: "#FBD5D2"
                        font.pixelSize: 32
                    }
                    MouseArea {
                        anchors.fill: parent
                        cursorShape: Qt.PointingHandCursor
                    }
                }
                Rectangle {
                    width: 48; height: 48; radius: 12
                    color: "#1A1111"
                    Layout.alignment: Qt.AlignTop
                    Text {
                        anchors.centerIn: parent
                        font.pixelSize: 32
                        font.family: boxicons.font.family
                        text: "\uEDCB"
                        color: "#D8C1C0"
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
                    color: "#1A1111"
                    Layout.alignment: Qt.AlignBottom
                    Text {
                        anchors.centerIn: parent
                        font.family: boxicons.font.family
                        text: "\uED52"
                        color: "#D8C1C0"
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
            color: "#231919"
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
                    color: "#231919"
                    gradient: Gradient {
                        orientation: Gradient.Horizontal
                        GradientStop { position: 0.0; color: "#00FFFFFF" }
                        GradientStop { position: 1.0; color: "#13FFA09D" }
                    }

                    Text {
                        anchors.centerIn: parent
                        text: "waveform or something"
                        color: "#F0DDDC"
                    }
                }

                RowLayout {
                    Text {
                        font.family: boxicons.font.family
                        text: "\uED70"
                        font.pixelSize: 32
                        font.bold: true
                        color: "#F0DDDC"
                    }

                    Text {
                        text: " Synthesize"
                        font.family: "Inter"
                        font.bold: true
                        font.pixelSize: 24
                        color: "#F0DDDC"
                    }
                }

                RowLayout {
                    spacing: 8

                    ButtonGroup { id: waveGroup }

                    // Sine
                    Rectangle {
                        id: sineButton
                        radius: 16
                        height: 32
                        color: waveGroup.checkedButton === sineButton ? "#ECAAA7" : "#4C3837"
                        Layout.preferredWidth: 80

                        RowLayout {
                            VectorImage {
                                preferredRendererType: VectorImage.CurveRenderer
                                source: "assets/svg/sine.svg"
                            }

                            Text {
                                anchors.centerIn: parent
                                text: "    Sine"   // Unicode sine-like symbol
                                color: waveGroup.checkedButton === sineButton ? "#2A1818" : "#F8D6D5"
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

                    // Square
                    Rectangle {
                        id: squareButton
                        radius: 16
                        height: 32
                        color: waveGroup.checkedButton === squareButton ? "#ECAAA7" : "#4C3837"
                        Layout.preferredWidth: 90

                        Text {
                            anchors.centerIn: parent
                            text: "\u23F9  Square"  // Unicode stop-square symbol
                            color: waveGroup.checkedButton === squareButton ? "#2A1818" : "#F8D6D5"
                            font.pixelSize: 16
                            font.family: "Inter"
                        }

                        MouseArea {
                            anchors.fill: parent
                            onClicked: waveGroup.checkedButton = squareButton
                            cursorShape: Qt.PointingHandCursor
                        }
                    }

                    // Saw
                    Rectangle {
                        id: sawButton
                        radius: 16
                        height: 32
                        color: waveGroup.checkedButton === sawButton ? "#ECAAA7" : "#4C3837"
                        Layout.preferredWidth: 80

                        Text {
                            anchors.centerIn: parent
                            text: "\u2A57  Saw"  // Unicode angled wave
                            color: waveGroup.checkedButton === sawButton ? "#2A1818" : "#F8D6D5"
                            font.pixelSize: 16
                            font.family: "Inter"
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
                    color: "#1A1111"
                    border.color: "#3a2a2a"
                    border.width: 1

                    Text {
                        anchors.centerIn: parent
                        text: "insert sine waves here"
                        color: "#f8d0d0"
                    }
                }
            }
        }
    }
}
