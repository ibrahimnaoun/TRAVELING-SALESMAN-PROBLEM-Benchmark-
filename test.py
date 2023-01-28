# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'untitled.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QTimer
from PyQt5.QtWidgets import QPushButton


class Ui_SecWindow(object):
    def setupUi(self, SecWindow):
        SecWindow.setObjectName("MainWindow")
        SecWindow.resize(853, 540)
        self.centralwidget = QtWidgets.QWidget(SecWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.CEntralLabel = QtWidgets.QLabel(self.centralwidget)
        self.CEntralLabel.setGeometry(QtCore.QRect(0, 0, 853, 540))
        movie = QtGui.QMovie("test.gif")
        movie.start()
        self.CEntralLabel.setMovie(movie)
        self.CEntralLabel.setObjectName("CEntralLabel")
        self.Title = QtWidgets.QLabel(self.centralwidget)
        self.Title.setGeometry(QtCore.QRect(170, 0, 521, 71))
        self.Title.setObjectName("Title")
        self.GA_Frame = QtWidgets.QFrame(self.centralwidget)
        self.GA_Frame.setGeometry(QtCore.QRect(69, 110, 191, 271))
        self.GA_Frame.setFrameShape(QtWidgets.QFrame.Box)
        self.GA_Frame.setFrameShadow(QtWidgets.QFrame.Plain)
        self.GA_Frame.setObjectName("GA_Frame")
        self.GA_title = QtWidgets.QLabel(self.GA_Frame)
        self.GA_title.setGeometry(QtCore.QRect(0, 0, 191, 51))
        self.GA_title.setObjectName("GA_title")
        self.GA_res = QtWidgets.QLabel(self.GA_Frame)
        self.GA_res.setGeometry(QtCore.QRect(6, 50, 181, 171))
        self.GA_res.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.GA_res.setObjectName("GA_res")
        self.GA_progressBar = QtWidgets.QProgressBar(self.GA_Frame)
        self.GA_progressBar.setGeometry(QtCore.QRect(20, 230, 171, 31))
        self.GA_progressBar.setProperty("value", 0)
        self.GA_progressBar.setObjectName("GA_progressBar")
        self.ABC_Frame = QtWidgets.QFrame(self.centralwidget)
        self.ABC_Frame.setGeometry(QtCore.QRect(340, 110, 191, 271))
        self.ABC_Frame.setFrameShape(QtWidgets.QFrame.Box)
        self.ABC_Frame.setFrameShadow(QtWidgets.QFrame.Plain)
        self.ABC_Frame.setObjectName("ABC_Frame")
        self.ABC_title = QtWidgets.QLabel(self.ABC_Frame)
        self.ABC_title.setGeometry(QtCore.QRect(0, 0, 191, 51))
        self.ABC_title.setObjectName("ABC_title")
        self.ABC_res = QtWidgets.QLabel(self.ABC_Frame)
        self.ABC_res.setGeometry(QtCore.QRect(6, 50, 181, 171))
        self.ABC_res.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.ABC_res.setObjectName("ABC_res")
        self.ABC_progressBar = QtWidgets.QProgressBar(self.ABC_Frame)
        self.ABC_progressBar.setGeometry(QtCore.QRect(20, 230, 171, 31))
        self.ABC_progressBar.setProperty("value", 0)
        self.ABC_progressBar.setObjectName("ABC_progressBar")
        self.DSMO_frame = QtWidgets.QFrame(self.centralwidget)
        self.DSMO_frame.setGeometry(QtCore.QRect(620, 110, 191, 271))
        self.DSMO_frame.setFrameShape(QtWidgets.QFrame.Box)
        self.DSMO_frame.setFrameShadow(QtWidgets.QFrame.Plain)
        self.DSMO_frame.setObjectName("DSMO_frame")
        self.DSMO_title = QtWidgets.QLabel(self.DSMO_frame)
        self.DSMO_title.setGeometry(QtCore.QRect(0, 0, 191, 51))
        self.DSMO_title.setObjectName("DSMO_title")
        self.DSMO_res = QtWidgets.QLabel(self.DSMO_frame)
        self.DSMO_res.setGeometry(QtCore.QRect(6, 50, 181, 171))
        self.DSMO_res.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.DSMO_res.setObjectName("DSMO_res")
        self.DSMO_progressBar = QtWidgets.QProgressBar(self.DSMO_frame)
        self.DSMO_progressBar.setGeometry(QtCore.QRect(20, 230, 171, 31))
        self.DSMO_progressBar.setProperty("value", 0)
        self.DSMO_progressBar.setObjectName("DSMO_progressBar")
        self.loading = QtWidgets.QLabel("Loading...",self.centralwidget)
        self.loading.setGeometry(QtCore.QRect(320, 420, 240, 60))
        self.loading.setStyleSheet("font-size:36pt; font-weight:600; align=\"center\"; color:black")
        self.loading.setObjectName("Loading")
        self.timer = QTimer()
        self.timer.timeout.connect(lambda: self.loading.setText("Done"))
        self.timer.start(5000)
        SecWindow.setCentralWidget(self.centralwidget)
        self.retranslateUi(SecWindow)
        QtCore.QMetaObject.connectSlotsByName(SecWindow)

    def retranslateUi(self, SecWindow):
        _translate = QtCore.QCoreApplication.translate
        SecWindow.setWindowTitle(_translate("MainWindow", "TSP Benchmark"))
        self.Title.setText(_translate("MainWindow", "<html><head/><body><p align=\"center\"><span style=\" font-size:48pt; font-weight:600; color:#d61cd3;\">TSP Benchmark</span></p></body></html>"))
        self.GA_title.setText(_translate("MainWindow", "<html><head/><body><p align=\"center\"><span style=\" font-size:20pt; font-weight:600; color:#fc07fc;\">GA</span></p></body></html>"))
        self.GA_res.setText(_translate("MainWindow", "<html><head/><body><p><br/></p></body></html>"))
        self.ABC_title.setText(_translate("MainWindow", "<html><head/><body><p align=\"center\"><span style=\" font-size:20pt; font-weight:600; color:#fc07fc;\">ABC</span></p></body></html>"))
        self.ABC_res.setText(_translate("MainWindow", "<html><head/><body><p><br/></p></body></html>"))
        self.DSMO_title.setText(_translate("MainWindow", "<html><head/><body><p align=\"center\"><span style=\" font-size:20pt; font-weight:600; color:#fc07fc;\">DSMO</span></p></body></html>"))
        self.DSMO_res.setText(_translate("MainWindow", "<html><head/><body><p><br/></p></body></html>"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    SecWindow = QtWidgets.QMainWindow()
    ui = Ui_SecWindow()
    ui.setupUi(SecWindow)
    SecWindow.show()

    sys.exit(app.exec_())