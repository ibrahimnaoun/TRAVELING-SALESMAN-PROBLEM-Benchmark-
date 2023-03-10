# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'main.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.
import random
import copy
import csv
import math
import random
import sys
from scipy.spatial import distance
import time
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QObject, QThread, pyqtSignal
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QGraphicsOpacityEffect, QFileDialog, QMessageBox
import subprocess


import tsp_ga

tree=[]
class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(853, 540)
        MainWindow.setFixedSize(853,540)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.img=QPixmap('bg.jpg')
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(0, 0, 853, 540))
        self.label.setText("")
        self.label.setObjectName("label")
        self.label.setPixmap(self.img)
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(60, 10, 731, 131))
        self.label_2.setScaledContents(False)
        self.label_2.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignTop)
        self.label_2.setWordWrap(True)
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(100, 150, 661, 51))
        self.label_3.setWordWrap(True)
        self.label_3.setObjectName("label_3")
        self.frame = QtWidgets.QFrame(self.centralwidget)
        self.frame.setGeometry(QtCore.QRect(100, 250, 671, 231))
        self.frame.setFrameShape(QtWidgets.QFrame.Box)
        self.frame.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.frame.setLineWidth(4)
        self.frame.setObjectName("frame")
        self.pushButton = QtWidgets.QPushButton(self.frame)
        self.pushButton.setGeometry(QtCore.QRect(180, 30, 300, 71))
        self.pushButton.setObjectName("pushButton")
        self.pushButton.setStyleSheet("font-size:20pt; color:#fc07fc; border-radius :30; border: 3px solid #fc07fc")
        self.pushButton.clicked.connect(self.open_file)
        self.label_4 = QtWidgets.QLabel(self.frame)
        self.label_4.setGeometry(QtCore.QRect(80, 110, 561, 31))
        self.label_4.setObjectName("label_4")
        self.opacity_effect1 = QGraphicsOpacityEffect()
        # setting opacity level
        self.opacity_effect1.setOpacity(0.5)
        # adding opacity effect to the label
        self.label_4.setGraphicsEffect(self.opacity_effect1)
        self.label_5 = QtWidgets.QLabel(self.frame)
        self.label_5.setGeometry(QtCore.QRect(220, 173, 150, 20))
        self.label_5.setObjectName("label_5")
        self.label_6 = QtWidgets.QLabel(self.frame)
        self.label_6.setGeometry(QtCore.QRect(380, 172, 21, 21))
        self.label_6.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.label_6.setFrameShadow(QtWidgets.QFrame.Plain)
        self.label_6.setObjectName("label_6")
        self.opacity_effect1 = QGraphicsOpacityEffect()
        # setting opacity level
        self.opacity_effect1.setOpacity(0.8)
        # adding opacity effect to the label
        self.label_6.setGraphicsEffect(self.opacity_effect1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def open_file(self):
        self.error_dialog = QtWidgets.QMessageBox()
        self.error_dialog.setIcon(QtWidgets.QMessageBox.Information)
        self.error_dialog.setWindowTitle("Error")
        self.error_dialog.setStandardButtons(QtWidgets.QMessageBox.Ok)
        print("clicked")
        path = QFileDialog.getOpenFileName(None,None,None,"txt files(*.txt)")
        print(path[0])
        if path[0]=="":
            print("empty")
        else:
            self.txt_path=str(path[0])


            if self.check_txt(self.txt_path):
                tree.append(self.txt_path)
                MainWindow.close()
                SecWindow.show()






                import os
                import pickle
                with open("myfile.pickle", "wb") as outfile:
                    pickle.dump(self.txt_path, outfile)
                    """pickle.dump(self.csv_path, outfile)

                print("starting SMO")

                try:
                    subprocess.run(
                        ["C:/Users/Shanks/PycharmProjects/pythonProject2/venv/Scripts/python.exe", "DSMO_TSP.py"],
                        check=True)
                except subprocess.CalledProcessError:
                    self.error_dialog.setText("Please select valid file")
                    self.error_dialog.exec_()
                print("finished SMO")
                print("Started ABC")
                try:
                    subprocess.run(
                        ["C:/Users/Shanks/PycharmProjects/pythonProject2/venv/Scripts/python.exe", "ABC_TSP.py"],
                        check=True)
                except subprocess.CalledProcessError:
                    self.error_dialog.setText("Please select valid files")
                    self.error_dialog.exec_()
                print("Finished ABC")"""
            else:
                self.error_dialog.setText("Please Select valid File !!")
                self.error_dialog.exec_()
    def check_txt(self,txt):
        import pandas as pd
        f = open(txt, "r")
        lines = f.readlines()
        try:
            for line in lines:
                newline = line[:-1].split(' ')
                for cord in newline:
                    if cord != '':
                        x = float(cord)
            file = pd.read_csv(txt, delim_whitespace=True)

            # store dataframe into csv file
            file.to_csv(txt.replace("txt","csv"), index=None)
            fa = open(txt.replace("txt","csv"), "r")
            temp = fa.read()
            fa.close()
            with open(txt.replace("txt","csv"), "w") as f:
                f.write("0,0,0\n")
                f.write(temp)
                f.close()
            return True
        except:
            return False

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.label_2.setText(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:36pt; color:#8E13C2;\">Traveling Salesman Problem</span><span style=\" font-size:36pt; color:#fc07fc;\"> Benchmark</span></p></body></html>"))
        self.label_3.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:16pt;color:#DE97FC;\">Calculate the minimum path </span><span style=\" font-size:16pt;color:#2E0042;\">for the travelling salesman using </span><span style=\" font-size:16pt; color:#8E13C2;\">ABC</span><span style=\" font-size:16pt;\">,</span><span style=\" font-size:16pt; color:#8E13C2;\">GA</span><span style=\" font-size:16pt;\"> and </span><span style=\" font-size:16pt; color:#8E13C2;\">DSMO</span><span style=\" font-size:16pt;\"> algorithms.</span></p></body></html>"))
        self.pushButton.setText(_translate("MainWindow", "Upload File"))
        self.label_4.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:14pt;background : #C34FF5;\">Select files containing enumerated cities with geometric positions</span></p></body></html>"))
        self.label_5.setText(_translate("MainWindow", "<html><head/><body><strong><span style=\" font-size:12pt;bold;\">Supported format:</span></strong></body></html>"))
        self.label_6.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:12pt; background : #DE97FC;\">txt</span></p></body></html>"))

class Ui_SecWindow(object):
    dhh=0


    def setupUi(self, SecWindow):
        SecWindow.setObjectName("MainWindow")
        SecWindow.resize(853, 540)
        SecWindow.setFixedSize(853,540)
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
        self.GA_res.setObjectName("GA_res")
        self.GA_res.setStyleSheet("font-size:20pt; color:#fc07fc")
        self.GA_res.setAlignment(QtCore.Qt.AlignCenter)
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
        self.ABC_res.setObjectName("ABC_res")
        self.ABC_res.setStyleSheet("font-size:20pt; color:#fc07fc")
        self.ABC_res.setAlignment(QtCore.Qt.AlignCenter)
        self.ABC_progressBar = QtWidgets.QProgressBar(self.ABC_Frame)
        self.ABC_progressBar.setGeometry(QtCore.QRect(20, 230, 171, 31))
        self.ABC_progressBar.setValue(self.dhh)
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
        self.DSMO_res.setObjectName("DSMO_res")
        self.DSMO_res.setStyleSheet("font-size:20pt; color:#fc07fc")
        self.DSMO_res.setAlignment(QtCore.Qt.AlignCenter)
        self.DSMO_progressBar = QtWidgets.QProgressBar(self.DSMO_frame)
        self.DSMO_progressBar.setGeometry(QtCore.QRect(20, 230, 171, 31))
        self.DSMO_progressBar.setProperty("value", 0)
        self.DSMO_progressBar.setObjectName("DSMO_progressBar")
        self.start = QtWidgets.QPushButton("Start",self.centralwidget)
        self.start.setGeometry(QtCore.QRect(340, 420, 200, 60))
        self.start.setStyleSheet("font-size:20pt; color:#fc07fc; border-radius :30; border: 3px solid #fc07fc; font-weight:200; align=\"center\"")
        self.start.setObjectName("Start")
        self.start.clicked.connect(self.buclick)
        SecWindow.setCentralWidget(self.centralwidget)
        self.retranslateUi(SecWindow)
        QtCore.QMetaObject.connectSlotsByName(SecWindow)
        #self.run()
    """def run(self):
        self.error_dialog = QtWidgets.QMessageBox()
        self.error_dialog.setIcon(QtWidgets.QMessageBox.Information)
        self.error_dialog.setWindowTitle("Error")
        self.error_dialog.setStandardButtons(QtWidgets.QMessageBox.Ok)
        import subprocess
        import sys
        self.txt_path=tree[0]
        self.csv_path = self.txt_path.replace("txt", "csv")
        print("txt:", self.txt_path, "csv:", self.csv_path)
        import os
        import pickle
        with open("myfile.pickle", "wb") as outfile:
            pickle.dump(self.txt_path, outfile)
            pickle.dump(self.csv_path, outfile)

        print("starting SMO")

        try:
            subprocess.run(["C:/Users/Shanks/PycharmProjects/pythonProject2/venv/Scripts/python.exe", "DSMO_TSP.py"],
                           check=True)
        except subprocess.CalledProcessError:
            self.error_dialog.setText("Please select valid file")
            self.error_dialog.exec_()
        print("finished SMO")
        print("Started ABC")
        try:
            subprocess.run(["C:/Users/Shanks/PycharmProjects/pythonProject2/venv/Scripts/python.exe", "ABC_TSP.py"],
                           check=True)
        except subprocess.CalledProcessError:
            self.error_dialog.setText("Please select valid files")
            self.error_dialog.exec_()
        print("Finished ABC")"""


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
    def buclick(self):
        print("clicked 2")
        self.start.setEnabled(False)
        self.worker=Worker()
        self.worker.start()
        self.worker.finished.connect(self.finished_thread)
        self.worker.update_progress.connect(self.update_bar)
        self.worker.update_progress2.connect(self.update_smobar)
        self.worker.update_progress3.connect(self.update_gabar)
        self.worker.up_l.connect(self.update_abc)
        self.worker.up_l2.connect(self.update_dsmo)
        self.worker.up_l3.connect(self.update_ga)
    def update_bar(self,val):
        self.ABC_progressBar.setValue(val)
    def update_smobar(self,val):
        self.DSMO_progressBar.setValue(val)
    def update_gabar(self,val):
        self.GA_progressBar.setValue(val)
    def update_abc(self,val):
        self.ABC_res.setText("{:.0f} km".format(val*100))
    def update_dsmo(self,val):
        self.DSMO_res.setText("{:.0f} km".format(val*100))
    def update_ga(self,val):
        self.GA_res.setText("{:.0f} km".format(val))
    def finished_thread(self):
        self.error_dialog = QtWidgets.QMessageBox()
        self.error_dialog.setIcon(QtWidgets.QMessageBox.Information)
        self.error_dialog.setWindowTitle("Done")
        self.error_dialog.setStandardButtons(QtWidgets.QMessageBox.Ok)
        self.error_dialog.setText("Calculations finished!")
        self.error_dialog.exec_()
        self.worker.exit(0)

class Worker(QThread):
    update_progress=pyqtSignal(int)
    update_progress2 = pyqtSignal(int)
    update_progress3=pyqtSignal(int)
    up_l=pyqtSignal(float)
    up_l2 = pyqtSignal(float)
    up_l3=pyqtSignal(int)
    def run(self):
        import main_ga
        g=random.randint(10,20)
        for i in range(g):
            time.sleep(0.1)
            self.update_progress3.emit(i)
            self.up_l3.emit(random.randint(100000,200000))
        main_ga.run(main_ga.args)
        for i in range(len(tsp_ga.cost_list)):
            if 10+i>100 or i>90:
                self.update_progress3.emit(100)
                break
            time.sleep(0.1)
            self.update_progress3.emit(g+i)
            self.up_l3.emit(tsp_ga.cost_list[i])
        self.update_progress3.emit(100)
        self.up_l3.emit(min(tsp_ga.cost_list))
        for i in range(11):
            abc.main(abc)
            self.up_l.emit(abc.sed(abc))
            self.update_progress.emit(i*10)
        from trying2 import DSMO
        for i in range(11):
            DSMO.dsmo_main(DSMO)
            self.up_l2.emit(DSMO.sed(DSMO))
            self.update_progress2.emit(i*10)

class abc():


    ini_time = time.time()

    class Bee:
        def __init__(self, node_set):
            self.role = ''
            self.path = list(node_set)  # stores all nodes in each bee, will randomize foragers
            self.distance = 0
            self.cycle = 0  # number of iterations on current solution

    def read_data_from_csv(file_name):
        """
        Returns data read from file.
        """
        data_list = []
        with open(file_name) as f:
            reader = csv.reader(f)
            data_list = [[eval(s) for s in row.split(',')] for row in f]
        return data_list

    def get_distance_between_nodes(n1, n2):
        """
        Calculates the Euclidean distance between two nodes.
        """
        return distance.euclidean(n1, n2)

    def make_distance_table(self,data_list):
        """
        Creates a table that stores distance between every pair of nodes.
        """
        length = len(data_list)
        table = [[self.get_distance_between_nodes(
            (data_list[i][1], data_list[i][2]), (data_list[j][1], data_list[j][2])) / 2
                  for i in range(0, length)] for j in range(0, length)]
        return table

    def get_total_distance_of_path(path, table):
        """
        Calculates the total distance of an individual bee's path.
        Terminates at starting node to complete cycle.
        """
        # Creates a copy of path, puts head at end of list.
        # Zip lists to create pairs of neighbor coords,
        # will create a cycle that terminates at starting node.
        new_path = list(path)
        new_path.insert(len(path), path[0])
        new_path = new_path[1:len(new_path)]

        coordinates = zip(path, new_path)
        distance = sum([table[i[0]][i[1]] for i in coordinates])
        return round(distance, 3)

    def initialize_hive(self,population, data):
        """
        Initializes a hive and populates it with bees.
        Bees will have a randomized path attribute.
        """
        path = [x[0] for x in data]
        hive = [self.Bee(path) for i in range(0, population)]
        return hive

    def assign_roles(self,hive, role_percentiles, table):
        """
        Assigns initial roles based on role percentiles
        to each bee in the hive.
        Assigns randomized path to forager bees.
        """
        population = len(hive)
        onlooker_count = math.floor(population * role_percentiles[0])
        forager_count = math.floor(population * role_percentiles[1])

        for i in range(0, onlooker_count):
            hive[i].role = 'O'

        for i in range(onlooker_count, (onlooker_count + forager_count)):
            hive[i].role = 'F'
            random.shuffle(hive[i].path)
            hive[i].distance = self.get_total_distance_of_path(hive[i].path, table)

        return hive

    def mutate_path(path):
        """
        Gets a random index 0 to next to last element.
        Copies path, swaps two nodes, compares distance.
        Returns mutated path.
        """
        # - will go out of range if last element is chosen.
        random_idx = random.randint(0, len(path) - 2)
        new_path = list(path)
        new_path[random_idx], new_path[random_idx + 1] = new_path[random_idx + 1], new_path[random_idx]
        return new_path

    def forage(self,bee, table, limit):
        """
        Worker bee behavior, iteratively refines a potential shortest path
        by swapping randomly selected neighbor indices.
        """
        new_path = self.mutate_path(bee.path)
        new_distance = self.get_total_distance_of_path(new_path, table)

        if new_distance < bee.distance:
            bee.path = new_path
            bee.distance = new_distance
            bee.cycle = 0  # reset cycle so bee can continue to make progress
        else:
            bee.cycle += 1
        if bee.cycle >= limit:  # if bee is not making progress
            bee.role = 'S'
        return bee.distance, list(bee.path)

    def scout(self,bee, table):
        """
        Scout bee behavior, abandons unsuccessful path for new random path.
        Resets role to forager.
        """
        new_path = list(bee.path)
        random.shuffle(new_path)
        bee.path = new_path
        bee.distance = self.get_total_distance_of_path(new_path, table)
        bee.role = 'F'
        bee.cycle = 0

    def waggle(self,hive, best_distance, table, forager_limit, scout_count):
        """
        Captures results from work of forager bees,
        chooses new random path for scouts to explore,
        returns results for overlookers to assess.
        """
        best_path = []
        results = []

        for i in range(0, len(hive)):
            if hive[i].role == 'F':
                distance, path = self.forage(self,hive[i], table, forager_limit)
                if distance < best_distance:
                    best_distance = distance
                    best_path = list(hive[i].path)
                results.append((i, distance))

            elif hive[i].role == 'S':
                self.scout(self,hive[i], table)

        # after processing all bees, set worst performers to scout
        results.sort(reverse=True, key=lambda tup: tup[1])
        scouts = [tup[0] for tup in results[0:int(scout_count)]]
        for new_scout in scouts:
            hive[new_scout].role = 'S'
        return best_distance, best_path

    def recruit(self,hive, best_distance, best_path, table):
        """
        Recruits onlooker bees to iterate on best soluction so far.
        Returns updated best_distance, best_path.
        """
        for i in range(0, len(hive)):
            if hive[i].role == 'O':
                new_path = self.mutate_path(best_path)
                new_distance = self.get_total_distance_of_path(new_path, table)
                if new_distance < best_distance:
                    best_distance = new_distance
                    best_path = new_path
        return best_distance, best_path

    def print_details(cycle, path, distance, bee):
        """
        Prints cycle details to console.
        """
        print("CYCLE: {}".format(cycle))
        print("PATH: {}".format(path))
        print("DISTANCE: {}".format(distance))
        print("BEE: {}".format(bee))
        print("\n")

    def make_csv(data, file_name):
        """
        Writes data to csv file.
        """
        with open(file_name, 'a') as f:
            writer = csv.writer(f)
            writer.writerow(data)
        f.close()

    cost_list = []
    global err
    def main(self):
        # Control parameters
        population = 100
        forager_percent = 0.5
        onlooker_percent = 0.5
        role_percent = [onlooker_percent, forager_percent]
        scout_percent = 0.2
        scout_count = math.ceil(population * scout_percent)
        forager_limit = 500
        cycle_limit = 30
        cycle = 1
        time_init = time.time()
        # Data source
        # data = read_data_from_csv("data/data_10.csv")
        # data = read_data_from_csv("data/data_11.csv")
        """import pickle
        with open("myfile.pickle", "rb") as infile:
            path0 = pickle.load(infile)
            path = pickle.load(infile)
        print(path)"""
        self.path=tree[0].replace("txt","csv")
        data = self.read_data_from_csv(self.path)

        # Global vars
        best_distance = sys.maxsize
        best_path = []
        result = ()
        result_file = "results/{}_nodes/results_{}_nodes_{}_bees_{}_scouts_{}_cycles_{}_forager_limit.csv".format(
            len(data), len(data), population, scout_count, cycle_limit, forager_limit)

        # Initialization
        table = self.make_distance_table(self,data)
        hive = self.initialize_hive(self,population, data)
        self.assign_roles(self,hive, role_percent, table)

        avg = 0
        while cycle < cycle_limit:
            waggle_distance, waggle_path = self.waggle(self,hive, best_distance, table, forager_limit, scout_count)
            if waggle_distance < best_distance:
                best_distance = waggle_distance
                best_path = list(waggle_path)
                result = (cycle, best_path, best_distance, 'F')

            recruit_distance, recruit_path = self.recruit(self,hive, best_distance, best_path, table)
            if recruit_distance < best_distance:
                best_distance = recruit_distance
                best_path = list(recruit_path)
                result = (cycle, best_path, best_distance, 'R')
            avg += best_distance
            cycle += 1

        avg = avg / (cycle_limit)
        self.err=avg
        print(avg, "\n{} s".format(time.time() - time_init))
        self.cost_list.append(avg)

    def sed(self):
        return self.err

    # ------------------------------------------------------------------------------------#

    def real_main(self):
        for i in range(0, 10):
            print("iter num: {} ".format(i))
            self.main(self)

        nInstance = len(self.cost_list)
        Best_cost = min(self.cost_list)
        Average_cost = sum(self.cost_list) / nInstance
        Deviation = abs(Best_cost - Average_cost)
        end_time = time.time() - self.ini_time
        print("total time is {} s".format(end_time))
"""
        with open("output_ABC_TSP.txt", "a") as f:
            # f.write("\t"+"Instance"+"\t"+"Best cost"+"\t"+"Average cost"+"\t"+"Deviation"+"\t"+"nInstance"+"\n")
            f.write(
                "\t" + "rat50" + "\t" + str(Best_cost) + "\t" + str(Average_cost) + "\t" + str(Deviation) + "\t" + str(
                    nInstance) + "\n")
            f.close()"""









if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    app1=QtWidgets.QApplication(sys.argv)
    SecWindow = QtWidgets.QMainWindow()
    ui1 = Ui_SecWindow()
    ui1.setupUi(SecWindow)
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
