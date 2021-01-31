from PyQt5 import QtGui, QtWidgets, QtCore
from functools import partial
import sys

class Window(QtWidgets.QMainWindow):
    def __init__(self, **kwargs):
        super(Window, self).__init__(**kwargs)
                
        self.setWindowTitle("Kalkukačka")
        self.setFixedSize(150, 150)
        self.center()
        self.init_gui()
        self.show()
    
    def center(self):
        """ Umístí okno do středu obrazovky """
        stredObrazu = QtWidgets.QDesktopWidget().availableGeometry().center()
        okno = self.frameGeometry()
        okno.moveCenter(stredObrazu)
        self.move(okno.topLeft())
        
    def init_gui(self):
    
        self.input = ""
        self.vysledek = 0
        
        grid = QtWidgets.QGridLayout()
        formular = QtWidgets.QWidget()
        formular.setLayout(grid)
        self.setCentralWidget(formular)
        
        self.entryBox = QtWidgets.QLabel("0", self)
        self.entryBox.setFont(QtGui.QFont("Arial", 14))
        self.entryBox.setAlignment(QtCore.Qt.AlignLeft)
        grid.addWidget(self.entryBox, 0, 0, 1, 5, QtCore.Qt.AlignRight)
        
        self.sevenButton = QtWidgets.QPushButton("7", self)
        grid.addWidget(self.sevenButton, 1, 0, QtCore.Qt.AlignLeft)
        self.eightButton = QtWidgets.QPushButton("8", self)
        grid.addWidget(self.eightButton, 1, 1, QtCore.Qt.AlignLeft)
        self.nineButton = QtWidgets.QPushButton("9", self)
        grid.addWidget(self.nineButton, 1, 2, QtCore.Qt.AlignLeft)
        self.CButton = QtWidgets.QPushButton("C", self)
        grid.addWidget(self.CButton, 1, 3, QtCore.Qt.AlignLeft)
        self.plusButton = QtWidgets.QPushButton("+", self)
        grid.addWidget(self.plusButton, 1, 4, QtCore.Qt.AlignLeft)
               
        self.fourButton = QtWidgets.QPushButton("4", self)
        grid.addWidget(self.fourButton, 2, 0, QtCore.Qt.AlignLeft)
        self.fiveButton = QtWidgets.QPushButton("5", self)
        grid.addWidget(self.fiveButton, 2, 1, QtCore.Qt.AlignLeft)
        self.sixButton = QtWidgets.QPushButton("6", self)
        grid.addWidget(self.sixButton, 2, 2, QtCore.Qt.AlignLeft)
        self.CEButton = QtWidgets.QPushButton("CE", self)
        grid.addWidget(self.CEButton, 2, 3, QtCore.Qt.AlignLeft)
        self.minusButton = QtWidgets.QPushButton("-", self)
        grid.addWidget(self.minusButton, 2, 4, QtCore.Qt.AlignLeft)
               
        self.oneButton = QtWidgets.QPushButton("1", self)
        grid.addWidget(self.oneButton, 3, 0, QtCore.Qt.AlignLeft)
        self.twoButton = QtWidgets.QPushButton("2", self)
        grid.addWidget(self.twoButton, 3, 1, QtCore.Qt.AlignLeft)
        self.threeButton = QtWidgets.QPushButton("3", self)
        grid.addWidget(self.threeButton, 3, 2, QtCore.Qt.AlignLeft)
        self.nullButton = QtWidgets.QPushButton("0", self)
        grid.addWidget(self.nullButton, 3, 3, QtCore.Qt.AlignLeft)
        self.multiplyButton = QtWidgets.QPushButton("*", self)
        grid.addWidget(self.multiplyButton, 3, 4, QtCore.Qt.AlignLeft)
        
        self.equalButton = QtWidgets.QPushButton("=", self)
        grid.addWidget(self.equalButton, 4, 0, 1, 3, QtCore.Qt.AlignLeft)
        self.pointButton = QtWidgets.QPushButton(",", self)
        grid.addWidget(self.pointButton, 4, 3, QtCore.Qt.AlignLeft)
        self.divideButton = QtWidgets.QPushButton("/", self)
        grid.addWidget(self.divideButton, 4, 4, QtCore.Qt.AlignLeft)
                
        self.sevenButton.clicked.connect(partial(self.click, "7"))
        self.eightButton.clicked.connect(partial(self.click, "8"))
        self.nineButton.clicked.connect(partial(self.click, "9"))
        self.CButton.clicked.connect(partial(self.click, "C"))
        self.plusButton.clicked.connect(partial(self.click, "+"))
        
        self.fourButton.clicked.connect(partial(self.click, "4"))
        self.fiveButton.clicked.connect(partial(self.click, "5"))
        self.sixButton.clicked.connect(partial(self.click, "6"))
        self.CEButton.clicked.connect(partial(self.click, "CE"))
        self.minusButton.clicked.connect(partial(self.click, "-"))
        
        self.oneButton.clicked.connect(partial(self.click, "1"))
        self.twoButton.clicked.connect(partial(self.click, "2"))
        self.threeButton.clicked.connect(partial(self.click, "3"))
        self.nullButton.clicked.connect(partial(self.click, "0"))
        self.multiplyButton.clicked.connect(partial(self.click, "*"))
        
        self.equalButton.clicked.connect(partial(self.click, "="))
        self.pointButton.clicked.connect(partial(self.click, ","))
        self.divideButton.clicked.connect(partial(self.click, "/"))
        
    def click(self, parametr):
        if parametr == "C" or parametr == "CE":
            self.input = ""
            self.entryBox.setText("")
        elif parametr == "=":
            self.vysledek = eval(self.input)
            self.entryBox.setText(str(self.vysledek))
        else:
            self.input = self.input + parametr
            self.entryBox.setText(self.input)
            
aplikace = QtWidgets.QApplication(sys.argv)
okno = Window()
sys.exit(aplikace.exec_())
