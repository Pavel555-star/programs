import sys
from OpenGL.GL import *
from PyQt5.QtOpenGL import *
from PyQt5 import QtCore, QtWidgets, QtGui

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__()

        self.widget = OpenGLWidget()
        mainLayout = QtWidgets.QVBoxLayout()
        mainLayout.addWidget(self.widget)
        
        central_widget = QtWidgets.QWidget()
        central_widget.setLayout(mainLayout)
        self.setCentralWidget(central_widget)
        
        timer = QtCore.QTimer(self)
        timer.setInterval(20)
        timer.timeout.connect(self.widget.updateGL)
        timer.start()

class OpenGLWidget(QGLWidget):
    def __init__(self, parent=None):
        QGLWidget.__init__(self, parent)
        QGLWidget.__init__(self, parent)
        self.rotY = 2.0

    def initializeGL(self):
        self.qglClearColor(QtGui.QColor(0, 0, 0))
        gl.glEnable(gl.GL_DEPTH_TEST)
        
    def resizeGL(self, width, height):
        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        aspect = width / float(height)

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glRotate(self.rotY, 0.0, 1.0, 0.0)
        glBegin(GL_TRIANGLES)
        glColor3f(1,0,0)
        glVertex3f(0,0.707,0)
        glColor3f(0,1,0)
        glVertex3f(-0.5,-0.5,0)
        glColor3f(0,0,1)
        glVertex3f(0.5,-0.5,0)
        glEnd()

app = QtWidgets.QApplication(sys.argv)    
Form = QtWidgets.QMainWindow()
window = MainWindow(Form)    
window.show()    
sys.exit(app.exec_())
glFlush()
