from PyQt5 import QtCore
from PyQt5 import QtGui
from PyQt5 import QtWidgets
from PyQt5 import QtOpenGL

import OpenGL.GL as GL
from OpenGL import GLU    
from OpenGL.arrays import vbo

import sys
import numpy as np

class GLWidget(QtOpenGL.QGLWidget):
    def __init__(self, parent=None):
        self.parent = parent
        QtOpenGL.QGLWidget.__init__(self, parent)
            
    def initializeGL(self):
        self.qglClearColor(QtGui.QColor(0, 0, 0))    # inicializace černé obrzovky
        GL.glEnable(GL.GL_DEPTH_TEST)                  # zapnutí testu hloubky

        self.initGeometry()
        self.rotationX = 0.0
        self.rotationY = 0.0
        self.rotationZ = 0.0
         
    def resizeGL(self, width, height):
        GL.glViewport(0, 0, width, height)
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()
        aspect = width / float(height)

        GLU.gluPerspective(45.0, aspect, 1.0, 100.0)
        GL.glMatrixMode(GL.GL_MODELVIEW)

    def paintGL(self):
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

        GL.glPushMatrix()    # kopírování aktuální matice do zásobníku

        GL.glTranslate(0.0, 0.0, -50.0)    # posunutí kostky dozadu
        GL.glScale(20.0, 20.0, 20.0)       # Nastavení měřítka kostky
        GL.glRotate(self.rotationX, 1.0, 0.0, 0.0) # Nastavení rotace kostky
        GL.glRotate(self.rotationY, 0.0, 1.0, 0.0)
        GL.glRotate(self.rotationZ, 0.0, 0.0, 1.0)
        GL.glTranslate(-0.5, -0.5, -0.5)   # posunutí kostky na střed

        GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
        GL.glEnableClientState(GL.GL_COLOR_ARRAY)

        GL.glVertexPointer(3, GL.GL_FLOAT, 0, self.vertVBO)
        GL.glColorPointer(3, GL.GL_FLOAT, 0, self.colorVBO)

        GL.glDrawElements(GL.GL_QUADS, len(self.cubeIndexArray), GL.GL_UNSIGNED_INT, self.cubeIndexArray)

        GL.glDisableClientState(GL.GL_VERTEX_ARRAY)
        GL.glDisableClientState(GL.GL_COLOR_ARRAY)

        GL.glPopMatrix()    # obnova předchozí matice v zásobníku
        
    def initGeometry(self):
        self.cubeVertexArray = np.array(
                [[0.0, 0.0, 0.0],
                 [1.0, 0.0, 0.0],
                 [1.0, 1.0, 0.0],
                 [0.0, 1.0, 0.0],
                 [0.0, 0.0, 1.0],
                 [1.0, 0.0, 1.0],
                 [1.0, 1.0, 1.0],
                 [0.0, 1.0, 1.0]])
        self.vertVBO = vbo.VBO(np.reshape(self.cubeVertexArray,
                                          (1, -1)).astype(np.float32))
        self.vertVBO.bind()
        
        self.cubeColorArray = np.array(
                [[0.0, 0.0, 0.0],
                 [1.0, 0.0, 0.0],
                 [1.0, 1.0, 0.0],
                 [0.0, 1.0, 0.0],
                 [0.0, 0.0, 1.0],
                 [1.0, 0.0, 1.0],
                 [1.0, 1.0, 1.0],
                 [0.0, 1.0, 1.0 ]])
        self.colorVBO = vbo.VBO(np.reshape(self.cubeColorArray,
                                           (1, -1)).astype(np.float32))
        self.colorVBO.bind()

        self.cubeIndexArray = np.array(
                [0, 1, 2, 3,
                 3, 2, 6, 7,
                 1, 0, 4, 5,
                 2, 1, 5, 6,
                 0, 3, 7, 4,
                 7, 6, 5, 4 ])

    def setRotationX(self, value):
        self.rotationX = value

    def setRotationY(self, value):
        self.rotationY = value

    def setRotationZ(self, value):
        self.rotationZ = value
        
class MainWindow(QtWidgets.QMainWindow):

    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)    # volání nadřazené třídy
        
        self.resize(300, 300)
        self.setWindowTitle('OpenGL Cube')

        self.glWidget = GLWidget(self)
        self.initGUI()
        
        timer = QtCore.QTimer(self)
        timer.setInterval(50)   # perioda obnovy v milisekundách
        timer.timeout.connect(self.glWidget.updateGL)
        timer.start()
        
    def initGUI(self):
        central_widget = QtWidgets.QWidget()
        gui_layout = QtWidgets.QVBoxLayout()
        central_widget.setLayout(gui_layout)

        self.setCentralWidget(central_widget)

        gui_layout.addWidget(self.glWidget)

        slider_x = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        slider_x.setMaximum(360)
        slider_x.valueChanged.connect(lambda value: self.glWidget.setRotationX(value))

        slider_y = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        slider_y.setMaximum(360)
        slider_y.valueChanged.connect(lambda value: self.glWidget.setRotationY(value))

        slider_z = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        slider_z.setMaximum(360)
        slider_z.valueChanged.connect(lambda value: self.glWidget.setRotationZ(value))
        
        gui_layout.addWidget(slider_x)
        gui_layout.addWidget(slider_y)
        gui_layout.addWidget(slider_z)

app = QtWidgets.QApplication(sys.argv)
window = MainWindow()
window.show()

sys.exit(app.exec_())
