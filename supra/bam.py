import os, sys

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *

from supra.GUI.Dialogs.SolutionGUI import SolutionGUI

app = QApplication(sys.argv)

splash_pix = QPixmap(os.path.join('supra', 'Docs', '_images', 'wmpl.png'))
splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
splash.setMask(splash_pix.mask())
splash.show()

app.processEvents()

gui = SolutionGUI()

gui.showFullScreen()
gui.showMaximized()
gui.show()

splash.finish(gui)

sys.exit(app.exec_())