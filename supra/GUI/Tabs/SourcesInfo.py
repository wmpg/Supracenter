import os

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *

from supra.GUI.Tools.CustomWidgets import SourceEx
from supra.GUI.Tools.GUITools import clearLayout
from supra.GUI.Dialogs.AddSource import SourceWindow
from supra.Files.SaveLoad import save, loadSourcesIntoBam

def addSource(obj):
    obj.p = SourceWindow(obj.bam)
    obj.p.setGeometry(QRect(500, 400, 500, 400))
    obj.p.show()

    saveSource(obj)
    loadSource(obj)

def saveSource(obj):

    save(obj)
    loadSourcesIntoBam(obj.bam)

def loadSource(obj):

    clearSource(obj)

    for src in obj.bam.source_list:
        print(src.source)
        widget = SourceEx(src)
        obj.sources_table_layout.addWidget(widget)

def clearSource(obj):
    clearLayout(obj.sources_table_layout)

def delSource(obj):

    new_list = []

    for i in reversed(range(obj.sources_table_layout.count())):     
        if obj.sources_table_layout.itemAt(i).widget().toggle.getState():
            new_list.append(obj.sources_table_layout.itemAt(i).widget().source)


    obj.bam.source_list = new_list

    clearSource(obj)
    loadSource(obj)
    saveSource(obj)