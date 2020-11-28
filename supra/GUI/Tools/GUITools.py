import os

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from pyqtgraph.Qt import QtGui, QtCore

from functools import partial

from supra.Utils.Formatting import checkExt

import numpy as np

def clearLayout(layout):
    """ Clears a layout by orphaning its children
    """

    for i in reversed(range(layout.count())): 
        layout.itemAt(i).widget().setParent(None)

def errorMessage(message, level, info='', title='Yikes!', detail=''):
    """
    message: base text
    level: 0 - info, 1 - warning, 2 - critical
    info: more detail
    detail: even more detail
    """
    msg = QMessageBox()
    
    app_icon = QtGui.QIcon()


    if level == 0:
        msg.setIcon(QMessageBox.Information)
        app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'Info.png'), QtCore.QSize(16, 16))
        msg.setWindowIcon(app_icon)
    elif level == 1:
        msg.setIcon(QMessageBox.Warning)
        app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'Warning.png'), QtCore.QSize(16, 16))
        msg.setWindowIcon(app_icon)
    else:
        msg.setIcon(QMessageBox.Critical)
        app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'Critical.png'), QtCore.QSize(16, 16))
        msg.setWindowIcon(app_icon)

    msg.setText(message)
    
    msg.setInformativeText(info)
    msg.setWindowTitle(title)
    msg.setDetailedText(detail)

    msg.setStandardButtons(QMessageBox.Ok)

    msg.exec_()

def toolTime(var):

    tool_time_dict = {}

    with open(os.path.join('supra', 'Misc', 'tool_tips.csv')) as f:
        for line in f:
            line = line.split(':')
            tool_time_dict[line[0]] = line[1]

    return tool_time_dict[var]

def createTab(parent, name):

    tab = QWidget()
    tab_layout = QGridLayout()
    tab.setLayout(tab_layout)
    parent.addTab(tab, name)

    return tab, tab_layout

def createLabelEditObj(label_name, parent, row, width=1, h_shift=0, tool_tip='', validate='', default_txt=''):
    """ Creates a label and line edit object beside each other

    Arguments:
    label_name [String]: string to be displayed in the label
    parent [obj]: the object to hold this item. Must be in QGridLayout
    row [int]: the row the object will sit on
    
    Keyword Arguments:
    width [int]: width of the lineedit box
    h_shift [int]: horizontal position to place the combination
    tool_tip [String]: identifier name in tool_tips.csv that the tool tip is under

    Returns:
    label_obj, edits_obj [obj]: the label and lineedit objects
    """
    
    label_obj = QLabel(label_name)
    edits_obj = QLineEdit(default_txt)
    parent.addWidget(label_obj, row, 1 + h_shift)
    parent.addWidget(edits_obj, row, 2 + h_shift, 1, width)

    if tool_tip != '':
        label_obj.setToolTip(toolTime(tool_tip))

    if validate == 'int':
        edits_obj.setValidator(QIntValidator())
    elif validate == 'float':
        edits_obj.setValidator(QDoubleValidator())
    elif validate == 'regexp':
        edits_obj.setValidator(QRegExpValidator())

    return label_obj, edits_obj

def createButton(text, parent, r, c, button_func, args=[]):
    button = QPushButton(text)
    parent.addWidget(button, r, c)
    button.clicked.connect(partial(button_func, *args))
    return button

def createFileSearchObj(label_name, parent, row, width=1, h_shift=0, tool_tip=''):

    label_obj, edits_obj = createLabelEditObj(label_name, parent, row, width=width, h_shift=h_shift, tool_tip=tool_tip)

    buton_obj = QPushButton('Browse')
    parent.addWidget(buton_obj, row, 3 + h_shift)

    return label_obj, edits_obj, buton_obj

def createComboBoxObj(label_name, parent, row, items=[], width=1, h_shift=0, tool_tip=''):

    label_obj = QLabel(label_name)
    combo_obj = QComboBox()
    parent.addWidget(label_obj, row, 1 + h_shift)
    parent.addWidget(combo_obj, row, 2 + h_shift, 1, width)

    for item in items: combo_obj.addItem(item)

    if tool_tip != '':
        label_obj.setToolTip(toolTime(tool_tip))

    return label_obj, combo_obj

def createToggle(label_name, parent, row, width=1, h_shift=0, tool_tip=''):

    toggle_obj = QCheckBox(label_name)
    parent.addWidget(toggle_obj, row, 1 + h_shift, 1, width)

    return toggle_obj

def createLabelDateEditObj(label_name, parent, row, width=1, h_shift=0, tool_tip='', popup=True):
    
    label_obj = QLabel(label_name)
    dedit_obj = QDateTimeEdit()
    parent.addWidget(label_obj, row, 1 + h_shift)
    parent.addWidget(dedit_obj, row, 2 + h_shift, 1, width)
    dedit_obj.setCalendarPopup(popup)
    dedit_obj.setDisplayFormat("yyyy-MM-dd HH:mm:ss.zzzzzz")

    if tool_tip != '':
        label_obj.setToolTip(toolTime(tool_tip))

    return label_obj, dedit_obj

def comboSet(obj, text):
    
    index = obj.findText(str(text).lower(), Qt.MatchFixedString)
    if index >= 0:
        obj.setCurrentIndex(index)

def toTable(obj, table):
    
    if len(table) > 0:
        X = len(table)
        Y = len(table[0])

        obj.setRowCount(X)
        obj.setColumnCount(Y)

        for x in range(X):
            for y in range(Y):
                if QTableWidgetItem(str(table[x][y])) == None:
                    obj.setItem(x, y, QTableWidgetItem(''))
                else:
                    obj.setItem(x, y, QTableWidgetItem(str(table[x][y])))

    else:
        print("Warning: Table has no length")

def fromTable(obj):

    X = obj.rowCount()
    Y = obj.columnCount()

    table = []

    for x in range(X):
        line = [None]*Y
        for y in range(Y):
            try:
                line[y] = float(obj.item(x, y).text())
            except ValueError:
                line[y] = obj.item(x, y).text()
            except AttributeError:
                line[y] = None
        table.append(line)
    return table

def toTableFromStn(obj, table):
    
    X = len(table)

    obj.setRowCount(X)

    for x in range(X):
        stn = table[x]

        obj.setItem(x, 0, QTableWidgetItem(str(stn.network)))
        obj.setItem(x, 1, QTableWidgetItem(str(stn.code)))
        obj.setItem(x, 2, QTableWidgetItem(str(stn.position.lat)))
        obj.setItem(x, 3, QTableWidgetItem(str(stn.position.lon)))
        obj.setItem(x, 4, QTableWidgetItem(str(stn.position.elev)))
        obj.setItem(x, 5, QTableWidgetItem(str(stn.name)))
        obj.setItem(x, 7, QTableWidgetItem(str(stn.file_name)))

def defTable(obj, h, w, headers=None):

    obj.setRowCount(h)
    obj.setColumnCount(w)

    if len(headers) == w:
        obj.setHorizontalHeaderLabels(headers)
        header = obj.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Stretch)
    else:
        print("ERROR: Number of headers is not the same as number of columns!")

def setTableRow(obj, row, terms=[]):

    for i in range(len(terms)):
        obj.setItem(row, i, QTableWidgetItem(str(terms[i])))

def createScatter(obj, x, y, title='', xlabel='', ylabel='', xunit='', yunit='', pen='r', update=True):
    obj.scatterPlot(x=x, y=y, pen=pen, update=update)
    obj.setTitle(title)
    obj.setLabel('bottom', xlabel, units=xunit)
    obj.setLabel('left', ylabel, units=yunit)

def createGrad(resid):

    max_error = max(resid)
    min_error = min(resid)

    error_range = max_error - min_error 

    resid = (resid - min_error)/error_range
    if np.isnan(resid):
        resid = max_error

    return resid

def sphereData(a, b, c, X, Y, Z):
    # Make data
    u = np.linspace(0, 2 * np.pi, 10)
    v = np.linspace(0, np.pi, 10)
    x = a * np.outer(np.cos(u), np.sin(v)) + X
    y = b * np.outer(np.sin(u), np.sin(v)) + Y
    z = c * np.outer(np.ones(np.size(u)), np.cos(v)) + Z

    return x, y, z

def changeRows(obj, change):
    n_rows = obj.rowCount()
    if change == 1:
        obj.insertRow(n_rows)
    else:
        obj.removeRow(n_rows-1)

def fileSearch(filters, obj):
    # obj - where the text goes
    # filter example: ['NetCDF (*.nc)', 'HDF (*.HDF)', 'CSV (*.csv)', 'TXT (*.txt)']
    dlg = QFileDialog()
    dlg.setFileMode(QFileDialog.AnyFile)
    dlg.setNameFilters(filters)
    dlg.selectNameFilter(filters[0])
    #filenames = QStringList()

    dlg.exec_()

    filename = dlg.selectedFiles()


    if obj != None:
        try:
            obj.setText(filename[0])
        except:
            errorMessage("File not Found!", 1)
    else:
        try:
            return filename[0]
        except IndexError:
            return filename

def saveFile(ext):
    dlg = QFileDialog.getSaveFileName(None, 'Save File')
    file_name = checkExt(dlg[0], ext)
    return file_name

def folderSearch(obj):
    # dlg = QFileDialog()
    # dlg.getExistingDirectory(self, "Select Directory")

    # obj - where the text goes
    dlg = QFileDialog()
    dlg.setFileMode(QFileDialog.Directory)
    #dlg.setNameFilters()
    # dlg.selectNameFilter(filters[0])
    #filenames = QStringList()

    dlg.exec_()

    filename = dlg.selectedFiles()

    obj.setText(filename[0])

def canvasFix(canvas):
    pass
