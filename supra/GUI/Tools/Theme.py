
def theme(obj):

    stylesheet = """ 

    QMainWindow{background: black;}
    QPushButton{color: white; background-color: rgb(0, 100, 200);}
    QDateTimeEdit {
    background: rgb(64, 64, 64);
    selection-background-color: rgb(0, 100, 200);
    selection-color: white;
    color: white;
    }
    QCalendarWidget QWidget { alternate-background-color: rgb(128, 128, 128); }
    QComboBox, QAbstractItemView{
    background: rgb(64, 64, 64);
    selection-background-color: rgb(0, 100, 200);
    selection-color: white;
    color: white;
    }
    QCalendarWidget QAbstractItemView:disabled 
    { 
    color: rgb(128, 128, 128); 
    }

    QCalendarWidget QWidget{
    color: white;
    background: rgb(64, 64, 64);
    }
    QCalendarWidget QAbstractItemView:enabled 
    { 
    color: white;  
    background-color: black;  
    selection-background-color: rgb(0, 100, 200); 
    selection-color: white; 
    }
    
    QCheckBox{color: white;}
    QCheckBox::indicator:checked{
    color: black;
    }
    QTableWidget{color: white; background: rgb(64, 64, 64);}
    QMenuBar{background-color: black;}
    QMenuBar::item {color: white;}
    QMenuBar::item:selected {
    background: rgb(0, 100, 200);
    }
    QMenu{
    background-color: black;
    color: white;
    }
    QMenu::item:selected{
    background: rgb(0, 100, 200);
    }

    QTabWidget>QWidget>QWidget{
    background: black;}
    QTabWidget::pane{
    background: black;
    }
    QTabWidget::tab-bar{
    background: black;
    }
    QTabBar::tab{
    background: black;
    color: white;
    }
    QTabBar::tab:selected {
    background: rgb(0, 100, 200);
    }
    QTabBar::tab:hover {
    background: rgb(64, 64, 64);
    }
    QLabel{color: white;}
    QLineEdit {
    background: rgb(64, 64, 64);
    selection-background-color: rgb(0, 100, 200);
    selection-color: white;
    color: white;
    }
    QPlainTextEdit{
    background: rgb(64, 64, 64);
    selection-background-color: rgb(0, 100, 200);
    selection-color: white;
    color: white;
    }
    QProgressBar {
    border: 2px solid grey;
    border-radius: 5px;
    text-align: center;
    }
    QProgressBar::chunk{background-color: rgb(0, 100, 200);}
    QDockWidget{color: white; background: black;}
    QGroupBox{color: white; background: rgb(25, 25, 25);}
    QGroupBox{ 
    border: 2px white; 
    border-radius: 0px; }
    QMessageBox{color: white; background: black;} 

    QWindow{background: black;}
    QScrollArea{color: white; background: black;}

    QCheckBox::indicator:checked {
    image: url(../Supracenter/supra/GUI/Images/Fin.png);}

    QCheckBox::indicator:unchecked {
    image: url(../Supracenter/supra/GUI/Images/NoFin.png);}

    QScrollArea{background: black;}
    QHeaderView::section{Background-color:black;}
    QTableView QTableCornerButton::section {
    background: black;

}
    """

    obj.setStyleSheet(stylesheet)