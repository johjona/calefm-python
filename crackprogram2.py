#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Main branch

"""


import calefm_base as un
import sys
import numpy as np
from PyQt5.QtCore import pyqtSlot, pyqtSignal, QThread, QSize
from PyQt5.QtWidgets import QApplication, QDialog, QWidget, QMainWindow, QFileDialog, QMessageBox
from PyQt5.uic import loadUi
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5.QtWidgets import QApplication, QMainWindow, QMenu, QVBoxLayout, QSizePolicy, QMessageBox, QWidget, QPushButton
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout
import matplotlib.pyplot as plt
import random
import calfem.vis_mpl as cfv
from math import sqrt
              

class SolverThread(QThread):
    """ Solver class 
    """
    def __init__(self, solver):
        QThread.__init__(self)
        self.solver = solver
        
    def __del__(self):
        self.wait()
        
    def run(self):
        self.solver.execute()

        
class MainWindow(QMainWindow):
    """ Main window class for GUI
    """
    
    def __init__(self):
        super(QMainWindow, self).__init__()
        
        # Initialise input and output data variables
        self.savePoints = True
        self.crackAddedFlag = False
        
        # Instansiate input and output data objects
        self.input_data = un.InputData()
        self.output_data = un.OutputData()
        
        # Instansiate solver objects
        self.solver = un.Solver(self.input_data, self.output_data)
        
        # Instansiate visualisation object
        self.vis = un.Visualisation(self.input_data, self.output_data)
        
        # Show UI
        self.app = app
        self.ui = loadUi('crack_program.ui', self)
        self.ui.show()
        self.ui.raise_()
        
        #  Initialise model on program start 
        self.initModel()
        self.linearElasticCalc.setChecked(True)
        self.IsoButton.setChecked(True)
        self.ui.elTypeAEdit.setChecked(True)
        self.enableComplianceEdits()
        self.enableMaterialEdits()
                
        # Create canvas and add as a layout to the QWidget
        vbox = QVBoxLayout()
        self.plotWidget.setLayout(vbox) 
        
        self.figure = plt.figure(dpi=150)
        self.canvas = FigureCanvas(self.figure)
        
        self.figure_toolbar = NavigationToolbar(self.canvas, self)
        
        self.plotWidget.layout().addWidget(self.figure_toolbar) 
        self.plotWidget.layout().addWidget(self.canvas)        
        
        self.myToolBar.setIconSize(QSize(16,16))
        
        # Connect qlines to update model
        # 1 - Crack parameters
        self.ui.crackLengthEdit.editingFinished.connect(self.updateModel)
        self.ui.crackDirectionEditX.editingFinished.connect(self.updateModel)
        self.ui.crackDirectionEditY.editingFinished.connect(self.updateModel)
        self.ui.crackLocationEdit.editingFinished.connect(self.updateModel)
        self.ui.crackIncrementEdit.editingFinished.connect(self.updateModel)
        self.ui.crackLengthEndEdit.editingFinished.connect(self.updateModel)
        
        # 2 - Element size parameter
        self.ui.elSizeEdit.editingFinished.connect(self.updateModel)
        
        # 3 - Load and boundary condition parameters
        self.ui.loadMarkerEdit.editingFinished.connect(self.updateModel)
        self.ui.loadEdit.editingFinished.connect(self.updateModel)
        self.ui.bcEdit.editingFinished.connect(self.updateModel)
        self.ui.bcMarkerEdit.editingFinished.connect(self.updateModel)
        self.ui.bcDirectionEdit.editingFinished.connect(self.updateModel)
        self.ui.loadDirectionEdit.editingFinished.connect(self.updateModel)
        
        # 4 Otrhotropic parameters
        self.ui.E1_edit.editingFinished.connect(self.updateModel)
        self.ui.E2_edit.editingFinished.connect(self.updateModel)
        self.ui.E3_edit.editingFinished.connect(self.updateModel)
        
        self.ui.v12_edit.editingFinished.connect(self.updateModel)
        self.ui.v13_edit.editingFinished.connect(self.updateModel)
        self.ui.v23_edit.editingFinished.connect(self.updateModel)
        
        self.ui.G12_edit.editingFinished.connect(self.updateModel)
        self.ui.G13_edit.editingFinished.connect(self.updateModel)
        self.ui.G23_edit.editingFinished.connect(self.updateModel)
                
        #  Connect buttons to methods
        self.ui.addPointButton.clicked.connect(self.onAddPoint)
        self.ui.removePointButton.clicked.connect(self.onRemovePoint)
        self.ui.showGeometryButton.clicked.connect(self.onShowGeometry)        
        self.ui.showMeshButton.clicked.connect(self.onShowMesh)
        self.ui.showGeometryButton2.clicked.connect(self.onShowGeometry)
        self.ui.addCrackButton.clicked.connect(self.onAddCrack)
        self.ui.runCalcButton.clicked.connect(self.onActionRunCalc)
        self.ui.showDisplacementsButton.clicked.connect(self.onShowDisplacements)
        self.ui.elTypeAEdit.toggled.connect(self.defineElementType)
        self.ui.complianceCalcButton.toggled.connect(self.enableComplianceEdits)
        self.ui.plotEnergRateButton.clicked.connect(self.showEnergyReleaseRate)
        self.ui.addLoadButton.clicked.connect(self.onCreateLoad)
        self.ui.addBcButton.clicked.connect(self.onCreateBC)
        self.ui.IsoButton.clicked.connect(self.enableMaterialEdits)
        self.ui.OrtoButton.clicked.connect(self.enableMaterialEdits)
        self.ui.clearBCButton.clicked.connect(self.onClearBC)
        self.ui.clearLoadButton.clicked.connect(self.onClearLoad)
        
        # Connect actions
        self.actionNew.triggered.connect(self.onActionNew)
        
    def onActionNew(self):
        """ Initiating a new model
        """
        self.initModel()
        self.figure.clear()
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()
        self.figure.canvas.update()
        
        
    def onCreateLoad(self):
        """Function for appending the user input load to the load list, later to be used for the FE-simulation
        """
        self.input_data.load_marker_list.append(int(self.loadMarkerEdit.text()))
        self.input_data.load_list.append(float(self.loadEdit.text()))
        self.input_data.load_dir_list.append(int(self.loadDirectionEdit.text()))
        self.outputTextBrowser.append("Load created")
        
        
    def onClearLoad(self):
        """ Function for clearing all current loads
        """
        self.input_data.load_marker_list = []
        self.input_data.load_list = []
        self.input_data.load_dir_list = []
        self.outputTextBrowser.append("Loads cleared.")
        self.input_data.load_flag = False
        
        
    def onCreateBC(self):
        """Function for appending the user input boundary condition to the load list, later to be used for the FE-simulation
        """
        self.input_data.bc_marker_list.append(int(self.bcMarkerEdit.text()))
        self.input_data.bc_list.append(float(self.bcEdit.text()))
        self.input_data.bc_dir_list.append(int(self.bcDirectionEdit.text()))
        self.outputTextBrowser.append("BC created")
        
        
    def onClearBC(self):
        """ Function for clearing the current boundary conditions
        """
        self.input_data.bc_marker_list = []
        self.input_data.bc_list = []
        self.input_data.bc_dir_list = []
        self.outputTextBrowser.append("Boundary conditions cleared.")
        self.input_data.bc_flag = False


    def defineElementType(self):
        """Function for switching element type based on which radiobutton user clicks
        """
        if self.ui.elTypeAEdit.isChecked():
            self.input_data.el_type = 2
        elif self.ui.elTypeBEdit.isChecked():
            self.input_data.el_type = 3
            
        
    def initModel(self):
        """Function for initalising all 
        """
        self.input_data.points = []      
        self.ui.addPointEditX.setText("0")
        self.ui.addPointEditY.setText("0")
        self.ui.crackDirectionEditX.setText("0")
        self.ui.crackDirectionEditY.setText("1")
        self.ui.crackLengthEdit.setText("0.1")
        self.ui.crackLocationEdit.setText("0")
        self.ui.elSizeEdit.setText("0.02")
        self.ui.crackIncrementEdit.setText("0.03")
        self.ui.crackLengthEndEdit.setText("0.5")
        self.ui.loadMarkerEdit.setText("1")
        self.ui.loadEdit.setText("1")
        self.ui.bcMarkerEdit.setText("1")
        self.ui.bcEdit.setText("0")
        self.ui.EEdit.setText("210e9")
        self.ui.vEdit.setText("0.3")
        self.ui.tEdit.setText("0.01")
        self.ui.loadDirectionEdit.setText("2")
        self.ui.bcDirectionEdit.setText("0")
        
        self.ui.E1_edit.setText("11000e6")
        self.ui.E2_edit.setText("690e6")
        self.ui.E3_edit.setText("690e6")
        
        self.ui.v12_edit.setText("0.3")
        self.ui.v13_edit.setText("0.3")
        self.ui.v23_edit.setText("0.3")
        
        self.ui.G12_edit.setText("440e6")
        self.ui.G13_edit.setText("390e6")
        self.ui.G23_edit.setText("49e6")
        
        self.ui.fractureEnergyEdit.setText("2500")
        
    
    def updateModel(self):
        """ Function for updating model parameters in the input data class when the user changes the input parameters
        """
        self.input_data.crack_length = float(self.ui.crackLengthEdit.text())
        self.input_data.crack_direction = [float(self.ui.crackDirectionEditX.text()), float(self.ui.crackDirectionEditY.text())]
        self.input_data.crack_location = int(self.ui.crackLocationEdit.text())
        self.input_data.crack_increment = float(self.ui.crackIncrementEdit.text())
        self.input_data.crack_length_end = float(self.ui.crackLengthEndEdit.text())
        
        self.input_data.el_size_factor = float(self.ui.elSizeEdit.text())
        
        self.input_data.load_marker = int(self.ui.loadMarkerEdit.text())
        self.input_data.load = float(self.ui.loadEdit.text())
        self.input_data.load_dir = int(self.ui.loadDirectionEdit.text())
        
        self.input_data.bc_marker = int(self.ui.bcMarkerEdit.text())
        self.input_data.bc = float(self.ui.bcEdit.text())
        self.input_data.bc_dir = int(self.ui.bcDirectionEdit.text())
    
        self.input_data.v = float(self.ui.vEdit.text())
        self.input_data.E = float(self.ui.EEdit.text())
        self.input_data.t = float(self.ui.tEdit.text())
                
        self.input_data.engineering_constants[0] = float(self.ui.E1_edit.text())
        self.input_data.engineering_constants[1] = float(self.ui.E2_edit.text())
        self.input_data.engineering_constants[2] = float(self.ui.E3_edit.text())
        
        self.input_data.engineering_constants[3] = float(self.ui.v12_edit.text())
        self.input_data.engineering_constants[4] = float(self.ui.v13_edit.text())
        self.input_data.engineering_constants[5] = float(self.ui.v23_edit.text())
        
        self.input_data.engineering_constants[6] = float(self.ui.G12_edit.text())
        self.input_data.engineering_constants[7] = float(self.ui.G13_edit.text())
        self.input_data.engineering_constants[8] = float(self.ui.G23_edit.text())
        
        self.input_data.fracture_energy = float(self.ui.fractureEnergyEdit.text())

              
    def onAddPoint(self):
        """ Method for adding point to the geometry
        """
        # If crack hasn't been added to the geometry, adding points is still possible
        if not self.crackAddedFlag:
            
            point = [float(self.ui.addPointEditX.text()), float(self.ui.addPointEditY.text())]
            
            if not point in self.input_data.points:
                
                # Append point from user to point list by calling append point function
                self.input_data.append_point_ui(float(self.ui.addPointEditX.text()), float(self.ui.addPointEditY.text()))
                output_string = "Added point: (" + str(float(self.ui.addPointEditX.text())) + ',' + str(float(self.ui.addPointEditY.text())) + ")"
                
                # Print command data in output text browser
                self.outputTextBrowser.append(output_string)
                
                # Update geometry plot with the new point
                self.onShowGeometry()
                
            else: 
                # Show warning when user tries to add duplicate point
                msgbox = QMessageBox()
                msgbox.setText('It is not possiblele to add duplicate points.')
                msgbox.setStandardButtons(QMessageBox.Ok)
                msgbox.exec()
                
        # If crack already exists in the geometry a warning is displayed that no more points can be added
        else:
            msgbox = QMessageBox()
            msgbox.setText('It is not possiblele to add/remove points after adding a crack to the geometry')
            msgbox.setStandardButtons(QMessageBox.Ok)
            msgbox.exec()
            
     
    def onRemovePoint(self):
        """ Function for undoing the last created point the in the geometry
        """
        # Show error message if no points have been defined yet
        if not self.input_data.points:
            msgbox = QMessageBox()
            msgbox.setText('No points to undo.')
            msgbox.setStandardButtons(QMessageBox.Ok)
            msgbox.exec()
        else:
            point = self.input_data.points[-1]
            del self.input_data.points[-1]
            output_string = "Removed point: (" + str(point[0]) + ',' + str(point[1]) + ")"
            self.outputTextBrowser.append(output_string)
            self.onShowGeometry()
        
        
    def createGeometry(self):
        """ Function for creating CALFEM geometry based on the points and crack geometry supplied by the user
        """
        self.input_data.g = self.input_data.geometry()


    def createMesh(self):
        """ Function for calling the CALFEM mesh creation based on the geometry
        """
        self.solver.create_mesh()

    
    def onShowGeometry(self):  
        """ Function for plotting the geometry when the user wants to display it
        """
        # Create/Update geometry with the current points defined by the user 
        self.createGeometry()
        self.figure.clear()
        
        # If the list of points is empty pass, otherwise plot the figure
        if not self.input_data.points:
            pass
        else:
            self.vis.showGeometry()
        
        self.canvas.draw()
        self.canvas.flush_events()
        self.canvas.update()            


    def onShowMesh(self):
        """ Function for plotting the mesh when the user wants to display it
        """
        output_string = "Creating mesh with element size " + str(self.input_data.el_size_factor) + "..."
        self.outputTextBrowser.append(output_string)
        self.createMesh()
        self.figure.clear()
        self.vis.showMesh()
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()
        self.figure.canvas.update()
        self.outputTextBrowser.append("...Mesh created.")
        

    def onShowDisplacements(self):
        """ Function for plotting the displacements when the user wants to display them
        """
        self.figure.clear()
        self.vis.showDisplacements()
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()
        self.figure.canvas.update()
        
        
    def showEnergyReleaseRate(self):  
        """ Function for plotting the energy release rate or the critical load as a function of crack length
        """
        self.G, self.cracks = self.compliance_object.calculateEnergyReleaseRate()
        self.figure.clear()
        axes = self.figure.gca()
        
        # Try to plot the critical load as a function of crack length
        # and if negative values of G has been acquired and the square
        # not can't be calculated instead show error message
        
        try:
            load = [0.001*sqrt((self.input_data.t*self.input_data.fracture_energy)/G) for G in self.G]
            axes.plot(self.cracks, load)
            axes.set_xlabel('Crack length [m]')
            axes.set_ylabel('Critical load [kN]')
            self.figure.tight_layout()
            self.figure.canvas.draw()
            self.figure.canvas.flush_events()
            self.figure.canvas.update()
        except ValueError:
            msgbox = QMessageBox()
            msgbox.setText('Negative values of crack driving force acquired - check loading and boundary conditions.')
            msgbox.setStandardButtons(QMessageBox.Ok)
            msgbox.exec()
    
    
    def onAddCrack(self):
        """ Function for adding a crack geometry to the list of points
        """
        self.input_data.add_crack(self.input_data.crack_location, self.input_data.crack_length, self.input_data.crack_direction)
        
        output_string = "Crack added to the geometry. Press display geometry to see results."
        self.outputTextBrowser.append(output_string)
        
        
    def onActionRunCalc(self):
        """ Function for running a CALFEM calculation
        """
        # Disable UI and update model parameters
        self.ui.setEnabled(False)
        self.updateModel()
        
        # Single linear elastic calculation
        if self.ui.linearElasticCalc.isChecked():
            self.solverThread = SolverThread(self.solver)
            self.solverThread.finished.connect(self.onSolverFinished)
            self.solverThread.start()
            
       # Linear compliance calculation
        elif self.ui.complianceCalcButton.isChecked(): 
            
            # Instansiate a compliance calculation object
            self.compliance_object = un.ComplianceCalc(self.input_data, self.output_data)
            
            # Reset compliance calculation lists
            self.G = []
            self.input_data_list = []
            self.output_data_list = []
            self.solver_list = []
            self.vis_list = []
            
            # Display warning for small increment size in relation to mesh size
            if (self.input_data.el_size_factor > self.input_data.crack_increment):
                msgbox = QMessageBox()
                msgbox.setText('Using a increment size smaller then the element size might yield inaccurate results.')
                msgbox.setStandardButtons(QMessageBox.Ok)
                msgbox.exec()
            
            # Call a function which creates a list of tuples with input/output/visuaisation/solver objects.
            self.solver_list = self.compliance_object.generateComplianceList()
            
            # Solve all solver objects in the tuple list
            for idx1, solver_pair in enumerate(self.solver_list):
                for solver in solver_pair:
                    self.solverThread = SolverThread(solver)
                    self.solverThread.finished.connect(self.onSolverFinished)
                    self.solverThread.start()
                    
                           
    def onSolverFinished(self):
        """ Function for enabling the UI after the calculation is finished
        """
        self.ui.setEnabled(True)
        
    
    def enableComplianceEdits(self):
        """ Function for enabling buttons depending on radiobuttons clicked
        """
        if self.ui.complianceCalcButton.isChecked():
            self.ui.crackIncrementEdit.setEnabled(True)
            self.ui.crackLengthEndEdit.setEnabled(True)
            self.ui.plotEnergRateButton.setEnabled(True)
        else:
            self.ui.crackIncrementEdit.setEnabled(False)
            self.ui.crackLengthEndEdit.setEnabled(False)
            self.ui.plotEnergRateButton.setEnabled(False)
 
            
    def enableMaterialEdits(self):
        """ Function for enabling buttons for material properties depending on radiobuttons clicked
        """
        if self.ui.IsoButton.isChecked():
            self.input_data.isotropic_flag = True
            self.input_data.orthotropic_flag = False
            
            self.ui.E1_edit.setEnabled(False)
            self.ui.E2_edit.setEnabled(False)
            self.ui.E3_edit.setEnabled(False)
            self.ui.v12_edit.setEnabled(False)
            self.ui.v13_edit.setEnabled(False)
            self.ui.v23_edit.setEnabled(False)
            self.ui.G12_edit.setEnabled(False)
            self.ui.G13_edit.setEnabled(False)
            self.ui.G23_edit.setEnabled(False)
            
            self.ui.vEdit.setEnabled(True)
            self.ui.EEdit.setEnabled(True)
            
        elif self.ui.OrtoButton.isChecked():
            self.input_data.isotropic_flag = False
            self.input_data.orthotropic_flag = True
            
            self.ui.E1_edit.setEnabled(True)
            self.ui.E2_edit.setEnabled(True)
            self.ui.E3_edit.setEnabled(True)
            self.ui.v12_edit.setEnabled(True)
            self.ui.v13_edit.setEnabled(True)
            self.ui.v23_edit.setEnabled(True)
            self.ui.G12_edit.setEnabled(True)
            self.ui.G13_edit.setEnabled(True)
            self.ui.G23_edit.setEnabled(True)
            
            self.ui.vEdit.setEnabled(False)
            self.ui.EEdit.setEnabled(False)
        
    

if __name__ == '__main__':
    
    app = QApplication(sys.argv)
    
    widget = MainWindow()
    widget.show()
    
    sys.exit(app.exec_())