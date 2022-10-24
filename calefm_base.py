"""
Created on Mon Sep 12 08:26:11 2022
Written by Johannes Jonasson
"""

import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv
import calfem.utils as cfu
import calfem.editor as cfe
import numpy as np
from math import sqrt


class InputData:
    def __init__(self):
        # Crack parameters
        self.crack_length = None
        self.crack_location = None
        self.crack_direction = None
        self.crack_length_end = None
        self.crack_increment = None
 
        # Mesh generation
        self.el_type = 2
        self.dofs_per_node = 2
        self.el_size_factor = 0.05
        
        # Geometry
        self.points = []
        self.original_points = None
        self.splines = []
        self.g = None
        
        # Material parameters
        self.t = 0.01
        self.ptype = 1
        self.ep = [self.ptype, self.t]
        self.E = 210e9
        self.v = 0.3
        self.engineering_constants = [0, 0, 0, 0, 0, 0, 0, 0, 0]
     
        # I/O
        self.filename = 'newfile.json'
        
        # flags
        self.bc_flag = False
        self.load_flag = False
        self.isotropic_flag = False
        self.orthotropic_flag = False
        
        self.bc_marker_list = []
        self.bc_list = []
        self.bc_dir_list = []
        
        self.load_marker_list = []
        self.load_list = []
        self.load_dir_list = []
        
        self.fracture_energy = 0
        
   
    def geometry(self):
        """
        Function for returning a calfem geometry based on the users input data 

        """
        # Create markers for each point and line
        self.point_identifiers = list(range(1,len(self.points)+1))
        self.identifier = list(range(10,10*len(self.points)+10, 10))
        
        # 1 - create empty geometry canvas 
        g = cfg.Geometry()
        
        self.g = g
        
        # 3 - create geometry and points for geometry
        self.create_points(self.points, g)
        
        # 4 - create splines for geometry -> call create_splines()
        self.create_splines(g, self.identifier)
        
        # 5 - create surface with CALFEM function
        g.surface(self.splines2surface)
             
        # 6 - return geometry
        return g


    def create_points(self, points, geometry):
        """
        Function for creating calfem geometry points based on the empty geometry object and the users input points.
        In addition, markers are added for each point to give the user possibility to later define loads for 
        individual points.

        """
        g = geometry
        for idx, position in enumerate(points):
            g.point([position[0], position[1]], marker = self.point_identifiers[idx])
            
            
    def create_splines(self, geometry, identifier):
        """
        Function for creating calfem geometry splines based on the geometry object consisting of points. In addition, 
        it creates markers for the lines to give the user possibility to later define loads/BCs for individual lines.

        """
        self.splines = []
        g = geometry
        
        for idx, point in enumerate(self.points):
            if idx == len(self.points)-1:
                self.splines.append([idx, 0])
            else:
                self.splines.append([idx, idx+1])
      
        for idx, s in enumerate(self.splines):
            g.spline(s, marker = identifier[idx])
            
        self.splines2surface = [i for i in range(len(self.splines))]
    
    
    def add_crack(self, crack_location, crack_length, crack_direction):
        """
        A functino for adding a crack to an arbitrary points in the geometry. The crack consists out of three points,
        2 points at the border and one "inside" of the geometry. All three together creates a sharp crack, with 
        essentially a triangular form (very sharp triangular form)

        """
        # Localize the crack locations index
        idx = crack_location
        
        # Define the point that the crack is going to start in
        crack_begin = [self.points[idx - 1][0], self.points[idx - 1][1]]
  
        # Define new coordinates for the "middle" crack point which is located "inside" the geometry
        end_x = crack_begin[0] + crack_length*crack_direction[0] 
        end_y = crack_begin[1] + crack_length*crack_direction[1]

    
        bcoord = np.array([self.points[idx - 2][0], self.points[idx - 2][1]])
        coord = np.array([self.points[idx - 1][0], self.points[idx - 1][1]])
        
        # Try to set the coordinates of the node in front, and if it doesn't exist
        # use the first point since the only way that it doesn't exist is if the user 
        # is trying to create a crack in the last point
        
        try:
            fcoord = np.array([self.points[idx][0], self.points[idx][1]])
        except IndexError: 
            fcoord = np.array([self.points[0][0], self.points[0][1]]) 
            
        # Calculate direction vectors for the forward and backward point respectively
        v_back = (bcoord - coord)
        dir_back = v_back/(sqrt(sum(v_back**2)))
        dxy_back = 0.00001*dir_back
        
        v_front = (fcoord - coord)
        dir_front = v_front/(sqrt(sum(v_front**2)))
        dxy_front = 0.00001*dir_front
        
        # Replace the pre-existing point with the first out of three crack points
        self.replace_point(idx-1, self.points[idx-1][0] + dxy_back[0], self.points[idx-1][1] + dxy_back[1])
        
        # Add the last out of three points to the geometry
        self.add_point(end_x, end_y, idx)
        
        # Add the crack point that is "inside" the geometry
        self.add_point(crack_begin[0] + dxy_front[0], crack_begin[1] + dxy_front[1], idx+1)

        
    def add_point(self, x, y, idx):
        self.points.insert(idx, [x, y]) 
        
    def replace_point(self, idx, x, y):
        self.points[idx] = [x, y]
        
    def append_point_ui(self, x, y):
        self.points.append([x, y])

        
class OutputData:
    def __init__(self):
        self.a = None
        self.r = None
        self.K = None
        self.coords = None
        self.edof = None
        self.U = None
        self.U_list = []    
        self.f = None
    
class Visualisation:
    def __init__(self, input_data, output_data):
        self.input_data = input_data
        self.output_data = output_data
        
    def showGeometry(self):
        g = self.input_data.g
        
        cfv.draw_geometry(g, draw_points = True, label_points = True, label_curves = True)   
        
        
    def showMesh(self):
        coords = self.output_data.coords
        edof = self.output_data.edof
        dofs_per_node = self.input_data.dofs_per_node
        el_type = self.input_data.el_type
    
        cfv.draw_mesh(coords, edof, dofs_per_node, el_type, filled = True)
        
    def showDisplacements(self):
        a = self.output_data.a
        coords = self.output_data.coords
        edof = self.output_data.edof
        dofs_per_node = self.input_data.dofs_per_node
        el_type = self.input_data.el_type

        cfv.draw_displacements(a, coords, edof, dofs_per_node, el_type, draw_undisplaced_mesh=False, title="Displacements")
        
        
class ComplianceCalc:
    def __init__ (self, input_data, output_data):
        self.input_data = input_data
        self.output_data = output_data
        self.crack_length = self.input_data.crack_length # Start length
        self.crack_increment = self.input_data.crack_increment 
        self.crack_length_end = self.input_data.crack_length_end # End length
    
    def calculateEnergyReleaseRate(self):
        """ Function for calculating the energy release rate for each step in the compliance calculation
        """
        G = []
    
        for output_pair in self.output_data_list:
            print(output_pair[1].U, output_pair[0].U)
            G.append(((output_pair[1].U - output_pair[0].U)/self.input_data.crack_increment))
            
        return G, self.cracks
        

    def generateComplianceList(self):
        """ Function creating a list of tuples of input/output/visualisation/solver objects
        """
        # Initialise crack propagation features
        self.input_data_list = []
        self.output_data_list = []
        self.geometry_list = []
        self.solver_list = []
        self.vis_list = []
        
        # Create cracks in between the user start and end arguments
        self.cracks = np.linspace(self.crack_length, self.crack_length_end, 25) 
        
        crack_list = []
        
        # Create a list of tuples where all cracks and respective cracks incremented crack are a tuple
        for crack in self.cracks:    
            crack_list.append((crack, crack + self.crack_increment))
            
        for crack_tuple in crack_list:
            # For each crack tuple
            temp_list_solver = []
            temp_list_output = []
            temp_list_input = []
            temp_list_vis = []
            
            for crack in crack_tuple:
                # For each crack and incremented crack in crack tuple
                # instantsiate input data object
                input_data = InputData() 
                
                # Copy input data from the user defined object
                input_data.points = self.input_data.points.copy()
                input_data.crack_length = crack
                
                input_data.crack_direction = self.input_data.crack_direction 
                input_data.crack_location = self.input_data.crack_location
                input_data.el_size_factor = self.input_data.el_size_factor 
                
                input_data.orthotropic_flag = self.input_data.orthotropic_flag
                input_data.isotropic_flag = self.input_data.isotropic_flag
                
                input_data.engineering_constants = self.input_data.engineering_constants
                
                input_data.bc_marker_list = self.input_data.bc_marker_list
                input_data.bc_list = self.input_data.bc_list
                input_data.bc_dir_list = self.input_data.bc_dir_list
                input_data.load_marker_list = self.input_data.load_marker_list
                input_data.load_list = self.input_data.load_list
                input_data.load_dir_list = self.input_data.load_dir_list
                
                input_data.add_crack(input_data.crack_location, input_data.crack_length, input_data.crack_direction)
                input_data.g = input_data.geometry()
                
                # Instansiate output data object and solver object
                output_data = OutputData()
                
                solver = Solver(input_data, output_data)
                
                # Create mesh for each object
                solver.create_mesh()
                
                # Create visulatisation object
                vis = Visualisation(input_data, output_data)
                
                temp_list_solver.append(solver)
                temp_list_output.append(output_data)
                temp_list_input.append(input_data)
                temp_list_vis.append(vis)
                
            # Append objects to global list
            self.solver_list.append(temp_list_solver)
            self.output_data_list.append(temp_list_output)
            self.input_data_list.append(temp_list_input)
        
        return self.solver_list


class Solver:
    def __init__(self, input_data, output_data):
        self.input_data = input_data
        self.output_data = output_data
        
    def create_mesh(self):
        """ Function for creating CALFEM mesh
        """
        
        g = self.input_data.g
            
        # Define mesh
        
        mesh = cfm.GmshMeshGenerator(g)
        
        # Define mesh attributes
        
        mesh.el_type = self.input_data.el_type
        mesh.dofs_per_node = self.input_data.dofs_per_node
        mesh.el_size_factor = self.input_data.el_size_factor
    
        coords, edof, dofs, bdofs, element_markers = mesh.create()
        
        self.output_data.coords = coords
        self.output_data.edof = edof
        self.output_data.dofs = dofs
        self.output_data.bdofs = bdofs
        self.output_data.element_markers = element_markers

        
    def createLoad(self):
        """ Function for creating boundary conditions based on the user input boundary conditions
        """
        bdofs = self.output_data.bdofs
        nDofs = np.size(self.output_data.dofs)
  
        load_marker_list = self.input_data.load_marker_list
        load_list = self.input_data.load_list
        load_dir_list = self.input_data.load_dir_list
        
        # Not sure this is needed anymore
        if not self.input_data.load_flag:
            self.output_data.f = np.zeros([nDofs, 1])
            self.input_data.load_flag = True
        
        for load_marker, load, load_dir in zip(load_marker_list, load_list, load_dir_list):
            cfu.applyforce(bdofs, self.output_data.f, load_marker, load, load_dir)
                
    def createBC(self):
        """ Function for creating boundary conditions based on the user input boundary conditions
        """
        bdofs = self.output_data.bdofs
        bc_marker_list = self.input_data.bc_marker_list
        bc_list = self.input_data.bc_list
        bc_dir_list = self.input_data.bc_dir_list
        
        # Not sure this if statement is needed anymore
        if not self.input_data.bc_flag:
            self.bc = np.array([], int)
            self.bcVal = np.array([], float)
            self.input_data.bc_flag = True
            
            
        # Create boundary conditions
        for bc_marker, bc, bc_dir in zip(bc_marker_list, bc_list, bc_dir_list):
            self.bc, self.bcVal = cfu.applybc(bdofs, self.bc, self.bcVal, bc_marker, bc, bc_dir)

 
    def execute(self):
        """ Function that performs a finite element calculation and also calculates the strain energy
            for the simulation
        """
        
        self.createBC()
        self.createLoad()
        
        f = self.output_data.f
        bc = self.bc
        bcVal = self.bcVal
        
        ep = self.input_data.ep
        E = self.input_data.E
        ptype = self.input_data.ptype
        v = self.input_data.v
               
        nDofs = np.size(self.output_data.dofs)        
        
        if self.input_data.isotropic_flag:
 
            D = cfc.hooke(ptype, E, v)
            
        elif self.input_data.orthotropic_flag:
            
            D = self.plan2eortho(self.input_data.engineering_constants)

        ex, ey = cfc.coordxtr(self.output_data.edof, self.output_data.coords, self.output_data.dofs)
                
        K = np.zeros([nDofs,nDofs])

        for eltopo, elx, ely in zip(self.output_data.edof, ex, ey):
            Ke = cfc.plante(elx, ely, ep, D)
            cfc.assem(eltopo, K, Ke)
                
        a,r = cfc.solveq(K,f,bc,bcVal)
                
        # Store results in output data class
        
        self.output_data.a = a
        self.output_data.r = r
        self.output_data.K = K 
      
        a = self.output_data.a
        coords = self.output_data.coords
        edof = self.output_data.edof
        dofs_per_node = self.input_data.dofs_per_node
        el_type = self.input_data.el_type
                
        U = np.transpose(a)*K*a
        self.output_data.U = U.item()
        
    def plan2eortho(self, engineering_constants):
        """
        Parameters
        ----------
        engineering_constants : constants describing the constitutive behaviour of ortotropic material

        Returns
        -------
        D : Returns the constitutive D-matrix for an orthotropic material.

        """
        
        E1 = engineering_constants[0]
        E2 = engineering_constants[1]
        E3 = engineering_constants[2]
        
        v12 = engineering_constants[3]
        v13 = engineering_constants[4]
        v23 = engineering_constants[5]
        
        G12 = engineering_constants[6]
        G13 = engineering_constants[7]
        G23 = engineering_constants[8]
        
        # Assemble compliance matrix
        
        C = np.array([[1/E1, v12/E2, v13/E3, 0, 0, 0],
                      [-v12/E1, 1/E2, -v23/E3, 0, 0, 0],
                      [-v13/E1, -v23/E2, 1/E3, 0, 0, 0],
                      [0, 0, 0, 1/(2*G23), 0, 0],
                      [0, 0, 0, 0, 1/(2*G13), 0],
                      [0, 0, 0, 0, 0, 1/(2*G12)]])
        
        # Inverse compliance matrix
        
        D =  np.linalg.inv(C)
        
        return D
        
    
################### MAIN STARTS HERE #####################   

# Initialize input and output data classes


#rack_length_list = np.linspace(0.1, 0.15, 1)

# if __name__ == "__main__":
    
#     U = []
    
#     #for idx, crack_length in enumerate(crack_length_list):   
#     output_data = OutputData()
#     input_data = InputData()
#     input_data.crack_length = 0.25
#     input_data.crack_location = 2
#     input_data.crack_increment = 0.01
#    # input_data.points = [[0,0], [0.5,0], [1,0], [1,0.25], [0.5,0.25], [0,0.25]]
#     input_data.points = [[0,0], [0.3, 0], [0.3,-0.1], [1.5,-0.1], [1.5, 0], [1.5,0.3], [0.5, 0.3], [0,0.3]]
    
#     # BOUNDARY CONDITIONS 
    
#     bc_identifier = [10]*len(input_data.points)
    
#     bcmarker = 20
#     bcmarker2 = 40
#     loadmarker = 30
    
#     bc_points = [5,7]
#     bc_identifier[bc_points[0]:bc_points[1]] = [bcmarker]*(bc_points[1]-bc_points[0])
    
#     bc2_points = [1,0]
#     bc_identifier[bc2_points[0]:bc2_points[1]] = [bcmarker2]*(bc2_points[1]-bc2_points[0])
    
#     # Load conditions
    
#     load_points = [7,9]
    
#     bc_identifier[load_points[0]:load_points[1]] = [loadmarker]*(load_points[1]-load_points[0])
    
#     # Initialize solver
    
#     solver = Solver(input_data, output_data)
    
#     # Execute solver
    
#     solver.execute()
    
#     vis = Visualisation(input_data, output_data)
#     vis.showGeometry()
#     vis.showMesh()
#     vis.showDeformedMesh()

#     # U.append(np.transpose(output_data.a)*output_data.K*output_data.a)
    
#     # print('The strain energy is: ', U[idx])
    
#     # G = (U[1] - U[0])/0.001
    
#     # mu = input_data.E/(2*(1+input_data.v))
#     # kappa = (3-input_data.v)/(1+input_data.v)
#     # K_1 = sqrt((8*mu*G)/(1 + kappa))
    
#     # print('The energy release rate is', G)
#     # print('The stress intensity factor is', K_1)