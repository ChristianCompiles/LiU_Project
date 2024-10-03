# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 12:49:43 2024

@author: joseph.iannelli
"""
#-------------------------------------------------------
"""

    This program solves an arbitrary set of first-order
    ordinary differential equations

    All that the user would need to do is to implement the right-hand sides
    of the differential equations in the 'ode_sys' function and the
    right-hand side jacobians in the 'jacob' function, both in the
    My_Differential_System() class

    The parameters of the differential equations are 
    stored in a 'self.param[]' array

    The remaining pieces of needed information are requested of the user at run time
"""
#-------------------------------------------------------
#-------------------------------------------------------
import numpy
import math # it may be needed for some differential equations
from   scipy.integrate import solve_ivp

import matplotlib.pyplot as plt


#-------------------------------------------------------
#-------------------------------------------------------
def main():
      
    # Class objects
    diff_sys        = My_Differential_System()   
    solve_system    = Solve_Differential_System()
      
    # object utilization
    solve_system.solve_sys( diff_sys )
      
    solve_system.print_to_file()
    #solve_system.print_to_screen()
    
    my_plot = Plot()
    my_plot.set_data( diff_sys.labels )
    
    my_plot.chart_it( solve_system.data )

    return
#-------------------------------------------------------
#-------------------------------------------------------
class Define_Differential_System( object ):

    #---------------------------------------
    def __init__( self ):
        
        self.labels = [ "Differential System Solution", "time", "values"]
        
        self.number_of_equations()
        self.specify_parameters()
    
        self.time_set_form()

        return
    #---------------------------------------
    #---------------------------------------
    def number_of_equations( self ):

        self.numb_of_eq         = 50
        
        
        return
    #---------------------------------------
    #---------------------------------------
    def specify_parameters( self ):
        
       #---------------------------------------
       # non-dimensional parameters and keys
        
        #lambd = 2.0
        #lambd = np.array(self.numb_of_eq)
     
        #values = [  lambd  ]
        #keys   = [ 'lambd' ]
        
        #---------------------------------------        
        #self.param = dict( zip( keys, values ) )
        #---------------------------------------
        
        self.param = [float(i) for i in range (self.numb_of_eq)]
        self.dx = 1.0/float(self.numb_of_eq - 1)
        
        self.k = 2.0
        self.source = 1.5
        
        return
    #---------------------------------------
    #---------------------------------------
    def time_set_form( self ):
    
        print()
        self.t_begin     =  float( input( 'Enter initial time: ' ) )
        self.t_end       =  float( input( 'Enter final time: ' ) )
        self.t_stations  =  int(   float( input( 'Enter number of time stations: ') ))
        
        self.delta_t = ( self.t_end - self.t_begin ) / self.t_stations

        self.time = numpy.linspace( self.t_begin, self.t_end, self.t_stations )

        self.time_range = numpy.zeros( 2, dtype = float )
        self.time_range[ 0 ] = self.t_begin
        self.time_range[ 1 ] = self.t_end

        return
    #---------------------------------------
    #---------------------------------------
    def init_cond( self ):
    
        return( self.y_in )

    #---------------------------------------
    #---------------------------------------
    def time_set( self ):
    
        return( self.time_range, self.time )
    #---------------------------------------
    
#-------------------------------------------------------
#-------------------------------------------------------
class My_Differential_System( Define_Differential_System ):
    
    #--------------------------------------- 
    def __init__( self ):
        
        #Define_Differential_System.__init__( self )
        super( My_Differential_System, self ).__init__()        
        
        self.init_cond_form()
        
        
    #--------------------------------------- 
    #---------------------------------------
    def init_cond_form( self ):
    
        self.y_in      =  numpy.zeros(  self.numb_of_eq, dtype = float )
        y = self.y_in 
        
        #----------------------------------------------------------------------
        # set initial conditions
        
        y = [(float(i)  / float(self.numb_of_eq - 1.0))**2 for i in range (self.numb_of_eq)]
        
        #----------------------------------------------------------------------
        #----------------------------------------------------------------------              
        self.y_in = y
        #----------------------------------------------------------------------
         
        return
    #---------------------------------------
    #---------------------------------------
    def ode_sys( self, t, y ):  # implement the differential equations in this function
                                # parameters may be specified through the self.param[] array
        
        #-----------------------------------------------------
        
        dy_dt = numpy.zeros(  self.numb_of_eq, dtype = float )
        par = self.param
        dx = self.dx
        
        # lambd = par[ 'lambd' ]
 
        self.dx_2 = dx**2
        dx_2 = self.dx_2
        
        k = self.k
        source = self.source
        
        #-----------------------------------------------------
       
        dy_dt[0] = 0.0
        
        for i in range(1, self.numb_of_eq-1):
            #dy_dt[ i ] =  - par[i] * y[ i ]
            #dy_dt[ i ] =  -(y[ i ] - y[i-1]) / dx
            dy_dt[i] = k * (y[i-1] - 2.0 * y[i] + y[i+1]) / dx_2  - source

        dy_dt[self.numb_of_eq-1] = 0.0
        #-----------------------------------------------------
        
        print('-> ', end = '')
        
        return( dy_dt )
    #---------------------------------------
    #---------------------------------------
    def jacob( self, t, y ): # implement the differential-equation jacobian in this function 
                             # parameters may be specified through the self.param[] array
        #-----------------------------------------------------
                             
        jac_mtrx = numpy.zeros( ( self.numb_of_eq, self.numb_of_eq ), dtype = float )
        par = self.param
        dx = self.dx
        dx_2 = self.dx_2

        #lambd = par[ 'lambd' ]
        
        k = self.k
        source = self.source
        #-----------------------------------------------------
        #jac_mtrx[0][0] = 0.0
        
        #dy_dt[ i ] =  -(y[ i ] - y[i-1])
        #for i in range(1, self.numb_of_eq):
            #jac_mtrx[i][i-1] = 1.0 / dx
           # jac_mtrx[i][i]  = - (1.0) / dx

        #jac_mtrx[ 0 ][ 0 ] = - lambd
        
        #dy_dt[i] = k * (y[i-1] - 2.0 * y[i] + y[i+1]) / dx_2  - source
        
        #dy_dt[0] = 0.0
        jac_mtrx[0][0] = 0.0
        for i in range(1, self.numb_of_eq-1):
            jac_mtrx[i][i-1] = k * (1.0 / dx_2)
            jac_mtrx[i][i]  = - 2*k / dx_2
            jac_mtrx[i][i+1] = k *  (1.0) / dx_2
            
        jac_mtrx[self.numb_of_eq-1][self.numb_of_eq -1] = 0.0
        
        #-----------------------------------------------------
        
        return( jac_mtrx )
    #---------------------------------------
#-------------------------------------------------------
#-------------------------------------------------------
class Solve_Differential_System( object ):

    #---------------------------------------
    def __init__( self ):
        
        self.comp_sol = []
         
        return
    #---------------------------------------
    #---------------------------------------
    def solve_sys( self, diff_sys ):

        time_range, time = diff_sys.time_set()
        
        y_in = diff_sys.init_cond()
 
        print('...working...')
        self.comp_sol = solve_ivp( diff_sys.ode_sys, time_range, y_in, method = 'Radau', t_eval = time, dense_output = True, atol = 1.0e-13, rtol = 1.0e-13, jac = diff_sys.jacob )
        
        
        self.data_matrix()
       
        return
    #---------------------------------------   
    #---------------------------------------
    def data_matrix( self ):    

        self.data = numpy.zeros( ( len( self.comp_sol.t) , len( self.comp_sol.y ) + 1 ), dtype = float )
        
        number_of_rows = len( self.data )
        
        for i in range( number_of_rows ):
            self.data[ i ][ 0 ] = self.comp_sol.t[ i ]
            j = 1
            for var in ( self.comp_sol.y ):
                self.data[ i ][ j ] = var[ i ]
                j = j + 1
                
            

        return
    #---------------------------------------
    #---------------------------------------
    def print_to_file( self ):
        
        x = self.comp_sol
        
        number_of_rows      = len( x.t )
        
        fs = open( 'solution.txt', 'wt')

        for i in range( number_of_rows ):
            string = str( x.t[ i ] ) + ' '  
            for y in ( x.y ):
                string = string + str( y[ i ] ) + ' '  
                
           

            
            fs.write( string )
      
            fs.write( '\n' )

        fs.close()

        return
    #---------------------------------------    
    #---------------------------------------
    def print_to_screen( self ):

        x = self.comp_sol
   
        number_of_rows      = len( x.t )
        print()
        
       
        for i in range( number_of_rows ):
            string = str( x.t[ i ] ) + ' '  
            for y in ( x.y ):
                string = string + str( y[ i ] ) + ' '  
	
            print( string )

        return
    #---------------------------------------
    
#-------------------------------------------------------
#-------------------------------------------------------- 
class Plot( object ):
       
    #----------------------------------------------------
    def __init__( self ):
        
        self.plt = plt
    
        self.data = []  
        
        self.labels = ["Title", "h-label", "v-label"]
        
        return
    #---------------------------------------------------- 
    #----------------------------------------------------    
    def set_data( self, data ):
        
        self.labels = data
        
        return

    #---------------------------------------------------- 
    #----------------------------------------------------    
    def chart_it( self, data ):
        
        self.plt = plt
    
        self.data = []  
            
        self.i_plot = True
        
        while ( self.i_plot == True ):
            
            self.chart_it_now( data )
            
            
        
        return
    #---------------------------------------------------- 
    #----------------------------------------------------    
    def chart_it_now( self, data ):
        

        number_of_rows    = len( data )
        number_of_columns = len( data[ 0 ] )
        
        #-------------------------------------------       
        def function_to_plot():
                
            i   = number_of_columns - 1
                
            if ( i > 1 ):
                
                i_msg_1 = "There are " + str( i ) + " functions in this data set \n"
                i_msg_2 = "Enter the number of the function to plot, between 1 and " + str(i) +"\n" 
                i_msg_3 = "Instead, enter " + str( i + 1 ) + " to plot all \n"
                i_msg_4 = "Otherwise, enter '0' to exit:  "
                i_msg   = i_msg_1 + i_msg_2  + i_msg_3 + i_msg_4

                x_sel = 1

            else:
                
                i_msg_1 = "There is 1 function in this data set \n"
                i_msg_2 = "Enter 1 to plot it \n"
                i_msg_4 = "Otherwise, enter '0' to exit: "
                i_msg   = i_msg_1 + i_msg_2 + i_msg_4
                
                x_sel = 0             
            
            print()
            self.ifunct = int( input( i_msg  ) )
            print()

            if ( ( self.ifunct == 0 ) or ( self.ifunct < 1 ) or ( self.ifunct > i + x_sel ) ):
                    
                self.i_plot = False
                    
            else:
                    
                if ( self.ifunct == i + x_sel ):
                        
                    self.plot_all = True
                    
                else:
                        
                    self.plot_all = False
             
            return
        #-------------------------------------------
        #-------------------------------------------
        def data_preparation():
                
            self.data = [ [ data[ i ][ j ] for i in range( number_of_rows ) ] for j in range( number_of_columns ) ]
            
            self.number_of_functions = number_of_columns
            
            return
        #-------------------------------------------
            
        function_to_plot()
        
        if ( self.i_plot == True ):
        
            data_preparation()
           
            self.init_plot()
            self.plot_it()
            
        return
    
    #----------------------------------------------------
    #----------------------------------------------------
    def init_plot( self ):
        
        # https://matplotlib.org/stable/api/markers_api.html
       
        """
        ls = '-' 	solid    line style
        ls = '--' 	dashed   line style
        ls = '-.' 	dash-dot line style
        ls = ':' 	dotted   line style
        """
        
        """
        'r' 	Red 	
        'g' 	Green 	
        'b' 	Blue 	
        'c' 	Cyan 	
        'm' 	Magenta 	
        'y' 	Yellow 	
        'k' 	Black 	
        'w' 	White
        """
        
        self.fig, self.ax = self.plt.subplots()
        
        self.ax.grid( True )
        self.fig.suptitle(  self.labels[ 0 ] )
        self.ax.set_xlabel( self.labels[ 1 ] ) ; self.ax.set_ylabel( self.labels[ 2 ] )
        
        if ( self.plot_all == True ):
            
            for self.ifunct in range( 1, self.number_of_functions ):
                
                self.ax.plot( self.data[ 0 ], self.data[ self.ifunct ],  ls = "-", linewidth = 1, color = 'g', marker = 'o', ms = 4, mec = 'b', mfc = 'b'  )
        
        else :
            
            self.ax.plot( self.data[ 0 ], self.data[ self.ifunct ],  ls = "-", linewidth = 1, color = 'g', marker = 'o', ms = 4, mec = 'b', mfc = 'b'  )        
        
        return
    #----------------------------------------------------         
    #----------------------------------------------------  
    def plot_it( self ):
        
        self.plt.show()
        
        return
    #----------------------------------------------------      
#-------------------------------------------------------
#-------------------------------------------------------
main()
#-------------------------------------------------------



