# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 15:40:16 2024

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
    
    my_plot.set_data( diff_sys.labels_t )
    my_plot.chart_it( solve_system.data )
    
    my_plot.set_data( diff_sys.labels_x ) 
    my_plot.chart_it( solve_system.solution_data )

    return
#-------------------------------------------------------

#-------------------------------------------------------
class Define_Differential_System( object ):

    #---------------------------------------
    def __init__( self ):
        
        self.labels_t = [ "Temperature Distribution", "time", "values"]
        self.labels_x = [ "Temperature Distribution", "x-coordinate", "values"]   
        
        self.number_of_equations_and_parameters()
        self.specify_parameters()
    
        self.time_set_form()

        return
    #---------------------------------------
    #---------------------------------------
    def number_of_equations_and_parameters( self ):

        self.numb_of_eq         = 101
        
        return
    #---------------------------------------
    #---------------------------------------
    def specify_parameters( self ):
        
       #---------------------------------------
       # non-dimensional physical, boundary condition, temperature parameters
        
        self.alpha         = 1.0    
        self.h             = 1.0
        self.per_a         = 0.0
        self.hg            = -5.0      
        self.convection    = 0.0

        """
            boundary conditions codes:
            1 = fixed temperature
            2 = fixed temperature gradient ( equivalent to fixed heat flux )
            3 = convection boundary condition ( heat flux = Newton's cooling law )
       """

        self.b_cond_left   = 1
        self.b_cond_right  = 3
               
        self.T_left        = 10.0
        self.T_right       = 20.0
        
        self.T_free_stream = 1.0
        
        self.dT_dx_left    =  10.0
        self.dT_dx_right   =  10.0
        #---------------------------------------
        
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
        
        self.set_parameters()
        self.init_cond_form()
        
        
    #--------------------------------------- 
    #--------------------------------------- 
    def set_parameters( self ):

        #---------------------------------------
        self.b_con = numpy.zeros( 2, dtype = int )
        
        self.b_con[ 0 ] = self.b_cond_left
        self.b_con[ 1 ] = self.b_cond_right     
        
        self.i_a = 1
        self.i_b = self.numb_of_eq - 1
        
        Dx   = 1.0 / float( self.i_b ) 
        Dx_2 = Dx**2.0
                    
        keys   = [ 'h', 'alpha', 'hg', 'per_a', 'u', 'T_left', 'T_right', 'T_free_stream', 'dT_dx_left', 'dT_dx_right', 'Dx', 'Dx_2']
        values = [ self.h, self.alpha, self.hg, self.per_a, self.convection, self.T_left, self.T_right, self.T_free_stream, self.dT_dx_left, self.dT_dx_right, Dx, Dx_2 ]
        
        self.param = dict( zip( keys, values ) )
        
        print( self.param )
        #---------------------------------------
        
        return
    #---------------------------------------
    #---------------------------------------
    def get_parameters( self ):
        
        par = self.param
        
        h     = par[ 'h' ]
        alpha = par[ 'alpha' ]
        hg    = par[ 'hg' ] 
        
        per_a = par[ 'per_a' ]
        
        u      = par[ 'u' ]
        
        T_left        = par[  'T_left' ]
        T_right       = par[ 'T_right' ]                
        T_free_stream = par[ 'T_free_stream' ]
        
        dT_dx_left    = par[ 'dT_dx_left'  ]
        dT_dx_right   = par[ 'dT_dx_right' ]                
                
        Dx    = par[ 'Dx' ]
        Dx_2  = par[ 'Dx_2' ]
        
        return( h, alpha, hg, per_a, u, T_left, T_right, T_free_stream, dT_dx_left, dT_dx_right, Dx, Dx_2 )
    #---------------------------------------
    #---------------------------------------
    def init_cond_form( self ):
    
        self.y_in      =  numpy.zeros(  self.numb_of_eq, dtype = float )
        par = self.param
        
        for i in range( self.numb_of_eq ):
            
            T_x = par[ 'T_left' ] + ( par[ 'T_right' ] - par[ 'T_left' ] ) * float( i ) / float( self.numb_of_eq - 1 )
               
            self.y_in[ i ] = T_x
         
        return
    #---------------------------------------
    #---------------------------------------
    def ode_sys( self, t, y ):  # implement the differential equations in this function
                                
        
        #-----------------------------------------------------
        
        dy_dt = numpy.zeros(  self.numb_of_eq, dtype = float )
        
        h, alpha, hg, per_a, u, T_left, T_right, T_free_stream, dT_dx_left_p, dT_dx_right_p, Dx, Dx_2 = self.get_parameters()        
        
        #-----------------------------------------------------

        for i in range( self.i_a, self.i_b ):
            
            dy_dt[ i ] = alpha * ( y[ i - 1 ] - 2.0 * y[ i ] + y[ i + 1 ] ) / Dx_2 - alpha * u * ( y[ i + 1 ] - y[ i - 1 ] ) / ( 2.0 *  Dx ) + alpha * h * per_a * ( T_free_stream - y[ i ] ) + alpha * hg
        
        
        dy_dt = self.boundary_conditions_ode_sys( t, y, dy_dt )
        #-----------------------------------------------------
        
        print('-> ', end = '')
        
        return( dy_dt )
    #---------------------------------------
    #---------------------------------------
    def boundary_conditions_ode_sys( self, t, y, dy_dt ):
        
        h, alpha, hg, per_a, u, T_left, T_right, T_free_stream, dT_dx_left_p, dT_dx_right_p, Dx, Dx_2 = self.get_parameters()        

        if ( self.b_con[ 0 ] == 1 ):
            dy_dt[ 0 ] = 0.0
        
        elif ( ( self.b_con[ 0 ] == 2 ) or ( self.b_con[ 0 ] == 3 ) ): 
            
            if ( self.b_con[ 0 ] == 2 ) :                
                
                dT_dx_left = dT_dx_left_p           
            else :               
                dT_dx_left = - h * ( T_free_stream - y[ 0 ] )
                
            dy_dt[ 0 ] = alpha * ( ( y[ 1 ] - y[ 0 ] ) / Dx - dT_dx_left ) / ( 0.5 * Dx ) - alpha * u * dT_dx_left + alpha * h * per_a * ( T_free_stream - y[ 0 ] ) + alpha * hg

        

        if ( self.b_con[ 1 ] == 1 ):
            dy_dt[ self.numb_of_eq - 1 ] = 0.0
            
        elif ( ( self.b_con[ 1 ] == 2 ) or ( self.b_con[ 1 ] == 3 ) ):
            
            if ( self.b_con[ 1 ] == 2 ) :                
                
                dT_dx_right = dT_dx_right_p           
            else :                
                dT_dx_right = h * ( T_free_stream - y[ self.numb_of_eq - 1 ] )
                       
            dy_dt[ self.numb_of_eq - 1 ] = alpha * ( dT_dx_right - ( y[ self.numb_of_eq - 1 ] - y[ self.numb_of_eq - 2 ] ) / Dx ) / ( 0.5 * Dx ) - alpha * u * dT_dx_right + alpha * h * per_a * ( T_free_stream - y[ self.numb_of_eq - 1 ] ) + alpha * hg
     
        return( dy_dt )
    #---------------------------------------     
    #---------------------------------------
    def jacob( self, t, y ): # implement the differential-equation jacobian in this function 
                             
        #-----------------------------------------------------
                             
        jac_mtrx = numpy.zeros( ( self.numb_of_eq, self.numb_of_eq ), dtype = float )

        h, alpha, hg, per_a, u, T_left, T_right, T_free_stream, dT_dx_left_p, dT_dx_right_p, Dx, Dx_2 = self.get_parameters() 

        #-----------------------------------------------------

        for i in range( self.i_a, self.i_b ):
            
            #dy_dt[ i ] = alpha * ( y[ i - 1 ] - 2.0 * y[ i ] + y[ i + 1 ] ) / Dx_2 - alpha * u * ( y[ i + 1 ] - y[ i - 1 ] ) / ( 2.0 * Dx ) + alpha * h * per_a * ( T_free_stream - y[ i ] ) + alpha * hg
            jac_mtrx[ i ][ i - 1 ] =         alpha / Dx_2 + alpha * u / ( 2.0 * Dx )
            jac_mtrx[ i ][   i   ] = - 2.0 * alpha / Dx_2 - alpha * h * per_a
            jac_mtrx[ i ][ i + 1 ] =         alpha / Dx_2 - alpha * u / ( 2.0 * Dx )
            
        self.boundary_conditions_jacob( t, y, jac_mtrx )
        #-----------------------------------------------------
        
        return( jac_mtrx )
    #---------------------------------------
    #---------------------------------------
    def boundary_conditions_jacob( self, t, y, jac_mtrx ):
        
        h, alpha, hg, per_a, u, T_left, T_right, T_free_stream, dT_dx_left, dT_dx_right, Dx, Dx_2 = self.get_parameters()        
        
        if ( ( self.b_con[ 0 ] == 2 ) or ( self.b_con[ 0 ] == 3 ) ): 
            
            if ( self.b_con[ 0 ] == 2 ) :
                
                # dT_dx_left = dT_dx_left_p
                hc = 0.0
            else:
                # dT_dx_left = - h * ( T_free_stream - y[ 0 ] )
                hc = h


            #dy_dt[ 0 ] = alpha * ( ( y[ 1 ] - y[ 0 ] ) / Dx - dT_dx_left ) / ( 0.5 * Dx ) - alpha * u * dT_dx_left + alpha * h * per_a * ( T_free_stream - y[ 0 ] ) + alpha * hg
            
            jac_mtrx[ 0 ][ 0 ] =   - alpha / ( 0.5 * Dx_2 ) - alpha * h * per_a  - hc
            jac_mtrx[ 0 ][ 1 ] =     alpha / ( 0.5 * Dx_2 )    
        
      
        if ( ( self.b_con[ 1 ] == 2 ) or ( self.b_con[ 1 ] == 3 ) ) : 

            if ( self.b_con[ 1 ] == 2 ) :
                
                # dT_dx_right = dT_dx_right_p
                hc = 0.0
            else:
                # dT_dx_right = h * ( T_free_stream - y[ self.numb_of_eq - 1 ] )
                hc = h
            

            #dy_dt[ self.numb_of_eq - 1 ] = alpha * ( dT_dx_right - ( y[ self.numb_of_eq - 1 ] - y[ self.numb_of_eq - 2 ] ) / Dx ) / ( 0.5 * Dx ) - alpha * u * dT_dx_right + alpha * h * per_a * ( T_free_stream - y[ self.numb_of_eq - 1 ] ) + alpha * hg
            
            jac_mtrx[ self.numb_of_eq - 1 ][ self.numb_of_eq - 2 ] =     alpha / ( 0.5 * Dx_2 )  
            jac_mtrx[ self.numb_of_eq - 1 ][ self.numb_of_eq - 1 ] =  -  alpha / ( 0.5 * Dx_2 ) - alpha * h * per_a - hc
        
    
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
        
        solution_data = []
        
        fs = open( 'solution.txt', 'wt')
        
        number_of_rows    = len( x.t )       
        number_of_columns = len( x.y )
        
        row = []
        string = ' '
        for i in range( number_of_columns ):
            coord  = float( i ) / float( number_of_columns - 1 )
            string = string + str( coord ) + ' '
            
            row.append( coord )
        
        solution_data.append( row )
        
        
        fs.write( string )
  
        fs.write( '\n' )
        
        coord_string = string
        
            
        for i in range( number_of_rows ):
            string = str( x.t[ i ] ) + ' ' 
            
            row = []
            for y in ( x.y ):
                string = string + str( y[ i ] ) + ' ' 
                
                row.append( y[ i ] )
                
            solution_data.append( row )    
            

            fs.write( string )
      
            fs.write( '\n' )
            
        fs.write( coord_string )
  
        fs.write( '\n' )
        
        for i in range( number_of_rows - 1, number_of_rows ):
            string = ' '  
            for y in ( x.y ):
                string = string + str( y[ i ] ) + ' '  

            fs.write( string )
      
            fs.write( '\n' )

        fs.close()
        
        solution_data = [ [ solution_data[ i ][ j ] for i in range( len( solution_data ) ) ] for j in range ( len( solution_data[ 0 ] ) ) ]
        
        self.solution_data = solution_data
        
      
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



