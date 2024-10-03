# -*- coding: utf-8 -*-

#-----------------------------------------------------------------------------


class Print_Department(object):
    def __init__(self):
        self.printer()
        self.variable2 = 5
    
    def printer(self):
        
        variable = 'hello world'
        print(variable)
        
        variable = 45.67 
        self.fun2(variable)
        
        variable = self.fun2
        
        #variable("again")
    
    #-----------------------------------------------------------------------------

    def fun2(self, variable):
        
        print(variable)
        
        return
        
#-----------------------------------------------------------------------------
def main():
    
    Print_Department()
    
    
    return
    
#-----------------------------------------------------------------------------
main()

#-----------------------------------------------------------------------------