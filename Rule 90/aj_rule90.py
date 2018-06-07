#------------------------------------------------------------------------------#
#-------------------------------- aj_rule90.py --------------------------------#
#------------------------------------------------------------------------------#

"""
Alex's attempt at coding 'Rule 90'.
"""

#---- Imports -----------------------------------------------------------------#

import numpy as np                          #for making a random array
import time

#---- Function ----------------------------------------------------------------#

def rule90(num_val, steps, delta=1, min_val=0, max_val=1):
    """
    Implementation of 'Rule 90'. Returns an array of values that vary over
    a number of time steps based on their neighbouring values.

    num_val: the number of values in the array
    steps:   the number of steps to iterate over
    delta:   the index distance between neigbours
    min_val: the minimum value possible in the initial array
    max_val: the maximum value possible in the initial array
    
    Example:
        arr = rule90(10, 1)
    Returns:
        Initial array: [0 1 0 0 0 1 0 0 0 0]
        Final array:   [1 0 1 0 1 0 1 0 0 0]
    """

    #Create an array of values that are random between the given min and max
    arr = np.random.randint(min_val, max_val+1, size=num_val)

    print ('Initial array:', arr)               #show the initial array

    for t in range(steps):                      #iterate over a number of steps:

        #Create an array of enough zeroes to hold the values for the next step
        next_arr = np.full(num_val, 0)
        
        for v, val in enumerate(arr):           #for each value in array:

            neighbour1 = arr[v-delta]           #get previous neighbour value
            if v > (num_val-1)-delta:           #if near the end of the array:
                neighbour2 = arr[(v+delta)-num_val] #get the next neighbour value (from the beginning of the array)
            else:                               #for the rest of the array:
                neighbour2 = arr[v+delta]       #get next neighbour value
                
            if neighbour1==neighbour2:      #if neighbour values are the same:
                next_arr[v] = 0             #value will become 0 for next step
            else:                           #if neighbour values are different:
                next_arr[v] = 1             #value will become 1 for next step

        arr = next_arr                      #update the array for next step

    print ('Final array:', arr)             #show the final array

    return arr                              #output the final array

#---- Inputs ------------------------------------------------------------------#

num_val = 100                               #number of values
steps = 100                                 #time steps to iterate over

#---- Timing ------------------------------------------------------------------#
start_time = time.clock() 

#---- Output ------------------------------------------------------------------#

#Execute 'Rule 90' and return the final array
arr = rule90(num_val, steps)

#---- Timing ------------------------------------------------------------------#
print("--- %s seconds ---" % (time.clock() - start_time)) 

#---- End ---------------------------------------------------------------------#
