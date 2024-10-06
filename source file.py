import numpy as np
import random
from scipy.linalg import expm


def dotproduct(state):
    product=state[0]
    for counter in range(1,len(state)):
        product=np.dot(product,state[counter])
    return(product[0])


hamiltonian=[]
for i in range (100):
    row=[]
    for j in range (100):
        if i==j:
            element=np.random.normal(0,2/np.sqrt(3))
        elif j==i+1 or j==i-1 :
            element=1
        else:
            element=0
        row.append(element)   
    hamiltonian.append(row)        

def unitary_operator(hamiltonian):
    operator=expm(complex(0,-10)*np.array(hamiltonian))
    return operator
unitary=unitary_operator(hamiltonian)

def initial_state(n):
    initial_state=[]
    
    for counter in range (n):
        if counter==n/2:
            initial_state.append(1)
        else:
            initial_state.append(0)
    return np.array(initial_state)
initial_state=initial_state(100)


def projection_operator_present(position_measurement):
    operator=[]
    for i in range (100):
        row=[]
        for j in range (100):
            if i==j==position_measurement:
                element=1
            else:
                element=0
            row.append(element)
        operator.append(row)   
    return np.array(operator)

def projection_operator_absent(position_measurement):
    operator=[]
    for i in range(100):
        row=[]
        for j in range(100):
            if i==j!=position_measurement:
                element=1
            else:
                element=0
            row.append(element)
        operator.append(row)      
    return np.array(operator)


# finding the sqrt of probability amplitude of wave function at the measurement site
def probabilty_of_projection(position,state):
    projection_operator=projection_operator_present(position)
    element_1=np.dot(np.matrix.getH(state),np.matrix.getH(projection_operator))
    element_2=np.dot(projection_operator,state)
    probability=np.sqrt(np.real(np.dot(element_1,element_2)))
    return (probability)

# determines the outcome of measurement, particle present or absent, based on a special random number generator
def projection_operator(position_measurement,state):
    measurement_trajectory=[]
    probability_threshold=random.uniform(0,1)
    comparison_probability=probabilty_of_projection(position_measurement, state)
    if (probability_threshold<comparison_probability**2):
        operator=projection_operator_present(position_measurement)
        measurement_trajectory.append([position_measurement,'present'])
    else:
        operator=projection_operator_absent(position_measurement)  
        measurement_trajectory.append([position_measurement,'absent'])
    print(measurement_trajectory)
    return operator


# probability of the particle being present at each site, for all generation
def probabilty_site(state_list):
    matrix=[]
    for counter_time in range(len(state_list)):
        row=[]
        for counter_position in range(100):
           probability=(probabilty_of_projection(counter_position, state_list[counter_time])) 
           if probability<10**(-5):
               element=-5
           else:
               element=np.log10(probability)
           row.append(element)     
        matrix.append(row)
    return matrix


def state_list_with_measurements(initial_state):
    state_list=[initial_state]
    for counter_measurement in range(1000):       
        for counter_evolution in range (10):
            counter=len(state_list)-1
            evolved_state=np.dot(unitary,state_list[counter])
            state_list.append(evolved_state)
        position_measurement=random.randint(0, 99)
        counter=len(state_list)-1
        state=state_list[counter]
        measurement_operator=projection_operator(position_measurement, state)
        measured_state=np.dot(measurement_operator,state_list[counter])
        norm=np.sqrt(np.dot(np.matrix.getH(measured_state),measured_state))
        measured_state_normalised=measured_state/norm
        state_list.append(measured_state_normalised)
    return state_list

def r(state_list,power):
    average_distance_list=[]
    for counter_time in range(len(state_list)):
        average_distance_each_gen=0
        for counter_position in range(100):
            moment=counter_position**power
            probability=probabilty_of_projection(counter_position, state_list[counter_time])
            element=moment*probability*probability
            average_distance_each_gen=average_distance_each_gen+element
        if power==1:
            average_distance_list.append(round(average_distance_each_gen))
        else:
            average_distance_list.append(average_distance_each_gen)
    return average_distance_list
      