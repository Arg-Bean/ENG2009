# -*- coding: utf-8 -*-
#Created on Fri Dec 23 10:36:46 2023

#@author: Tom Giglio

#Importing the necessary modules for the necessary arrays, simulations and plots.
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

#defining the unit impulse input, f(t) with an arbitrary amplitude of 10, duration 0.1s from t=0.
def f(t):
    return np.where((0 < t) & (t <= 0.1), 10.0, 0.0)

#Defining the mass-spring damper system, using eqauations 6 and 7. To be used for Euler, Heun and ODEint (exact) Simulations:
def System(y, t):
    #Stating the variables for our system
    x1, x1_prime, x2, x2_prime = y
    #Setting dx1(t)/dt as x1 prime, x1'. And dx2(t)/dt as x2 prime, x2'.
    #So the d^2x1(t)/dt^2 and d^2x2(t)/dt^2 have an easier notation of dx1dt prime. Essentially double prime: x1'' 
    dx1dt = x1_prime
    dx2dt = x2_prime
    dx1dt_prime = -x1_prime - 2 - 4 * x1 + x2_prime #Eq. 6 where d^2x1(t)/dt^2 is the subject
    dx2dt_prime = -x2_prime + f(t) #Eq7. where d^2x2(t)/dt^2 is the subject
    #Then extracting the variables we want
    return [dx1dt, dx1dt_prime, dx2dt, dx2dt_prime]

##Euler Simulation:
    
#Defining the length of t for the Euler Simulation & the amount of data points.
t_Euler = np.linspace(0, 6,60, endpoint=False) #time start at 0s, ends at 6s, with a step size of 6/60=0.1 
#Formula for finding the step size h:
h = t_Euler[1] - t_Euler[0] 
print ('Euler Step size:',h) #Priniting the chosen step size

#For Exact Final Value of x2(t), from Q7.
x2_Exact = np.ones((len(t_Euler)))

#Step size of 0.1 chosen as it provides an accurate simulation, as it is very similar to the ODEint simulation, but also does not take too many data points, so the program can run quickly. smaller step size took too long to run, especialy when running both simulations at once, larger step sizes resulted in inaccurate simulations.
#6s is selected for the end of the end time as the the simulation reaches 1% error value to the steady state value by 5s so one more second is given to show its plateau on the graph.

#Defining Euler Method:
def Euler_Method(System, y0, t_Euler):
    #Creating an empty array for values of y to be stored. 
    y = np.zeros((len(t_Euler), len(y0)))
    # Euler's Method: y[i+1]=yi + h*f(xi,yi)
    for i in range (1, len(t_Euler)): #Setting the time bounds for the system to run
        y[i] = y[i - 1] + h * np.array(System (y[i - 1], t_Euler[i - 1])) 
    return y

#Setting the inital conditions for x1, x2 (Displacements) and x1_prime, x2_prime (velocities). For both simulations
y0 = [0, 0, 0, 0] 

#Running the Euler Method code and storing it in <Euler_Solution>
Euler_Solution = Euler_Method(System, y0, t_Euler)

#From here the values of x1, x1_prime, x2, x2_prime, need to extracted from the array: Euler_Solution.
#At the moment the array has the x1, x1_prime, x2 and the x2_prime values in the columns and the respective t values in the rows. i.e. the array is (4,60)
#The variables and t values need to be swapped so that the array is [x1] 'respective values of x1':[..][..][..], then on the next row: [x2] 'respective values of x2':[..][..][..]
#This way it is in a readable format to plot each variable.
#The variables are flipped by using .T (transpose array), so the array is now (60,4), as needed.
#Extracting the Displacements (x1, x2) and Velocities (x1_prime, x2prime) of x1 and x2. In a useable array
x1_Euler, x1_prime_Euler, x2_Euler, x2_prime_Euler = Euler_Solution.T

#Defining the Error function to be used for the Euler Simulation:
def Error_Euler(x2_Euler):
    #% difference from 1, being the steady state value calculated in Q7.
    ErrE_Percent = (1-x2_Euler)/1 *100
    return ErrE_Percent
#Running the Error_Euler function and storing it in an array: Error_EulerArray
Error_EulerArray = np.array(Error_Euler(x2_Euler))
#Finding when the %Error of the Euler Simulation reahces l%, compared to the steady state value of Q7: 1
onePC_Error_Euler = t_Euler[Error_EulerArray <= 1][0]
#Prining the time at wich the %Error of the Euler method reaches 1%
print('Time at which Euler Simulation %Error drops to 1% =', onePC_Error_Euler,'s')

#Plot 1: Euler Simulation with x2(t)=1 also plotted (Steady State value form Q7)
plt.plot(t_Euler, x2_Euler,label='x2(t)_Euler')
plt.plot(t_Euler, x2_Exact , label= 'Exact Solution', color ='green')#steady state value from Q7.
plt.title('Euler Simulation Graph')
plt.grid(True)
plt.xlabel('Time')
plt.ylabel('Displacement')
plt.legend()
plt.show()

#Plot 2: %Error from steady state value graph for Euler Simulation
plt.plot(t_Euler, Error_Euler(x2_Euler), label = 'Euler Error% plot')
plt.title('Euler Error Graph')
plt.xlabel('Time')
plt.ylabel('Error %')
plt.grid(True)
plt.legend()
plt.show()

#Comments on System Behaviour in Euler Method:
    #x2(t) is overdamped in the Euler simulation as it takes longer than necessary (i.e. critical damping) to reach it's steady state value. It is not overdamped as it never exceeds its final value. 

###Runge-Kutta 2nd Order, Heun Method::
    
#defining the length of t and the amount of points in t, i.e. the step size
t_Heun = np.linspace(0, 6,60, endpoint=False) 
#Formula for finding the step size h:
h  = t_Heun[1] - t_Heun[0] 
print ('Heun Step size:',h)
#time start at 0s, ends at 6s, as the system is within 1% of the steady state value within 5seconds ,so allow for one extra second for visualising on the graph.
#step size is 6/60 =0.1, this is when the simulation gives its most accurate value, viaually comparing to the ODEint graph(green line) without taking excessive data points and overwhelming the program.

#Function for Heun Method:
def Heun_Method(System, y0, t_Heun):
    #Creating an empty array for values of y to be stored. 
    y = np.zeros((len(t_Heun), len(y0)))
    #Heun's Method,:y(i+1) = yi +(0.5k1 + 0.5k2)h. Where k1=f(ti,yi) & k2=f(ti+h, yi+k1*h)
    for i in range(1,len(t_Heun)):
        k1 = np.array(System(y[i-1],t_Heun[i-1]))
        k2 =np.array(System(y[i- 1] + k1*h, t_Heun[i-1] + h))
        y [i] = y [i -1] + h* 0.5* ( k1 +  k2)
    return y
#Initial conditions are already set to 0:@ line 52

#For Exact Final Value of x2(t), from Q7.
x2_Exact = np.ones((len(t_Heun)))

#Running the Heun Method code and storing it in <Heun_Solution>
Heun_Solution = Heun_Method(System, y0, t_Heun)

#extract variables: and put in array form: (60,4) Explained in more detail in Euler simulation lines 57-62
x1_Heun,x1_prime_Heun, x2_Heun, x2_prime_Heun = Heun_Solution.T #Need to unack all 4 variables defined in System

#Defining the Error function for the Heun simulation.
def Error_Heun(x2_Heun):
    ErrH_Percent = (1-x2_Heun)/1 *100
    return ErrH_Percent
#Running the Error_Heun function and storing it in an array: Error_HeunArray
Error_HeunArray = np.array(Error_Heun(x2_Heun))
#Finding when the %Error of the Heun Simulation reahces l%, compared to the steady state value of Q7: 1
onePC_Error_Heun = t_Heun[Error_HeunArray <= 1][0]
#Printing the value of when the %Error reaches 1
print('Time at which Heun Simulation %Error drops to 1% =', onePC_Error_Heun,'s')
#This gives a value of 4.8s this is the first time the %Error of the Heun Simulation dropps below 1%

#Plot 3: Heun Simulation with x2(t)=1 also plotted (Steady State value form Q7)
plt.plot(t_Heun, x2_Heun,label='x2(t)_Heun',color='orange')
plt.plot(t_Heun, x2_Exact , label= 'Exact Solution', color ='green')#steady state value from Q7.
plt.grid(True)
plt.title('Heun Simulation Graph')
plt.xlabel('Time')
plt.ylabel('Displacement')
plt.legend()
plt.show()

#Plot 4: %Error from steady state value graph for Heun Simulation
plt.plot(t_Heun, Error_Heun(x2_Heun),label='Heun Error% plot')
plt.title('Heun Error Graph')
plt.xlabel('Time')
plt.ylabel('Error %')
plt.grid(True)
plt.legend()
plt.show()

#Comments on System Behaviour in Heun Simulaytion:
    #x2(t) is overdamped in the simulation as it takes longer than necessary (i.e. critical damping) to reach it's steady state value. It is not overdamped as it never exceeds its final value. 
    #However it does have a slightly slower response time than the Euler method as it takes 0.2s longer to reach the 1% Error difference. Which is counter intuitive as the Heun method is considered to be more accurate than the Euler method.

#Plotting the System using ODEint, for an comparison of efficacy of simulations.
ODE = integrate.odeint(System, y0 , t_Euler)
#Creating the ODEint simulation using the System, 0 initial conditions, and the same length of time used in the Euler simulation (t_Heun caould also be used)
plt.plot(t_Euler,ODE, label = 'ODEint')
plt.grid(True)
plt.xlabel('Time')
plt.ylabel('Displacement')
plt.title('ODEint Graph')
plt.show()