import numpy as np
import matplotlib.pyplot as plt

#STEP 1

def thetalist(N=50,bottom=-np.pi/2,top=np.pi/2):
    thetalist = np.linspace(bottom,top,N)
    return thetalist

def h(theta0=0, thetalist=thetalist(), c=3.0, epsilon=0.1):
    h = c*((1-epsilon)+epsilon*np.cos(2*(thetalist - theta0)))
    return h

def g(h_i=h(),T=0,beta=0.1):
    g =[]
    for v in h_i:
        if v<=T:
            g.append(0)
        elif T<v<=T+1/beta:
            g.append(beta*(v-T))
        elif v>T+1/beta:
            g.append(1)            
    return np.array(g)

def plot1():
    plt.close()
    plt.figure("h(theta)")
    plt.xlabel('theta')
    plt.ylabel('h')
    plt.plot(thetalist(),h(),'b')
    
    hrange= np.linspace(-15,15,300)
    plt.figure("g(h)")
    plt.xlabel("h")
    plt.ylabel("g(h)")
    plt.plot(hrange,g(h_i=hrange),'b')
    plt.show()

#STEP 2

def ratebasednm(h_i=h(),m0=np.zeros(len(h())),deltat=1.0,tau=5.0,N_i=30):
    m_list = [m0]
    m=m0
    g_h_i=g(h_i=h_i)
    
    #Explicit Euler method to solve
    for i in range(N_i):
        m = m + deltat/tau*(-m+g_h_i)
        m_list = np.append(m_list,m)
        
    twoD = m_list.reshape((N_i+1,50))
    
    #plot out the solution for all
    fig = plt.figure()
    axcolmap = fig.add_subplot(111)
    axcolmap.set_title('Neuron Activity', y=1.02)
    plt.imshow(twoD,extent=[-np.pi/2,np.pi/2,30,0],interpolation='none')
    axcolmap.set_aspect('auto')
    axcolmap.set_ylabel('Time t (steps)')
    axcolmap.set_xlabel('Position x')
    plt.colorbar(orientation='vertical')
    return m
    
def thalamus(theta0=0, m0=np.zeros(len(h())), N=50, N_i=30, c=1.5, epsilon=0.9):
    h_i = h(theta0=theta0, thetalist=thetalist(N=N), c=c, epsilon=epsilon)
    m = ratebasednm(h_i=h_i,m0=m0,N_i=N_i)     
    return m
    
#STEP 3

def J_ij(thetalist=thetalist(), J0=86.0, J2=112.0):
    xx,yy = np.meshgrid(thetalist,thetalist,indexing='ij')
    J_ij = -J0+J2*np.cos(2*(xx-yy))
    return J_ij

def plotJ(J_ij=J_ij()):
    fig = plt.figure()
    axcolmap = fig.add_subplot(111)
    axcolmap.set_title('Neuron connections', y=1.02)
    plt.imshow(J_ij,extent=[-np.pi/2,np.pi/2,np.pi/2,-np.pi/2])
    axcolmap.set_aspect('auto')
    axcolmap.set_ylabel('Theta_i')
    axcolmap.set_xlabel('Theta_j')
    plt.colorbar(orientation='vertical')
    
def thalamusfull(theta0=0, m0=np.zeros(len(h())), N=50, N_i=30, c=1.2, epsilon=0.9, deltat=1.0, tau=5.0, hflag=True):
    m_list = [m0]
    m=m0
    if hflag == True:
        h_ext=h(theta0=theta0, thetalist=thetalist(N=N), c=c, epsilon=epsilon)
    elif hflag ==False:
        h_ext=np.zeros(N)
    
    #Explicit Euler method to solve
    for i in range(N_i):
        interactions = J_ij()*m
        h_i = np.sum(interactions,axis=1)+h_ext
        g_h_i =g(h_i)
        m = m + deltat/tau*(-m+g_h_i)
        m_list = np.append(m_list,m)
        
    twoD = m_list.reshape((N_i+1,N))
    
    #plot out the solution for all
    fig = plt.figure()
    axcolmap = fig.add_subplot(111)
    axcolmap.set_title('Neuron Activity', y=1.02)
    plt.imshow(twoD,extent=[-np.pi/2,np.pi/2,30,0],interpolation=None)
    axcolmap.set_aspect('auto')
    axcolmap.set_ylabel('Time t (steps)')
    axcolmap.set_xlabel('Position x')
    plt.colorbar(orientation='vertical')
    plt.show()
    return m     

#STEP 4 - NEED TO DISCUSS COMMENTS!
def change(theta0_1=0,theta0_2=2.0*np.pi/3.0,N_i1=30,N_i2=500,c=100,epsilon=0.8):
    out = thalamusfull(theta0=theta0_1, m0=np.zeros(len(h())), N=50, N_i=N_i1, c=c, epsilon=epsilon, deltat=1.0, tau=5.0)
    thalamusfull(theta0=theta0_2,m0=out, N=50, N_i=N_i2, c=c, epsilon=epsilon, deltat=1.0, tau=5.0)

def changenoconnectivity(theta0_1=0,theta0_2=2.0*np.pi/3.0,N_i1=30,N_i2=500,c=100,epsilon=0.8):
    out = thalamus(theta0=theta0_1, m0=np.zeros(len(h())), N=50, N_i=N_i1, c=c, epsilon=epsilon)
    thalamus(theta0=theta0_2,m0=out, N=50, N_i=N_i2, c=c, epsilon=epsilon)

def remove(c=1.2,epsilon=0.1,N_i1=30,N_i2=30):
    out = thalamusfull(theta0=0, m0=np.zeros(len(h())), N=50, N_i=N_i1, c=c, epsilon=epsilon, deltat=1.0, tau=5.0)
    thalamusfull(theta0=0,m0=out, N=50, N_i=N_i2, c=c, epsilon=epsilon, deltat=1.0, tau=5.0, hflag=False)
