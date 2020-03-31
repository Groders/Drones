import numpy as np
from scipy.stats import norm
from matplotlib import pyplot as plt

class UAVCluster:
    def __init__(self, numUAVs):
        self.uavs = generate_uavs(numUAVs)
        self.T = 0

    def update(self):
        # TODO

class UAV:
    def __init__(self, r, v, T=0):
        '''
        R -> position vector of R^2
        V -> velecity vector of R^2
        '''
        self.global_states = []
        self.local_states = [np.array([r, v]).reshape((4,1))]
        
        self.T = T
        # self.start = 

    def get_global_state(self, T=-1):
        return self.global_state[T]
        
    def get_local_state(self, T=-1):
        return self.global_state[T]

    def get_local_position(self, T=-1):
        return self.local_states[T][0:2, :]

    def get_local_velocity(self, T=-1):
        return self.local_states[T][2:4, :]

    def get_nearby_uavs(self, uavs, d, T):
        nearby_uavs = []
        for uav in uavs:
            if np.linalg.norm(uav.get_local_position(T) - self.get_local_position(T)) < d and (uav is not self):
                nearby_uavs.append(uav)
        return nearby_uavs

    def delta_v(self, a, c0, v0, dW, T, debug=False):
        '''
        We assume that dt = 1 unit of time
        Temporal state dynamics for velocity.
        a -> accel
        c0 -> some positive constant
        v0 -> wind velocity
        W -> standard wiener process i.i.d across UAVs
        '''
        # covariance of wind velocity
        V0 = np.eye(2) * v0
        delta_W = dW.reshape((2,1))
        if debug:
            print("a:{}\nc0:{}\nvel:{}\nv0:{}\nV0:{}\ndW:{}".format(a, c0, self.get_local_velocity(T), v0, V0, delta_W))
        dv = a - c0 * (self.get_local_velocity(T) - v0) + V0.dot(delta_W)
        return dv

    def delta_r(self, T):
        '''
        Since we assume that dt=1, then the  change in position is simply the amount of 
        distance we achieved in the current velocity
        '''
        return self.get_local_velocity(T)

    def average_cost(self, uavs, d):
        nearby_states = []
        x = [uav.s for uav in self.get_nearby_uavs(uavs, d)]
        print(x)
        print(x.shape)
        global_state = [self.s].concatenate(np.array([uav.s for uav in self.get_nearby_uavs(uavs, d)]))


    def __str__(self):
        return "{} {}".format(self.get_position(), self.get_velocity())

    def __repr__(self):
        return "p:{} v:{}".format(self.get_position().T, self.get_velocity().T)#{'pos': self.get_position(), 'vel': self.get_velocity()}
        




def wind_velocity():
    t_0 = 0 # define model parameters
    t_end = 2
    length = 1000
    theta = 1.1
    mu = 0.8
    sigma = 0.3
    t = np.linspace(t_0,t_end,length) # define time axis
    dt = np.mean(np.diff(t))
    y = np.zeros(length)
    y0 = np.random.normal(loc=0.0,scale=1.0) # initial condition
    drift = lambda y,t: theta*(mu-y) # define drift term, google to learn about lambda
    diffusion = lambda y,t: sigma # define diffusion term
    noise = np.random.normal(loc=0.0,scale=1.0,size=length)*np.sqrt(dt) #define noise process
    # solve SDE
    for i in range(1,length):
        y[i] = y[i-1] + drift(y[i-1],i*dt)*dt + diffusion(y[i-1],i*dt)*noise[i]

    plt.plot(t,y)
    plt.show()

def generate_uavs(n):
    uavs = []
    for i in range(n):
        r = [np.random.randint(0, 10), np.random.randint(0, 10)]
        v = [np.random.randint(0, 10), np.random.randint(0, 10)]
        uavs.append(UAV(r,v))
    return uavs

def brownian(T=1, N=100):
    '''
    Implementation of a standard wiener process 
    source: https://sites.me.ucsb.edu/~moehlis/APC591/tutorials/tutorial7/node2.html

    returns dW and W
    '''
    dt = T/N 
    dW = np.zeros((N,2))
    W = np.zeros((N,2))
    dW[0] = [np.sqrt(dt)*np.random.normal(loc=0.0, scale=1.0), np.sqrt(dt)*np.random.normal(loc=0.0, scale=1.0)]
    W[0] = dW[0]
    for i in range(1, N):
        dW[i] = [np.sqrt(dt)*np.random.normal(loc=0.0, scale=1.0), np.sqrt(dt)*np.random.normal(loc=0.0, scale=1.0)]
        W[i] = dW[i-1] + dW[i]
    return dW, W

def main():
    uavs = generate_uavs(5)
    print(uavs)

    dW, W = brownian(T=1, N=1000)
    wind_vel = np.array([100,100])[np.newaxis, :].T
    a = np.array([0.2,0.2]).reshape((2,1))
    c0 = 0.1
    change = uavs[0].delta_v(a, c0, wind_vel, dW[0], True)
    print(change)

    
    nearby = uavs[0].get_nearby_uavs(uavs, 100)
    print(len(nearby))

    uavs[0].average_cost(uavs, 100)


if __name__ == "__main__":
    main()