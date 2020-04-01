import numpy as np
from scipy.stats import norm
from matplotlib import pyplot as plt



class UAVCluster:
    def __init__(self, numUAVs, d=100):
        self.uavs = generate_uavs(numUAVs)
        self.T = 0
        self.d = d

    def update(self):
        '''
        Update all uav states, call this after each iteration
        '''
        #TODO figure out how to update local state, maybe do it after we calc dv dr 
        for uav in self.uavs:
            uav.update(self.uavs, self.d)
    #self, a, v0, dW, debug=False, c0 = 0.1
    c0 = 0.1

    def sim (self, T= 2):
        dW, W = brownian(T=1, N=1000)
        wind_vel = np.array([100,100])[np.newaxis, :].T
        for t in range(T):
            for i, uav in enumerate (self.uavs):
                accel = np.array([np.random.uniform(0, 1),np.random.uniform(0, 1)]).reshape((2,1))
                delta_v = uav.delta_v(accel,wind_vel,dW[t])
                # print("delta_v", delta_v)
                uav.apply_dv(delta_v)
            
            self.update()
        for t in range(1):
            for i, uav in enumerate (self.uavs):
                uav.hjb_control(self.uavs, self.d, t)
                break
class UAV:
    #constants
    C_0 = 0.1
    C_1 = 100
    C_2 = 1.5
    C_3 = 1.5
    C_4 = 0.5
    BETA = 1
    E = 0.001

    A_mat = np.array([[0,0, 1, 0],
                  [0,0, 0,  1],
                  [0,0, -C_0, 0],
                  [0,0, 0,  -C_0]])
    B_mat = np.array([[0, 0],
                    [0, 0],
                  [1, 0],
                  [0, 1]])
    
    # Wind Velocity
    v0 = np.array([5,5])[np.newaxis, :].T

    V0 = np.eye(2) * v0
    
    G_mat = np.array([[0,0],
                        [0,0],
                  V0[0, 0], V0[0,1],
                  V0[1, 0], V0[1,1]])


    # a = []
    def __init__(self, r, v, T=0):
        '''
        R -> position vector of R^2
        V -> velecity vector of R^2
        '''
        self.global_states = []
        self.local_states = [np.array([r, v]).reshape((4,1))]
        self.a = []
        self.T = T
        # self.start = 

    def update(self,swarm,d):
        self.global_states.append(self.get_nearby_uavs_states(swarm,d,-1))
        self.T += 1

    def get_global_state(self, T=-1):
        return self.global_states[T]
        
    def get_local_state(self, T=-1):
        return self.local_states[T]

    def get_local_position(self, T=-1):
        return self.local_states[T][0:2, :]

    def get_local_velocity(self, T=-1):
        return self.local_states[T][2:4, :]

    def apply_dv(self, a, T=-1):
        new_state = np.zeros((4,1))
        new_state[2:4] = self.local_states[T][2:4] + a
        new_state[0:2] = self.local_states[T][0:2] + new_state[2:4]
        self.local_states.append(new_state)

    def get_nearby_uavs_states(self, uavs, d, T):
        nearby_uavs = []
        for uav in uavs:
            if np.linalg.norm(uav.get_local_position(T) - self.get_local_position(T)) < d and (uav is not self):
                nearby_uavs.append(uav.get_local_state())
        return nearby_uavs

    def delta_v(self, a, v0, dW, debug=False, c0 = 0.1):
        '''
        We assume that dt = 1 unit of time
        Temporal state dynamics for velocity.
        a -> accel
        c0 -> some positive constant
        v0 -> wind velocity
        W -> standard wiener process i.i.d across UAVs
        '''
        # covariance of wind velocity
        self.a.append(a)
        # V0 = np.eye(2) * v0
        delta_W = dW.reshape((2,1))
        if debug:
            print("a:{}\nc0:{}\nvel:{}\nv0:{}\nV0:{}\ndW:{}".format(a, c0, self.get_local_velocity(), UAV.v0, UAV.V0, delta_W))
        dv = a - c0 * (self.get_local_velocity() - UAV.v0) + UAV.V0.dot(delta_W)
        return dv

    def delta_r(self, T):
        '''
        Since we assume that dt=1, then the  change in position is simply the amount of 
        distance we achieved in the current velocity
        '''
        return self.get_local_velocity(T)

    def average_cost(self, T):
        total_cost = 0
        for t in range(T):
            total_cost += UAV.C_4*self.collision_cost(T) + self.kinetic_cost(T)
        
        return total_cost

    def collision_cost(self, t):
        global_state = self.global_states[t]
        prefix = 1/(len(global_state)+0.000001) # avoid 0 division

        local_state = self.local_states[t]

        v = local_state[2:4]
        r = local_state[0:2]
        s = 0
        for state in global_state:
            v_i = state[2:4]
            r_i = state[0:2]
            s += (np.linalg.norm(v - v_i)**2) / ((UAV.E + (np.linalg.norm(r - r_i)**2))**UAV.BETA)
        return s*prefix


    def kinetic_cost(self, t):
        #position velocity
        state = self.local_states[t]

        v = state[2:4]
        r = state[0:2]
        # print(np.linalg.norm(r))
        return (np.dot(v.T, r)/np.linalg.norm(r)) + UAV.C_1*(np.linalg.norm(r)**2) + UAV.C_2*(np.linalg.norm(v)**2) + UAV.C_3*(np.linalg.norm(self.a[t])**2)

    def hjb_control(self, uavs, d, T):
        '''
        uavs - > the total UAVS
        d -> distance of communication
        T -> time instance, an integer
        '''
        # print(type(uavs))
        # print(uavs)
        nearby_uavs = self.get_nearby_uavs_states(uavs, d, T)
        min_cost = self._hjb(T)

    def _hjb(self, T):
        d_avg_cost = self.C_4*self.collision_cost(T) + self.kinetic_cost(T)
        term2 = UAV.A_mat.dot(self.local_states[T])
        term3 = 1/(4*UAV.C_3) * UAV.B_mat.dot(UAV.B_mat.T)

        # print(UAV.A_mat)
        # print(self.local_states[T])
        print(term2)
        print(term3)

    def __str__(self):
        return "{} {}".format(self.get_local_position(), self.get_local_velocity())

    def __repr__(self):
        return "p:{} v:{}".format(self.get_local_position().T, self.get_local_velocity().T)#{'pos': self.get_position(), 'vel': self.get_velocity()}
        




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
        r = [np.random.uniform(0, 10), np.random.uniform(0, 10)]
        v = [np.random.uniform(0, 10), np.random.uniform(0, 10)]
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


    swarm = UAVCluster(2, d=5)
    swarm.sim()

    print(swarm.uavs[0].average_cost(1))

    '''
    uavs = generate_uavs(5)
    print(uavs)

    dW, W = brownian(T=1, N=1000)
    wind_vel = np.array([100,100])[np.newaxis, :].T
    a = np.array([0.2,0.2]).reshape((2,1))
    c0 = 0.1
    change = uavs[0].delta_v(a, wind_vel, dW[0], True)
    print(change)

    nearby = uavs[0].get_nearby_uavs_states(uavs, 100)
    print(len(nearby))

    uavs[0].average_cost(uavs, 100)
    '''

if __name__ == "__main__":
    main()