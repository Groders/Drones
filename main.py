import numpy as np
from scipy.stats import norm
from matplotlib import pyplot as plt
import sys
# from sympy import Symbol, lambdify, Matrix, Sum, MatrixSymbol, Derivative
import sympy
from sympy import Symbol, lambdify, Matrix, Sum, MatrixSymbol, Derivative
# from sympy import *
"""
"""

class UAVCluster:
    def __init__(self, numUAVs, d=5):
        self.uavs = generate_uavs(numUAVs)
        self.T = 0
        self.d = d
        # intialize each UAV's nearest neighors for state 0
        # for uav in self.uavs:
        #     uav.global_states.append(uav.get_nearby_uavs_states(self.uavs,self.d,-1))

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
                  [V0[0, 0], V0[0,1]],
                  [V0[1, 0], V0[1,1]]])


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


        """define our symbols and functions"""
        px, py, vx, vy = Symbol('px') ,Symbol('py'), Symbol('vx'), Symbol('vy') 
        opx, opy, ovx, ovy = Symbol('opx') ,Symbol('opy'), Symbol('ovx'), Symbol('ovy') 
        oposition = Matrix([opx, opy])
        ovelocity = Matrix([ovx, ovy])
        ostate = Matrix([oposition, ovelocity])
        ax, ay = Symbol('ax'), Symbol('ay')
        n = Symbol('n')

        position = Matrix([px, py])
        velocity = Matrix([vx, vy])
        state = Matrix([position, velocity])
        acceleration = Matrix([ax, ay])

        kinetic_cost_expression = (velocity.dot(position)/((px**2 + py**2)**0.5) + UAV.C_1 * (px**2 + py**2) + UAV.C_2*(vx**2 + vy**2) + UAV.C_3*(ax**2 + ay**2))
        delta1_kinetic_cost = kinetic_cost_expression.diff(state)
        delta2_kinetic_cost = delta1_kinetic_cost.diff(state)
        # ||v-vi||^2 -> 
        collision_avoidance_cost = ((vx-ovx)**2 + (ovy-vy)**2)/((UAV.E + ((px-opx)**2 + (py - opy)**2))**UAV.BETA)
        delta_collision_avoidance_cost = collision_avoidance_cost.diff(state)
        second_delta_collision_avoidance_cost = delta_collision_avoidance_cost.diff(state)
        # collision_avoidance_cost = 1/n * (ovelocity-velocity).norm()**2/(UAV.E + (oposition-position).norm()**2)**UAV.BETA
        # collision_avoidance_cost = (ovelocity-velocity).T*(ovelocity-velocity)/((UAV.E + (oposition-position).T*(oposition-position))**UAV.BETA)
        self.collision_cost_lambda = lambdify([px,py,vx,vy,opx,opy,ovx,ovy], collision_avoidance_cost, "numpy")
        self.delta1_collision_cost_lambda = lambdify([px,py,vx,vy,opx,opy,ovx,ovy], delta_collision_avoidance_cost, "numpy") # this is the first derivative with respect to state
        self.delta2_collision_cost_lambda = lambdify([px,py,vx,vy,opx,opy,ovx,ovy], second_delta_collision_avoidance_cost, "numpy") # this is the first derivative with respect to state
        # print("asdfasdfasdf")
        print(delta_collision_avoidance_cost)
        # print("col", collision_avoidance_cost)
        # print(delta_kinetic_cost)
        self.kinetic_cost_lambda = lambdify([px, py, vx, vy, ax, ay], kinetic_cost_expression, "numpy" )
        self.delta1_kinetic_cost_lambda = lambdify([px, py, vx, vy, ax, ay], delta1_kinetic_cost, "numpy" )
        self.delta2_kinetic_cost_lambda = lambdify([px, py, vx, vy, ax, ay], delta2_kinetic_cost, "numpy" )



    def total_cost(self, t):
        total_cost = 0
        for t in range(T):
            px, py, vx, vy, ax, ay = self.local_states[T][0,0], self.local_states[T][1,0], self.local_states[T][2,0], self.local_states[T][3,0], self.a[T][0], self.a[T][1]
            total_cost += UAV.C_4*self.collision_avoidance_cost(T) + self.kinetic_cost_lambda(px,py,vx,vy, ax,ay)
        return total_cost

    def delta_total_cost(self, t, delta=1):
        # print("kinetic cost: ", self.delta_kinetic_cost(t, delta))
        # print("collision cost:", self.delta_collision_cost(t, delta))
        # print(self.delta_kinetic_cost(t, delta) + self.delta_collision_cost(t, delta))
        total_cost = self.delta_kinetic_cost(t, delta) + self.delta_collision_cost(t, delta)
        return total_cost

        # print(self.delta1_collision_cost_lambda())
    def delta_kinetic_cost(self, T, delta=1):
        px, py, vx, vy, ax, ay = self.local_states[T][0,0], self.local_states[T][1,0], self.local_states[T][2,0], self.local_states[T][3,0], self.a[T][0], self.a[T][1]
        if delta == 1:
            return self.delta1_kinetic_cost_lambda(px,py, vx,  vy, ax, ay)
        else:
            delta_2 = self.delta2_kinetic_cost_lambda(px,py, vx,  vy, ax, ay)
            delta_2 = np.array(delta_2).reshape((4,4))
            # print("kinetic_:", delta_2)
            # print(delta_2.shape)
            # sys.exit(0)
            return delta_2
       

    def delta_collision_cost(self, t, delta=1):
        global_state = self.global_states[t]
        n = 1/(len(global_state)+0.000001) # avoid 0 division
        local_state = self.local_states[t]
        r = local_state[0:2]
        v = local_state[2:4]
        if delta==1:
            s = np.zeros((4,1))
            for state in global_state:
                v_i = state[2:4]
                r_i = state[0:2]
                cost = self.delta1_collision_cost_lambda(r[0], r[1], v[0],v[1], r_i[0], r_i[1], v_i[0], v_i[1])
                cost = cost.reshape((4,1))
                s = np.add(s, cost)
            return s*n   
        else:
            s = np.zeros((4,4))
            for state in global_state:
                v_i = state[2:4]
                r_i = state[0:2]
                cost = self.delta2_collision_cost_lambda(r[0], r[1], v[0],v[1], r_i[0], r_i[1], v_i[0], v_i[1])
                cost = np.array(cost).flatten().tolist()
                new_cost = np.zeros((16))
                for i, x in enumerate(cost):
                    if type(x) == list:
                        new_cost[i] = x[0]
                    else:
                        new_cost[i] = x
                cost = new_cost.reshape((4,4))
                # print("new cost", new_cost)
                # print("flattened_cost", cost)
                # print(np.concatenate(cost).ravel())
                # print("delta2", cost.shape)
                # sys.exit(0)
                s = np.add(s, cost)

            return s*n
           


    def collision_cost(self, t):
        global_state = self.global_states[t]
        prefix = 1/(len(global_state)+0.000001) # avoid 0 division

        local_state = self.local_states[t]

        r = local_state[0:2]
        v = local_state[2:4]
        s = 0
        for state in global_state:
            v_i = state[2:4]
            r_i = state[0:2]
            s += (np.linalg.norm(v - v_i)**2) / ((UAV.E + (np.linalg.norm(r - r_i)**2))**UAV.BETA)
        print("orig cost:", s)
        return s*prefix        
    def collision_avoidance_cost(self, t):
        global_state = self.global_states[t]
        n = 1/(len(global_state)+0.000001) # avoid 0 division
        local_state = self.local_states[t]
        r = local_state[0:2]
        v = local_state[2:4]
        s = 0

        for state in global_state:
            r_i = state[0:2]
            v_i = state[2:4]
            # print("ostate {}, lstate {}".format(state, local_state))
            cost = self.collision_cost_lambda(r[0], r[1], v[0],v[1], r_i[0], r_i[1], v_i[0], v_i[1])
            # print("lambda cost: ", cost)
            s += cost
        print("col avoid", s)
        return s*n
    def update(self,swarm,d):
        print(swarm)
        print(d)

        print(self.get_nearby_uavs_states(swarm,d,-1))
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
        # print("inside get nearby", uavs)
        for uav in uavs:
            # print("is {} smaller than {}: {}".format(np.linalg.norm(uav.get_local_position(T) - self.get_local_position(T)), d, np.linalg.norm(uav.get_local_position(T) - self.get_local_position(T)) < d) )
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
            print("orig",self.kinetic_cost(T))
            px, py, vx, vy, ax, ay = self.local_states[T][0,0], self.local_states[T][1,0], self.local_states[T][2,0], self.local_states[T][3,0], self.a[T][0], self.a[T][1]
            print("lambda", self.kinetic_cost_lambda(px,py, vx,  vy, ax, ay))

            print("col_orig", self.collision_cost(T))
            print("lambda_col", self.collision_avoidance_cost(T))
            total_cost += UAV.C_4*self.collision_cost(T) + self.kinetic_cost(T)
        
        return total_cost


    def kinetic_cost(self, t):
        #position velocity
        state = self.local_states[t]
        r = state[0:2]
        v = state[2:4]
        px, py, vx, vy, ax, ay = r[0,0], r[1,0], v[0,0], v[1,0], self.a[t][0], self.a[t][1]
        return self.kinetic_cost_lambda(px,py, vx,  vy, ax, ay)#(np.dot(v.T, r)/np.linalg.norm(r)) + UAV.C_1*(np.linalg.norm(r)**2) + UAV.C_2*(np.linalg.norm(v)**2) + UAV.C_3*(np.linalg.norm(self.a[t])**2)

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
        print("hjb_cost ", min_cost)

    def _hjb(self, T):
        #the derivative of the cost function is simply the sum of the two energy functions (since derivative of integral cancels out)
        delta1_cost = self.delta_total_cost(T, delta=1)
        delta2_cost = self.delta_total_cost(T, delta=2)
        kin_cost = self.kinetic_cost(T) 
        col_cost = self.collision_avoidance_cost(T) 
        term1 = kin_cost + col_cost
        term2 = (UAV.A_mat.dot(self.local_states[T]) + 1/(4.0 * UAV.C_3) * (UAV.B_mat.dot(UAV.B_mat.T)).dot(delta1_cost) + UAV.B_mat.dot(UAV.v0) * UAV.C_0).T.dot(delta1_cost)
        term3 = 0.5 * (UAV.G_mat.dot(UAV.G_mat.T).dot(delta2_cost)).trace()
        hjb_total = term1 + term2 + term3 + kin_cost + col_cost
        return hjb_total

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


    # a = np.array([1,1])
    # b = np.array([2,2])

    # print(np.linalg.norm((a-b))**2 )
    # print(a-b)
    # print(abs((a[0]-b[0])**2 + (a[1]-b[1]**2)))
    # sys.exit(0)

    swarm = UAVCluster(3, d=20)
    swarm.sim()
    print(swarm.uavs[0].average_cost(1))
    # delta_total_cost = swarm.uavs[0].delta_total_cost(1)
    # print("delta_total_cost", delta_total_cost)


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