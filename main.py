import numpy as np
from scipy.stats import norm
from matplotlib import pyplot as plt
import sys
# from sympy import Symbol, lambdify, Matrix, Sum, MatrixSymbol, Derivative
import sympy
from sympy import Symbol, lambdify, Matrix, Sum, MatrixSymbol, Derivative
# from sympy import *
from apm import *

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

    def sim (self, T=1):
        dW, W = brownian(T=1, N=1000)
        wind_vel = np.array([0,0])[np.newaxis, :].T
        a = []
        for t in range(T):
            for i, uav in enumerate (self.uavs):
                accel = np.array([np.random.uniform(0, 1),np.random.uniform(0, 1)]).reshape((2,1))
                delta_v = uav.delta_v(accel,wind_vel,dW[t])
                # print("delta_v", delta_v)
                uav.apply_dv(delta_v)
                a.append(accel)
            
            self.update()
        hjb_cost = [0] * len(self.uavs)
        for t in range(T):
            for i, uav in enumerate (self.uavs):
                hjb_cost[i] = uav.hjb_control(t)
        print("hjb cost:")
        print(hjb_cost)
        print("acceleration: {}".format(a))

    def hjb_control(self, uav):
        accels = generate_mesh_grid_points((-10,10), 11)
        print("accelerations: ",accels)
        hjb_values = []
        print("length: {}".format(len(accels)))
        for i in range(len(accels)):
            accel = np.array(accels[i]).reshape((2,1))
            # get the last state this uav was in
            local_state = uav.local_states[-1]
            global_state = uav.global_states[-1]

            updated_local_state = uav.get_new_local_state(local_state, accel)
            hjb_output = uav._hjb_control(updated_local_state, global_state, accel)
            # hjb_output = uav._hjb_control(local_state, global_state, accel)
            hjb_values.append(hjb_output)

        hjb_values = np.array(hjb_values)
        min_value = hjb_values.min()
        index = hjb_values.argmin()
        print(hjb_values)
        print(hjb_values.min())
        print(hjb_values.argmin())
        print(uav.local_states[-1])
        print("accel: {}".format(accels[index]))
        opt_accel = np.array(accels[index]).reshape((2,1))

        # delta_v = uav.delta_v(opt_accel,wind_vel,brown)
        # print("deltaV = {}".format(delta_v))
        updated_local_state = uav.get_new_local_state(uav.local_states[-1], opt_accel)
        # print("updated_local_state: {}".format(updated_local_state))
        # out = uav.optimal_accel(updated_local_state, uav.global_states[-1], opt_accel)
        out = uav.optimal_accel(uav.local_states[-1], uav.global_states[-1], opt_accel)
        print(out)
        
    
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
    v0 = np.array([2,2])[np.newaxis, :].T
    




    # a = []
    def __init__(self, r, v, brownian_noise, T=0):
        '''
        R -> position vector of R^2
        V -> velecity vector of R^2
        '''
        self.global_states = []
        self.local_states = [np.array([r, v]).reshape((4,1))]
        self.a = []
        self.T = T
        # self.acceleration_points = generate_mesh_grid_points((-10,10), 111)
        self.brownian_noise = brownian_noise #dW and w

        self.V0 = calculate_wind_covariance(UAV.v0, brownian_noise[1])
        self.G_mat = np.array([[0,0],
                        [0,0],
                  [self.V0[0, 0], self.V0[0,1]],
                  [self.V0[1, 0], self.V0[1,1]]])

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

###############################################################################
    ## Updated stuff here 
    def single_div (self, local_state, global_state, accel):
        #return row vector
        diff = 0.00000001
        diffs = np.array(np.zeros(4))
        # initial = self.average_cost_interval (0, T, state = state_1)
        initial = self.local_state_cost(local_state, accel) + self.global_state_cost(local_state, global_state)
        for i in range (4):
            state = np.copy(local_state)
            state[i,0] += diff
            diffs[i] = (self.local_state_cost(state, accel) + self.global_state_cost(state, global_state) - initial)/diff
        print(diffs)
        return diffs


    def optimal_accel(self, local_state, global_state, accel):
        # deriv = self.partial_psi_wrt_state(local_state, global_state, accel)
        deriv = self.single_div(local_state, global_state, accel)
        print("derivative: {}".format(deriv))
        print("1/2c3: {}".format(1.0/(2*UAV.C_3)))
        print("b_mat.t: {}".format(UAV.B_mat.T))
        a = 1.0/(2*UAV.C_3) * UAV.B_mat.T.dot(deriv)
        return a

    def get_new_local_state(self, old_state, accel):
        new_state = np.zeros((4,1))
        new_state[2:4] = old_state[2:4] + accel
        new_state[0:2] = old_state[0:2] + old_state[2:4]
        return new_state

    def d_psi_wrt_all(self, local_state, global_state, accel):
        '''
        Calculates the derivative of the cost function with respect to all variables (in this case it's the addition of the 2 cost functions)
        '''
        total = self.global_state_cost(local_state, global_state) + self.local_state_cost(local_state, accel)
        return total

    def partial_psi_wrt_state(self, local_state, global_state, accel, delta=1):
        total = self.partial_global_state_cost(local_state, global_state, delta) + self.partial_local_state_cost(local_state, accel, delta)
        return total

    def local_state_cost(self, local_state, accel):
        #position velocity
        state = local_state
        r = state[0:2]
        v = state[2:4]
        px, py, vx, vy, ax, ay = r[0,0], r[1,0], v[0,0], v[1,0], accel[0,0], accel[1,0]
        return self.kinetic_cost_lambda(px,py, vx,  vy, ax, ay)     

    def global_state_cost(self, local_state, global_state):
        n = 1/(len(global_state)+0.000001) # avoid 0 division
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
        return s*n
    

    def partial_local_state_cost(self, local_state, accel, delta=1):
        px, py, vx, vy, ax, ay = local_state[0,0], local_state[1,0], local_state[2,0], local_state[3,0], accel[0], accel[1]
        if delta == 1:
            return self.delta1_kinetic_cost_lambda(px,py,vx,vy,ax,ay)
        else:
            delta_2 = self.delta2_kinetic_cost_lambda(px,py, vx,vy,ax,ay)
            delta_2 = np.array(delta_2).reshape((4,4))
            return delta_2


    def partial_global_state_cost(self, local_state, global_state, delta=1):
        n = 1/(len(global_state)+0.000001) # avoid 0 division
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
                s = np.add(s, cost)
            return s*n

    def _hjb_control(self, local_state, global_state, accel):
        d_cost = self.d_psi_wrt_all(local_state, global_state, accel)

        delta1_cost = self.partial_psi_wrt_state(local_state, global_state, accel, delta=1)
        delta2_cost = self.partial_psi_wrt_state(local_state, global_state, accel, delta=2)
        global_cost = self.global_state_cost(local_state, global_state)
        local_cost  = self.local_state_cost(local_state, accel)

        term1 = d_cost
        term2 = (UAV.A_mat.dot(local_state) + 1/(4.0 * UAV.C_3) * (UAV.B_mat.dot(UAV.B_mat.T)).dot(delta1_cost) + UAV.B_mat.dot(UAV.v0) * UAV.C_0).T.dot(delta1_cost)
        term3 = 0.5 * (self.G_mat.dot(self.G_mat.T).dot(delta2_cost)).trace()
        hjb_total = term1 + term2 + term3 + local_cost + global_cost
        hjb_total = hjb_total.flatten()[0]

        return hjb_total

    def psi_cost(self, t):
        total_cost = 0
        T = self.T
        # list_of_costs = []
        for t in range(T):
            px, py, vx, vy, ax, ay = self.local_states[t][0,0], self.local_states[t][1,0], self.local_states[t][2,0], self.local_states[t][3,0], self.a[t][0], self.a[t][1]
            total_cost += UAV.C_4*self.collision_avoidance_cost(T) + self.kinetic_cost_lambda(px,py,vx,vy, ax,ay)
        return total_cost

    ############################################################################################
    def total_cost(self, t):
        total_cost = 0
        T = self.T
        # list_of_costs = []
        for t in range(T):
            px, py, vx, vy, ax, ay = self.local_states[t][0,0], self.local_states[t][1,0], self.local_states[t][2,0], self.local_states[t][3,0], self.a[t][0], self.a[t][1]
            total_cost += UAV.C_4*self.collision_avoidance_cost(T) + self.kinetic_cost_lambda(px,py,vx,vy, ax,ay)
        return total_cost

    def delta_total_cost(self, t, delta=1):
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

    def delta_kinetic_cost2(self, state, accel, delta=1):
        px, py, vx, vy, ax, ay = state[0,0], state[1,0], state[2,0], state[3,0], accel[0][0], accel[1][0]
        if delta == 1:
            return self.delta1_kinetic_cost_lambda(px,py, vx,  vy, ax, ay)
        else:
            delta_2 = self.delta2_kinetic_cost_lambda(px,py, vx,  vy, ax, ay)
            delta_2 = np.array(delta_2).reshape((4,4))
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
            print("a:{}\nc0:{}\nvel:{}\nv0:{}\nV0:{}\ndW:{}".format(a, c0, self.get_local_velocity(), UAV.v0, self.V0, delta_W))
        dv = a - c0 * (self.get_local_velocity() - UAV.v0) + self.V0.dot(delta_W)
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

    # def hjb_control(self, T):
    #     '''
    #     uavs - > the total UAVS
    #     d -> distance of communication
    #     T -> time instance, an integer
    #     '''
    #     # print(type(uavs))
    #     # print(uavs)
    #     nearby_uavs = self.get_nearby_uavs_states(uavs, d, T)
    #     min_cost = self._hjb(T)
    #     print("hjb_cost ", min_cost)
    #     return 
    def hjb_control(self, T):
        #the derivative of the cost function is simply the sum of the two energy functions (since derivative of integral cancels out)
    

        delta1_cost = self.delta_total_cost(T, delta=1)
        delta2_cost = self.delta_total_cost(T, delta=2)
        kin_cost = self.kinetic_cost(T) 
        col_cost = self.collision_avoidance_cost(T) 
        term1 = kin_cost + col_cost
        term2 = (UAV.A_mat.dot(self.local_states[T]) + 1/(4.0 * UAV.C_3) * (UAV.B_mat.dot(UAV.B_mat.T)).dot(delta1_cost) + UAV.B_mat.dot(UAV.v0) * UAV.C_0).T.dot(delta1_cost)
        term3 = 0.5 * (self.G_mat.dot(self.G_mat.T).dot(delta2_cost)).trace()
        hjb_total = term1 + term2 + term3 + kin_cost + col_cost
        hjb_total = hjb_total.flatten()[0]
        # print("flattened",hjb_total)
        return hjb_total


    def __str__(self):
        return "{} {}".format(self.get_local_position(), self.get_local_velocity())

    def __repr__(self):
        return "p:{} v:{}".format(self.get_local_position().T, self.get_local_velocity().T)#{'pos': self.get_position(), 'vel': self.get_velocity()}
        

def apply_acceleration(accel, old_state):
    new_state = np.zeros((4,1))
    new_state[2:4] = old_state[2:4] + accel
    new_state[0:2] = old_state[0:2] + new_state[2:4]
    return new_state


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

        r = [np.random.uniform(-10, 10), np.random.uniform(-10, 10)]
        v = [np.random.uniform(-5, 5), np.random.uniform(-5, 5)]
        uavs.append(UAV(r,v,brownian()))
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


def generate_mesh_grid_points(interval=(-5,5), samples=11):
    nx = np.linspace(interval[0], interval[1], samples)
    x, y = np.meshgrid(nx, nx, indexing='xy')
    point = []
    for i in range(len(nx)):
        for j in range(len(nx)):
            point.append([x[i,j], y[i,j]])
    points = np.array(point)
    return points

def calculate_wind_covariance(expected_wind_velocity, brownian_noise):
    wind_samples = np.zeros((100,2))
    for i, noise in enumerate(brownian_noise):
        real_wind_vel = expected_wind_velocity + np.array(noise).reshape(2,1)
        wind_samples[i] = real_wind_vel.reshape(2,)   
    cov = np.cov(wind_samples.T)
    return cov
    
def main():
    swarm = UAVCluster(1, d=20)
    swarm.sim()
    uav1 = swarm.uavs[0]
    
    swarm.hjb_control(uav1)
    # print(uav1._hjb_control(uav1.local_states[0], uav1.global_states[0], uav1.a[0]))
    # print(uav1.hjb_control(0))
    # print(uav1.a)
    

    dw, w = brownian()
    # print(brown[0][0])
    # sys.exit(0)
    # wind_vel = np.array([2,2])[np.newaxis, :].T
    # swarm.hjb_control(uav1, wind_vel, brown[0][0])
    # points = generate_mesh_grid_points((-5,5), 11)



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