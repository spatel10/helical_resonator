import math
import matplotlib.pyplot as plt
import numpy as np

#constants
pi = math.pi
mhz = pow(10,6) #megahertz to hertz
pf = pow(10,-12) #picofarad to farad
mm = pow(10,-3) #millimeter to meter

#define the input parameters
#using typical values from the silverns paper
w0 = 2*pi*20*mhz                #w0: resonant freq
c_trap = 15*pf                  #c_trap: trap capacitance
c_wire = 15*pf                  #c_wire: capacitance due to connecting wires
r_trap = 1                      #r_trap: trap resistance (in ohms)
r_junct = 0.5                   #r_junct: helical coil to shield junction resistance
d0 = 5*mm                       #d0: coil wire diameter
tau = 10*mm                      #tau: coil winding pitch
c_sum = c_trap + c_wire

#---------------------------------------------------------
def coil_height(d,D):
    k_cd = 35 * d * pf  # function of d
    k_cb = 11.26 * pow(10, -12)  # constant (farad/meter)
    k_cs = 39.37 * 0.75 * pow(10, -12) / math.log(D / d, 10)  # function of d,D
    k_lc = 39.37 * 0.025 * pow(d, 2) * (1 - pow(d / D, 2)) * pow(10, -6) / pow(tau,2)  # function of d,D,tau

    b = ((c_sum + k_cd) / (k_cs + k_cb)) * (
                math.sqrt(((k_cs + k_cb) / (pow(c_sum + k_cd, 2) * k_lc * pow(w0, 2))) + (1 / 4)) - (1 / 2))

    return k_cd, k_cb, k_cs, k_lc, b;

def cap_shield(b,k_cs):
    cs = b * k_cs
    return cs;

def cap_coil(b,d):
    h = 11.26*(b/d) + 8 + 27/(math.sqrt(b/d))
    cc = h*d*pow(10,-12)
    return cc;

def res_shield_coil(b,d,D):
    res_copper = 0.0168*pow(10,-6)
    permiability = 0.999991*4*pi*pow(10,-7)
    skin_depth = math.sqrt((2*res_copper)/(permiability*w0))
    l_turn = math.sqrt(pow(tau,2)+pow(pi*d,2))
    num_turns = b/tau
    lc = l_turn*num_turns
    num_shield = (b*lc)/(4*pi*pow(D-d,2))
    ls = num_shield*math.sqrt(pow(pi*D,2)+pow(b/num_shield,2))
    rs = (num_shield*res_copper*ls)/(b*skin_depth)
    rc = (res_copper*lc)/(d0*pi*skin_depth)
    return rs, rc;

def res_esr(rs, rc, cs):
    r_esr = r_junct + rc + rs + r_trap*pow(c_trap/(c_trap+cs+c_wire),2)
    return r_esr;

def indc_coil(b, k_lc):
    lc = b*k_lc
    return lc;

def quality_factor(lc, r_esr):
    q = w0*lc/r_esr
    return q;

def check_res_freq(cs, cc, lc):
    check_w0 = 1/math.sqrt((cs+c_trap+c_wire+cc)*lc)
    return check_w0;

#--------------------------------------------------------------
def main():
    d_list = []
    d_D_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    q_list = []
    D_list = []
    constrain = []

    count_i = 0

    for d in np.arange(0.05,0.2,0.005):
        d_list.append(d)
        q_list.append([])
        constrain.append([])
        for d_D in d_D_list:
            D = d/d_D
            k_cd, k_cb, k_cs, k_lc, b = coil_height(d,D)
            cs = cap_shield(b,k_cs)
            if b >= 0:
                cc = cap_coil(b,d)
                rs, rc = res_shield_coil(b,d,D)
                r_esr = res_esr(rs,rc,cs)
                lc = indc_coil(b, k_lc)
                q = quality_factor(lc, r_esr)
                check_w0 = check_res_freq(cc, cs, lc)
                q_list[count_i].append(q)

                if b/d < 1:
                    constrain[count_i].append(1)
                else:
                    constrain[count_i].append(0)
        count_i += 1

    #check the dimensions of the lists are consistent
    if len(d_list) != len(q_list):
        print("number of elements in d_list do not match number of rows of Q")

    for i in range(len(q_list)):
        if len(q_list[i]) != len(d_D_list):
            print("Number of elements in ", i, "th row of Q does not match number of elements in d_D list")


    #plotting the contour using mathplotlib
    d_D_array = np.array(d_D_list)
    d_array = np.array(d_list)
    q_array = np.array(q_list)
    constrain_array = np.array(constrain)

    fig,ax = plt.subplots(1,1)
    cp = ax.contourf(d_D_array, d_array, constrain_array,  levels = [0,0.5,1], colors = ['white','black','black'])
    cp = ax.contourf(d_D_array, d_array, q_array, alpha=0.9)
    ax.set_ylabel('d [m]')
    ax.set_xlabel('d/D')
    ax.set_title('Q')
    fig.colorbar(cp)
    plt.show()

#-------------------------------------------------------
if __name__ == "__main__":
    main()
