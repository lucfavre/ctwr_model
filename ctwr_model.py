import numpy as np

# Bissection method to find the zero of a function
def dicho(f,a1,b1,eps):
    while (b1 - a1) > 2*eps :
        m = (a1 + b1)/2
        if f(a1)*f(m) < 0 :
            b1 = m
        else :
            a1 = m
    return (a1+b1)/2

#Function that returns the average Reynolds number of a falling droplet of a
#given diameter d_drop in an air flow of velocity u_h
#WARNING : air velocity must be negative since z axis is oriented downwards
def avg_rey(d_drop, u_h, H):

    #Physical constant
    g = 9.81 #Gravity
    rho_h = 1.2 #Air density
    rho_l = 1000 #Water density
    mu_h = 1.8e-5 #Air viscosity

    #Droplet mass, it is assumed that it remains constant through its fall
    md=(4/3) * np.pi * (d_drop/2)**3 * rho_l

    v=0 #Initial drop velocity
    dt=0.1 #Timestep
    z=H #Initial height

    #Droplet Reynolds number
    Re_d = rho_h * abs(v-u_h) * d_drop / mu_h

    Relist = [Re_d]

    while z>0:
        if Re_d > 0:
            CD = (24 / Re_d) * (1 + 0.15 * Re_d**0.687)
        else:
            CD = 0

        #Newton's 2nd law : velocity increase md * dv/dt = md * g + drag
        dv = dt * (g + 0.5 * rho_h * CD * (u_h - v) * abs(u_h - v) * np.pi * (d_drop/2)**2 / md)

        v += dv
        z = z - v*dt

        Re_d = rho_h * abs(v-u_h) * d_drop / mu_h
        Relist.append(Re_d)

    return np.average(Relist)

#1D solving of the air/water exchange equations for cooling towers
def solve_ctwr(h_cold, h_pack, h_hot, S_pack, P_0, T_h_i, T_w_i, mpt_d_i, mpt_w_i, x_rel_i, lamb_evap=0.5, n_evap=0.5, n_cells=3000):

    # Geometrical parameters
    H = h_pack # Packing height
    h2 = h_hot # Hot rain height
    h1 = h_cold # Cold rain height
    S = S_pack  # Packing cross section

    # Packing evaporation coefficients
    lamb = lamb_evap
    n = n_evap

    # Numerical parameters
    #dz = 5e-3 # Cell size
    #N2 = int((H+h1+h2)/dz) # Number of vertical cells
    N2 = n_cells
    dz = (H+h1+h2)/N2 # Cell size

    # Physical constants and properties
    cp_d = 1005 #Dry air heat capacity
    cp_w = 4185 #Liquid water heat capacity
    cp_v = 2010 #Water vapour heat capacity
    l_0 = 2501e3 #Evaporation latent heat
    l_100 = 2.3e6 #Idem
    g = 9.81  #Gravity constant
    rho_w = 998 #Liquid water density
    mu_h = 1.85e-5 #Humid air dynamic viscosity
    k_h = 0.0262   #Humid air thermal conductivity
    M_d = 28.9647e-3 #Dry air molar mass
    M_w = 18.01528e-3
    delta = M_w/M_d #Ratio of molar masses between water and dry air
    R = 8.31446261815324 # Perfect gases constant

    def P_sat(T): #Saturation pressure
        if T > 350 :
            T = 350
        T_c = T - 273.15
        return 6.112*np.exp(17.67*T_c/(T_c+243.5))*1e2

    def x_sat(T): #Humidity
        return delta*(P_sat(T)/(P_0 - P_sat(T)))

    def cp_h(x, T_h): #Humid air heat capacity
        if x < x_sat(T_h):
            resu = cp_d + x*cp_v
        else:
            resu = (cp_d + x_sat(T_h)*cp_v + (x - x_sat(T_h))*cp_w)
        return resu


    #Quantities of interest
    Fa_sur_Fe = mpt_d_i/mpt_w_i
    x_i = x_rel_i * x_sat(T_h_i) #Initial humidity
    rho_h = (1 + x_i) * delta / (delta + x_i) * P_0 / ((R/M_d) * T_h_i)  #Humid air density

    mpt_h_i = mpt_d_i*(1+x_i) #Humid air initial mass flow
    Eg = (mpt_w_i/rho_w)/((mpt_h_i/rho_h)+(mpt_w_i/rho_w)) #Estimated volume fraction of water in rain zones
    u_h = mpt_h_i / (rho_h * S) #Humid air inlet velocity

    def Lef(x, T_h, T_w): #Bosjnakovic formula
        Le = 0.866
        if x < x_sat(T_h):
            x_Lef = x
        else:
            x_Lef = x_sat(T_h)

        ksi = (x_sat(T_w)+delta)/(x_Lef+delta)
        if (ksi-1) < 1e-15:
            resu = Le**(2.0/3.0)
        else :
            resu = (Le**(2.0/3.0))*((ksi-1)/np.log(ksi))
        return resu


    #Droplet Reynolds number : cold rain zone
    Dg_cold = 5.0e-3
    Re_cold = avg_rey(Dg_cold, -1.0 * u_h, h1)

    Dg_hot = 3.0e-3
    Re_hot = avg_rey(Dg_hot, -1.0 * u_h, h2)

    def beta_ai(altiz,x,T_h,T_w):
        Pr = mu_h * cp_h(x,T_h)/k_h
        #Cold rain zone
        if altiz < h1:
            Nu = 2 + 0.6 * (Re_cold**(1/2)) * (Pr**(1/3))
            alpha = k_h * Nu/Dg_cold
            a_i = 6 * Eg/Dg_cold
            resu = a_i * alpha * Lef(x, T_h, T_w) / cp_h(x, T_h)

        #Hot rain zone
        elif altiz > h1+H:
            Nu = 2 + 0.6*(Re_hot**(1/2))*(Pr**(1/3))
            alpha = k_h*Nu/Dg_hot
            a_i = 6*Eg/Dg_hot
            resu = a_i*alpha*Lef(x, T_h, T_w)/cp_h(x, T_h)

        #Packing zone
        else:
            resu = lamb*(mpt_w_i*(1/S))*(mpt_d_i/mpt_w_i)**n
        return resu

    ### Resolution

    def dTh_dz(x,T_w,T_h,altiz):
        if x < x_sat(T_h):
            resu = (beta_ai(altiz,x,T_h,T_w)*S/(mpt_d_i*cp_h(x, T_h)))*(Lef(x,T_h,T_w)*cp_h(x, T_h)*(T_w-T_h) + (x_sat(T_w)-x)*(cp_v*T_w-cp_v*T_h))
        else :
            resu = (beta_ai(altiz,x,T_h,T_w)*S/(mpt_d_i*cp_h(x, T_h)))*(Lef(x,T_h,T_w)*cp_h(x, T_h)*(T_w-T_h) + (x_sat(T_w)-x_sat(T_h))*(cp_v*T_w-cp_v*T_h))
        return resu

    def dx_dz(x,T_w,T_h,altiz):
        if x < x_sat(T_h):
            resu = (beta_ai(altiz,x,T_h,T_w)*S/mpt_d_i)*(x_sat(T_w)-x)
        else :
            resu = (beta_ai(altiz,x,T_h,T_w)*S/mpt_d_i)*(x_sat(T_w)-x_sat(T_h))
        return resu

    def dTw_dz(x, T_w, T_h,altiz):
        if x < x_sat(T_h):
            resu = beta_ai(altiz, x, T_h, T_w)*S*(Lef(x,T_h,T_w)*cp_h(x,T_h)*(T_w-T_h)+(x_sat(T_w)-x)*(cp_v*T_w-cp_w*T_w+l_0))*(1/(cp_w*(mpt_w_i-(x-x_i)*mpt_d_i)))
        else:
            resu = beta_ai(altiz, x, T_h, T_w)*S*(Lef(x,T_h,T_w)*cp_h(x,T_h)*(T_w-T_h)+(x_sat(T_w)-x_sat(T_h))*(cp_v*T_w-cp_w*T_w+l_0))*(1/(cp_w*(mpt_w_i-(x-x_i)*mpt_d_i)))
        return resu

    def euler_asc_final(T_w_z0):
        T_wz = T_w_z0
        T_hz = T_h_i
        xz = x_i
        altiz = 0
        T_w = [T_w_z0]
        T_h = [T_h_i]
        x = [x_i]
        x_sat_list = [x_sat(T_h_i)]
        alti = [0]
        for i in range(1,N2):
            T_wz = T_wz + dz*dTw_dz(xz,T_wz,T_hz,altiz)
            T_hz = T_hz + dz*dTh_dz(xz,T_wz,T_hz,altiz)
            xz = xz + dz*dx_dz(xz,T_wz,T_hz,altiz)
            altiz = altiz + dz
            alti.append(altiz)
            T_w.append(T_wz)
            T_h.append(T_hz)
            x.append(xz)
            x_sat_list.append(x_sat(T_hz))

        res = {'T_w':np.array(T_w), 'T_h':np.array(T_h), 'Tw_in':T_w[-1],
               'z':np.array(alti), 'x':np.array(x), 'x_sat':np.array(x_sat_list)}
        return res


    def comp_asc_final(T_w_z0):
        Tw_in = euler_asc_final(T_w_z0)['Tw_in']
        return (Tw_in - T_w_i)

    ### Dichotomie
    a = 1 + 273.15
    b = 50 + 273.15

    Tw_out = dicho(comp_asc_final, a, b, 0.1)
    res_final = euler_asc_final(Tw_out)

    return res_final
