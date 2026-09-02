import numpy as np
import pandas as pd
from scipy.optimize import root_scalar
from scipy.interpolate import interp1d

# --- PVT Correlations ---

def standing_rs(p, t, api, sg_gas):
    """Standing correlation for solution gas-oil ratio (Rs) in scf/STB."""
    x = 0.0125 * api - 0.00091 * t
    rs = sg_gas * (((p / 18.2) + 1.4) * 10**x)**1.2048
    return np.maximum(rs, 0)

def standing_bo(rs, t, api, sg_gas):
    """Standing correlation for oil formation volume factor (Bo) in rb/STB."""
    sg_o = 141.5 / (131.5 + api)
    f = rs * np.sqrt(sg_gas / sg_o) + 1.25 * t
    bo = 0.9759 + 0.00012 * f**1.2
    return bo

def vasquez_beggs_rs(p, t, api, sg_gas, sep_p=100, sep_t=100):
    """Vasquez-Beggs correlation for solution gas-oil ratio (Rs) in scf/STB."""
    sg_gas_adj = sg_gas * (1.0 + 5.912e-5 * api * sep_t * np.log10(sep_p / 114.7))
    if api <= 30:
        c1, c2, c3 = 0.0362, 1.0937, 25.724
    else:
        c1, c2, c3 = 0.0178, 1.1870, 23.931
    rs = c1 * sg_gas_adj * p**c2 * np.exp(c3 * api / (t + 460))
    return np.maximum(rs, 0)

def vasquez_beggs_bo(rs, t, api, sg_gas, sep_p=100, sep_t=100):
    """Vasquez-Beggs correlation for oil formation volume factor (Bo) in rb/STB."""
    sg_gas_adj = sg_gas * (1.0 + 5.912e-5 * api * sep_t * np.log10(sep_p / 114.7))
    if api <= 30:
        c1, c2, c3 = 4.677e-4, 1.751e-5, -1.811e-8
    else:
        c1, c2, c3 = 4.67e-4, 1.1e-5, 1.337e-9
    bo = 1.0 + c1 * rs + c2 * (t - 60) * (api / sg_gas_adj) + c3 * rs * (t - 60) * (api / sg_gas_adj)
    return bo

def calc_bg(p, t, sg_gas):
    """Gas formation volume factor (Bg) in rb/scf (Simplified using Ideal Gas Law)."""
    z = 0.85 
    t_R = t + 460
    bg = 0.02827 * z * t_R / p # cu ft/scf
    bg_rb = bg / 5.615 # rb/scf
    return bg_rb

def calc_visc_dead_oil(t, api):
    """Beggs-Robinson correlation for dead oil viscosity (cp)."""
    x = t**-1.163 * np.exp(6.9824 - 0.04658 * api)
    mu_od = 10**x - 1
    return np.maximum(mu_od, 0.1)

def calc_visc_sat_oil(mu_od, rs):
    """Chew-Connally correlation for saturated oil viscosity (cp)."""
    a = 10.715 * (rs + 100)**-0.515
    b = 5.44 * (rs + 150)**-0.338
    mu_os = a * mu_od**b
    return mu_os

def calc_visc_gas(p, t, sg_gas):
    """Simplified gas viscosity (cp)."""
    return 0.015

def calc_pvt(p, t, api, sg_gas, water_cut, pvt_method='Vasquez-Beggs'):
    """Calculates all necessary PVT properties at given P and T."""
    # Ensure pressure is strictly positive to prevent division by zero or log(0) errors
    p = max(p, 14.7)
    
    if pvt_method == 'Vasquez-Beggs':
        rs = vasquez_beggs_rs(p, t, api, sg_gas)
        bo = vasquez_beggs_bo(rs, t, api, sg_gas)
    else:
        rs = standing_rs(p, t, api, sg_gas)
        bo = standing_bo(rs, t, api, sg_gas)
        
    bg = calc_bg(p, t, sg_gas)
    
    mu_od = calc_visc_dead_oil(t, api)
    mu_o = calc_visc_sat_oil(mu_od, rs)
    mu_g = calc_visc_gas(p, t, sg_gas)
    
    bw = 1.0
    mu_w = 0.5
    
    return rs, bo, bg, mu_o, mu_g, bw, mu_w

# --- IPR Models ---

def vogel_ipr(pr, pb, q_test, pwf_test):
    """Calculates Vogel's two-phase IPR curve."""
    pwf_arr = np.linspace(pr, 14.7, 50)
    q_arr = np.zeros_like(pwf_arr)
    
    if pwf_test >= pb:
        j = q_test / (pr - pwf_test)
        qb = j * (pr - pb)
        qmax = qb + j * pb / 1.8
    else:
        if pr <= pb:
            qmax = q_test / (1 - 0.2*(pwf_test/pr) - 0.8*(pwf_test/pr)**2)
            qb = 0
            j = qmax / (pr / 1.8)
        else:
            j = q_test / ((pr - pb) + (pb / 1.8) * (1 - 0.2*(pwf_test/pb) - 0.8*(pwf_test/pb)**2))
            qb = j * (pr - pb)
            qmax = qb + j * pb / 1.8
            
    for i, pwf in enumerate(pwf_arr):
        if pwf >= pb:
            q_arr[i] = j * (pr - pwf)
        else:
            q_arr[i] = qb + (qmax - qb) * (1 - 0.2*(pwf/pb) - 0.8*(pwf/pb)**2)
            
    return q_arr, pwf_arr, qmax

def fetkovich_ipr(pr, q_test, pwf_test, n):
    """Calculates Fetkovich IPR curve."""
    pwf_arr = np.linspace(pr, 14.7, 50)
    c = q_test / (pr**2 - pwf_test**2)**n
    qmax = c * (pr**2)**n
    q_arr = c * (pr**2 - pwf_arr**2)**n
    q_arr[np.isnan(q_arr)] = 0
    return q_arr, pwf_arr, qmax

# --- VLP Model (Beggs & Brill) ---

def beggs_brill_gradient(p, t, vsl, vsg, rho_l, rho_g, mu_l, mu_g, d_in, roughness, angle=90):
    """Calculates the pressure gradient dp/dz using Beggs & Brill correlation."""
    d = d_in / 12.0 # ft
    area = np.pi * d**2 / 4.0
    vm = vsl + vsg
    lambda_l = vsl / vm if vm > 0 else 0
    
    if vm == 0:
        return rho_l * np.sin(np.radians(angle)) / 144.0
    
    n_fr = vm**2 / (32.174 * d)
    # Prevent extremely small values from causing errors
    n_fr = max(n_fr, 1e-6)
    lambda_l = max(lambda_l, 1e-6)
    
    l1 = 316 * lambda_l**0.302
    l2 = 0.0009252 * lambda_l**-2.4684
    l3 = 0.1 * lambda_l**-1.4516
    l4 = 0.5 * lambda_l**-6.738
    
    regime = 'transition'
    if (lambda_l < 0.01 and n_fr < l1) or (lambda_l >= 0.01 and n_fr < l2):
        regime = 'segregated'
    elif (0.01 <= lambda_l < 0.4 and l3 < n_fr <= l1) or (lambda_l >= 0.4 and l3 < n_fr <= l4):
        regime = 'intermittent'
    elif (lambda_l < 0.4 and n_fr >= l1) or (lambda_l >= 0.4 and n_fr > l4):
        regime = 'distributed'
    
    if regime == 'segregated':
        a, b, c = 0.98, 0.4846, 0.0868
    elif regime == 'intermittent':
        a, b, c = 0.845, 0.5351, 0.0173
    elif regime == 'distributed':
        a, b, c = 1.065, 0.5824, 0.0609
    else: 
        a, b, c = 0.98, 0.4846, 0.0868
        
    hl_0 = a * lambda_l**b / n_fr**c
    if hl_0 < lambda_l:
        hl_0 = lambda_l
        
    theta = np.radians(angle)
    # Beggs & Brill correction for inclination
    if regime == 'segregated':
        e, f, g, h = 0.011, -3.768, 3.539, -1.614
    elif regime == 'intermittent':
        e, f, g, h = 2.96, 0.305, -0.4473, 0.0978
    elif regime == 'distributed':
        e, f, g, h = 0.0, 0.0, 0.0, 0.0 # No correction for distributed
    else:
        e, f, g, h = 0.011, -3.768, 3.539, -1.614
        
    n_lv = 1.938 * vsl * (rho_l / 0.00204)**0.25 # using surface tension approx
    n_lv = max(n_lv, 1e-6)
    c_inc = (1 - lambda_l) * np.log(max(1 + e * np.sin(1.8*theta) * n_lv**f * n_fr**g * lambda_l**h, 0.1))
    if np.isnan(c_inc): c_inc = 0
    
    hl = hl_0 * (1 + c_inc)
    hl = max(min(hl, 1.0), 0.0) 
    
    rho_m = rho_l * hl + rho_g * (1 - hl)
    mu_m = mu_l * hl + mu_g * (1 - hl)
    
    n_re = 1488 * rho_m * vm * d / max(mu_m, 1e-6)
    if n_re > 0:
        f_n = (-1.8 * np.log10((roughness/d)/3.7 + 6.9/n_re))**-2
    else:
        f_n = 0
        
    y = lambda_l / (hl**2) if hl > 0 else 0
    y = max(y, 1e-6)
    if 1 < y < 1.2:
        s = np.log(y) / (-0.0523 + 3.182 * np.log(y) - 0.8725 * np.log(y)**2 + 0.01853 * np.log(y)**4)
    else:
        s = np.log(max(2.2*y - 1.2, 1e-6))
        
    f_tp = f_n * np.exp(s) if n_re > 0 else 0
    
    dp_dz_elev = rho_m * np.sin(theta) / 144.0
    dp_dz_fric = 2 * f_tp * rho_m * vm**2 / (32.174 * d * 144.0)
    
    dp_dz = dp_dz_elev + dp_dz_fric
    return dp_dz

def calculate_vlp_curve(p_wh, t_wh, t_bh, depth, q_liqs, water_cut, gor, api, sg_gas, d_in, roughness, pvt_method='Vasquez-Beggs'):
    """Calculates VLP curve for a range of liquid rates."""
    pwf_arr = []
    dz = 500 # step size in ft
    steps = int(depth / dz)
    
    sg_o = 141.5 / (131.5 + api)
    rho_w_std = 62.4
    rho_o_std = 62.4 * sg_o
    rho_g_std = 0.0764 * sg_gas
    
    for q_liq in q_liqs:
        if q_liq == 0:
            pwf_arr.append(p_wh + depth * (rho_o_std*(1-water_cut) + rho_w_std*water_cut) / 144.0)
            continue
            
        q_o = q_liq * (1 - water_cut)
        q_w = q_liq * water_cut
        q_g = q_o * gor
        
        p = p_wh
        t_grad = (t_bh - t_wh) / depth
        
        for i in range(steps):
            z = i * dz
            t = t_wh + t_grad * z
            
            rs, bo, bg, mu_o, mu_g, bw, mu_w = calc_pvt(p, t, api, sg_gas, water_cut, pvt_method)
            
            q_o_in_situ = q_o * bo * 5.615 / 86400
            q_w_in_situ = q_w * bw * 5.615 / 86400
            
            free_gas = q_g - q_o * rs
            if free_gas < 0: free_gas = 0
            q_g_in_situ = free_gas * bg * 5.615 / 86400
            
            area = np.pi * (d_in/24.0)**2
            vsl = (q_o_in_situ + q_w_in_situ) / area
            vsg = q_g_in_situ / area
            
            rho_o = (rho_o_std + rs * rho_g_std / 5.615) / bo
            rho_w = rho_w_std / bw
            rho_l = rho_o * (1 - water_cut) + rho_w * water_cut
            
            rho_g = rho_g_std * (p / 14.7) * (520 / (t + 460)) / 0.85
            
            mu_l = mu_o * (1 - water_cut) + mu_w * water_cut
            
            dp_dz = beggs_brill_gradient(p, t, vsl, vsg, rho_l, rho_g, mu_l, mu_g, d_in, roughness)
            p += dp_dz * dz
            
        pwf_arr.append(p)
        
    return np.array(pwf_arr)

def find_intersection(q_ipr, p_ipr, q_vlp, p_vlp):
    """Finds intersection point between IPR and VLP."""
    sort_ipr = np.argsort(q_ipr)
    q_ipr_sorted = q_ipr[sort_ipr]
    p_ipr_sorted = p_ipr[sort_ipr]
    
    sort_vlp = np.argsort(q_vlp)
    q_vlp_sorted = q_vlp[sort_vlp]
    p_vlp_sorted = p_vlp[sort_vlp]
    
    q_ipr_u, ind = np.unique(q_ipr_sorted, return_index=True)
    p_ipr_u = p_ipr_sorted[ind]
    
    q_vlp_u, ind = np.unique(q_vlp_sorted, return_index=True)
    p_vlp_u = p_vlp_sorted[ind]
    
    if len(q_ipr_u) < 2 or len(q_vlp_u) < 2:
        return None, None
        
    f_ipr = interp1d(q_ipr_u, p_ipr_u, bounds_error=False, fill_value="extrapolate")
    f_vlp = interp1d(q_vlp_u, p_vlp_u, bounds_error=False, fill_value="extrapolate")
    
    def diff(q):
        return f_ipr(q) - f_vlp(q)
        
    try:
        diff_arr = f_ipr(q_vlp_u) - p_vlp_u
        sign_changes = np.where(np.diff(np.sign(diff_arr)))[0]
        if len(sign_changes) == 0:
            # Check if there is any point where VLP is below IPR, meaning it can flow.
            # If not, it means there is no intersection and the well is dead.
            return None, None
            
        idx = sign_changes[-1]
        q_guess1 = q_vlp_u[idx]
        q_guess2 = q_vlp_u[idx+1]
        
        sol = root_scalar(diff, bracket=[q_guess1, q_guess2], method='brentq')
        if sol.converged:
            q_int = sol.root
            p_int = f_ipr(q_int)
            if q_int >= 0 and q_int <= max(q_ipr_u):
                return q_int, p_int
    except Exception as e:
        pass
        
    return None, None
