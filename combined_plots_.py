import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from CoolProp.CoolProp import PropsSI

# ==========================================
# 1. Plot Styling (Publication Standard)
# ==========================================
sns.set_style("whitegrid")

# Update rcParams: force Times New Roman and set all base sizes to 16
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
    'font.size': 16,
    'font.weight': 'normal',         # Normal font weight overall
    'axes.labelweight': 'bold',      # Bold axis labels
    'axes.titleweight': 'bold',      # Bold title
    'axes.labelsize': 16,          
    'axes.titlesize': 16,
    'xtick.labelsize': 16,         
    'ytick.labelsize': 16,         
    'legend.fontsize': 16,
    'legend.edgecolor': 'black',
    'legend.fancybox': False
})

refrigerant = 'R1234yf'
Q_target = 3517
eta_n = 0.85
eta_isen = 0.85
T0 = 298.15 

# ==============================================================================
# PART 1: Evaporator Temperature Variation
# ==============================================================================
T_evap_array = np.array([275.15, 277.15, 279.15, 281.15, 283.15])
T_cond_fixed = 308.15
P_cond = PropsSI('P', 'T', T_cond_fixed, 'Q', 0, refrigerant)
n1 = len(T_evap_array)

EDR_sys_vcrs_1 = np.zeros(n1)
EDR_sys_ej_1 = np.zeros(n1)
Ex_eff_vcrs_1 = np.zeros(n1)
Ex_eff_ej_1 = np.zeros(n1)

for i in range(n1):
    T_evap = T_evap_array[i]
    P_evap = PropsSI('P', 'T', T_evap, 'Q', 1, refrigerant)
    Ep = Q_target * (T0 / T_evap - 1)
    
    # Simple VCRS
    h1_vcrs = PropsSI('H', 'P', P_evap, 'Q', 1, refrigerant)
    s1_vcrs = PropsSI('S', 'P', P_evap, 'Q', 1, refrigerant)
    h2s_vcrs = PropsSI('H', 'P', P_cond, 'S', s1_vcrs, refrigerant)
    h2_vcrs = h1_vcrs + (h2s_vcrs - h1_vcrs) / eta_isen
    h3 = PropsSI('H', 'P', P_cond, 'Q', 0, refrigerant)
    h4_vcrs = h3
    
    m_comp_vcrs = Q_target / (h1_vcrs - h4_vcrs)
    W_comp_vcrs = m_comp_vcrs * (h2_vcrs - h1_vcrs)
    
    s2_vcrs = PropsSI('S', 'P', P_cond, 'H', h2_vcrs, refrigerant)
    s3_vcrs = PropsSI('S', 'P', P_cond, 'Q', 0, refrigerant)
    s4_vcrs = PropsSI('S', 'P', P_evap, 'H', h4_vcrs, refrigerant)
    
    ED_comp_vcrs = T0 * m_comp_vcrs * (s2_vcrs - s1_vcrs)
    ED_cond_vcrs = m_comp_vcrs * ((h2_vcrs - h3) - T0 * (s2_vcrs - s3_vcrs))
    ED_exp_vcrs  = T0 * m_comp_vcrs * (s4_vcrs - s3_vcrs)
    ED_evap_vcrs = T0 * (m_comp_vcrs * (s1_vcrs - s4_vcrs) - Q_target / T_evap)
    ED_total_vcrs = ED_comp_vcrs + ED_cond_vcrs + ED_exp_vcrs + ED_evap_vcrs
    
    EDR_sys_vcrs_1[i] = ED_total_vcrs / Ep
    Ex_eff_vcrs_1[i] = (Ep / W_comp_vcrs) * 100
    
    # Ejector VCRS
    h4s = PropsSI('H', 'P', P_evap, 'S', s3_vcrs, refrigerant)
    h4 = h3 - eta_n * (h3 - h4s)
    u4 = np.sqrt(2 * (h3 - h4))
    h10 = PropsSI('H', 'T', T_evap, 'Q', 1, refrigerant)
    s10 = PropsSI('S', 'T', T_evap, 'Q', 1, refrigerant)
    
    mu_low, mu_high = 0.01, 1.0
    for _ in range(100):
        mu = (mu_low + mu_high) / 2
        u6 = u4 / (1 + mu)
        h6 = ((h4 + 0.5 * u4**2) + mu * h10) / (1 + mu) - 0.5 * u6**2
        s6 = PropsSI('S', 'P', P_evap, 'H', h6, refrigerant)
        h7 = (h3 + mu * h10) / (1 + mu)
        h7s = h6 + eta_n * (h7 - h6)
        P_c = PropsSI('P', 'S', s6, 'H', h7s, refrigerant)
        try: x7 = max(0, min(1, PropsSI('Q', 'P', P_c, 'H', h7, refrigerant)))
        except: x7 = 1.0 if h7 >= PropsSI('H', 'P', P_c, 'Q', 1, refrigerant) else 0.0
        residual = (1 + mu) * x7 - 1
        if residual > 0: mu_high = mu
        else: mu_low = mu
        if abs(residual) < 1e-4: break
            
    h8 = PropsSI('H', 'P', P_c, 'Q', 0, refrigerant)
    m_evap = Q_target / (h10 - h8)
    m_comp_ej = m_evap / mu
    m_total_ej = m_comp_ej + m_evap
    
    h1_ej = PropsSI('H', 'P', P_c, 'Q', 1, refrigerant)
    s1_ej = PropsSI('S', 'P', P_c, 'Q', 1, refrigerant)
    h2s_ej = PropsSI('H', 'P', P_cond, 'S', s1_ej, refrigerant)
    h2_ej = h1_ej + (h2s_ej - h1_ej) / eta_isen
    W_comp_ej = m_comp_ej * (h2_ej - h1_ej)
    
    s2_ej = PropsSI('S', 'P', P_cond, 'H', h2_ej, refrigerant)
    s7 = PropsSI('S', 'P', P_c, 'H', h7, refrigerant)
    s8 = PropsSI('S', 'P', P_c, 'Q', 0, refrigerant)
    s9 = PropsSI('S', 'P', P_evap, 'H', h8, refrigerant)
    
    ED_comp_ej = T0 * m_comp_ej * (s2_ej - s1_ej)
    ED_cond_ej = m_comp_ej * ((h2_ej - h3) - T0 * (s2_ej - s3_vcrs))
    ED_ej      = T0 * (m_total_ej * s7 - m_comp_ej * s3_vcrs - m_evap * s10)
    ED_sep_ej  = T0 * (m_comp_ej * s1_ej + m_evap * s8 - m_total_ej * s7)
    ED_exp_ej  = T0 * m_evap * (s9 - s8)
    ED_evap_ej = T0 * (m_evap * (s10 - s9) - Q_target / T_evap)
    ED_total_ej= ED_comp_ej + ED_cond_ej + ED_ej + ED_sep_ej + ED_exp_ej + ED_evap_ej
    
    EDR_sys_ej_1[i] = ED_total_ej / Ep
    Ex_eff_ej_1[i] = (Ep / W_comp_ej) * 100

# ==============================================================================
# PART 2: Gas Cooler Exit Temperature Variation (Transcritical)
# ==============================================================================
T_cond_array = np.array([308.15, 310.65, 313.15, 315.65, 318.15])
T_evap_fixed = 275.15
P_gc = 110e5 
n2 = len(T_cond_array)

EDR_sys_vcrs_2 = np.zeros(n2)
EDR_sys_ej_2 = np.zeros(n2)
Ex_eff_vcrs_2 = np.zeros(n2)
Ex_eff_ej_2 = np.zeros(n2)

P_evap_2 = PropsSI('P', 'T', T_evap_fixed, 'Q', 1, refrigerant)
h10_2 = PropsSI('H', 'T', T_evap_fixed, 'Q', 1, refrigerant)
s10_2 = PropsSI('S', 'T', T_evap_fixed, 'Q', 1, refrigerant)
h1_vcrs_2 = PropsSI('H', 'P', P_evap_2, 'Q', 1, refrigerant)
s1_vcrs_2 = PropsSI('S', 'P', P_evap_2, 'Q', 1, refrigerant)
h2s_vcrs_2 = PropsSI('H', 'P', P_gc, 'S', s1_vcrs_2, refrigerant)
h2_vcrs_2 = h1_vcrs_2 + (h2s_vcrs_2 - h1_vcrs_2) / eta_isen
s2_vcrs_2 = PropsSI('S', 'P', P_gc, 'H', h2_vcrs_2, refrigerant)
Ep_2 = Q_target * (T0 / T_evap_fixed - 1)

for i in range(n2):
    T_gc_exit = T_cond_array[i]
    h3 = PropsSI('H', 'P', P_gc, 'T', T_gc_exit, refrigerant)
    s3 = PropsSI('S', 'P', P_gc, 'T', T_gc_exit, refrigerant)
    
    # Simple VCRS
    m_vcrs = Q_target / (h1_vcrs_2 - h3)
    W_comp_vcrs = m_vcrs * (h2_vcrs_2 - h1_vcrs_2)
    s4_vcrs = PropsSI('S', 'P', P_evap_2, 'H', h3, refrigerant)
    
    ED_comp_vcrs = T0 * m_vcrs * (s2_vcrs_2 - s1_vcrs_2)
    ED_gc_vcrs   = m_vcrs * ((h2_vcrs_2 - h3) - T0 * (s2_vcrs_2 - s3))
    ED_exp_vcrs  = T0 * m_vcrs * (s4_vcrs - s3)
    ED_evap_vcrs = T0 * (m_vcrs * (s1_vcrs_2 - s4_vcrs) - Q_target / T_evap_fixed)
    ED_total_vcrs = ED_comp_vcrs + ED_gc_vcrs + ED_exp_vcrs + ED_evap_vcrs
    
    EDR_sys_vcrs_2[i] = ED_total_vcrs / Ep_2
    Ex_eff_vcrs_2[i] = (Ep_2 / W_comp_vcrs) * 100
    
    # Ejector Cycle
    h4s = PropsSI('H', 'P', P_evap_2, 'S', s3, refrigerant)
    h4 = h3 - eta_n * (h3 - h4s)
    u4 = np.sqrt(2 * (h3 - h4))
    
    mu_low, mu_high = 0.01, 1.0
    for _ in range(100):
        mu = (mu_low + mu_high) / 2
        u6 = u4 / (1 + mu)
        h6 = ((h4 + 0.5 * u4**2) + mu * h10_2) / (1 + mu) - 0.5 * u6**2
        s6 = PropsSI('S', 'P', P_evap_2, 'H', h6, refrigerant)
        h7 = (h3 + mu * h10_2) / (1 + mu)
        h7s = h6 + eta_n * (h7 - h6)
        P_c = PropsSI('P', 'S', s6, 'H', h7s, refrigerant)
        try: x7 = max(0, min(1, PropsSI('Q', 'P', P_c, 'H', h7, refrigerant)))
        except: x7 = 1.0 if h7 >= PropsSI('H', 'P', P_c, 'Q', 1, refrigerant) else 0.0
        residual = (1 + mu) * x7 - 1
        if residual > 0: mu_high = mu
        else: mu_low = mu
        if abs(residual) < 1e-4: break
            
    h8 = PropsSI('H', 'P', P_c, 'Q', 0, refrigerant)
    m_evap = Q_target / (h10_2 - h8)
    m_comp_ej = m_evap / mu
    m_total_ej = m_comp_ej + m_evap
    
    h1_ej = PropsSI('H', 'P', P_c, 'Q', 1, refrigerant)
    s1_ej = PropsSI('S', 'P', P_c, 'Q', 1, refrigerant)
    h2s_ej = PropsSI('H', 'P', P_gc, 'S', s1_ej, refrigerant)
    h2_ej = h1_ej + (h2s_ej - h1_ej) / eta_isen
    W_comp_ej = m_comp_ej * (h2_ej - h1_ej)
    
    s2_ej = PropsSI('S', 'P', P_gc, 'H', h2_ej, refrigerant)
    s7 = PropsSI('S', 'P', P_c, 'H', h7, refrigerant)
    s8 = PropsSI('S', 'P', P_c, 'Q', 0, refrigerant)
    s9 = PropsSI('S', 'P', P_evap_2, 'H', h8, refrigerant) 

    ED_comp_ej = T0 * m_comp_ej * (s2_ej - s1_ej)
    ED_gc_ej   = m_comp_ej * ((h2_ej - h3) - T0 * (s2_ej - s3))
    ED_ej      = T0 * (m_total_ej * s7 - m_comp_ej * s3 - m_evap * s10_2)
    ED_sep_ej  = T0 * (m_comp_ej * s1_ej + m_evap * s8 - m_total_ej * s7)
    ED_exp_ej  = T0 * m_evap * (s9 - s8)
    ED_evap_ej = T0 * (m_evap * (s10_2 - s9) - Q_target / T_evap_fixed)
    ED_total_ej= ED_comp_ej + ED_gc_ej + ED_ej + ED_sep_ej + ED_exp_ej + ED_evap_ej
    
    EDR_sys_ej_2[i] = ED_total_ej / Ep_2
    Ex_eff_ej_2[i] = (Ep_2 / W_comp_ej) * 100


# ==============================================================================
# PLOTTING ROUTINES (Times New Roman, 16pt Formatting)
# ==============================================================================

# Note: Removed the "title" parameter here to match your function calls
def plot_comparison(x, y1, y2, x_label, y_label):
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Dropped linewidth to 2 to keep lines proportionate
    ax.plot(x, y1, marker='o', markersize=8, linewidth=2, color='#1f77b4', label='Simple VCRS', zorder=3)
    ax.plot(x, y2, marker='s', markersize=8, linewidth=2, color='#ff7f0e', label='Ejector VCRS', zorder=3)
    
    # Explicitly enforce Times New Roman, 16pt, bold labels
    ax.set_xlabel(x_label, fontname='Times New Roman', fontweight='bold', fontsize=16, labelpad=10)
    ax.set_ylabel(y_label, fontname='Times New Roman', fontweight='bold', fontsize=16, labelpad=10)
    
    # Enforce Times New Roman, 16pt, normal weight for ticks
    for label in ax.get_xticklabels():
        label.set_fontname('Times New Roman')
        label.set_fontsize(16)
        label.set_fontweight('normal')
        
    for label in ax.get_yticklabels():
        label.set_fontname('Times New Roman')
        label.set_fontsize(16)
        label.set_fontweight('normal')
        
    # Enforce Times New Roman, 16pt, normal weight for legend
    ax.legend(frameon=True, edgecolor='black', fancybox=False, loc='best', 
              prop={'family': 'Times New Roman', 'size': 16, 'weight': 'normal'})
              
    plt.tight_layout()
    plt.show()

# 1. System EDR vs Evaporator Temperature
plot_comparison(T_evap_array, EDR_sys_vcrs_1, EDR_sys_ej_1, 
                'Evaporator Temperature (K)', 'System Exergy Destruction Ratio')

# 2. System EDR vs Gas Cooler Exit Temperature
plot_comparison(T_cond_array, EDR_sys_vcrs_2, EDR_sys_ej_2, 
                'Gas Cooler Exit Temperature (K)', 'System Exergy Destruction Ratio')

# 3. Exergetic Efficiency vs Evaporator Temperature
plot_comparison(T_evap_array, Ex_eff_vcrs_1, Ex_eff_ej_1, 
                'Evaporator Temperature (K)', 'Exergetic Efficiency (%)')

# 4. Exergetic Efficiency vs Gas Cooler Exit Temperature
plot_comparison(T_cond_array, Ex_eff_vcrs_2, Ex_eff_ej_2, 
                'Gas Cooler Exit Temperature (K)', 'Exergetic Efficiency (%)')