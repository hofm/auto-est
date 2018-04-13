#!/usr/bin/env python3

def model_setup(f_core, f_uncore):
    # Machine
    global n_cores, f_core_nom, p0_nom, alpha
    n_cores     = 8
    f_core_nom  = 2.7
    p0_nom      = 13
    alpha       = 0.35
    sust_mem_bw = { 1.2:28.2, 1.3:29.9, 1.4:31.3, 1.5:32.5, 1.6:33.5,
                    1.7:34.2, 1.8:34.8, 1.9:35.3, 2.0:35.8, 2.1:36.2,
                    2.2:36.5, 2.3:36.8, 2.4:37.1, 2.5:37.4, 2.7:37.8 }
    

    # Power model
    global pi_base_0, pi_base_1, pi_base_2, pi_core_0, pi_core_1, pi_core_2 
    pi_base_0 = 14.62
    pi_base_1 = 1.07
    pi_base_2 = 1.02
    pi_core_0 = 0.68
    pi_core_1 = 1.16
    pi_core_2 = 1.08

    # ECM model
    global T_L3Mem, T_ECM
    T_OL    = 52
    T_nOL   = 56
    T_L1L2  = (20*64)/32
    T_L2L3  = (12*64)/(32 * min(f_uncore/f_core, 1))
    T_L3Mem = (12*64)/(sust_mem_bw[f_uncore]/f_core)
    T_ECM   = max(T_OL, T_nOL+T_L1L2+T_L2L3+T_L3Mem)

def memory_bus_utilization(n_core, T_L3Mem, T_ECM, f_core):
    if n_core == 1:
        return T_L3Mem / T_ECM
    else:
        p0 = p0_nom * f_core/f_core_nom
        utilization = (n_core * T_L3Mem) / (T_ECM + (n_core-1) * p0 * memory_bus_utilization(n_core-1, T_L3Mem, T_ECM, f_core))
        return min(utilization, 1)

def ecm_est(n_core, f_core):
    ECM_perf_sat = 16 / T_L3Mem * (f_core * 1000**3)
    return memory_bus_utilization(n_core, T_L3Mem, T_ECM, f_core) * ECM_perf_sat

def pow_est(f_uncore, f_core, n_core, pareff):
    baseline_est = pi_base_0 + pi_base_1 * float(f_uncore) + pi_base_2 * float(f_uncore)**2
    core_est = pi_core_0 + pi_core_1 * float(f_core) + pi_core_2 * float(f_core)**2
    dampening = pareff**alpha
    return baseline_est + n_core * core_est*dampening

def generate_estimates(f_core, f_uncore):
    f = open('est-f_core-{}-f_uncore-{}.dat'.format(f_core, f_uncore), 'w+')
    model_setup(f_core, f_uncore)
    for n_core in range(1, n_cores+1):
        perf_est = ecm_est(n_core, f_core)
        pareff = (perf_est/ecm_est(1, f_core))/n_core
        energy_est = pow_est(f_uncore, f_core, n_core, pareff)/perf_est
        f.write('{} {}\n'.format(perf_est, energy_est))

freqs = (1.2, 1.8, 2.7)
for f in freqs:
    generate_estimates(f, f)
