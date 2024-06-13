from sympy import sqrt

def schwarzschild_event_horizon(M, **kwargs):
    return [2 * M]

def kerr_event_horizon(M, a, **kwargs):
    r_plus = M + sqrt(M**2 - a**2)
    r_minus = M - sqrt(M**2 - a**2)
    return [r_plus, r_minus]

def reissner_nordstrom_event_horizon(M, Q, **kwargs):
    r_plus = M + sqrt(M**2 - Q**2)
    r_minus = M - sqrt(M**2 - Q**2)
    return [r_plus, r_minus]

def kerr_newman_event_horizon(M, a, Q, **kwargs):
    r_plus = M + sqrt(M**2 - a**2 - Q**2)
    r_minus = M - sqrt(M**2 - a**2 - Q**2)
    return [r_plus, r_minus]

def flrw_event_horizon(a, k, **kwargs):
    # FLRW metric's event horizon depends on the scale factor and curvature
    return []

def de_sitter_event_horizon(Lambda, **kwargs):
    return [sqrt(3 / Lambda)]

def anti_de_sitter_event_horizon(Lambda, **kwargs):
    # Anti-de Sitter has no event horizon in this context
    return []

def minkowski_event_horizon(**kwargs):
    # Minkowski space has no event horizon
    return []

def vaidya_event_horizon(M, **kwargs):
    return [2 * M]

def gottingen_event_horizon(**kwargs):
    # Göttingen metric's event horizon needs a specific definition, typically none
    return []

def bertotti_robinson_event_horizon(r0, **kwargs):
    # Bertotti-Robinson metric is typically horizon-less in common contexts
    return []

def detect_event_horizon(metric_name, M=0, a=0, Q=0, r0=0, Lambda=0, k=0):
    if metric_name == "schwarzschild":
        return schwarzschild_event_horizon(M)
    elif metric_name == "kerr":
        return kerr_event_horizon(M, a)
    elif metric_name == "reissner_nordstrom":
        return reissner_nordstrom_event_horizon(M, Q)
    elif metric_name == "kerr_newman":
        return kerr_newman_event_horizon(M, a, Q)
    elif metric_name == "flrw":
        return flrw_event_horizon(a, k)
    elif metric_name == "de_sitter":
        return de_sitter_event_horizon(Lambda)
    elif metric_name == "anti_de_sitter":
        return anti_de_sitter_event_horizon(Lambda)
    elif metric_name == "minkowski":
        return minkowski_event_horizon()
    elif metric_name == "vaidya":
        return vaidya_event_horizon(M)
    elif metric_name == "gottingen":
        return gottingen_event_horizon()
    elif metric_name == "bertotti_robinson":
        return bertotti_robinson_event_horizon(r0)
    else:
        raise ValueError("Unknown metric")
