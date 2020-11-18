import numpy as np

def calculate_emittances(x,px,y,py,sigma,delta,beta_gamma=1):
    delta_square = delta*delta
    mean_delta_square = np.mean(delta_square)
    def emittance(u,up):
        mean_u_delta = np.mean(u*delta)
        mean_up_delta = np.mean(up*delta)
        uu = np.mean(u*u) - mean_u_delta**2/mean_delta_square
        uup = np.mean(u*up) - mean_u_delta*mean_up_delta/mean_delta_square
        upup = np.mean(up*up) - mean_up_delta**2/mean_delta_square
        return np.sqrt(uu*upup - uup**2)*beta_gamma
    return emittance(x,px), emittance(y,py)