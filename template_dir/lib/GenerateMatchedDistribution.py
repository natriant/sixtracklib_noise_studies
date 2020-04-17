import math
import sys
from itertools import chain
import numpy as np
import csv
import random

from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy.constants import speed_of_light

class DoubleRF(object):
   
    def __init__(self,R,eta,beta,E,phi,h,V):
        cos = np.cos
        sin = np.sin
        pi = np.pi
        sqrt = np.sqrt
        c = speed_of_light

        phi = np.atleast_1d(phi)
        h = np.atleast_1d(h)
        V = np.atleast_1d(V)
        self.beta = beta
        self.eta = eta
        self.T = lambda dp: 0.5*beta*c*eta*dp**2
        
        def v(V,h,phi):
            return lambda z: V/h*(cos(h*z/R+phi) - cos(phi))     
        self.V = lambda z: 1./(E*beta/c)/2/pi * sum([v(*i)(z) for i in zip(V,h,phi)])
    
    def get_H(self,z,dp):
        return self.T(dp) + self.V(z)
    
    def get_dp(self,z,H):
        c = speed_of_light
        return np.abs(np.sqrt((H-self.V(z))/(0.5*self.beta*c*self.eta)))


def _Gauss(x,x0,a,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))


def _GaussianFit(x, y):
    mean = sum(x*y)/sum(y)
    sigma = np.sqrt(sum(y*(x-mean)**2)/sum(y))
    amplitude = max(y)
    popt,pcov = curve_fit(_Gauss,x,y,p0=[mean,amplitude,sigma])
    amplitude_norm = popt[1]*np.sqrt(2*np.pi)/(x[1]-x[0]) * popt[2] / np.float(sum(y))
    return popt, amplitude_norm


def _Gaussian_sigma_from_FWHM(x,y):
    from scipy.interpolate import UnivariateSpline
    spline = UnivariateSpline(x, y-np.max(y)/2, s=0)
    r1, r2 = spline.roots() 
    return (r2-r1)/2.3548



class LongitudinalBinomialDistribution():

    def __init__(self, RF, z_max, m):
        self.RF = RF
        self.z_max = z_max
        self.H_max = RF.get_H(z_max,0.0)
        self.m = m
        self.dist = lambda z,m: (1-np.clip(z,0,1)**2)**(m-1) # REPRESENTATION OF BEAM ELLIPSES FOR TRANSPORT CALCULATIONS, W. JOHO, SIN-REPORT TM-11-14

    def getCoordinates(self, n_mp=1):
        dist = self.dist
        RF = self.RF
        z_max = self.z_max
        H_max = self.H_max
        m = self.m

        z = np.linspace(-z_max,z_max,100)
        dp_max = 1.2*np.nanmax(RF.get_dp(z[1:-1], H_max))
        U_ = []
        V_ = []
        W_ = []		
        while len(U_)<n_mp:
            u = np.random.uniform(-z_max,z_max,n_mp)
            v = np.random.uniform(-dp_max,dp_max,n_mp)
            w = np.random.uniform(0,1,n_mp)
            d = dist(np.sqrt(RF.get_H(u,v)/H_max), m)
            mask = np.where(w < d)[0]
            U_.extend(u[mask])
            V_.extend(v[mask])
            W_.extend(w[mask])
            # print len(U_)
        z_rand = np.array(U_[:n_mp])
        dp_rand = np.array(V_[:n_mp])
        return z_rand, dp_rand

    def getBunchProfile(self, n_steps=100):
        dist = self.dist
        RF = self.RF
        z_max = self.z_max
        H_max = self.H_max
        m = self.m

        z = np.linspace(-z_max,z_max,n_steps)
        dp_max = 1.2*np.nanmax(RF.get_dp(z[1:-1], H_max))
        dp = np.linspace(-dp_max,dp_max,n_steps)
        xx, yy = np.meshgrid(z, dp, sparse=False)
        hh = dist(RF.get_H(xx,yy)/H_max, m)
        hh_ysum = np.sum(hh,axis=0)
        z_step = np.mean(np.diff(z))
        z_profile = hh_ysum/np.sum(hh_ysum*z_step)
        z_mean = sum(z*z_profile)/sum(z_profile)    
        z_rms = np.sqrt( sum(z_profile * (z-z_mean)**2)/sum(z_profile) )
        hh_xsum = np.sum(hh,axis=1)
        dp_step = np.mean(np.diff(dp))
        dp_profile = hh_xsum/np.sum(hh_xsum*dp_step)
        dp_mean = sum(dp*dp_profile)/sum(dp_profile)
        dp_rms = np.sqrt( sum(dp_profile*(dp-dp_mean)**2)/sum(dp_profile) )
        return z, z_profile, z_rms, dp, dp_profile, dp_rms



def generate_longitudinal_distribution(partCO, gamma_transition,
                                        circumference, harmonic_number, 
                                        rf_voltage, rf_phase, z_max, 
                                        JohoParameter, n_macroparticles):

    eta = 1./gamma_transition**2 - 1./partCO.gamma0**2
    Radius = circumference/2/np.pi
    h_main = np.atleast_1d(harmonic_number)[0]
    RF = DoubleRF(Radius, eta, partCO.beta, partCO.energy0, rf_phase, harmonic_number, rf_voltage)
    Longitudinal_distribution = LongitudinalBinomialDistribution(RF, z_max, JohoParameter)
    sigma, delta = Longitudinal_distribution.getCoordinates(n_macroparticles)
    # sigma_arr, sigma_profile, sigma_rms, \
    # 	delta_arr, delta_profile, delta_rms = Longitudinal_distribution.getBunchProfile()
    # output_dictionary = {}
    # output_dictionary['dpp_sigma'] = _GaussianFit(dp, dp_profile)[0][2]
    # output_dictionary['dpp_sigma_from_FWHM'] = _Gaussian_sigma_from_FWHM(dp, dp_profile)
    # output_dictionary['dpp_profile'] = np.array([dp, dp_profile])
    # output_dictionary['dpp_rms'] = dpp_rms
    # output_dictionary['linedensity_profile'] = np.array([z_arr, z_profile])
    return sigma, delta


def generate_transverse_distribution(partCO,twiss_at_start,neps_x,neps_y,
                                        delta,n_macroparticles):
    assert(len(delta) == n_macroparticles)
    alfa_x = twiss_at_start['alfx']
    alfa_y = twiss_at_start['alfy']
    beta_x = twiss_at_start['betx']
    beta_y = twiss_at_start['bety']
    dx = twiss_at_start['dx']
    dy = twiss_at_start['dy']
    dpx = twiss_at_start['dpx']
    dpy = twiss_at_start['dpy']

    sigma_x  = np.sqrt(neps_x * beta_x / partCO.beta0 / partCO.gamma0)
    sigma_y  = np.sqrt(neps_y * beta_y / partCO.beta0 / partCO.gamma0)
    sigma_px = np.sqrt(neps_x / beta_x / partCO.beta0 / partCO.gamma0)
    sigma_py = np.sqrt(neps_y / beta_y / partCO.beta0 / partCO.gamma0)

    x_wrt_CO = np.random.normal(loc=0.0, scale=sigma_x, size=n_macroparticles)
    y_wrt_CO = np.random.normal(loc=0.0, scale=sigma_y, size=n_macroparticles)
    px_wrt_CO = np.random.normal(loc=0.0, scale=sigma_px, 
                            size=n_macroparticles) - alfa_x/beta_x * x_wrt_CO
    py_wrt_CO = np.random.normal(loc=0.0, scale=sigma_py, 
                            size=n_macroparticles) - alfa_y/beta_y * y_wrt_CO
    # account for dispersion
    x_wrt_CO += delta * dx
    y_wrt_CO += delta * dy
    px_wrt_CO += delta * dpx
    py_wrt_CO += delta * dpy

    return x_wrt_CO, y_wrt_CO, px_wrt_CO, py_wrt_CO
	



    # x,xp,y,yp = [],[],[],[]
    # for epsn_x, epsn_y, intensity in zip(np.atleast_1d(p['epsn_x']), 
    # 		np.atleast_1d(p['epsn_y']), np.atleast_1d(p['intensity'])):

    # 	Transverse_distribution = GaussDist2D(twissX, twissY, cut_off=p['TransverseCut'])
    # 	n_macroparticles_tmp = int(p['n_macroparticles']*(intensity/np.sum(p['intensity'])))
    # 	Transverse_coords = np.array(map(lambda i: Transverse_distribution.getCoordinates(), xrange(n_macroparticles_tmp)))
    # 	x.extend(Transverse_coords[:,0].tolist())
    # 	xp.extend(Transverse_coords[:,1].tolist())
    # 	y.extend(Transverse_coords[:,2].tolist())
    # 	yp.extend(Transverse_coords[:,3].tolist())
    # # in case x has not yet a length of n_macroparticles
    # while len(x)<p['n_macroparticles']:
    # 	Transverse_coords = Transverse_distribution.getCoordinates()
    # 	x.append(Transverse_coords[0])
    # 	xp.append(Transverse_coords[1])
    # 	y.append(Transverse_coords[2])
    # 	yp.append(Transverse_coords[3])
	# x = np.array(x) + p['x0']  + dpp * p['etax0']
	# xp = np.array(xp) + p['xp0'] + dpp * p['etapx0']
	# y = np.array(y) + p['y0']  + dpp * p['etay0']
	# yp = np.array(yp) + p['yp0'] + dpp * p['etapy0']
    
	# return output_file
