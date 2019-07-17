import numpy as np
from scipy import constants as const
from astropy import constants

from einsteinpy.utils import kerr_utils
from einsteinpy.coordinates import BoyerLindquistDifferential, Cartesian
from kerrmotionconstants import MotionConstants , LNRF


class KerrNullGeodesic:

    """
    Class for defining Kerr Null Geodesic if impact parameters given
    """

    @u.quantity_input(M=u.kg)
    def __init__(self, bl_coords, M ):
        self.input_coords = bl_coords
        self.a = self.input_coords.a.to(u.m)
        self.M = M

        pos_vec, vel_vec = (
            self.input_coords.si_values()[:3],
            self.input_coords.si_values()[3:],
        )

    if alpha, beta

        lambda_const = MotionConstants.lambda_const_impact_para(alpha, theta_obs)
        q_const = MotionConstants.q_const_impact_para(alpha, beta, theta_obs, a)

        def R_motion(r):
            return r ** 4 - (q_const + lambda_const * lambda_const - a * a)(r * r)
                    + 2 * (q_const + ((lambda_const - a) ** 2)) * r - (q * a * a)

        def theta_motion(theta):
            return q_const + (a * np.cos(theta)) ** 2
                            - (lambda_const * np.cot(theta)) ** 2

        def T_motion(r, a):
            return r ** 2 + a ** 2 - a * lambda_const

#covariant components of the four-momentum of a photon in B–L coordinates
        def four_momentum(E, r, theta):
            R , theta_m,= R_motion(r) , theta_motion(theta)
            dl, sg = LNRF.delta(r, a) , LNRF.sigma(r, theta, a)

            p_momentum = np.zeros(shape=(4, 0), dtype=float)
            p_momentum[0] =  - (E)
            p_momentum[1] = ((R ** 0.5) * E ) / dl
            p_momentum[2] = (theta_m ** 0.5)  * E
            p_momentum[3] = lambda_const  * E
            return p_momentum

#components of the four-momentum of a photon in LNRF are given
    if p_phi_bar, p_t_bar

        lambda_const = MotionConstants.lambda_const_from_momentum
                                        (p_phi_bar, p_t_bar, theta, r, a)
        q_const = MotionConstants.q_const_from_momentum
                                        (p_phi_bar, p_t_bar, theta, r, a)

        def R_motion(r):
            return r ** 4 - (q_const + lambda_const * lambda_const - a * a)(r * r)
                    + 2 * (q_const + ((lambda_const - a) ** 2)) * r - (q_const  * a * a)

        def theta_motion(theta):
            return q_const + (a * np.cos(theta)) ** 2 - (lambda_const * np.cot(theta)) ** 2

        def T_motion(r, a):
            return r ** 2 + a ** 2 - a * lambda_const
