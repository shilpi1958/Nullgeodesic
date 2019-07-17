
import numpy as np
from scipy import constants as const
from astropy import constants

from einsteinpy.utils import kerr_utils
from einsteinpy.coordinates import BoyerLindquistDifferential,

_h = const.h.value
_hbar = constants.hbar.value
_c = constants.c.value
Rs = 2

# LNRF frame
#$\bar{p}_{(a)}$, which are the components of the four-momentummeasured in the
#locally nonrotating reference frame (LNRF) and have been specified by the user.
class LNRF:

    @u.quantity_input(M=u.kg)
    def __init__(self, bl_coords, M):
        self.input_coords = bl_coords
        self.a = self.input_coords.a.to(u.m)
        self.M = M

        )

    def sigma(r, theta, a):
        """
        Returns the value r^2 + a^2 * cos^2(theta)
        Specific to Boyer-Lindquist coordinates
        Parameters
        ----------
        r : float
            Component r in vector
        theta : float
            Component theta in vector
        a : float
            Spin factor
        Returns
        -------
        float
            The value r^2 + a^2 * cos^2(theta)

        """
        return (r ** 2) + ((a * np.cos(theta)) ** 2)



    def delta(r,a):
        """
        Returns the value r^2 - 2 * r + a^2
        Specific to Boyer-Lindquist coordinates
        Parameters
        ----------
        r : float
            Component r in vector
        M : float
            Mass of massive body
        a : float
            Spin factor
        Returns
        -------
        float
            The value r^2 - 2 * r + a^2

        """
        return (r ** 2) - (2 * r) + (a ** 2)



    def AA(r, theta, a):
        """
        Returns the value r^2 - 2 * r + a^2
        Specific to Boyer-Lindquist coordinates
        Parameters
        ----------
        r : float
            Component r in vector
        M : float
            Mass of massive body
        a : float
            Spin factor
        Returns
        -------
        float
            The value r^2 - 2 * r + a^2

        """
        return ((r ** 2) + (a ** 2)) - delta(r, M, a)(a ** 2)(np.sin(theta) ** 2)




    def observer_metric(r, theta, M, a):

        """
        Returns the basis vectors of the orthonormal tetrad of the observers
        Specific to Boyer-Lindquist coordinates

        Parameters
        ----------
        theta : float
            Angle from z-axis
        r : float
            Component r in vector
        M : float
            Mass of massive body
        a : float
            Spin factor
        Returns
        -------
        ~numpy.array
            Numpy array of shape (4,4)

        """

        sg, dl ,AA = sigma(r, theta, a), delta(r, M, a), AA(r, theta, a)
        # set the diagonal/off-diagonal terms of metric
        m[0, 0] = (AA / (sg * dl)) ** 0.5
        m[1, 1] = (dl / sg) ** 0.5
        m[2, 2] = (1 / sg) ** 0.5
        m[3, 3] = (sg / (AA * (np.sin(theta) ** 2))) ** 0.5
        m[0, 3] = ((AA / (sg * dl)) ** 0.5) * ((2 * a * r) / AA)
        return m



class MotionConstants:

    """
    Class for defining Kerr Motion constant
    """

    @u.quantity_input(M=u.kg)
    def __init__(self, bl_coords, M):
        self.input_coords = bl_coords
        self.a = self.input_coords.a.to(u.m)
        self.M = M

    def carter_const(m , a , E, theta L_theta, Lz):

        """
        Returns the value of Carter Constant
        Specific to Boyer-Lindquist coordinates

        Parameters
        ----------
        theta : float
            Component theta in pos_vector
        M : float
            Mass of the particle
        a : float
            Black hole Spin factor
        E : float
            Energy measured by aa distant observer
        L_z : float
            Particle's axial angular momentum
        L_theta: float
            Latitudinal component of the particle's angular momentum

        Returns
        -------
        float
        The value p_theta^2 + cos(theta)^2(a^2(M^2 - E^2) + (J/sin(theta)^2))

        """

        Q = (L_theta ** 2) +
            ((np.cos(theta) ** 2)((a ** 2)(M ** 2 - E * *2) +
            ((L_z/sin(theta)) ** 2)))
        return Q


#Calculation of Motion Constants from Impact Parameters
#impact parameters α(alpha), β(beta) coordinates
#theta_obs is the inclination angle of the observer.
#four momentum not given by user

    def lambda_const_impact_para(alpha, theta_obs):

        """
        Returns the value of lambda constant
        Given the impact parameters

        Parameters
        ----------
        alpha : float
            coordinate
        theta_obs: float
            Mass of the particle

        Returns
        -------
        float
            The value - (alpha * sin(theta_obs))

        """
        return  = - (alpha * np.sin(theta_obs))

    def q_const_impact_para(alpha, beta, theta_obs, a):

        """
        Returns the value of q constant
        Given the impact parameters

        Parameters
        ----------
        alpha : float
            coordinate
        beta : float
            coordinate
        theta_obs: float
            Mass of the particle
        a : float
            Spin factor

        Returns
        -------
        float
            The value of
             (beta ** 2) + ((alpha ** 2) - (a ** 2))(cos(theta_obs) ** 2)
        """
        return (beta ** 2) + ((alpha ** 2) -
                (a ** 2))(np.cos(theta_obs) ** 2)




#Calculation of Motion Constants from four momentum
#compute λ and q and pr, pθ from the components of the four-momentum
#measured in the (LNRF) and have been specified by the user

    def lambda_const_from_momentum(p_phi_bar, p_t_bar, theta, r, a):

        """
        Returns the value of q constant
        Given the impact parameters

        Parameters
        ----------
        theta : float
            Angle from z-axis
        r : float
            Component r in vector
        p_phi: float
            phi component of four momrntum in LNRF frame
        p_t : float
            t component of four momrntum in LNRF frame
        a : float
            Spin factor

        Returns
        -------
        float
            The value of (sin(theta)) * x) /
                        ((w * (sin(theta)) * x) - (((dl ** 0.5) * sg) / aa))

        """

        dl, sg, aa  = LNRF.delta(r, a) , LNRF.sigma(r, theta, a) , LNRF.AA(r, theta, a)
        w =((2 * a * r) / AA)
        x = (p_phi_bar/p_t_bar)

        return ((np.sin(theta)) * x) /
                        ((w * (np.sin(theta)) * x) - (((dl ** 0.5) * sg) / aa))


    def q_const_from_momentum(p_phi_bar, p_t_bar, theta, r, a):

        """
        Returns the value of q constant
        Given the impact parameters

        Parameters
        ----------
        theta : float
            Angle from z-axis
        r : float
            Component r in vector
        p_phi: float
            phi component of four momrntum in LNRF frame
        p_t : float
            t component of four momrntum in LNRF frame
        a : float
            Spin factor

        Returns
        -------
        float
            The value of ((((lambda_const/sin(theta)) ** 2)
                        - (a ** 2))(cos(theta) ** 2)) +
                        (((x * (1 - (lambda_const * w))) ** 2) * (aa / sg))
        """


        dl, sg, aa  = LNRF.delta(r, a) , LNRF.sigma(r, theta, a) , LNRF.AA(r, theta, a)
        w =((2 * a * r) / AA)
        x = (p_phi_bar/p_t_bar)
        lambda_const = lambda_const_from_momentum(p_phi_bar, p_t_bar, theta, r, a)

        return ((((lambda_const/np.sin(theta)) ** 2)
                    - (a ** 2))(np.cos(theta) ** 2)) +
                    (((x * (1 - (lambda_const * w))) ** 2) * (aa / sg))
