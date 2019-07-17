import numpy as np
from astropy import constants

from einsteinpy.coordinates import BoyerLindquistDifferential
from odessay_initialcondition import initial_velocity_photon , coord_photon

class MotionConstants:

    """
    Class for defining Kerr Motion constant
    """

    @u.quantity_input(M=u.kg,
            r_obs =u.km,
            theta_obs=u.rad,
            phi_obs=u.rad ,
            a=u.km)

    def __init__(self, r_obs ,theta_obs, phi_obs,  a ):
        self.a = a
        self.M = M
        self.r_obs= r_obs
        self.theta_obs= theta_obs
        self.phi_obs= r_obs


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



    def delta(r, a , M):
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
            The value r^2 - 2 * r * M+ a^2

        """
        return (r ** 2) - (2 * r * M) + (a ** 2)


    def energy(self, x , y , z):

        """
        Returns Energy (E) which is constant along the geodesic.
        --------------------------------------------------

        Parameters
        ----------
        r_obs : float
                The observer is located at a distance r_obs from
                the black hole center

        theta_obs: float
                The observer is located at an angle theta_obs
                from the positive black hole z'-axis
                (coinciding with the spin axis)

        phi_obs: float
                The observer is located at an angle phi_obs
                with respect to the black hole’s x′-axis

        M : float
            Mass of massive body

        a : float
            Spin factor

        x , y , z: float
            Coordinates in observer's coordinate system i.e. observer grid

        x_bh , y_bh , z_bh:
            Coordinates in black hole's  coordinate system

        """

        coord = coord_photon(self)
        r = coord[0]
        theta = coord[1]
        sg, dl ,AA = sigma(r, theta, a), delta(r, M, a), AA(r, theta, a)
        v_ini = initial_velocity_photon(self)

        return (((sg - (2 * r))/(sg * dl))((sg * v_ini[0] * v_ini[0]) +
                (sg * dl * v_ini[1] * v_ini[1])) - (dl * restmass))
                + (dl * (v_ini[2] * np.cos(theta)) ** 2 )


    def angular_momentum_z(self):

        """
        Returns Angular momentum (Lz) which is constant along the geodesic.
        """

        coord = coord_photon(self)
        r = coord[0]
        theta = coord[1]
        sg, dl ,AA = sigma(r, theta, a), delta(r, M, a), AA(r, theta, a)
        v_ini = initial_velocity_photon(self)
        E = energy(self)

        return (((dl * sg * v_ini[2]) - (2 * a * r * E)) * ((np.sin(theta)) ** 2))
                / (sg - (2 * r))


    def carter_const():

        """
        Returns Carter Constant (Q) which is constant along the geodesic.
        """
        lz = angular_momentum_z(self)
        E = energy(self)
        v_ini = initial_velocity_photon(self)
        p_theta =  v_ini[1] * E

        coord = coord_photon(self)
        theta = coord[1]

        return ((p_theta) ** 2) +
                (((lz * np.cosec(theta)) ** 2)
                    - ((a * a) 8 (E * E + restmass))) * (np.cos(theta) * np.cos(theta))
