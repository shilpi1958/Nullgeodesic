import astropy.units as u
import numpy as np

from einsteinpy.utils import kerr_utils
from initialcondition import InitialConditionsBHFrame


class MotionConstants(InitialConditionsBHFrame):

    """
    Class for defining Kerr Motion constant
    """

    @u.quantity_input(
        M=u.kg,
        r_obs=u.km,
        theta_obs=u.rad,
        phi_obs=u.rad,
        a=u.km,
        x=u.km,
        y=u.km,
        z=u.km,
    )
    def __init__(self, r_obs, theta_obs, phi_obs, a, M, x, y, z):

        super(MotionConstants, self).__init__(r_obs, theta_obs, phi_obs, a)
        self.M = M
        self.x = x
        self.y = y
        self.z = z

    def _delta(self, r, a):
        """
        Parameters
        ----------
        r : float
            Component r in vector
        a : float
            Any constant
        Returns
        -------
        float
            The value r^2 + r + a^2
        
        """
        return (r ** 2) - (r * u.km) + (a ** 2)

    def energy(self):

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

        coordinates = self.coord_photon(self.x, self.y, self.z)
        coords = coordinates.si_values()
        r = coords[0] * u.km
        theta = coords[1] * u.rad
        sg = kerr_utils.sigma(r, theta, self.a)
        dl = self._delta(r, self.a)
        v_ini = self.initial_velocity_photon(self.x, self.y, self.z)
        restmass = 0

        A = (sg - (2 * r * u.km)) / dl
        B = (v_ini[0] * v_ini[0]) + (dl * v_ini[1] * v_ini[1])
        return np.sqrt((A * B) - (dl * (v_ini[2] * np.cos(theta)) ** 2))

    def angular_momentum_z(self):

        """
        Returns Angular momentum (Lz) which is constant along the geodesic.
        """

        coordinates = self.coord_photon(self.x, self.y, self.z)
        coords = coordinates.si_values()
        r = coords[0] * u.km
        theta = coords[1] * u.rad
        sg = kerr_utils.sigma(r, theta, self.a)
        dl = kerr_utils.delta(r, self.M, self.a, c=constant.c.value, G=constant.G.value)
        v_ini = self.initial_velocity_photon(self.x, self.y, self.z)
        restmass = 0
        E = self.energy()

        A = (dl * sg * v_ini[2]) - (2 * a * r * E)
        B = sg - (2 * r)

        return (A * np.sin(theta) * np.sin(theta)) / B

    def carter_const(self):

        """
        Returns Carter Constant (Q) which is constant along the geodesic.
        """
        lz = self.angular_momentum_z()
        E = self.energy()
        v_ini = self.initial_velocity_photon(self.x, self.y, self.z)
        p_theta = v_ini[1] * E
        restmass = 0

        coordinates = self.coord_photon(self.x, self.y, self.z)
        coords = coordinates.si_values()
        theta = coords[1] * u.rad

        A = lz * np.cosec(theta)
        B = (a * a) * (E * E + restmass)

        return (p_theta * p_theta) + (A - B) * (np.cos(theta) * np.cos(theta))
