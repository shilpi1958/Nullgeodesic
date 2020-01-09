import numpy as np
from astropy import constants
import astropy.units as u

from einsteinpy.utils import kerr_utils
from einsteinpy.coordinates import Cartesian ,CartesianDifferential

class InitialConditionsBHFrame:

    """
    Class for defining initial conditions for photon rays in black hole frame
    """

    @u.quantity_input(
        r_obs =u.km,
        theta_obs=u.rad,
        phi_obs=u.rad ,
        a=u.km
        )

    def __init__(self, r_obs ,theta_obs, phi_obs,  a ):
        self.a = a
       #self.M = M
        self.r_obs= r_obs
        self.theta_obs= theta_obs
        self.phi_obs= phi_obs

    def coord_photon(self , x , y , z):

        """

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

        Returns
        -------
        array
            (r, θ, φ) conditions for a photon on the observer grid

        """

        D = (((np.sqrt((self.r_obs * self.r_obs) + (self.a * self.a))) - z) * np.sin(self.theta_obs))
        -(y * np.cos(self.theta_obs))


        x_bh = (D * np.cos(self.phi_obs)) - (x * np.sin(self.phi_obs))
        y_bh = (D * np.sin(self.phi_obs)) + (x * np.cos(self.phi_obs))
        z_bh = (self.r_obs - z) * np.cos(self.theta_obs) + (y * np.sin(self.theta_obs))

        blackhole_grid = Cartesian(x_bh , y_bh , z_bh)

        return blackhole_grid

    def coord_photon_bl(self , x , y , z):

        """

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

        Returns
        -------
        array
            (r, θ, φ) conditions for a photon on the observer grid

        """

        D = (((np.sqrt((self.r_obs * self.r_obs) + (self.a * self.a))) - z) * np.sin(self.theta_obs))
        -(y * np.cos(self.theta_obs))


        x_bh = (D * np.cos(self.phi_obs)) - (x * np.sin(self.phi_obs))
        y_bh = (D * np.sin(self.phi_obs)) + (x * np.cos(self.phi_obs))
        z_bh = (self.r_obs - z) * np.cos(self.theta_obs) + (y * np.sin(self.theta_obs))

        blackhole_grid = Cartesian(x_bh , y_bh , z_bh)

        return blackhole_grid.to_bl(self.a)

    def initial_velocity_photon(self, x , y ,z):

        """

        Each ray arrives perpendicular to the image plane, moving
        parallel to the z-axis, hence we set velocity
        in observer grid as = (0, 0, 1)

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

        x_vel , y_vel, z_vel: (0, 0, 1)
            Cartesian components of the ray’s velocity in observer grid

        x_vel_bh , y_vel_bh , z_vel_bh:
            Cartesian components of the ray’s velocity in black hole
            coordinates

        Returns
        -------
        array
             Initial Velocity  of photon in BL coordinates
             on the observer grid

            """

        coordinates = self.coord_photon(x , y , z)
        coords = coordinates.si_values()
        x_bh = coords[0]  * u.km
        y_bh = coords[1]  * u.km
        z_bh = coords[2]  * u.km


        # velocity in observer grid : (0, 0, 1)
        x_vel_bh = - (np.cos(self.phi_obs) * np.sin(self.theta_obs)) * u.km / u.s
        y_vel_bh = - (np.sin(self.phi_obs) * np.sin(self.theta_obs)) * u.km / u.s
        z_vel_bh = - (np.cos(self.theta_obs)) * u.km / u.s

        photon_velocity = CartesianDifferential(x_bh, y_bh, z_bh, x_vel_bh , y_vel_bh , z_vel_bh )

        return photon_velocity.bl_differential(self.a)
