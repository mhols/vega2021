import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.widgets import Slider
from scipy.interpolate import PPoly, make_interp_spline


def TPXYZ(theta, phi):
    return np.array(
        [
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta)
        ]
    )

def XYZTP(x,y,z):
    r = np.sqrt( x**2 + y**2 + z**2)
    return np.array(
        [
            np.arccos(z),
            np.arctan2(y,x)
        ]
    )

def YZTP(y,z):
    x = np.sqrt(1-y**2 -z**2)
    return XYZTP(x,y,z)

def rotation_matrix_Z(om):
    co = np.cos(om)
    so = np.sin(om)
    return np.array(
        [[co, -so, 0],[so, co, 0],[0,0,1]]
    )

def rotation_matrix_X(om):
    co = np.cos(om)
    so = np.sin(om)
    return np.array(
        [[1, 0, 0],[0, co, -so],[0,so,co]]
    )

def rotation_matrix_Y(om):
    co = np.cos(om)
    so = np.sin(om)
    return np.array(
        [[co, 0, so],[0, 1, 0],[-so, 0, co]]
    )

def rot_Y_Z(omx, omz, p):
    
    v = np.dot(rotation_matrix_Z(omz),p)
    return np.dot( rotation_matrix_Y(omx), v)

def vega_trace(theta, phi, phases=None):
    if phases is None:
        phases = np.linspace(0, 2*np.pi, 1024)
    
    resvrad, resvisi = [], []
    for phase in phases:
        p = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
        v = rot_Y_Z(83*np.pi/180, phase, p.T)
        resvrad.append(v[1])
        resvisi.append(v[0]<0)

    return np.array(resvrad), np.array(resvisi), phases


class StarSpotFinder:

    def __init__(self, **kwargs):
        self.kwargs = kwargs

        self.fig, (self.ax0, self.ax1) = plt.subplots(1,2, width_ratios=[1, 1])
        
        self.axamp = self.fig.add_axes([0.9, 0.1, 0.1, 0.8])
        
        # Make a vertically oriented slider to control the amplitude
        self.amp_slider = Slider(
        ax=self.axamp,
            label="p",
            valmin=0,
            valmax=2*np.pi,
            valinit=0,
            orientation="vertical"
        )



    

        # register the update function with each slider
        self.amp_slider.on_changed(self._change_phase)

        self.velocity = np.loadtxt(self.kwargs['vbins'])

        self.velocity = velocity

        # inernal state of possible spots (for the moment only one)
        self.phase = 0
        self.theta = np.pi/4
        self.phi = 0 # longitude on star in its restframe

        self.clickvelocity = 0
        self.clickphase = 0

        self._plot_background()
        self._plot_star_disk()

        self._plot_orbit_on_star_and_picture()

        self.fig.canvas.mpl_connect("button_press_event", self._switch_on)
        self.fig.canvas.mpl_connect("motion_notify_event", self._move_the_point)
        self.fig.canvas.mpl_connect("button_release_event", self._switch_off)

        self._on = False
        self._mod = 'phase'


    def _switch_on(self, event):

        if event.inaxes is self.ax0:
            self.clickvelocity, self.clickphase = event.xdata, event.ydata
            print(self.clickvelocity, self.clickphase )
        else:
            self._on = True



    def _move_the_point(self, event):
        if event.inaxes is self.ax1 and self._on and self._mod == "phase":

            self.theta, self.phi = self.YZTP(event.xdata, event.ydata)
            

        """
        if event.inaxes is self.ax0:
            if self._mod == 'phase':
                self.phase = event.ydata
            else:
                roots = self._points_with_v_on_orbit( self.clickvelocity)
                if len(roots) == 2:
                    self.phase = self.clickphase-roots[0]

            for art in list(self.ax1.lines)+list(self.ax0.lines):
                art.remove()
            self._plot_orbit_on_star_and_picture()

        self.fig.canvas.draw()
        """

        for art in list(self.ax1.lines)+list(self.ax0.lines):
            art.remove()
        self._plot_orbit_on_star_and_picture()
        self.fig.canvas.draw()

    def _switch_off(self, event):
        self._on = False

    # The function to be called anytime a slider's value changes
    def _change_phase(self, val):
        self.phase = self.amp_slider.val
        for art in list(self.ax1.lines)+list(self.ax0.lines):
            art.remove()
         
        self._plot_orbit_on_star_and_picture()

        self.fig.canvas.draw()

    def velocity(self, y):
        return y * self.kwargs['vmax']
    
    def y(self, v):
        return v / self.kwargs['vmax']

    @property
    def sampled_phases(self):
        return np.linspace(0, 2*np.pi, 1024, endpoint=False)


    def points_on_star_of_velocity(self,v):
        """the xyz coordinates in star coordinates or points of velocity v

        :param v: velocity
        :type v: float
        """
        y = self.y(v)

        phases = self.sampled_phases

        orbit_o = np.row_stack(
            (
                np.sin(phases),
                y,
                np.cos(phases)
            )
        )
        orbit = np.dot(self.rotation_earth_to_star(), orbit_o)

        return orbit
    
    def points_on_star_passing_through_velocity_phase(self, v, phase):
        """the point in xyz star coordinates that pass through the point v, phase

        :param v: velocity
        :type v: float
        :param phase: the phase
        :type phase: float (mod 2)
        """

        pass


    def _points_with_v_on_orbit(self, v):
        """returns the normlized y coordinate and
            points ([phasef, zb], [y], or []) with prescribed v

        :param v: radial velocity
        :type v: float
        """

        orbit, visi, phases = self.YZOrbit()

        yv = orbit[1,:]

        y = self.y(v)

        p = PPoly.from_spline(make_interp_spline(phases, yv-y, k=1 )).roots()
        
        return [pp for pp in p if pp*(2*np.pi - pp) > 0]


    def spots_compatible_with_v_at_phase(self, v, phase):
        """locus of points at phase=0 that pass through v, phase        

        :param v: _description_
        :type v: _type_
        :param phase: _description_
        :type phase: _type_
        """

        pass

    def rotation_to_observer(self):
        """turns the star to the observer
        """
        
        return rotation_matrix_Y(np.pi/2 - self.kwargs['angle'])
    
    def rotation_star_to_earth(self):
        return np.dot(self.rotation_to_observer(), rotation_matrix_Z(self.phase))
    
    def rotation_earth_to_star(self):
        return self.rotation_star_to_earth().T
    
    def orbit_in_star_coord_XYZ(self, y=None, z=None):
        pass
    

    def YZTP(self, y, z):
        t, p = YZTP(y, z)
        v = TPXYZ(t,p)
        v = np.dot(self.rotation_to_observer().T , v)
        t, p = XYZTP(*v)
        return t, p+self.phase



    def XYZ(self):
        """the projected star disk coordinates of the  spot and a boolean visibility
        """
        v = TPXYZ(self.theta, self.phi-self.phase)
        v = np.dot( self.rotation_to_observer(self.phase) , v)
        return v
    
    def YZOrbit(self):
        n = 1024

        phases = self.sampled_phases

        xyz = []

        for ph in phases:
            v = TPXYZ(self.theta, self.phi-self.phase+ph)
            v = np.dot( self.rotation_to_observer() , v)
            xyz.append( v )

        orbit = np.column_stack(xyz)
        visi = orbit[0,:] > 0

        return orbit, visi, phases

    def _plot_background(self):
        res = np.loadtxt(self.kwargs['image'])
        velocity = self.velocity
        v0 = self.kwargs['v0']
        self.ax0.set_xlim(velocity[0]-v0, velocity[-1]-v0)
        self.ax0.imshow(-res, cmap=plt.cm.gray_r, 
                       aspect='auto',interpolation='bicubic', 
                       origin='lower',extent=[velocity[0]-v0, velocity[-1]-v0,0,2*np.pi],
                       vmin=-0.0002,vmax=0.0002)
        

    def _plot_star_disk(self):

        
        self.ax1.set_aspect('equal')
        self.ax1.set_xlim(-1.1, 1.1)
        self.ax1.set_ylim(-1.1, 1.1)

        stardisk =  Circle((0, 0), 1)

        north = Circle((0, np.sin(self.kwargs['angle'])), 0.02, color='red')

        
        self.ax1.add_artist(stardisk)
        self.ax1.add_artist(north)


    def _plot_orbit_on_star_and_picture(self):
        """_summary_

        :param y: y coordinate on projected unit star disk
        :type y: float in [-1, 1]
        :param z: z coordinate on projected unit star disk
        :type z: float in [-1, 1]
        """

        orbit, visi, phases = self.YZOrbit()

        self.ax1.plot( orbit[1,visi], orbit[2, visi], '.k', markersize=0.2)
        self.ax0.plot( orbit[1, visi]*self.kwargs['vmax'], phases[visi], '.r', markersize=0.2)
        self.ax1.plot( [orbit[1,0]], [orbit[2,0]], 'ok', markersize=10)

        self.ax0.plot( [self.clickvelocity], [self.clickphase], 'or', markersize=5)









    


if __name__ == '__main__':

    import matplotlib.pyplot as plt

    image = "moving_simple_per_night.txt"
    vbins = "velocitybins.txt"

    res = np.loadtxt(image)
    velocity = np.loadtxt(vbins)




            


    #for theta in [np.pi/2, np.pi/3, np.pi/4, np.pi/8, np.pi/16]:
    #    vr, vi, phases = vega_trace(theta, 0)
    #    plt.plot(vr[vi], phases[vi], 'o')


    star = StarSpotFinder(
        image = "moving_simple_per_night.txt",
        vbins = "velocitybins.txt",
        v0=-13.01,
        vmax = 22,
        angle= 7 * np.pi/180
        )
    
    


    plt.show()