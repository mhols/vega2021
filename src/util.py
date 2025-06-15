import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.widgets import Slider
from scipy.interpolate import PPoly, make_interp_spline
from mpl_toolkits.basemap import Basemap

def TPXYZ(theta, phi):
    return np.row_stack(
        [
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta)
        ]
    )

def XYZTP(x,y,z):
    r = np.sqrt( x**2 + y**2 + z**2)
    return np.row_stack(
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

    def __init__(self, velocity_bins=None, movingpeak=None, interactive=True, **kwargs):
        self.kwargs = kwargs

        self.phase = 0
        self.theta = np.pi/4
        self.phi = 0 # longitude on star in its restframe
        if velocity_bins is None:
            self.velocity_bins = np.loadtxt(self.kwargs['vbins'])
        else:
            self.velocity_bins=velocity_bins
        if movingpeak is None:
            self.movingpeak = np.loadtxt(self.kwargs['image'])
        else:
            self.movingpeak = movingpeak
 
        if interactive:
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



            # inernal state of possible spots (for the moment only one)
            self.clickvelocity = 0
            self.clickphase = 0
            self._slide_phase = 0

            self.line_image = []

            self._on = False
            self._mod = ''

            self.line_star_earth_phi = []  # XYZ earth coordinates of orbit under different phases
            self.line_star_earth_vel = []  # XYZ earth coordinates of orbit under constant velocity
            self.line_velocities = []      # velocities 
        
            self.line_star_star_phi = []   # XYZ star coordinates of orbit under different phases 
            self.line_star_star_vel = []   # XYZ star coordinates of orbit under different phases 


            self._update_line_image()
            self._update_line_star_earth_phi()
            self._update_line_star_earth_vel()
            self._update_line_star_star()


            self._plot_background()
            self._plot_star_disk()

            self._plot_orbit_on_star_and_picture()

            self.fig.canvas.mpl_connect("button_press_event", self._switch_on)
            self.fig.canvas.mpl_connect("motion_notify_event", self._move_the_point)
            self.fig.canvas.mpl_connect("button_release_event", self._switch_off)

 

    def _switch_on(self, event):

        self._on = True
        self._move_the_point(event)



    def _move_the_point(self, event):
        if event.inaxes is self.ax0 and self._on:
            self.clickvelocity, self.clickphase = event.xdata, event.ydata

            #r = self.points_on_star_of_velocity(self.clickvelocity, self.clickphase)
            #self.theta, self.phi = XYZTP(*r)
            self._update_line_star_earth_vel()
        
        if event.inaxes is self.ax1 and self._on: # and self._mod == "phase":

            return
            #self.theta, self.phi = self.YZTP(event.xdata, event.ydata)
            #self.clickvelocity = self.velocity(event.xdata)
            #self.clickphase = self.phase
            #self._update_line_star_earth_phi() 
        
        


                
        self._clear_and_redraw()

    def _clear_and_redraw(self):
        for art in list(self.ax1.lines)+list(self.ax0.lines):
            art.remove()
        self._plot_orbit_on_star_and_picture()
        self.fig.canvas.draw()

    def _switch_off(self, event):
        self._on = False

    # The function to be called anytime a slider's value changes

    def _change_phase(self, val):
        self._slide_phase = val
        #self.phi = self.amp_slider.val
        i = int(self.sampled_phases.shape[0] * val / (2 * np.pi) )
        self.theta, self.phi = XYZTP( *self.line_star_star_vel[:,i])

        self._update_line_star_earth_phi()
        

        #self._clear_and_redraw()

    def _update_line_image(self):
        pass

    def _update_line_star_earth_vel(self):
        if not self._on:
            return
        
        self.line_star_earth_vel, self.line_star_star_vel = \
            self.spots_compatible_with_v_at_phase(self.clickvelocity, self.clickphase)
        
        self.theta, self.phi = XYZTP( *self.line_star_star_vel[:, 0])

        self.line_star_earth_phi, self.line_star_star_phi = self.spots_constant_latitude_on_star(self.theta)
        self.line_velocities = self.kwargs['vmax'] * self.line_star_earth_phi[1, :]

    def _update_line_star_earth_phi(self):
        
        self.line_star_earth_phi, self.line_star_star_phi = self.spots_constant_latitude_on_star(self.theta)
        self.line_velocities = self.kwargs['vmax'] * self.line_star_earth_phi[1, :]



    def _update_line_star_star(self):
        pass

    def velocity(self, y):
        return y * self.kwargs['vmax']
    
    def y(self, v):
        return v / self.kwargs['vmax']

    @property
    def sampled_phases(self):
        return np.linspace(0, 2*np.pi, 1024, endpoint=False)

    

    def points_on_star_of_velocity(self, v, phases=None):
        """the xyz coordinates in star coordinates or points of velocity v

        :param v: velocity
        :type v: float
        """

        if phases is None:
            phases = self.sampled_phases
        else:
            phases = np.atleast_1d(phases)
        y = np.full(phases.shape[0], self.y(v))

        r = np.sqrt(1-y**2)
        orbit_o = np.row_stack(
            (
                r * np.sin(phases),
                y,
                r * np.cos(phases)
            )
        )
        
        return orbit_o
    

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
        raise Exception()
        
        return [pp for pp in p if pp*(2*np.pi - pp) > 0]


    def spots_compatible_with_v_at_phase(self, v, phase):
        """locus of points at phase=0 that pass through v, phase        

        :param v: _description_
        :type v: _type_
        :param phase: _description_
        :type phase: _type_
        """

        xyz = self.points_on_star_of_velocity(v)
        #t,p = XYZTP(*xyz)

        #p -= phase

        #v = TPXYZ(t, p)

        r = np.dot(rotation_matrix_Z(-phase), xyz)

        return np.dot(self.rotation_star_to_earth(), r), r
    
    def spots_constant_latitude_on_star(self, theta):
        phases = self.sampled_phases
        ct = np.cos(theta)
        st = np.sin(theta)

        v = np.row_stack(
            (
               st * np.cos(phases),
               st * np.sin(phases),
               np.full( phases.shape[0], ct) 
            )
        )

        return np.dot( self.rotation_star_to_earth(), v), v
    
    def likelihood_spot(self, theta, phi, phase=0):
        """the likelihood to have a starspot at theta, phi (at phase)

        :param theta: _description_
        :type theta: _type_
        :param phi: _description_
        :type phi: _type_
        :param phase: _description_, defaults to 0
        :type phase: int, optional
        """
        phases = np.linspace(0, 2*np.pi, self.movingpeak.shape[0])
        ct = np.cos(theta)
        st = np.sin(theta)

        v = np.vstack(
            (
               st * np.cos(phi+phase+phases),
               st * np.sin(phi+phase+phases),
               np.full( phases.shape[0], ct) 
            )
        )

        v0 = self.kwargs['v0']
        vv = np.dot( self.rotation_star_to_earth(), v) 
        vr = self.kwargs['vmax'] * vv[1, :] + self.kwargs['v0']
        visi = vv[0,:] > 0

        dv = self.velocity_bins[1] - self.velocity_bins[0]

        vel = np.zeros(self.velocity_bins.shape[0] + 1)
        vel[1:] = self.velocity_bins + dv/2
        vel[0] = self.velocity_bins[0] - dv/2

        j = np.digitize(vr, vel)
        # print(j)

        res = 0
        for i in range(phases.shape[0]):
            if visi[i]:
                res += self.movingpeak[i, j[i]] 

            #print(i, j[i], self.movingpeak[i, j[i]]) 
        return res




    def lh_spot(self):


        theta=np.linspace(0, np.pi, 128, endpoint=True)
        phi = np.linspace(0, 2*np.pi, 128)


        #p, t = map.makegrid(nx=128, ny=128)

        #print('couout', p, t)

        t, p = np.meshgrid(theta, phi)
        #x,y = map(t, p)

        res = np.zeros(len(theta)*len(phi))

        for i, tp in enumerate(zip(t.ravel(), p.ravel())):
            res[i] = self.likelihood_spot(tp[0], tp[1])
    
        res = np.reshape(res, (len(phi), len(theta) ))

        #print(res)
        return theta, phi, res



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
        """Y Z in earth coordinates, T P on star"""

        x = np.sqrt(1 - y**2 - z**2)
        v = np.row_stack( (x, y, z))

        v = np.dot(self.rotation_earth_to_star(), v)

        return XYZTP(*v)
    
    

    def XYZ(self):
        """the projected star disk coordinates of the  spot and a boolean visibility
        """
        v = TPXYZ(self.theta, self.phi)
        v = np.dot( self.rotation_star_to_earth() , v)
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

        raise Exception()

        return orbit, visi, phases

    def _plot_background(self):
        velocity = self.velocity_bins
        v0 = self.kwargs['v0']
        self.ax0.set_xlim(velocity[0]-v0, velocity[-1]-v0)
        res=self.movingpeak
        self.ax0.imshow(-res, cmap=plt.cm.gray_r,
                       aspect='auto',interpolation='none',
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

        #orbit, visi, phases = self.YZOrbit()

        #self.ax1.plot( orbit[1,visi], orbit[2, visi], '.k', markersize=0.2)
        #self.ax1.plot( [orbit[1,0]], [orbit[2,0]], 'ok', markersize=10)
        
        #print(self.likelihood_spot(self.theta, self.phi, self.phase))
        try:

            self.ax0.plot( self.line_velocities, np.mod( self.sampled_phases - self.phi, 2 * np.pi), '.r', markersize=1)

            self.ax1.plot( self.line_star_earth_phi[1,:], self.line_star_earth_phi[2,:], '.r', markersize=1)
            self.ax1.plot( self.line_star_earth_vel[1,:], self.line_star_earth_vel[2,:], '.y', markersize=1)

            self.ax0.plot( [self.clickvelocity], [self.clickphase], 'oy', markersize=5)

            x, y, z = self.XYZ()

            self.ax1.plot( [y], [z], '.y', markersize=10)

            # print(f'Spotposition (lat, lon) {180 * self.theta / np.pi},  {180 * self.phi / np.pi}')

            #print(self.likelihood_spot(self.theta, self.phi, self.phase))

        except:
            pass


    def plot_lh(self):

        theta, phi, res = self.lh_spot()

        map = Basemap(projection='moll', lon_0=0)
        plt.figure()
        #map.drawcoastlines()
        #map.imshow(res[:, ::-1].T, interpolation='none',vmin=-0.004,vmax=0.005)
        
        #### simultions showed that a bright line in the dynamic spectrum produces a bright
        # yellow spot in the map. A bright line in the dynamic spectrum corresponds to a dark spot
        # on the stellar surface. If we want to have bright yellow spot on map = bright spot on
        # star, we have to plot -res and not res!
        
        
        
        
        map.imshow(-res[:, ::-1].T, interpolation='none',vmin=-0.004,vmax=0.004)

        map.drawmeridians(np.linspace(-180, 180, 12))
        map.drawparallels(np.linspace(-90, 90, 6))
        #map.plot(180*phi/np.pi, 128*[180*self.kwargs['angle']/np.pi], '-b', latlon=True)
        map.plot(180*phi/np.pi, 128*[0], '-r', latlon=True)
        #map.plot(180*phi/np.pi, 128*[-180*self.kwargs['angle']/np.pi], '-b', latlon=True)
        plt.show()





    


if __name__ == '__main__':

    import matplotlib.pyplot as plt

    image = "moving_simple_per_night.txt"
    vbins = "velocitybins.txt"


    star = StarSpotFinder(
        image = image, #"test_spot.txt", "moving_simple_per_night.txt",
        vbins = "velocitybins.txt",
        v0=-13.01,
        vmax = 22,
        angle= 7 * np.pi/180,
        interactive = False
        )

    star.plot_lh()

    plt.show()
