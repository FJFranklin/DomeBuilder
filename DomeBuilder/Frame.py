import csv

import numpy as np

from PyNite import FEModel3D

from .MaterialLib import MaterialLib
from .SectionList import SectionList

class Frame(object):

    __fig = None
    __mpl = None
    __plt = None
    __p3d = None

    @staticmethod
    def figForPlotting():
        if Frame.__fig is None:
            from mpl_toolkits.mplot3d import Axes3D

            import matplotlib as mpl
            import matplotlib.pyplot as plt

            Frame.__mpl = mpl
            Frame.__plt = plt

            from mpl_toolkits.mplot3d.art3d import Poly3DCollection

            Frame.__p3d = Poly3DCollection

            xsize = 1500
            ysize = 1500
            dpi_osx = 192 # Something very illogical here.

            Frame.__fig = plt.figure(figsize=(xsize / dpi_osx, ysize / dpi_osx), dpi=(dpi_osx/2))

            Frame.__plt.ion()
            Frame.__plt.show()
        else:
            Frame.__plt.clf()

        return Frame.__fig

    @staticmethod
    def syncPlotting(dt=0.000001):
        if Frame.__fig is not None:
            Frame.__plt.draw()
            Frame.__plt.pause(dt)

    def __init__(self, tank_dimensions, dome_radius, rings, nodes_per_ring, coefficients, verbose=True):
        r, h = tank_dimensions
        self.Rt = r # radius of tank wall
        self.Ht = h # height of tank wall

        self.Nr = rings                  # number of rings of nodes
        self.Nn = nodes_per_ring         # number of nodes per ring
        self.Nt = 1 + self.Nr * self.Nn  # total number of nodes

        self.node = np.zeros((self.Nt, 6)) # x, y, z, Fx, Fy, Fz

        self.radii = np.zeros(self.Nr + 1) # includes ring 0, which has radius 0
        self.chord = np.zeros(self.Nr + 1) # includes ring 0, but meaningless
        self.ringz = np.zeros(self.Nr + 1) # includes ring 0, which is the top of the dome

        self.edge = [] # dict with name, associated nodes & other properties
        self.face = [] # dict with name, associated nodes, members & other properties

        self.dome_origin = None
        self.dome_area   = None
        self.dome_mass   = 0 # total structural mass

        self.CHS = SectionList.load_CHS()

        self.group = {}

        self.verbose = verbose

        # Central set of girders
        gname = self.__create_name('GX', 0)
        gdata = {}
        gdata['sec_type']   = self.CHS
        gdata['sec_index']  = 0
        gdata['material']   = MaterialLib.steel()
        gdata['max_stress'] = 0
        gdata['buckling']   = 0
        self.group[gname] = gdata

        for r in range(1, self.Nr):
            # Anticlockwise girders
            gname = self.__create_name('GA', r)
            gdata = {}
            gdata['sec_type']   = self.CHS
            gdata['sec_index']  = 0
            gdata['material']   = MaterialLib.steel()
            gdata['max_stress'] = 0
            gdata['buckling']   = 0
            self.group[gname] = gdata

            # Clockwise girders
            gname = self.__create_name('GC', r)
            gdata = {}
            gdata['sec_type']   = self.CHS
            gdata['sec_index']  = 0
            gdata['material']   = MaterialLib.steel()
            gdata['max_stress'] = 0
            gdata['buckling']   = 0
            self.group[gname] = gdata

            # Diagonal (radial/circumferential) girders
            gname = self.__create_name('GD', r)
            gdata = {}
            gdata['sec_type']   = self.CHS
            gdata['sec_index']  = 0
            gdata['material']   = MaterialLib.steel()
            gdata['max_stress'] = 0
            gdata['buckling']   = 0
            self.group[gname] = gdata

        self.frame = None
        self.max_stress = 0

        Ns, Nm = coefficients
        self.__position_rings(Ns)
        self.__determine_coordinates(dome_radius)
        self.__create_edges_and_faces(Nm)

    def get_group(self, region, ring=0): # region is one of 'anticlockwise', 'clockwise', 'central' and 'diagonal'
        group = None
        if region == 'central' and ring == 0:
            group = self.group[self.__create_name('GX', 0)]
        elif ring > 0 and ring < self.Nr:
            if region == 'anticlockwise':
                group = self.group[self.__create_name('GA', ring)]
            elif region == 'clockwise':
                group = self.group[self.__create_name('GC', ring)]
            elif region == 'diagonal':
                group = self.group[self.__create_name('GD', ring)]
        return group
        
    def RNI(self, r, ni): # wrap node & index
        if r == 0:
            no = 0
            index = 0
        else:
            if ni < 1:
                no = ni + self.Nn
            elif ni > self.Nn:
                no = ni - self.Nn
            else:
                no = ni
            index = (r - 1) * self.Nn + no
        return (r, no, index)

    def __position_rings(self, Ns): # Ns is a spacing coefficient
        self.radii[0] = 0
        self.chord[0] = 0

        for r in range(1, self.Nr + 1):
            radius = self.Rt * np.power(r / self.Nr, 20 / Ns)

            self.radii[r] = radius
            self.chord[r] = radius * np.sin(np.pi / self.Nn) * 2

    def __determine_ring_heights(self, Rd):
        rr = (Rd**2 - self.Rt**2)**0.5
        z0 = self.Ht - rr

        self.dome_origin = np.asarray([0,0,z0])
        self.dome_area   = 2 * np.pi * Rd * (Rd - rr)

        for r in range(0, self.Nr + 1):
            self.ringz[r] = z0 + (Rd**2 - self.radii[r]**2)**0.5

    def __nearest(self, x, y):
        radius = (x**2 + y**2)**0.5
        r_near = self.Nr
        for r in range(0, self.Nr):
            if radius < ((self.radii[r] + self.radii[r+1]) / 2):
                r_near = r
                break
        r = r_near

        theta = np.arctan2(y, x) * self.Nn / np.pi
        if r & 0x1:
            n = int(np.around(theta / 2))
        else:
            n = int(np.around((theta - 1) / 2))

        return self.RNI(r, n)

    def __determine_coordinate(self, rni):
        r, n, index = rni

        radius = self.radii[r]
        z      = self.ringz[r]

        if r & 0x1:
            theta = 2 * n * np.pi / self.Nn
        else:
            theta = (2 * n + 1) * np.pi / self.Nn

        x = radius * np.cos(theta)
        y = radius * np.sin(theta)

        self.node[index,0:3] = [x, y, z]

    def __determine_coordinates(self, dome_radius):
        self.__determine_ring_heights(dome_radius)

        for r in range(0, self.Nr + 1):
            if r == 0:
                self.__determine_coordinate(self.RNI(r, 0))
            else:
                for n in range(1, self.Nn + 1):
                    self.__determine_coordinate(self.RNI(r, n))

    def __create_name(self, prefix, value):
        if isinstance(value, tuple):
            v1, v2 = value
            if v2 < v1:
                v2, v1 = value
            name = '{p}_{a}_{b}'.format(p=prefix, a=v1, b=v2)
        else:
            name = '{p}_{a}'.format(p=prefix, a=value)
        return name

    def __create_edge(self, gname, rn0, rn1):
        r0, n0, i0 = rn0
        r1, n1, i1 = rn1
        mname = self.__create_name('M', (i0, i1))

        data = {}
        data['gname'] = gname
        data['mname'] = mname
        data['node0'] = rn0
        data['node1'] = rn1
        self.edge.append(data)

    def __create_face(self, rn0, rn1, rn2):
        r0, n0, i0 = rn0
        r1, n1, i1 = rn1
        r2, n2, i2 = rn2

        e10 = np.asarray([self.node[i0,0]-self.node[i1,0],self.node[i0,1]-self.node[i1,1],self.node[i0,2]-self.node[i1,2]])
        e12 = np.asarray([self.node[i2,0]-self.node[i1,0],self.node[i2,1]-self.node[i1,1],self.node[i2,2]-self.node[i1,2]])

        vec_cross = np.cross(e10, e12)
        magnitude = np.linalg.norm(vec_cross)

        data = {}
        data['node0'] = rn0
        data['node1'] = rn1
        data['node2'] = rn2
        data['norm']  = vec_cross / magnitude
        data['area']  = magnitude / 2
        self.face.append(data)

    def __create_edges_and_faces(self, Nm): # Nm is a coefficient that determines radial vs circumferential diagonals
        for r in range(0, self.Nr):
            if r == 0:
                gname = self.__create_name('GX', r)
                for n in range(1, self.Nn + 1):
                    self.__create_edge(gname, self.RNI(0, 0), self.RNI(1, n))
            else:
                cname = self.__create_name('GC', r)
                aname = self.__create_name('GA', r)
                dname = self.__create_name('GD', r)

                chord  = self.chord[r]
                radial = ((self.radii[r+1] - self.radii[r-1])**2 + (self.ringz[r+1] - self.ringz[r-1])**2)**0.5

                group = self.group[dname]

                if (radial / chord) > (Nm / 20):
                    bRadial = False
                    group['orientation'] = 'circumferential'
                else: # select radial
                    bRadial = True
                    group['orientation'] = 'radial'

                for n in range(1, self.Nn + 1):
                    if bRadial:
                        self.__create_edge(dname, self.RNI(r-1, n), self.RNI(r+1, n))
                        if r & 0x1:
                            self.__create_face(self.RNI(r-1, n), self.RNI(r+1, n), self.RNI(r, n))
                            self.__create_face(self.RNI(r+1, n), self.RNI(r-1, n), self.RNI(r, n+1))
                        else:
                            self.__create_face(self.RNI(r-1, n), self.RNI(r+1, n), self.RNI(r, n-1))
                            self.__create_face(self.RNI(r+1, n), self.RNI(r-1, n), self.RNI(r, n))
                    else:
                        self.__create_edge(dname, self.RNI(r, n), self.RNI(r, n+1))
                        if r & 0x1:
                            self.__create_face(self.RNI(r, n+1), self.RNI(r, n),   self.RNI(r-1, n))
                            self.__create_face(self.RNI(r, n),   self.RNI(r, n+1), self.RNI(r+1, n))
                        else:
                            self.__create_face(self.RNI(r, n+1), self.RNI(r, n),   self.RNI(r-1, n+1))
                            self.__create_face(self.RNI(r, n),   self.RNI(r, n+1), self.RNI(r+1, n+1))

                    if r & 0x1:
                        self.__create_edge(cname, self.RNI(r, n), self.RNI(r+1, n-1))
                        self.__create_edge(aname, self.RNI(r, n), self.RNI(r+1, n))
                    else:
                        self.__create_edge(cname, self.RNI(r, n), self.RNI(r+1, n))
                        self.__create_edge(aname, self.RNI(r, n), self.RNI(r+1, n+1))

    def __add_frame_dead_load(self):
        for e in self.edge:
            r0, n0, i0 = e['node0']
            r1, n1, i1 = e['node1']
            e10 = np.asarray([self.node[i0,0]-self.node[i1,0],self.node[i0,1]-self.node[i1,1],self.node[i0,2]-self.node[i1,2]])

            group = self.group[e['gname']]
            sectype  = group['sec_type']
            secindex = group['sec_index']
            material = group['material']

            length = np.linalg.norm(e10)
            area = sectype.get(secindex, 'area')
            density = material['density']

            mass = length * area * density
            weight = -9.81 * mass

            self.dome_mass = self.dome_mass + mass

            self.node[i0,5] = self.node[i0,5] + weight / 2
            self.node[i1,5] = self.node[i1,5] + weight / 2

    def __add_shell_dead_load(self, material, shell_thickness, snow_load): # snow load [+ve] per square metre (projected)
        mass_per_sqm = material['density'] * shell_thickness # shell weight per square metre

        for f in self.face:
            r0, n0, i0 = f['node0']
            r1, n1, i1 = f['node1']
            r2, n2, i2 = f['node2']

            norm = f['norm']
            area = f['area']

            mass = mass_per_sqm * area
            weight_total = -9.81 * mass - snow_load * area * norm[2]

            self.dome_mass = self.dome_mass + mass

            self.node[i0,5] = self.node[i0,5] + weight_total / 3
            self.node[i1,5] = self.node[i1,5] + weight_total / 3
            self.node[i2,5] = self.node[i2,5] + weight_total / 3

    def __apply_wind_pressure(self, file_name):
        with open(file_name, newline='') as csvfile:
            reader = csv.reader(csvfile)
            bHeader = True
            count = 0
            pdata = np.zeros(4)
            v     = np.zeros(3)
            for row in reader:
                if len(row) < 1:
                    break
                if len(row) < 2:
                    continue
                if bHeader:
                    bHeader = False
                    continue

                # x, z, y, pressure
                pdata[:] = list(map(float, row))

                radius = (pdata[0]**2 + pdata[2]**2)**0.5
                if radius > self.Rt - 0.01:
                    continue # tank wall pressure - skip

                # find the nearest node
                r, n, index = self.__nearest(pdata[0], pdata[2])

                # the direction the pressure is acting in:
                v[0] = self.dome_origin[0] - pdata[0]
                v[1] = self.dome_origin[1] - pdata[2]
                v[2] = self.dome_origin[2] - pdata[1]
                # add to nearest node, with magnitude equal to the pressure
                self.node[index,3:6] = self.node[index,3:6] + pdata[3] * v / np.linalg.norm(v)

                count = count + 1

            # multiply pressures times area, assuming each pressure point has equal area
            self.node[:,3:6] = self.node[:,3:6] * self.dome_area / count

    def apply_loads(self, wind_file_name, shell_material, shell_thickness, snow_load):
        if wind_file_name is not None:
            if self.verbose:
                print('Applying wind load to dome')
            self.__apply_wind_pressure(wind_file_name)

        if self.verbose:
            print('Adding frame member weights')
        self.__add_frame_dead_load()

        if self.verbose:
            print('Adding shell weight and snow load')
        self.__add_shell_dead_load(shell_material, shell_thickness, snow_load)

    def __bending_stress(self, member, material, sectype, secindex, moments=None): # ~~~~> TODO: Think about how to do this best... + check! <~~~~
        if moments is not None:
            My, Mz = moments
        else:
            My = max([abs(member.MinMoment('My')), abs(member.MaxMoment('My'))])
            Mz = max([abs(member.MinMoment('Mz')), abs(member.MaxMoment('Mz'))])

        buckling = 0

        FM = member.f()
        Fx = FM[0,0]
        Mx = FM[3,0]

        A    = sectype.get(secindex, 'area')
        Iyy  = sectype.get(secindex, 'Iyy')
        Izz  = sectype.get(secindex, 'Izz')
        J    = sectype.get(secindex, 'J')

        d_yy = sectype.get(secindex, 'd-yy')
        d_zz = sectype.get(secindex, 'd-zz')

        # axial stress
        sigma_yy = (d_yy / 2) * (My / Iyy) + abs(Fx) / A
        sigma_zz = (d_zz / 2) * (Mz / Izz) + abs(Fx) / A

        # torsion only - shear not yet included
        tau_yy = (d_yy / 2) * (Mx / J)
        tau_zz = (d_zz / 2) * (Mx / J)

        # von Mises stress
        stress_yy = (sigma_yy**2 + 3 * tau_yy**2)**0.5
        stress_zz = (sigma_zz**2 + 3 * tau_zz**2)**0.5

        if Fx < 0: # member is in compression
            E = material['Young']
            I = min([Iyy, Izz])

            # Critical force for bucklng
            Fc = 4 * E * I * (np.pi / member.L)**2

            # Buckling coefficient is then just the ratio
            buckling = -Fx / Fc

        return max([stress_yy, stress_zz]), buckling

    def __maximum_bending_stress(self):
        self.max_stress = 0

        for e in self.edge:
            member = self.frame.GetMember(e['mname'])

            group = self.group[e['gname']]
            sectype  = group['sec_type']
            secindex = group['sec_index']
            material = group['material']

            sigma, buckling = self.__bending_stress(member, material, sectype, secindex)

            if self.max_stress < sigma:
                self.max_stress = sigma

            # normalise w.r.t. material yield stress, and record group-maximum
            sigma = sigma / material['yield']

            if group['max_stress'] < sigma:
                group['max_stress'] = sigma
            if group['buckling'] < buckling:
                group['buckling'] = buckling

        return self.max_stress

    def plot_frame(self, bShowNodes, bShowFrame, bShowFaces, bShowResults): # ==== Plot the frame ====
        fig = Frame.figForPlotting()

        ax = fig.add_subplot(111, projection='3d')
        ax.set_position([0.07, 0.06, 0.90, 0.90])
        ax.set_facecolor('white')

        ax.set_xlim(-22, 22)
        ax.set_ylim(-22, 22)

        if bShowNodes:
            ax.scatter(self.node[:,0], self.node[:,1], self.node[:,2], marker='.', color='blue')

        if bShowFrame:
            for e in self.edge:
                r0, n0, i0 = e['node0']
                r1, n1, i1 = e['node1']

                x = [self.node[i0,0],self.node[i1,0]]
                y = [self.node[i0,1],self.node[i1,1]]
                z = [self.node[i0,2],self.node[i1,2]]

                ax.plot(x, y, z, color='black')

        if bShowFaces:
            for f in self.face:
                r0, n0, i0 = f['node0']
                r1, n1, i1 = f['node1']
                r2, n2, i2 = f['node2']

                x = [self.node[i0,0],self.node[i1,0],self.node[i2,0]]
                y = [self.node[i0,1],self.node[i1,1],self.node[i2,1]]
                z = [self.node[i0,2],self.node[i1,2],self.node[i2,2]]

                verts = [list(zip(x, y, z))]
                ax.add_collection3d(Frame.__p3d(verts), zs='z')

        if bShowResults:
            max_stress = self.max_stress / 1E6

            cax, _ = Frame.__mpl.colorbar.make_axes(Frame.__plt.gca(), shrink=0.5)
            cbar   = Frame.__mpl.colorbar.ColorbarBase(cax, cmap=Frame.__mpl.cm.coolwarm, norm=Frame.__mpl.colors.Normalize(vmin=0, vmax=max_stress))
            cbar.set_clim(0, max_stress)

            data = np.zeros((20,3))

            for e in self.edge:
                r0, n0, i0 = e['node0']
                r1, n1, i1 = e['node1']

                N0 = self.frame.GetNode(self.__create_name('N', i0))
                N1 = self.frame.GetNode(self.__create_name('N', i1))

                scale = 10

                xyz0 = self.node[i0,0:3] + np.asarray([N0.DX, N0.DY, N0.DZ]) * scale
                xyz1 = self.node[i1,0:3] + np.asarray([N1.DX, N1.DY, N1.DZ]) * scale

                member = self.frame.GetMember(e['mname'])

                group = self.group[e['gname']]
                sectype  = group['sec_type']
                secindex = group['sec_index']
                material = group['material']

                for ix in range(0, 20):
                    # initial beam positions along length
                    data[ix,0] = (xyz0[0] * (19 - ix) + xyz1[0] * ix) / 19
                    data[ix,1] = (xyz0[1] * (19 - ix) + xyz1[1] * ix) / 19
                    data[ix,2] = (xyz0[2] * (19 - ix) + xyz1[2] * ix) / 19
                for ix in range(0, 19):
                    mx = (ix + 0.5 ) * member.L / 19 # member-local x-coordinate
                    My = member.Moment('My', mx)
                    Mz = member.Moment('Mz', mx)
                    sigma, buckling = self.__bending_stress(member, material, sectype, secindex, (My, Mz))
                    sigma = sigma / 1E6

                    color = Frame.__mpl.cm.coolwarm(sigma / max_stress)

                    x = [data[ix,0],data[ix+1,0]]
                    y = [data[ix,1],data[ix+1,1]]
                    z = [data[ix,2],data[ix+1,2]]
                    ax.plot(x, y, z, color=color)

        Frame.syncPlotting()

    def analyse(self):
        self.frame = FEModel3D()

        # Add nodes
        for i in range(0, self.Nt):
            nname = self.__create_name('N', i)
            self.frame.AddNode(nname, self.node[i,0], self.node[i,1], self.node[i,2])
            self.frame.AddNodeLoad(nname, 'FX', self.node[i,3])
            self.frame.AddNodeLoad(nname, 'FY', self.node[i,4])
            self.frame.AddNodeLoad(nname, 'FZ', self.node[i,5])

        # Fix outer ring of nodes
        for n in range(1, self.Nn + 1):
            r, n, index = self.RNI(self.Nr, n)
            nname = self.__create_name('N', index)
            self.frame.DefineSupport(nname, True, True, True, True, True, True)

        for e in self.edge:
            r0, n0, i0 = e['node0']
            r1, n1, i1 = e['node1']

            mname  = e['mname']
            nname0 = self.__create_name('N', i0)
            nname1 = self.__create_name('N', i1)

            group = self.group[e['gname']]
            sectype  = group['sec_type']
            secindex = group['sec_index']
            material = group['material']

            area = sectype.get(secindex, 'area')
            Iyy  = sectype.get(secindex, 'Iyy')
            Izz  = sectype.get(secindex, 'Izz')
            J    = sectype.get(secindex, 'J')

            E = material['Young']
            G = material['Shear']

            self.frame.AddMember(mname, nname0, nname1, E, G, Iyy, Izz, J, area, self.dome_origin)

        if self.verbose:
            print('Setup complete. Analysing...')

        self.frame.Analyze()

        if self.verbose:
            print('Done.')

        max_deflection = 0
        for i in range(0, self.Nt):
            Ni = self.frame.GetNode(self.__create_name('N', i))
            deflection = np.linalg.norm(np.asarray([Ni.DX, Ni.DY, Ni.DZ]))
            if max_deflection < deflection:
                max_deflection = deflection

        return self.dome_mass, self.__maximum_bending_stress(), max_deflection
