import argparse
import csv
import sys

import numpy as np

from BeesEtAl.BA_Garden   import BA_Garden
from BeesEtAl.Base_Coster import Base_Coster

from DomeBuilder.MaterialLib import MaterialLib
from DomeBuilder.Frame       import Frame

# Defaults

Nr_min_max = ( 3,  8)
Nn_min_max = ( 5, 20)
Ns_min_max = (10, 30)
Nm_min_max = (10, 40)

bMESO         = True
bShowStress   = False
bVerbose      = False
maxsol_runs   = 1000
priorities    = [5,2,2,1]
prior_init    = None
fixedRings    = None
bSelectedOnly = False
bLoadCheats   = False
sourceFile    = None

parser = argparse.ArgumentParser(description="DomeBuilder dome optimisation tool.")

parser.add_argument('-v', '--verbose', help='Print additional commentary.', action='store_true')
parser.add_argument('--show-stress',   help='Plot each frame, showing stresses.', action='store_true')
parser.add_argument('--no-meso',       help='Use plain Bees Algorithm without MESO.', action='store_true')
parser.add_argument('--refinement',    help='Run a study with 72 elite patches, 5 bees each.', action='store_true')
parser.add_argument('--fix-rings',     help='Fix the number of node-rings.', dest='rings', default=0, type=int)
parser.add_argument('--pareto-out',    help='Specify output file name for Pareto-optimal set [pareto.csv].', type=str, default='pareto.csv')
parser.add_argument('--results-out',   help='Specify output file name for results [results.csv].', type=str, default='results.csv')
parser.add_argument('--selected-in',   help='Specify input file name: do not run the optimiser, just evaluate locations.', type=str)
parser.add_argument('--sources-in',    help='Specify input file name: pre-load scout locations.', type=str)

args = parser.parse_args()

resultsFile = args.results_out
paretoFile  = args.pareto_out

if args.verbose:
    bVerbose = True

if args.show_stress:
    bShowStress = True

if args.selected_in is not None:
    bSelectedOnly = True
    sourceFile    = args.selected_in
    print('Evaluating locations from: {f}'.format(f=sourceFile))
else:
    if args.no_meso:
        bMESO = False
        print('MESO disabled.')
    else:
        print('MESO enabled.')

    if args.rings > 0:
        fixedRings = args.rings

    if args.sources_in is not None:
        bLoadCheats = True
        sourceFile  = args.sources_in

    if args.refinement:
        if args.sources_in is None:
            print('A refinement study should be used with preloaded scout locations.')
            sys.exit()
        maxsol_runs = 3600
        priorities  = [5]*72 + [1]
        prior_init  = [1]*73
        print('Running a refinement study.')

class DomeCoster(Base_Coster):

    @staticmethod
    def create_garden(priorities, Nr_min_max, Nn_min_max, Ns_min_max, Nm_min_max):
        Nr_min, Nr_max = Nr_min_max
        Nn_min, Nn_max = Nn_min_max
        Ns_min, Ns_max = Ns_min_max
        Nm_min, Nm_max = Nm_min_max

        Nsec = (Nr_max - 1) * 3 + 1

        par = np.zeros((2, 4+Nsec))

        par[0,0:4] = [Nr_min, Nn_min, Ns_min, Nm_min]
        par[1,0:4] = [Nr_max, Nn_max, Ns_max, Nm_max]

        par[1,4:(4+Nsec)] = np.ones(Nsec) * 78

        G = BA_Garden(par[0,:], par[1,:], priorities)

        G.costfn = DomeCoster(G, par)

        return G

    def __init__(self, garden, parameter_matrix):
        Base_Coster.__init__(self, garden)

        self.par = parameter_matrix

        self.Ht =  8   # tank-wall height
        self.Rt = 21   # tank-wall radius
        self.Rd = 29.4 # dome radius

        self.wind_file_name  = 'resources/Pressures_29_4.axdt' # or None to skip
        self.shell_material  = MaterialLib.steel()
        self.shell_thickness = 0.007
        self.snow_load       = 1000  # snow load [+ve] per square metre (projected)

        self.dome = None

        if bVerbose:
            self.verbose = True # print patch info

        with open(resultsFile, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            minima = self.par[0,:]
            maxima = self.par[1,:]
            writer.writerow(['','','','','',*minima])
            writer.writerow(['','','','','',*maxima])

    def map_to_solution_space(self, X):
        # Let's restrict solution space to discrete values, and choose the nearest
        XA = np.copy(X)
        XA[0] = np.around(X[0], decimals=0)
        XA[1] = np.around(X[1], decimals=0)

        for s in range(4, len(X)):
            XA[s] = np.around(X[s], decimals=0)

        return XA

    def __determine_penalty(self, value, lower_limit, upper_limit, max_penalty=1000):
        penalty = 0

        if value > lower_limit:
            if value > upper_limit:
                penalty = 1
            else:
                penalty = (value - lower_limit) / (upper_limit - lower_limit)

        self.cost = self.cost + penalty * max_penalty

    def evaluate_cost(self):
        Nr = int(self.XA[0])   # no. rings
        Nn = int(self.XA[1])   # no. nodes per ring
        Ns = self.XA[2]        # ring spacing coefficient
        Nm = self.XA[3]        # lattice coefficient
        Nc = 0                 # number of circumferential rings

        Ng = 3 * (Nr - 1) + 1  # no. of distinct section-groups, so Nr=20 => 58 section-groups
        pg = 10000 / Ng        # adjusted maximum penalty

        if self.verbose:
            print('Creating dome')

        self.dome = Frame((self.Rt, self.Ht), self.Rd, Nr, Nn, (Ns, Nm), self.verbose)

        s = 4
        group = self.dome.get_group('central')
        group['sec_index'] = int(self.XA[s])
        for r in range(1, Nr):
            s = s + 1
            group = self.dome.get_group('anticlockwise', r)
            group['sec_index'] = int(self.XA[s])
            s = s + 1
            group = self.dome.get_group('clockwise', r)
            group['sec_index'] = int(self.XA[s])
            s = s + 1
            group = self.dome.get_group('diagonal', r)
            group['sec_index'] = int(self.XA[s])

        self.dome.apply_loads(self.wind_file_name, self.shell_material, self.shell_thickness, self.snow_load)

        if self.verbose:
            print('Setting up frame analysis')

        mass, sigma, deflection = self.dome.analyse()

        if self.verbose:
            print("mass = {:.2f} tonnes; max stress = {:.2f} MPa; deflection = {:.2f} mm.".format(mass/1E3, sigma/1E6, deflection*1E3))

        self.cost = mass / 1E3

        self.__determine_penalty(deflection, 0.0375, 0.05, pg)

        group = self.dome.get_group('central')
        self.__determine_penalty(group['max_stress'], 0.75, 1, pg)
        self.__determine_penalty(group['buckling'],   0.75, 1, pg)
        for r in range(1, Nr):
            group = self.dome.get_group('anticlockwise', r)
            self.__determine_penalty(group['max_stress'], 0.75, 1, pg)
            self.__determine_penalty(group['buckling'],   0.75, 1, pg)
            group = self.dome.get_group('clockwise', r)
            self.__determine_penalty(group['max_stress'], 0.75, 1, pg)
            self.__determine_penalty(group['buckling'],   0.75, 1, pg)
            group = self.dome.get_group('diagonal', r)
            self.__determine_penalty(group['max_stress'], 0.75, 1, pg)
            self.__determine_penalty(group['buckling'],   0.75, 1, pg)
            if group['orientation'] == 'circumferential':
                Nc = Nc + 1

        if self.verbose:
            print("cost = {:.2f} <~~~~~~~~~~~~~~~~".format(self.cost))

        if bShowStress:
            self.dome.plot_frame(False, True, False, True)

        with open(resultsFile, 'a') as csvfile:
            writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer.writerow([mass, sigma, deflection, self.cost, Nc, *self.XA])

        # let's make this multiobjective
        self.cost = [self.cost, (deflection * 1E3)]

    def __meso_check_group(self, s, group):
        secindex = group['sec_index']
        if secindex > 0:
            if group['max_stress'] > 0.75 or group['buckling'] > 0.75:
                self.XM[s] = secindex - 1 # suggest an increase in the second moment of area of the section
                return

        if secindex < 78:
            if self.strongest_group is None:
                self.strongest_group = group
                self.strongest_index = s
            elif group['max_stress'] < self.strongest_group['max_stress']:
                self.strongest_group = group
                self.strongest_index = s

    def meso(self):
        if bMESO:
            Nr = int(self.XA[0])   # no. rings - no. of distinct section-groups = 3 * (Nr - 1) + 1, so Nr=20 => 58 section-groups

            self.strongest_group = None
            self.strongest_index = None

            s = 4
            self.__meso_check_group(s, self.dome.get_group('central'))
            for r in range(1, Nr):
                s = s + 1
                self.__meso_check_group(s, self.dome.get_group('anticlockwise', r))
                s = s + 1
                self.__meso_check_group(s, self.dome.get_group('clockwise', r))
                s = s + 1
                self.__meso_check_group(s, self.dome.get_group('diagonal', r))

            if self.strongest_group is not None:
                secindex = self.strongest_group['sec_index']
                self.XM[self.strongest_index] = secindex + 1 # suggest a decrease in the second moment of area of the section

    def load_and_run(self, file_name):
        with open(file_name, newline='') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                values = list(map(float, row))
                X = np.asarray(values)
                cost = self.calculate_cost(X)

G = DomeCoster.create_garden(priorities, Nr_min_max, Nn_min_max, Ns_min_max, Nm_min_max)

if bSelectedOnly:
    G.costfn.load_and_run(sourceFile)
    G.pareto(paretoFile)

    print('==== Finished - quitting in 60 seconds ====')
    plt.pause(60)
else:
    method = 'gauss'
    Nfail  = 6      # i.e., stops at 6th failure
    rf     = 1 / 78 # distance between CHS sections in unit space
    r0, sf = G.initial_radius_and_shrinking(Nfail, rf, method)
    params = { 'radius': r0, 'shrink': sf, 'fail_at': Nfail, 'neighborhood': method, 'dynamic': True }

    G.set_search_params(**params)

    if fixedRings is not None:
        # let's fix the number of rings
        mask = np.ones(G.Ndim)
        defs = np.zeros(G.Ndim)
        mask[0] = 0
        defs[0] = fixedRings
        G.set_mask_and_defaults(mask, defs)
        print('Fixing number of rings = {n}.'.format(n=fixedRings))

    if sourceFile is not None:
        Ncheats = G.preload(sourceFile)
        print('Preloaded {n} solutions from {f}.'.format(n=Ncheats, f=sourceFile))

    solver_runs = 0
    it = 0
    while solver_runs < maxsol_runs:
        it = it + 1

        if it == 1 and prior_init is not None:
            solver_runs = G.iterate(maxsol_runs, override=prior_init)
        else:
            solver_runs = G.iterate(maxsol_runs)

        best_cost, best_X = G.best()

        print('Iteration {:4d}: {s}/{t} solver runs'.format(it, s=solver_runs, t=maxsol_runs))
        print('Global best = ' + str(best_cost) + ' { ' + str(best_X) + ' }')

    G.pareto(paretoFile)
