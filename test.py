from DomeBuilder.MaterialLib import MaterialLib

from DomeBuilder.Frame import Frame

Rt = 21   # radius of tank wall
Ht =  8   # height of tank wall
Rd = 29.4 # dome radius
Nr =  8   # no. rings
Nn = 20   # no. nodes per ring
Ns = 20   # ring spacing coefficient
Nm = 30   # lattice coefficient

shell_material  = MaterialLib.steel()
shell_thickness = 0.007
snow_load       = 1000  # snow load [+ve] per square metre (projected)

print('Creating dome')
dome = Frame((Rt, Ht), Rd, Nr, Nn, (Ns, Nm))

s = 4
group = dome.get_group('central')
group['sec_index'] = 0
for r in range(1, Nr):
    s = s + 1
    group = dome.get_group('anticlockwise', r)
    group['sec_index'] = 0
    s = s + 1
    group = dome.get_group('clockwise', r)
    group['sec_index'] = 0
    s = s + 1
    group = dome.get_group('diagonal', r)
    group['sec_index'] = 0

dome.apply_loads('resources/Pressures_29_4.axdt', shell_material, shell_thickness, snow_load)

print('Setting up frame analysis')
mass, sigma, deflection = dome.analyse()
print("mass = {:.2f} tonnes; max stress = {:.2f} MPa; deflection = {:.2f} mm.".format(mass/1E3, sigma/1E6, deflection*1E3))

cost = mass / 1E3

dome.plot_frame(False, True, False, True)
Frame.syncPlotting(10)
