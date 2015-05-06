from gearbox.transmission.gears import *
from gearbox.standards.iso import Pitting as isoPitting
from gearbox.standards.iso import Bending as isoBending
from gearbox.standards.agma import Pitting as agmaPitting
from gearbox.standards.agma import Bending as agmaBending
from gearbox.export.export import *

module = 2.5  # m
helix_angle = 0.0  # beta
pressure_angle = 20.0  # alpha

lubricant = Lubricant(
    name='Kiruna',
    v40=160
)

material = Material(
    name='AISI 2010',
    classification='NV(nitrocar)',
    sh_limit=1500.,
    sf_limit=460.,
    e=206000.,
    poisson=0.3,
    density=7.83e-6,
    brinell=286.6667
)

tool = Tool(
    ha_p=1,
    hf_p=1.25,
    rho_fp=0.38,
    x=0,
    rho_ao=0,
    delta_ao=0,
    nc=10.
)

pinion = Gear(
    profile=tool,
    material=material,
    z=22.,
    beta=helix_angle,
    alpha=pressure_angle,
    m=module,
    x=0.525,
    b=34.0,
    bs=34.0,
    sr=0.0,
    rz=3.67,
    precision_grade=6.0,
    shaft_diameter=35.0,
    schema=3.0,
    l=60.0,
    s=15.0,
    backlash=0.017,
    gear_crown=1,
    helix_modification=1,
    favorable_contact=True,
    gear_condition=1
)

gear = Gear(
    profile=tool,
    material=material,
    z=40.,
    m=module,
    beta=helix_angle,
    alpha=pressure_angle,
    x=-0.275,
    b=34.0,
    bs=34.0,
    sr=0.0,
    rz=3.67,
    precision_grade=6.0,
    shaft_diameter=50.0,
    schema=3.0,
    l=60.0,
    s=35.0,
    backlash=-0.017,
    gear_crown=1,
    helix_modification=1,
    favorable_contact=True,
    gear_condition=1
)

pair = [pinion, gear]

transmission = Transmission(
    gears=pair,
    lubricant=lubricant,
    rpm_in=1450.0,
    rpm_out=243.5,
    n=40.0,
    l=10000.0,
    gear_box_type=2,
    ka=1.3,
    sh_min=1,
    sf_min=1
)

print '========================================'
print 'ISO Pitting'
print isoPitting(transmission=transmission).calculate()
print '========================================'

print '========================================'
print 'ISO Bending'
print isoBending(transmission=transmission).calculate
print '========================================'

print '========================================'
print 'AGMA Pitting'
print agmaPitting(transmission=transmission).calculate()
print '========================================'

print '========================================'
print 'AGMA Bending'
print agmaBending(transmission=transmission).calculate()
print '========================================'

output_folder = os.path.join(os.path.dirname(__file__), 'output')
try:
    os.mkdir(output_folder)
except:
    pass

# ===============================================
# EXPORT TO MATLAB/COMSOL SCRIPT 
# ===============================================
# 2D model matlab-comsol export
# for 2D export type='2D' is optional because '2D' is the default output
comsol_transmission_model_name2d = 'transmission2d'
comsol_pinion_model_name2d = 'pinion2d'
comsol_wheel_model_name2d = 'wheel2d'
ExportGear(pinion).matlab_comsol(model_name=comsol_pinion_model_name2d, output_folder=output_folder, type='2D')
ExportGear(gear).matlab_comsol(model_name=comsol_wheel_model_name2d, output_folder=output_folder, type='2D')
ExportPair(pair).matlab_comsol(model_name=comsol_transmission_model_name2d, output_folder=output_folder, type='2D')

# 3D model matlab-comsol export
comsol_transmission_model_name3d = 'transmission3d'
comsol_pinion_model_name3d = 'pinion3d'
comsol_wheel_model_name3d = 'wheel3d'
ExportGear(pinion).matlab_comsol(model_name=comsol_pinion_model_name3d, output_folder=output_folder, type='3D')
ExportGear(gear).matlab_comsol(model_name=comsol_wheel_model_name3d, output_folder=output_folder, type='3D')
ExportPair(pair).matlab_comsol(model_name=comsol_transmission_model_name3d, output_folder=output_folder, type='3D')

# ===============================================
# EXPORT TO ABAQUS PYTHON 
# ===============================================
# 2D model abaqus export
# for 2D export type='2D' is optional because '2D' is the default output
abaqus_transmission_model_name2d = 'transmission2d'
abaqus_pinion_model_name2d = 'pinion2d'
abaqus_wheel_model_name2d = 'wheel2d'
ExportGear(pinion).abaqus(model_name=abaqus_pinion_model_name2d, output_folder=output_folder, type='2D')
ExportGear(gear).abaqus(model_name=abaqus_wheel_model_name2d, output_folder=output_folder, type='2D')
ExportPair(pair).abaqus(model_name=abaqus_transmission_model_name2d, output_folder=output_folder, type='2D')

# 3D model abaqus export
abaqus_transmission_model_name3d = 'transmission3d'
abaqus_pinion_model_name3d = 'pinion3d'
abaqus_wheel_model_name3d = 'wheel3d'
ExportGear(pinion).abaqus(model_name=abaqus_pinion_model_name3d, output_folder=output_folder, type='3D')
ExportGear(gear).abaqus(model_name=abaqus_wheel_model_name3d, output_folder=output_folder, type='3D')
ExportPair(pair).abaqus(model_name=abaqus_transmission_model_name3d, output_folder=output_folder, type='3D')

# ===============================================
# EXPORT TO ANSYS/WOKRBENCH Script
# ===============================================
# 2D model abaqus export
# for 2D export type='2D' is optional because '2D' is the default output
ansys_transmission_model_name2d = 'transmission2d'
ansys_pinion_model_name2d = 'pinion2d'
ansys_wheel_model_name2d = 'wheel2d'
ExportGear(pinion).ansys(model_name=ansys_pinion_model_name2d, output_folder=output_folder, type='2D')
ExportGear(gear).ansys(model_name=ansys_wheel_model_name2d, output_folder=output_folder, type='2D')
ExportPair(pair).ansys(model_name=ansys_transmission_model_name2d, output_folder=output_folder, type='2D')
#
# 3D model ansys export
ansys_transmission_model_name3d = 'transmission3d'
ansys_pinion_model_name3d = 'pinion3d'
ansys_wheel_model_name3d = 'wheel3d'
ExportGear(pinion).ansys(model_name=ansys_pinion_model_name3d, output_folder=output_folder, type='3D')
ExportGear(gear).ansys(model_name=ansys_wheel_model_name3d, output_folder=output_folder, type='3D')
ExportPair(pair).ansys(model_name=ansys_transmission_model_name3d, output_folder=output_folder, type='3D')
