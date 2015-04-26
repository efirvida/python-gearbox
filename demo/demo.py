from gearbox.transmition.gears import *
from gearbox.standards.iso import Pitting as iso_pitting
from gearbox.standards.iso import Bending as iso_bending
from gearbox.standards.agma import Pitting as agma_pitting
from gearbox.standards.agma import Bending as agma_bending
from gearbox.export.export import *


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
    beta=16.0,
    alpha=20.0,
    m=2.5,
    x=0.0,
    b=34.0,
    bs=34.0,
    sr=0.0,
    rz=3.67,
    precision_grade=6.0,
    shaft_diameter=35.0,
    schema=3.0,
    l=60.0,
    s=15.0,
    backlash=0.017
)

gear = Gear(
    profile=tool,
    material=material,
    z=40.,
    m=2.5,
    beta=16.0,
    alpha=20.0,
    x=0.0,
    b=34.0,
    bs=34.0,
    sr=0.0,
    rz=3.67,
    precision_grade=6.0,
    shaft_diameter=50.0,
    schema=3.0,
    l=60.0,
    s=35.0,
    backlash=-0.017
)

pair = [pinion, gear]

transmition = Transmition(
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
print iso_pitting(transmition=transmition).calculate()
print '========================================'

print '========================================'
print 'ISO Bending'
print iso_bending(transmition=transmition).calculate
print '========================================'

print '========================================'
print 'AGMA Pitting'
print agma_pitting(transmition=transmition).calculate()
print '========================================'

print '========================================'
print 'AGMA Bending'
print agma_bending(transmition=transmition).calculate()
print '========================================'

output_folder = os.path.join(os.path.dirname(__file__), 'comsol_output')
try:
    os.mkdir(output_folder)
except:
    pass

# 2D model matlab-comsol model export
# for 2D export type='2D' is optional because '2D' is the default output
transmition_model_name2d = 'transmition2d'
pinion_model_name2d = 'pinion2d'
wheel_model_name2d = 'wheel2d'
ExportGear(pinion).matlab_comsol(model_name=pinion_model_name2d, output_folder=output_folder, type='2D')
ExportGear(gear).matlab_comsol(model_name=wheel_model_name2d, output_folder=output_folder, type='2D')
ExportPair(pair).matlab_comsol(model_name=transmition_model_name2d, output_folder=output_folder, type='2D')

#3D model matlab-comsol model export
transmition_model_name3d = 'transmition3d'
pinion_model_name3d = 'pinion3d'
wheel_model_name3d = 'wheel3d'

ExportGear(pinion).matlab_comsol(model_name=pinion_model_name3d, output_folder=output_folder, type='3D')
ExportGear(gear).matlab_comsol(model_name=wheel_model_name3d, output_folder=output_folder, type='3D')
ExportPair(pair).matlab_comsol(model_name=transmition_model_name3d, output_folder=output_folder, type='3D')
