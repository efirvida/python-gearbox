from distutils.core import setup

setup(
    name='python_gearbox',
    version='0.1',
    packages=['', 'demo', 'libs', 'export', 'standards', 'transmition'],
    url='',
    license='GPL 3.0',
    author='efirvida',
    author_email='efirvida@gmail.com',
    description='Python library for gear transmission design', requires=['numpy', 'scipy'],
    url='https://github.com/efirvida/python-gearbox',
    download_url='https://github.com/efirvida/python-gearbox/archive/master.zip',
    keywords=['gearbox', 'gear', 'agma', 'iso', 'gear transmission', 'engineering'],
)
