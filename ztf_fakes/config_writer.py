from configparser import ConfigParser
import os

workdir = os.getcwd()

config = ConfigParser(allow_no_values=True)
config['file_structure'] = {
    'workdir'     : workdir,
    'logfile'     : os.path.join(workdir,'fakes.log'),
    'fakecatalog' : os.path.join(workdir,workdir+'fakes.cat'),
    'SEcatalog'   : os.path.join(workdir,workdir+'sex.cat'),
    'galaxydir'   : os.path.join(workdir,workdir+'galaxy.cat')
}

config['common_checks'] = {
    'min_num_obj'    : '40',
    'allowed_flag'   : '0',
    'edge_limit'     : '40'
}

config['star_checks'] = {
    'min_mag'         : '10.0',
    'max_mag'         : '20.0',
    'min_ellip'       : '0.2',
    'max_ellip'       : '',
    'min_elgn'        : '0.9' ,
    'max_elgn'        : '5.0',
    'min_class_star'  : '0.5',
    'max_class_star'  : '1.0',
    'min_FWHM'        : '2.0',
    'max_FWHM'        : '15.0'
}

config['galaxy_checks'] = {
    'min_mag'         : '10.0',
    'max_mag'         : '20.0',
    'min_ellip'       : '0.3',
    'max_ellip'       : '',
    'min_elgn'        : '0.3',
    'max_elgn'        : '3.0',
    'min_class_star'  : '0.5',
    'max_class_star'  : '1.0',
    'min_FWHM'        : '2.0',
    'max_FWHM'        : '15.0'
}

with open('etc/conf.ini','w') as configfile:
    config.write(configfile)
