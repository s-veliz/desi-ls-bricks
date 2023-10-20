# 

import os
import subprocess
import time
import urllib.request
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as numpy

def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def get_bricks(name_cl, ra_cl, dec_cl, r200_cl, n):
    pos_cl = SkyCoord(ra_cl, dec_cl, unit='deg', frame='icrs')
    pos_brick = SkyCoord(ra_brick, dec_brick, unit='deg', frame='icrs')
    offsets_cl = (pos_cl.separation(pos_brick).to('arcmin')).value
    brick_index = np.where(offsets_cl < (n*r200_cl + semi_diag))
    brick_list = survey_bricks[brick_index[0]]['BRICKNAME']
    print(f'{len(brick_list)} bricks found for {name_cl}')
    with open(f'{output_path}{name_cl}_bricks.info', 'w') as file:
        for item in brick_list:
            file.write('%s\n' % item)

def download_bricks(output_path, name_cl):
    nersc_url = 'https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr10/south/coadd/'
    print(f'Downloading bricks for {name_cl}...')
    brick_list = Table.read(f'{output_path}/{name_cl}_bricks.info', format='ascii.fast_no_header')
    N = len(brick_list)
    for i in range(N):
        nnn = brick_list[i][0][0:3]
        brick_name = brick_list[i][0]
        print(f'Brick {brick_name}: Downloading...')
        parent = f'{output_path}{name_cl}'
        path = os.path.join(parent, brick_name)
        os.makedirs(path)
        url_g = nersc_url+nnn+'/'+brick_name+f'/legacysurvey-{brick_name}-image-g.fits.fz'
        url_r = nersc_url+nnn+'/'+brick_name+f'/legacysurvey-{brick_name}-image-r.fits.fz'
        url_i = nersc_url+nnn+'/'+brick_name+f'/legacysurvey-{brick_name}-image-i.fits.fz'
        url_z = nersc_url+nnn+'/'+brick_name+f'/legacysurvey-{brick_name}-image-z.fits.fz'
        name_g = f'{path}/legacysurvey-{brick_name}-image-g.fits.fz'
        name_r = f'{path}/legacysurvey-{brick_name}-image-r.fits.fz'
        name_i = f'{path}/legacysurvey-{brick_name}-image-i.fits.fz'
        name_z = f'{path}/legacysurvey-{brick_name}-image-z.fits.fz'
        urllib.request.urlretrieve(url_g, name_g)
        urllib.request.urlretrieve(url_r, name_r)
        urllib.request.urlretrieve(url_i, name_i)
        urllib.request.urlretrieve(url_z, name_z)
        print(f'Brick {brick_name}: Done.')
    print(f'All bricks downloaded for {name_cl}.')

survey_bricks = Table.read('survey-bricks.fits.gz', format='fits')
ra_brick = survey_bricks['RA']
dec_brick = survey_bricks['DEC']
length_brick = 15 #arcmin
semi_diag = (np.sqrt(2) * length_brick)/2

output_path = input('Set the output path (full path): ')

q1 = input('How many clusters do you want to download bricks for? ')

while not q1.isdigit():
    print('The number entered is invalid. Enter an integer.')
    q1 = input('How many clusters do you want to download bricks for? ')

if q1 == '1':

    name_cl = str(input('Enter the name of the cluster: '))
    ra_cl = input('Enter the RA of the cluster (degrees): ')
    while not is_float(ra_cl):
        print('Invalid RA. Enter a float variable.')
        ra_cl = input('Enter the RA of the cluster (degrees): ')
    dec_cl = input('Enter the DEC of the cluster (degrees): ')
    while not is_float(dec_cl):
        print('Invalid DEC. Enter a float variable.')
        dec_cl = input('Enter the DEC of the cluster (degrees): ')
    r200_cl = input('Enter the R200 of the cluster (arcmin): ')
    while not is_float(r200_cl):
        print('Invalid R200. Enter a float variable.')
        r200_cl = input('Enter the R200 of the cluster (degrees): ')
    n = input('Enter the lenght of the side of the square corresponding to the area of your interest of the cluster in R200 units (ex: 1.5, 2, 5, etc.): ')
    while not is_float(n):
        print('Invalid R200 factor. Enter a float variable.')
        n = input('Enter the lenght of the side of the square corresponding to the area of your interest of the cluster in R200 units (ex: 1.5, 2, 5, etc.): ')
    ra_cl = float(ra_cl)
    dec_cl = float(dec_cl)
    r200_cl = float(r200_cl)
    n = float(n)

    get_bricks(name_cl, ra_cl, dec_cl, r200_cl, n)
    download_bricks(output_path, name_cl)

else:

    table_path = str(input('Enter the path of your table (full path): '))
    table_format = str(input('Enter the format of your table (csv, ascii, fits, etc.): '))
    name_table = str(input('Enter the name of the column with the cluster IDs or cluster names: '))
    ra_table = str(input('Enter the name of the column with the cluster RAs: '))
    dec_table = str(input('Enter the name of the column with the cluster DECs: '))
    r200_table = str(input('Enter the name of the column with the cluster R200s: '))
    n = input('Enter the lenght of the side of the square corresponding to the area of your interest of the cluster in R200 units (ex: 1.5, 2, 5, etc.): ')
    while not is_float(n):
        print('Invalid R200 factor. Enter a float variable.')
        n = input('Enter the lenght of the side of the square corresponding to the area of your interest of the cluster in R200 units (ex: 1.5, 2, 5, etc.): ')
    n = float(n)
    clusters = Table.read(table_path, format=table_format)

    M = len(clusters)
    for j in range(M):
        get_bricks(clusters[name_table][j], clusters[ra_table][j], clusters[dec_table][j], clusters[r200_table][j], n)
        download_bricks(output_path, clusters[name_table][j])
        subprocess.run(f'mv {output_path}{name_cl}_bricks.info {output_path}{name_cl}/{name_cl}_bricks.info')

    print('DONE!')
