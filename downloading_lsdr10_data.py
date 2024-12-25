####### LIBRARIES #######

# General libraries
import os
import shutil
import time
import subprocess
import urllib.request
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata import NDData, Cutout2D
from astropy.coordinates import SkyCoord
from astropy.io import fits
from getpass import getpass
from photutils.psf import extract_stars, EPSFBuilder
# Reproject libraries
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
# NOIRLab's libraries
from dl import authClient as ac, queryClient as qc
from dl.helpers.utils import convert

####### CONFIGS AND FIXED VARIABLES #######

brick_length = 15 # arcmin
brick_semi_diag = (np.sqrt(2) * brick_length) / 2
pixel_scale = 0.262 # arcsec/pixel

####### INPUT VARIABLES #######

ls_bricks = Table.read('survey-bricks.fits.gz', format='fits')
ra_bricks = ls_bricks['RA']
dec_bricks = ls_bricks['DEC']
initial_path = '/Users/sveliz/Documents/Projects/CHANCES_evolution/CHANCES_morphology/Downloading_DECaLS_data/'
output_path = '/Users/sveliz/Documents/Projects/CHANCES_evolution/CHANCES_morphology/Downloading_DECaLS_data/LS_DR10_data_local/'
galaxies = pd.read_csv('../ds18_gzd_balanced_local.csv')
ra_column = 'ra'
dec_column = 'dec'
name_column = 'iauname'

####### FUNCTIONS #######

def verify_url(url):
	try:
		request = urllib.request.Request(url, method='HEAD')
		response = urllib.request.urlopen(request)
		if response.status == 200:
			return True
		else:
			return False
	except:
		return False

def get_bricks(name, ra, dec, radius=2.5):

	pos_obj = SkyCoord(ra, dec, unit='deg', frame='icrs')
	pos_brick = SkyCoord(ra_bricks, dec_bricks, unit='deg', frame='icrs')
	offsets = (pos_obj.separation(pos_brick).to('arcmin')).value
	brick_indices = np.where(offsets < radius + brick_semi_diag)
	brick_list = ls_bricks[brick_indices[0]]['BRICKNAME'].value
	brick_list = [b.decode('utf-8') if isinstance(b, bytes) else b for b in brick_list]

	print(f'{len(brick_list)} bricks found for {name}')
	# Verify the available bands
	brick_g = []
	brick_r = []
	brick_i = []
	brick_z = []
	nersc_url = 'https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr10/south/coadd/'
	for i, brick_name in enumerate(brick_list, start=0):
		nnn = brick_name[0:3]
		url_g = nersc_url+nnn+'/'+brick_name+f'/legacysurvey-{brick_name}-image-g.fits.fz'
		url_r = nersc_url+nnn+'/'+brick_name+f'/legacysurvey-{brick_name}-image-r.fits.fz'
		url_i = nersc_url+nnn+'/'+brick_name+f'/legacysurvey-{brick_name}-image-i.fits.fz'
		url_z = nersc_url+nnn+'/'+brick_name+f'/legacysurvey-{brick_name}-image-z.fits.fz'
		if verify_url(url_g):
			print(f'g-band available for brick {brick_name}')
			brick_g.append(1)
		else:
			print(f'WARNING: g-band not available for brick {brick_name}')
			brick_g.append(0)
		if verify_url(url_r):
			print(f'r-band available for brick {brick_name}')
			brick_r.append(1)
		else:
			print(f'WARNING: r-band not available for brick {brick_name}')
			brick_r.append(0)
		if verify_url(url_i):
			print(f'i-band available for brick {brick_name}')
			brick_i.append(1)
		else:
			print(f'WARNING: i-band not available for brick {brick_name}')
			brick_i.append(0)
		if verify_url(url_z):
			print(f'z-band available for brick {brick_name}')
			brick_z.append(1)
		else:
			print(f'WARNING: z-band not available for brick {brick_name}')
			brick_z.append(0)
	bricks_dat = pd.DataFrame({'brick': brick_list, 'g': brick_g, 'r': brick_r, 'i': brick_i, 'z': brick_z})
	bricks_dat.to_csv(f'{output_path}{name}_bricks.csv', index=False)
	return brick_list

def download_bricks(output_path, name, mode):
	"""
	mode: image, weight or other
	"""
	nersc_url = 'https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr10/south/coadd/'
	print(f'Downloading bricks for {name}...')
	bricks_dat = pd.read_csv(f'{output_path}{name}_bricks.csv')
	brick_list = bricks_dat['brick'].values
	bricks_g = bricks_dat['g'].values
	bricks_r = bricks_dat['r'].values
	bricks_i = bricks_dat['i'].values
	bricks_z = bricks_dat['z'].values
	parent = f'{output_path}{name}'
	flag_g = np.sum(bricks_g)
	flag_r = np.sum(bricks_r)
	flag_i = np.sum(bricks_i)
	flag_z = np.sum(bricks_z)
	N = len(brick_list)
	for i, brick_name in enumerate(brick_list, start=0):
		nnn = brick_name[0:3]
		print(f'Brick {brick_name}: Downloading {mode}s...')
		path = os.path.join(parent, brick_name)
		os.makedirs(path, exist_ok=True)
		url_g = nersc_url+nnn+'/'+brick_name+f'/legacysurvey-{brick_name}-{mode}-g.fits.fz'
		url_r = nersc_url+nnn+'/'+brick_name+f'/legacysurvey-{brick_name}-{mode}-r.fits.fz'
		url_i = nersc_url+nnn+'/'+brick_name+f'/legacysurvey-{brick_name}-{mode}-i.fits.fz'
		url_z = nersc_url+nnn+'/'+brick_name+f'/legacysurvey-{brick_name}-{mode}-z.fits.fz'
		name_g = f'{path}/legacysurvey-{brick_name}-{mode}-g.fits.fz'
		name_r = f'{path}/legacysurvey-{brick_name}-{mode}-r.fits.fz'
		name_i = f'{path}/legacysurvey-{brick_name}-{mode}-i.fits.fz'
		name_z = f'{path}/legacysurvey-{brick_name}-{mode}-z.fits.fz'
		if flag_g == N:
			urllib.request.urlretrieve(url_g, name_g)
		if flag_r == N:
			urllib.request.urlretrieve(url_r, name_r)
		if flag_i == N:
			urllib.request.urlretrieve(url_i, name_i)
		if flag_z == N:
			urllib.request.urlretrieve(url_z, name_z)
		print(f'Brick {brick_name}: Done ({mode}s).')
	return N, flag_g, flag_r, flag_i, flag_z

def read_bricks(output_path, name, band, mode):
	bricks_dat = pd.read_csv(f'{output_path}{name}_bricks.csv')
	bricks = bricks_dat['brick'].values
	for brick in bricks:
		# Extract .fits.gz before reading
		path = f'{output_path}{name}/{brick}'
		print('Extracting files...')
		os.chdir(path)
		os.system(f'funpack {path}/*.fits.fz > /dev/null 2>&1')
	images = [output_path+name+f'/{brick}/legacysurvey-{brick}-{mode}-{band}.fits' for brick in bricks]
	print('All files extracted!')
	return images

def join_bricks(images, band, mode, output_path, name):
	print(f'Joining {band}-band {mode}...')
	hdus = [fits.open(image) for image in images]
	wcs_out, shape_out = find_optimal_celestial_wcs(hdus)
	array, footprint = reproject_and_coadd(hdus, wcs_out, shape_out=shape_out, reproject_function=reproject_interp)
	new_header = wcs_out.to_header()
	wcs = WCS(new_header)
	print(f'{band}-band {mode} joined!')
	if mode == 'image':
		joined = fits.PrimaryHDU(array, header=wcs.to_header())
		joined.writeto(f'{output_path}{name}_image_{band}_full.fits', overwrite=True)
	return array, wcs

def cutout_and_save(name, coord_sky, array, wcs, band, mode, size=5):
	print(f'Cutting {band}-band {mode}...')
	size_pix = (size*u.arcmin).to('arcsec').value
	x, y = wcs.world_to_pixel(coord_sky)
	coord_pix = (x, y)
	cutout = Cutout2D(array, position=coord_pix, size=size_pix, wcs=wcs, mode='partial')
	filename = output_path+f'{name}_{mode}_{band}.fits'
	header = cutout.wcs.to_header()
	hdu = fits.PrimaryHDU(cutout.data, header=header)
	hdu.writeto(filename, overwrite=True)
	print(f'{band}-band {mode} cut and saved!')

def has_neighbors(catalog, segmap, X='X_IMAGE', Y='Y_IMAGE', n='NUMBER', size=25):
	neighbors = []
	for i, candidate in catalog.iterrows():
		x = int(candidate[X])
		y = int(candidate[Y])
		position = (x, y)
		number = candidate[n]
		cutout = Cutout2D(data=segmap, position=position, size=size*u.pixel, mode='partial', fill_value=0)
		cutout_segmap = cutout.data
		mask = (cutout_segmap != 0) & (cutout_segmap != number)
		neighbor = np.any(mask)
		neighbors.append(neighbor)
	catalog['neighbors'] = neighbors
	return catalog

def make_psf(path, name, band, psf_fwhm, x='X_IMAGE', y='Y_IMAGE'):
	os.chdir(path)
	subprocess.run(f'sex {output_path}{name}_image_{band}_full.fits -c common.se', shell=True)
	image = fits.open(f'{output_path}{name}_image_{band}_full.fits')[0].data
	segmap = fits.open('segmentation.fits')[0].data
	sources = Table.read('se.dat', format='ascii')
	candidates = sources[(sources['FLUX_RADIUS'] < np.pi*psf_fwhm) & (sources['CLASS_STAR'] > 0.8)]
	candidates = candidates.to_pandas()
	stars = has_neighbors(candidates, segmap)
	stars_tbl = Table()
	stars_tbl['x'] = stars[x]
	stars_tbl['y'] = stars[y]
	nddata = NDData(data=image)
	stars = extract_stars(nddata, stars_tbl, size=25)
	epsf_builder = EPSFBuilder(oversampling=1, maxiters=3, progress_bar=False)
	epsf, fitted_stars = epsf_builder(stars)
	hdu = fits.PrimaryHDU(epsf.data)
	hdu.writeto(f'{output_path}{name}_psf_{band}.fits', overwrite=True)
	os.remove('segmentation.fits')
	os.remove('se.dat')
	os.remove(f'LS_DR10_data_local/{name}_image_{band}_full.fits')

def match_ls(ra, dec, radius=5):
    """
    ra: degrees
    dec: degrees
    radius: arcsec
    """
    radius = (radius*u.arcsec).to('deg').value
    query = f"""
            SELECT release, brickid, objid, brick_primary, type, ra, dec, psfsize_g, psfsize_r, psfsize_i, psfsize_z, ngood_g, ngood_r,
            ngood_i, ngood_z,
            q3c_dist(ra, dec, {ra}, {dec}) AS dist
            FROM ls_dr10.tractor
            WHERE q3c_radial_query(ra, dec, {ra}, {dec}, {radius})
            ORDER BY dist ASC
            LIMIT 1
            """
    result = qc.query(sql=query)
    catalog = convert(result, 'pandas')
    return catalog

####### TESTING #########

"""# Datos de la galaxia
name = 'J120354.19+015324.9'
ra = 180.97606246614552
dec = 1.8912349738097065

# Buscamos el PSF_FWHM en cada filtro
galaxy_data = match_ls(ra=ra, dec=dec)
psf_g = galaxy_data['psfsize_g'].values
psf_r = galaxy_data['psfsize_r'].values
psf_i = galaxy_data['psfsize_i'].values
psf_z = galaxy_data['psfsize_z'].values

# Creamos la psf en el filtro g
make_psf(initial_path, name, 'g', psf_g)"""

####### MAIN CODE #######

# Login DataLab
print('Log in DataLab:')
token = ac.login(input('Enter username: '),getpass('Enter password: '))
ac.whoAmI()

for i in range(11414, 11464):

	# Variable auxiliar para medir el tiempo
	t_i = time.time()

	# Input de datos básicos de la galaxia
	name = galaxies[name_column][i]
	ra = galaxies[ra_column][i]
	dec = galaxies[dec_column][i]
	coord_sky = SkyCoord(ra=ra, dec=dec, unit='deg', frame='icrs')

	# Mostrar por pantalla que comenzamos la iteración
	print(f'Initializing iteration {i+1}. Galaxy: {name}')

	# Identificar los bricks asociados a la galaxia
	get_bricks(name, ra, dec)

	# Descargar los bricks
	download_bricks(output_path=output_path, name=name, mode='image')

	# Verificar la cantidad de bandas disponibles y guardar la información en un archivo
	N, flag_g, flag_r, flag_i, flag_z = download_bricks(output_path=output_path, name=name, mode='invvar')
	os.chdir(initial_path)
	old_data = pd.read_csv('galaxies_data_local.csv')
	new_data = pd.DataFrame({'id': [name], 'ra': [ra], 'dec': [dec], 'N': [N],'g': [flag_g], 'r': [flag_r], 'i': [flag_i], 'z': [flag_z]})
	updated_data = pd.concat([old_data, new_data], ignore_index=True)
	updated_data.to_csv('galaxies_data_local.csv', index=False)

	# Buscamos el PSF_FWHM en cada filtro
	galaxy_data = match_ls(ra=ra, dec=dec)
	psf_g = galaxy_data['psfsize_g'].values
	psf_r = galaxy_data['psfsize_r'].values
	psf_i = galaxy_data['psfsize_i'].values
	psf_z = galaxy_data['psfsize_z'].values


	# Recorte en el filtro g y creación de PSF
	if flag_g == N:
		images_g = read_bricks(output_path=output_path, name=name, band='g', mode='image')
		array_image_g, wcs_image_g = join_bricks(images_g, band='g', mode='image', output_path=output_path, name=name)
		cutout_and_save(name=name, coord_sky=coord_sky, array=array_image_g, wcs=wcs_image_g, band='g', mode='image')
		invvars_g = read_bricks(output_path=output_path, name=name, band='g', mode='invvar')
		array_invvar_g, wcs_invvar_g = join_bricks(invvars_g, band='g', mode='invvar', output_path=output_path, name=name)
		cutout_and_save(name=name, coord_sky=coord_sky, array=array_invvar_g, wcs=wcs_image_g, band='g', mode='invvar')
		make_psf(initial_path, name, 'g', psf_g)
	else:
		print('g-band not available... Moving to the next band.')
	if flag_r == N:
		images_r = read_bricks(output_path=output_path, name=name, band='r', mode='image')
		array_image_r, wcs_image_r = join_bricks(images_r, band='r', mode='image', output_path=output_path, name=name)
		cutout_and_save(name=name, coord_sky=coord_sky, array=array_image_r, wcs=wcs_image_r, band='r', mode='image')
		invvars_r = read_bricks(output_path=output_path, name=name, band='r', mode='invvar')
		array_invvar_r, wcs_invvar_r = join_bricks(invvars_r, band='r', mode='invvar', output_path=output_path, name=name)
		cutout_and_save(name=name, coord_sky=coord_sky, array=array_invvar_r, wcs=wcs_image_r, band='r', mode='invvar')
		make_psf(initial_path, name, 'r', psf_r)
	else:
		print('r-band not available... Moving to the next band.')
	if flag_i == N:
					images_i = read_bricks(output_path=output_path, name=name, band='i', mode='image')
					array_image_i, wcs_image_i = join_bricks(images_i, band='i', mode='image', output_path=output_path, name=name)
					cutout_and_save(name=name, coord_sky=coord_sky, array=array_image_i, wcs=wcs_image_i, band='i', mode='image')
					invvars_i = read_bricks(output_path=output_path, name=name, band='i', mode='invvar')
					array_invvar_i, wcs_invvar_i = join_bricks(invvars_i, band='i', mode='invvar', output_path=output_path, name=name)
					cutout_and_save(name=name, coord_sky=coord_sky, array=array_invvar_i, wcs=wcs_image_i, band='i', mode='invvar')
					make_psf(initial_path, name, 'i', psf_i)
	else:
		print('i-band not available... Moving to the next band.')
	if flag_z == N:
		images_z = read_bricks(output_path=output_path, name=name, band='z', mode='image')
		array_image_z, wcs_image_z = join_bricks(images_z, band='z', mode='image', output_path=output_path, name=name)
		cutout_and_save(name=name, coord_sky=coord_sky, array=array_image_z, wcs=wcs_image_z, band='z', mode='image')
		invvars_z = read_bricks(output_path=output_path, name=name, band='z', mode='invvar')
		array_invvar_z, wcs_invvar_z = join_bricks(invvars_z, band='z', mode='invvar', output_path=output_path, name=name)
		cutout_and_save(name=name, coord_sky=coord_sky, array=array_invvar_z, wcs=wcs_image_z, band='z', mode='invvar')
		make_psf(initial_path, name, 'z', psf_z)
	else:
		print('z-band not available...')
	print(f'All data downloaded and processed for {name}... Cleaning trash...')
	os.remove(f'{output_path}{name}_bricks.csv')
	shutil.rmtree(f'{output_path}{name}')
	t_f = time.time()
	print(f'All done for galaxy {name}. Iteration {i+1} finished in {round(t_f-t_i)}s.')