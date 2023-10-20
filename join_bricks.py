import subprocess
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy import units as u
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd

pix_scale = 0.262 # arcsec/pixel

cl_name = input('Enter the cluster name: ')
cl_ra = float(input('Enter the cluster RA (degrees): '))
cl_dec = float(input('Enter the cluster DEC (degrees): '))
cl_path = input('Enter the path to the cluster bricks: ')
cl_r200_arcmin = float(input('Enter the cluster R200 (arcmin): '))

subprocess.run(f'mkdir {cl_path}{cl_name}_images', shell=True)

cl_r200_pixel = (cl_r200_arcmin * 60.) / pix_scale
cl_coord_sky = SkyCoord(cl_ra, cl_dec, frame='icrs', unit='deg')

cluster_bricks_table = Table.read(f'{cl_path}{cl_name}_bricks.info', format='ascii.fast_no_header')
cluster_bricks = np.array(cluster_bricks_table['col1'])

cl_g_images = [cl_path+cl_name+f'/{brick}/legacysurvey-{brick}-image-g.fits.fz' for brick in cluster_bricks]
cl_r_images = [cl_path+cl_name+f'/{brick}/legacysurvey-{brick}-image-r.fits.fz' for brick in cluster_bricks]
cl_i_images = [cl_path+cl_name+f'/{brick}/legacysurvey-{brick}-image-i.fits.fz' for brick in cluster_bricks]
cl_z_images = [cl_path+cl_name+f'/{brick}/legacysurvey-{brick}-image-z.fits.fz' for brick in cluster_bricks]

cl_hdus_g = [fits.open(image)[1] for image in cl_g_images]
cl_hdus_r = [fits.open(image)[1] for image in cl_r_images]
cl_hdus_i = [fits.open(image)[1] for image in cl_i_images]
cl_hdus_z = [fits.open(image)[1] for image in cl_z_images]

wcs_out_g, shape_out_g = find_optimal_celestial_wcs(cl_hdus_g)
wcs_out_r, shape_out_r = find_optimal_celestial_wcs(cl_hdus_r)
wcs_out_i, shape_out_i = find_optimal_celestial_wcs(cl_hdus_i)
wcs_out_z, shape_out_z = find_optimal_celestial_wcs(cl_hdus_z)
array_g, footprint_g = reproject_and_coadd(cl_hdus_g, wcs_out_g, shape_out=shape_out_g,
	reproject_function=reproject_interp)
array_r, footprint_r = reproject_and_coadd(cl_hdus_r, wcs_out_r, shape_out=shape_out_r,
	reproject_function=reproject_interp)
array_i, footprint_i = reproject_and_coadd(cl_hdus_i, wcs_out_i, shape_out=shape_out_i,
	reproject_function=reproject_interp)
array_z, footprint_z = reproject_and_coadd(cl_hdus_z, wcs_out_z, shape_out=shape_out_z,
	reproject_function=reproject_interp)

new_header_g = wcs_out_g.to_header()
new_header_r = wcs_out_r.to_header()
new_header_i = wcs_out_i.to_header()
new_header_z = wcs_out_z.to_header()

wcs_g = WCS(new_header_g)
wcs_r = WCS(new_header_r)
wcs_i = WCS(new_header_i)
wcs_z = WCS(new_header_z)

x_g, y_g = wcs_g.world_to_pixel(cl_coord_sky)
x_r, y_r = wcs_r.world_to_pixel(cl_coord_sky)
x_i, y_i = wcs_i.world_to_pixel(cl_coord_sky)
x_z, y_z = wcs_z.world_to_pixel(cl_coord_sky)

cl_coord_img_g = (x_g, y_g)
cl_coord_img_r = (x_r, y_r)
cl_coord_img_i = (x_i, y_i)
cl_coord_img_z = (x_z, y_z)

cl_cutout_g = Cutout2D(array_g, position=cl_coord_img_g, size=2*cl_r200_pixel*u.pixel,
	wcs=wcs_g, mode='partial')
cl_cutout_r = Cutout2D(array_r, position=cl_coord_img_r, size=2*cl_r200_pixel*u.pixel,
	wcs=wcs_r, mode='partial')
cl_cutout_i = Cutout2D(array_i, position=cl_coord_img_i, size=2*cl_r200_pixel*u.pixel,
	wcs=wcs_i, mode='partial')
cl_cutout_z = Cutout2D(array_z, position=cl_coord_img_z, size=2*cl_r200_pixel*u.pixel,
	wcs=wcs_z, mode='partial')

cl_filename_g = f'{cl_path}{cl_name}_images/{cl_name}_g.fits'
cl_filename_r = f'{cl_path}{cl_name}_images/{cl_name}_r.fits'
cl_filename_i = f'{cl_path}{cl_name}_images/{cl_name}_i.fits'
cl_filename_z = f'{cl_path}{cl_name}_images/{cl_name}_z.fits'
cl_header_g = cl_cutout_g.wcs.to_header()
cl_header_r = cl_cutout_r.wcs.to_header()
cl_header_i = cl_cutout_i.wcs.to_header()
cl_header_z = cl_cutout_z.wcs.to_header()
cl_hdu_g = fits.PrimaryHDU(cl_cutout_g.data, header=cl_header_g)
cl_hdu_r = fits.PrimaryHDU(cl_cutout_r.data, header=cl_header_r)
cl_hdu_i = fits.PrimaryHDU(cl_cutout_i.data, header=cl_header_i)
cl_hdu_z = fits.PrimaryHDU(cl_cutout_z.data, header=cl_header_z)
cl_hdu_g.writeto(cl_filename_g, overwrite=True)
cl_hdu_r.writeto(cl_filename_r, overwrite=True)
cl_hdu_i.writeto(cl_filename_i, overwrite=True)
cl_hdu_z.writeto(cl_filename_z, overwrite=True)