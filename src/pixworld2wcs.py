# Older version of code by Clare Shanahan
# Now part of astropy PR: https://github.com/astropy/astropy/pull/7884



from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
from astropy.wcs.utils import celestial_frame_to_wcs
from astropy.modeling import Fittable2DModel,Parameter, models, fitting
from astropy.modeling.models import custom_model
from astropy.modeling.fitting import LevMarLSQFitter

def _initial_guess_fit(x, y, cords, proj_point = None):

	if not proj_point:
		x0, y0 = x - np.mean(x),y - np.mean(y)
		ra0, dec0 = cords.ra.deg - np.mean(cords.ra.deg), \
					cords.dec.deg - np.mean(cords.dec.deg)
	else:
		(x_0, y_0), (ra_0, dec_0) = proj_point
		x0, y0 = x - proj_point[0][0],y - proj_point[0][1]
		ra0, dec0 = cords.ra.deg - proj_point[1][0], cords.dec.deg - proj_point[1][1]

	wcs_obj_0 = wcs.WCS()
		
	#lstsq fit, ra to x and y for PC11_0, PC12_0
	lstq_ra = np.linalg.lstsq(np.column_stack((x0,y0)),ra0)
	#lstsq fit, dec to x and y for PC11_0, PC12_0
	lstsq_dec = np.linalg.lstsq(np.column_stack((x0,y0)),dec0)

	wcs_obj_0.wcs.pc = ((lstq_ra[0]),(lstsq_dec[0]))
	wcs_obj_0.wcs.crpix = (np.mean(x),np.mean(y))
	wcs_obj_0.wcs.crval = (np.mean(cords.ra.deg),np.mean(cords.dec.deg))
		
	if not proj_point:
		return wcs_obj_0
	else:
		wcs_obj_0.wcs.crpix = proj_point[0]
		wcs_obj_0.wcs.crval = proj_point[1]
		return wcs_obj_0
	
class _projmodel2D_linear(Fittable2DModel):

	PC11 = Parameter()
	PC12 = Parameter()
	PC21 = Parameter()
	PC22 = Parameter()

	CRPIX0 = Parameter(fixed = True)
	CRPIX1 = Parameter(fixed = True)
	CRVAL0 = Parameter(fixed = True)
	CRVAL1 = Parameter(fixed = True)
	ax = Parameter(fixed = True)

	inputs = ('x','y')
	outputs = ('ax_cords',)
	
	@staticmethod
	def evaluate(x,y,PC11,PC12,PC21,PC22,CRPIX0,CRPIX1,CRVAL0,CRVAL1,ax):
	
		xy_cords = np.column_stack((x,y))
	
		wcs_obj = wcs.WCS()
	
		#need to figure out how to pass strings as parameters or make a mapping of 
		#frames and projection types to ints to generalize this...
		wcs_obj.wcs.ctype = 'RA---TAN',	 'DEC--TAN'
		wcs_obj.wcs.crpix = (CRPIX0[0],CRPIX1[0])
		wcs_obj.wcs.crval = (CRVAL0[0],CRVAL1[0])
	
		wcs_obj.wcs.pc = ((PC11[0],PC12[0]),(PC21[0],PC22[0]))

		world_cords = wcs_obj.wcs_pix2world(xy_cords,1)
		ax_cords = world_cords[:,int(ax) - 1] 
	
		return ax_cords
			
def wcs_pixworld2wcs(input_table, projection = 'TAN', proj_point = None):

	""" For input pixel positions and their corresponding sky positions, returns the best
		fitting WCS object. Only the linear terms are fit - No `SIP` or 
		`distortion paper` table lookup corrections are fit. The terms of the PC matrix
		are fit independently to RA and Dec, Levenberg-Marquardt algorithm and 
		least squares statistic.
		
	
		Parameters
		----------
		input_table: astropy table
			Astropy Table with three columns: X pixel positions, Y pixel positions,
			and a Skycoord object with the sky positions (which encodes the celestial 
			frame)
			
		projection: str
			Code to use in ctype, default is 'TAN'
			
		proj_point: list or tuple 
			((X,Y) and (RA, Dec)) of desired projection point. If None, the geometric
			center of input positions is chosen as the point for the projection.
	"""
		
	#unpack input table (x pixel positions, y pixel positions, skycoord object)
	x = input_table['x']
	y = input_table['y']
	skycords = input_table['coord']

	#empty WCS object with keywords populated based on frame and projection type
	wcs_obj = celestial_frame_to_wcs(frame = skycords.frame, projection = projection)

	#initial guess for PC matrix terms, and optionally for CRVAL
	wcs_initial_guess = _initial_guess_fit(x,y,skycords,proj_point = proj_point)
	(PC11_0,PC12_0),(PC21_0,PC22_0) = wcs_initial_guess.wcs.pc
	crpix = wcs_initial_guess.wcs.crpix
	crval = wcs_initial_guess.wcs.crval

	#fitter
	fit = LevMarLSQFitter()
									  		  
	#fit RA
	PC = _projmodel2D_linear(PC11_0,PC12_0,PC21_0,PC22_0,crpix[0],crpix[1],crval[0],
							crval[1],ax=1)	
	mPC = fit(PC, x, y, skycords.ra.deg)
	fit_pcA,fit_pcB  = mPC.PC11.value, mPC.PC12.value
	
	#fit Dec
	PC = _projmodel2D_linear(mPC.PC11.value, mPC.PC12.value, mPC.PC21.value, 
							mPC.PC22.value,crpix[0],crpix[1],crval[0],crval[1],ax=2)	
	mPC = fit(PC, x, y, skycords.dec.deg)
	fit_pcC,fit_pcD=mPC.PC21.value,mPC.PC22.value
		
	#add fit terms to wcs				 
	wcs_obj.wcs.crval = crval
	wcs_obj.wcs.crpix = crpix
	wcs_obj.wcs.pc = ((fit_pcA,fit_pcB),(fit_pcC,fit_pcD))
	
	return wcs_obj
	
