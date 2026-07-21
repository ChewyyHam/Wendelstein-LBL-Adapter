#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Wrap file for LBL recipes

Created on 2026-07-20 10:26:55.864

@author: Neil Cook, Etienne Artigau, Charles Cadieux, Thomas Vandal, Ryan Cloutier, Pierre Larue

user: liuyu@Beaglegeuse
lbl version: 0.67.008
lbl date: 2026-06-06
"""
from lbl import lbl_wrap

# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # set up parameters
    rparams = dict()
    # -------------------------------------------------------------------------
    # LBL parameters
    # -------------------------------------------------------------------------
    # You may also add any constant here to override the default value
    #     (see README for details) - this is NOT recommended for non developers
    #   Note this may have undesired affects as these parameters apply globally
    #     for all LBL recipes
    # -------------------------------------------------------------------------
    # This is the instrument name
    #   Currently supported instruments are 
	#		SPIROU
	#		HARPS
	#		ESPRESSO
	#		CARMENES
	#		NIRPS_HA
	#		NIRPS_HE
	#		HARPSN
	#		MAROONX
	#		SOPHIE
	#		CORALIE
	#		EXPRES
	#		NEID
	#		Generic
    rparams['INSTRUMENT'] = 'Generic'
    #   Data source must be as follows: 
	#		SPIROU: APERO or CADC
	#		NIRPS_HA: APERO or CADC or ESO
	#		NIRPS_HE: APERO or CADC or ESO
	#		HARPS: ORIG or ESO or ESSP
	#		CARMENES: None
	#		ESPRESSO: None
	#		HARPSN: ORIG or ESO or ESSP
	#		MAROONX: RED or BLUE
	#		SOPHIE: None
	#		CORALIE: None
	#		EXPRES: ORIG or ESSP
	#		NEID: ESSP
	#		Generic: None
    rparams['DATA_SOURCE'] = 'None'
    # The data directory where all data is stored under - this should be an
    #    absolute path
    rparams['DATA_DIR'] = r'E:\LBL_test\data'
    # The input file string (including wildcards) - if not set will use all
    #   files in the science directory (for this object name)
    # rparams['INPUT_FILE'] = '*'
    # The input science data are blaze corrected
    rparams['BLAZE_CORRECTED'] = True
    # Override the blaze filename
    #      (if not set will use the default for instrument)
    # rparams['BLAZE_FILE'] = 'blaze.fits'
    # -------------------------------------------------------------------------
    # science criteria
    # -------------------------------------------------------------------------
    # The data type (either SCIENCE or FP or LFC)
    rparams['DATA_TYPES'] = ['SCIENCE']
    # The object name (this is the directory name under the /science/
    #    sub-directory and thus does not have to be the name in the header
    rparams['OBJECT_SCIENCE'] = ['TOI5786']
    # This is the template that will be used or created (depending on what is
    #   run)
    rparams['OBJECT_COMPARISON'] = ['TOI5786']
    # This is the object temperature in K - used for getting a stellar model
    #   for the masks it only has to be good to a few 100 K
    rparams['OBJECT_TEFF'] = [5500]
    # -------------------------------------------------------------------------
    # what to run and skip if already on disk
    # -------------------------------------------------------------------------
    # Whether to reset all files before processing
    rparams['RUN_LBL_RESET'] = False
    # Whether to run the telluric cleaning process (NOT recommended for data
    #   that has better telluric cleaning i.e. SPIROU using APERO)
    rparams['RUN_LBL_TELLUCLEAN'] = False
    # Whether to create templates from the data in the science directory
    #   If a template has been supplied from elsewhere this set is NOT required
    rparams['RUN_LBL_TEMPLATE'] = True
    # Whether to create a mask using the template created or supplied
    rparams['RUN_LBL_MASK'] = True
    # Whether to run the LBL compute step - which computes the line by line
    #   for each observation
    rparams['RUN_LBL_COMPUTE'] = True
    # Whether to run the LBL compile step - which compiles the rdb file and
    #   deals with outlier rejection
    rparams['RUN_LBL_COMPILE'] = True
    # whether to skip observations if a file is already on disk (useful when
    #   adding a few new files) there is one for each RUN_XXX step
    #   - Note cannot skip tellu clean
    rparams['SKIP_LBL_TEMPLATE'] = True
    rparams['SKIP_LBL_MASK'] = True
    rparams['SKIP_LBL_COMPUTE'] = True
    rparams['SKIP_LBL_COMPILE'] = True
    # -------------------------------------------------------------------------
    # LBL settings
    # -------------------------------------------------------------------------
    # You can change any setting in parameters (or override those changed
    #   by specific instruments) here
    # -------------------------------------------------------------------------
    # Advanced settings
    #   Do not use without contacting the LBL developers
    # -------------------------------------------------------------------------
    # Dictionary of table name for the file used in the projection against the
    #     derivative. Key is to output column name that will propagate into the
    #     final RDB table and the value is the filename of the table. The table
    #     must follow a number of characteristics explained on the LBL website.
    # rparams['RESPROJ_TABLES'] = {'DTEMP3500': 'temperature_gradient_3500.fits'}

    # Rotational velocity parameters, should be a list of two values, one being
    #     the epsilon and the other one being the vsini in km/s as defined in the
    #     PyAstronomy.pyasl.rotBroad function
    # rparams['ROTBROAD'] = []

    # turn on plots
    rparams['PLOT'] = False

    # -------------------------------------------------------------------------
    # Other settings
    # -------------------------------------------------------------------------
    
    
    # =====================================================================
    # Generic instrument configuration for Wendelstein / MaHPS
    #
    # IMPORTANT:
    # All parameters in this section must contain explicit non-None values
    # when using the LBL Generic instrument.
    #
    # The wavelength limits must match the physical spectral-order range
    # selected during the Wendelstein-to-LBL conversion.
    # =====================================================================

    # ---------------------------------------------------------------------
    # Instrument identity
    # ---------------------------------------------------------------------

    rparams['GENERIC_INSTRUMENT'] = 'WEND'

    # Optional reduction-pipeline or observing-mode identifier.
    # This must be the string 'None', not the Python value None.
    rparams['GENERIC_DATA_SOURCE'] = 'None'

    # ---------------------------------------------------------------------
    # Instrument wavelength coverage
    # ---------------------------------------------------------------------

    # These values correspond to physical orders 84--114:
    #
    # order 84 lower wavelength boundary: 497.522 nm
    # order 114 upper wavelength boundary: 684.933 nm
    #
    # IMPORTANT:
    # Change these values if a different order range is selected in the
    # Wendelstein-to-LBL adapter. The values should be taken from the order
    # wavelength table in Appendix V of Hanna Kellermann's work.
    rparams['GENERIC_WAVEMIN'] = 497.522
    rparams['GENERIC_WAVEMAX'] = 684.933

    # ---------------------------------------------------------------------
    # Observatory and basic processing settings
    # ---------------------------------------------------------------------

    # Registered Astropy observatory-site name
    rparams['EARTH_LOCATION'] = 'wendelstein'

    # High-pass filtering width in km/s
    rparams['HP_WIDTH'] = 256

    # Minimum observation SNR accepted by LBL
    rparams['SNR_THRESHOLD'] = 5.0

    # Bands used when constructing the clean CCF
    rparams['CCF_CLEAN_BANDS'] = ['r']

    # Converted-array order index used for the diagnostic model plot.
    # This is an array index, not a physical echelle-order number.
    rparams['COMPUTE_MODEL_PLOT_ORDERS'] = [10]

    # ---------------------------------------------------------------------
    # Compilation settings
    # ---------------------------------------------------------------------

    # Keep compilation inside the selected instrument wavelength range
    rparams['COMPIL_WAVE_MIN'] = 497.522
    rparams['COMPIL_WAVE_MAX'] = 684.933

    rparams['COMPIL_MAX_PIXEL_WIDTH'] = 50
    rparams['COMPIL_CUT_PEARSONR'] = -1
    rparams['COMPIL_FP_EWID'] = 3.0

    rparams['COMPIL_ADD_UNIFORM_WAVEBIN'] = True
    rparams['COMPIL_NUM_UNIFORM_WAVEBIN'] = 15

    rparams['COMPILE_BINNED_BAND1'] = 'g'
    rparams['COMPILE_BINNED_BAND2'] = 'r'
    rparams['COMPILE_BINNED_BAND3'] = 'i'

    # Reference wavelength near the centre of the selected range
    rparams['COMPIL_SLOPE_REF_WAVE'] = 590.0

    # ---------------------------------------------------------------------
    # Detector, mask and blaze parameters
    # ---------------------------------------------------------------------

    rparams['READ_OUT_NOISE'] = 15.0
    rparams['MASK_SNR_MIN'] = 20.0

    rparams['BLAZE_SMOOTH_SIZE'] = 20.0
    rparams['BLAZE_THRESHOLD'] = 0.2

    rparams['BERVBIN_SIZE'] = 3000.0

    # ---------------------------------------------------------------------
    # Telluric-cleaning settings
    # ---------------------------------------------------------------------

    # Disable during initial compatibility testing.
    # The associated values must still be non-None.
    rparams['DO_TELLUCLEAN'] = False

    rparams['TELLUCLEAN_DV0'] = 0.0

    rparams['TELLUCLEAN_MASK_DOMAIN_LOWER'] = 550.0
    rparams['TELLUCLEAN_MASK_DOMAIN_UPPER'] = 670.0

    rparams['TELLUCLEAN_FORCE_AIRMASS'] = True
    rparams['TELLUCLEAN_CCF_SCAN_RANGE'] = 150.0
    rparams['TELLUCLEAN_MAX_ITERATIONS'] = 20

    rparams['TELLUCLEAN_KERNEL_WID'] = 1.4
    rparams['TELLUCLEAN_GAUSSIAN_SHAPE'] = 2.2

    # Telluric model grid restricted to the selected MaHPS order range
    rparams['TELLUCLEAN_WAVE_LOWER'] = 497.522
    rparams['TELLUCLEAN_WAVE_UPPER'] = 684.933

    rparams['TELLUCLEAN_TRANSMISSION_THRESHOLD'] = -1.0
    rparams['TELLUCLEAN_SIGMA_THRESHOLD'] = 10.0

    rparams['TELLUCLEAN_RECENTER_CCF'] = False
    rparams['TELLUCLEAN_RECENTER_CCF_FIT_OTHERS'] = True

    rparams['TELLUCLEAN_DEFAULT_WATER_ABSO'] = 5.0

    rparams['TELLUCLEAN_WATER_BOUNDS_LOWER'] = 0.05
    rparams['TELLUCLEAN_WATER_BOUNDS_UPPER'] = 15.0

    rparams['TELLUCLEAN_OTHERS_BOUNDS_LOWER'] = 0.05
    rparams['TELLUCLEAN_OTHERS_BOUNDS_UPPER'] = 15.0

    # ---------------------------------------------------------------------
    # Template-construction settings
    # ---------------------------------------------------------------------

    rparams['TEMPLATE_MEDBINMAX'] = 19
    rparams['MAX_CONVERGENCE_TEMPLATE_RV'] = 100.0
    
    # ---------------------------------------------------------------------
    # Parameters for the template construction
    # ---------------------------------------------------------------------
    # max number of bins for the median of the template. Avoids handling
    # too many spectra at once.
    # Example: 19
    # Type: INT
    rparams['TEMPLATE_MEDBINMAX'] = None

    # maximum RMS between the template and the median of the template
    # to accept the median of the template as a good template. If above
    # we iterate once more. Expressed in m/s
    # Example: 100
    # Type: FLOAT
    rparams['MAX_CONVERGENCE_TEMPLATE_RV'] = None
    
    # ---------------------------------------------------------------------
    # Keywords which must be in the header of science files
    # ---------------------------------------------------------------------
    
    # define the header key that gives the mid exposure time in MJD
    # Example: MJDMID
    rparams['KW_MID_EXP_TIME'] = "MJDMID"
    
    # define the header key that gives the start time of the observation
    # Example: MJSTART
    rparams['KW_MJDATE'] = "MJSTART"
    
    # define the header key that gives the barycentric julian date
    # Example: BJD
    rparams['KW_BJD'] = "BJD"
        
    # define the header key that gives human date of the observation 
    # yyyy-mm-dd HH:MM:SS
    # Example: DATE
    rparams['KW_DATE'] = "DATE"
    
    # define the header key that gives the berv (in km/s)
    # Example: BERV
    rparams['KW_BERV'] = "BERV"
    
    # define the header key that gives the original object name
    # Example: OBJECT
    rparams['KW_OBJNAME'] = "OBJECT"
    
    # define the header key that gives the exposure time of the observation
    # Example: EXPTIME
    rparams['KW_EXPTIME'] = "EXPTIME"
    
    # define the header key that gives the airmass of the observation 
    # (if value is unknown set to a sensible value)
    # Example: AIRMASS
    rparams['KW_AIRMASS'] = "AIRMASS"

    # define the header key that gives the snr
    # (if value is unknown set to a sensible value)
    # Example: SNR
    rparams['KW_SNR'] = "SNR"
    
    # define the header key that gives the SNR in chosen order 
    # (if value is unknown set to a sensible value)
    # Example: SNR
    rparams['KW_EXT_SNR'] = "EXT_SNR"
    
    # Approximate mean spectral resolving power, R = lambda / delta_lambda
    rparams['APPROX_RESOLUTION'] = 60000
    # -------------------------------------------------------------------------
    # Run the wrapper code using the above settings
    # -------------------------------------------------------------------------
    # run main
    lbl_wrap(rparams)

# =============================================================================
# End of code
# =============================================================================
