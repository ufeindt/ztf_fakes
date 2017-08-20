import numpy as np
import astropy
from astropy.io import fits
from astropy.table import Table
import logging
from ConfigParser import ConfigParser

class CalibratedImage:
    '''
    A calibrated ZTF/PTF image on the sky.
    '''
    __fits_cards__ = [
        'IMAGEZPT', 'GAIN', 'SEEING', 'CCDID', 'PTFFIELD',
        'READNOI', 'MOONILLF', 'AIRMASS', 'ELLIP',
        'MOONRA', 'MOONDEC', 'OBSMJD'
    ]
    ### new names of the SE table columns for ease of notation
    __working_params__ = [
        'id', 'flags', 'flux', 'mag', 'x', 'y', 'bkgnd', 'elgn',
        'ellip', 'A', 'B', 'theta', 'class_st'
    ]
    __working2SE__ = {
        'id'      : 'NUMBER',
        'flags'   : 'FLAGS',
        'flux'    : 'FLUX_AUTO',
        'mag'     : 'MAG_AUTO',
        'x'       : 'X_IMAGE',
        'y'       : 'Y_IMAGE',
        'bkgnd'   : 'BACKGROUND',
        'elgn'    : 'ELONGATION',
        'ellip'   : 'ELLIPTICITY',
        'A'       : 'A_IMAGE',
        'B'       : 'B_IMAGE',
        'theta'   : 'THETA_IMAGE',
        'class_st': 'CLASS_STAR',
        'zp'      : 'ZEROPOINT',
        'fwhm'    : 'FWHM_IMAGE'
        
    }
    ### SE parameters which are extracted from the SE catalog
    __SE_params__ = [
        'NUMBER', 'ZEROPOINT', 'FLAGS', 'FWHM_IMAGE', 'FLUX_AUTO', 'MAG_AUTO', 'X_IMAGE',
        'Y_IMAGE', 'BACKGROUND', 'ELONGATION', 'X2_IMAGE', 'Y2_IMAGE',
        'XY_IMAGE', 'THETA_IMAGE', 'A_IMAGE', 'B_IMAGE', 'CLASS_STAR'
    ]
    __SE2working__ = {
        'NUMBER'     : 'id',
        'ZEROPOINT'  : 'zp',
        'FLAGS'      : 'flags',
        'FWHM_IMAGE' : 'fwhm',
        'FLUX_AUTO'  : 'flux',
        'MAG_AUTO'   : 'mag',
        'X_IMAGE'    : 'x',
        'Y_IMAGE'    : 'y',
        'BACKGROUND' : 'bkgnd',
        'ELONGATION' : 'elgn',
        'X2_IMAGE'   : 'x2',
        'Y2_IMAGE'   : 'y2',
        'XY_IMAGE'   : 'xy',
        'THETA_IMAGE': 'theta',
        'A_IMAGE'    : 'A',
        'B_IMAGE'    : 'B',
        'CLASS_STAR' : 'class_st'
    }
    
    ###################### List of checks #########################
    _common_checks = [
        'flagCheck', 'edgeObjectCheck'
    ]

    _star_checks = [
        'starMagCheck', 'starEllipCheck', 'starElgnCheck', 
        'starClassStarCheck', 'starFWHMCheck'
    ]
    
    _gal_checks = [
        'galMagCheck', 'galEllipCheck', 'galElgnCheck',
        'galClassStarCheck', 'galFWHMCheck'
    ]
    
    ##################### Sanity checks ###########################
    assert len(__SE2working__) == len(__SE_params__)
    
    def __init__(self, fits_file, SE_cat=None, config=None, logger=None):
        self.fits_file    = fits_file # FITS filename
        self.SE_cat       = SE_cat    # SE catalog associated with fits file
        self.config       = config    # ConfigParser object
        self.logger       = logger    # Logger object
        self.io_tag       = None      # set to true of I/O succeeds
        self.img_quality  = None      # set to true if all checks succeed
        
        self.img_header   = None      # image header
        self.img_data     = None      # image data
        
        ### Load image using astropy fits ###
        self._loadFits()

        ### Extract metadata listed in __fits_cards__
        self._getMetaInfo()
        
        ### Parse SE catalog
        self._parseSECatalog()
    
    def _loadFits(self):
        '''
        Loads supplied fits image using astropy fits
        '''
        try:
            hdu_list = fits.open(self.fits_file)
        except IOError:
            self.io_tag = False
            if self.logger:
                self.logger.error('Error opening fits image %s'%(self.fits_file))
        else:
            self.io_tag = True
            self.img_header = hdu_list[0].header
            self.img_data   = hdu_list[0].data
            
    def _getMetaInfo(self):
        '''
        Extracts relevant meta information from image header and puts
        them in a meta dictionary
        '''
        assert self.io_tag
        
        try:
            self.meta = {key:self.img_header[key] for key in CalibratedImage.__fits_cards__}
        except KeyError:
            self.img_quality = False
            if self.logger:
                self.logger.error('Invalid fits card in %s'%(self.fits_file))
    
    def _parseSECatalog(self):
        '''
        Parse necessary SE catalog information in astropy table format
        '''
        try:
            tab = Table.read(self.SE_cat)
        except:
            if self.logger:
                self.logger.error('Error loading SE catalog for %s'%(self.fits_file))
        else:
            self.SE_info_table = tab[CalibratedImage.__SE_params__]
            ### Rename column names to ease usage
            for colname in self.SE_info_table.colnames:
                self.SE_info_table.rename_column(colname, CalibratedImage.__SE2working__[colname])
            
    def recordLog(foo):
        '''
        Decorator to add a log message
        '''
        check_name = foo.__name__
        def writeLog(self):
            if self.logger: ###FIXME: change logging formatting
                self.logger.info('Attempting %s in %s'%(check_name, self.fits_file))
            
            res = foo(self)
            if res is None and self.logger:
                self.logger.warn('Returned None for %s in %s'%(check_name, self.fits_file))
            return res
        return writeLog
    

    ########################################################
    ### Following methods perform checks to filter out 
    ### GALAXIES
    ########################################################
    @recordLog
    def galMagCheck(self):
        '''
        Returns boolean array after applying magnitude cuts
        '''
        assert self.SE_info_table is not None
        
        if not self.config.has_section('galaxy_checks') : return None
        
        # Lower value of the check
        minval = self.config.get('galaxy_checks','min_mag')
        # Upper value of the check
        maxval = self.config.get('galaxy_checks', 'max_mag')
        # Resulting boolean array
        res    = np.ones(len(self.SE_info_table)).astype(bool)
        # Apply filters
        if minval:
            res *= self.SE_info_table['mag'] + self.SE_info_table['zp'] >= float(minval)
        if maxval:
            res *= self.SE_info_table['mag'] + self.SE_info_table['zp'] <= float(maxval)
        
        return res
    
    @recordLog
    def galEllipCheck(self):
        '''
        Returns boolean array after ellipticity cuts
        '''
        assert self.SE_info_table is not None
        
        if not self.config.has_section('galaxy_checks') : return None
        
        # Lower value of the check
        minval = self.config.get('galaxy_checks','min_ellip')
        # Upper value of the check
        maxval = self.config.get('galaxy_checks', 'max_ellip')
        # Resulting boolean array
        # start with all True
        res    = np.ones(len(self.SE_info_table)).astype(bool)
        # Compute ellipticity if not provided
        ellip  = self.SE_info_table['ellip'] if 'ellip' in self.SE_info_table.colnames \
                    else 1.0 - self.SE_info_table['B']/self.SE_info_table['A']
        # Apply filters
        if minval:
            res *= ellip >= float(minval)
        if maxval:
            res *= ellip <= float(maxval)
        
        return res
    
    @recordLog
    def galElgnCheck(self):
        '''
        Returns boolean array after elongation cuts
        '''
        assert self.SE_info_table is not None
        
        if not self.config.has_section('galaxy_checks') : return None

        # Lower value of the check
        minval = self.config.get('galaxy_checks','min_ellip')
        # Upper value of the check
        maxval = self.config.get('galaxy_checks', 'max_ellip')
        # Resulting boolean array
        # star with all True
        res    = np.ones(len(self.SE_info_table)).astype(bool)
        # Compute ellipticity if not provided
        elgn   = self.SE_info_table['elgn'] if 'elgn' in self.SE_info_table.colnames \
                    else self.SE_info_table['A']/self.SE_info_table['B']
        # Apply filters
        if minval:
            res *= elgn >= float(minval)
        if maxval:
            res *= elgn <= float(maxval)
        
        return res
        
    @recordLog
    def galClassStarCheck(self):
        '''
        Returns boolean array with class star limit cut
        '''
        assert self.SE_info_table is not None
        
        if not self.config.has_section('galaxy_checks') : return None
        
        # Lower value of the check
        minval = self.config.get('galaxy_checks','min_class_star')
        # Upper value of the check
        maxval = self.config.get('galaxy_checks', 'max_class_star')
        # Resulting boolean array
        res    = np.ones(len(self.SE_info_table)).astype(bool)
        # Apply filters
        if minval:
            res *= self.SE_info_table['class_st'] >= float(minval)
        if maxval:
            res *= self.SE_info_table['class_st'] <= float(maxval)
        
        return res
    
    @recordLog
    def galFWHMCheck(self):
        '''
        Returns boolean array with FWHM cut
        '''
        assert self.SE_info_table is not None
        
        if not self.config.has_section('galaxy_checks') : return None
        
        # Lower value of the check
        minval = self.config.get('galaxy_checks','min_FWHM')
        # Upper value of the check
        maxval = self.config.get('galaxy_checks', 'max_FWHM')
        # Resulting boolean array
        res    = np.ones(len(self.SE_info_table)).astype(bool)
        # Apply filters
        if minval:
            res *= self.SE_info_table['fwhm'] >= float(minval)
        if maxval:
            res *= self.SE_info_table['fwhm'] <= float(maxval)
        
        return res
    
    ########################################################
    ### Following methods perform checks to filter out STARS
    ########################################################
    @recordLog
    def starMagCheck(self):
        '''
        Returns boolean array after applying magnitude cuts
        '''
        assert self.SE_info_table is not None
        
        if not self.config.has_section('star_checks') : return None
        
        # Lower value of the check
        minval = self.config.get('star_checks','min_mag')
        # Upper value of the check
        maxval = self.config.get('star_checks', 'max_mag')
        # Resulting boolean array
        res    = np.ones(len(self.SE_info_table)).astype(bool)
        # Apply filters
        if minval:
            res *= self.SE_info_table['mag'] + self.SE_info_table['zp'] >= float(minval)
        if maxval:
            res *= self.SE_info_table['mag'] + self.SE_info_table['zp'] <= float(maxval)
        
        return res
    
    @recordLog
    def starEllipCheck(self):
        '''
        Returns boolean array after ellipticity cuts
        '''
        assert self.SE_info_table is not None
        
        if not self.config.has_section('star_checks') : return None
        
        # Lower value of the check
        minval = self.config.get('star_checks','min_ellip')
        # Upper value of the check
        maxval = self.config.get('star_checks', 'max_ellip')
        # Resulting boolean array
        # start with all True
        res    = np.ones(len(self.SE_info_table)).astype(bool)
        # Compute ellipticity if not provided
        ellip  = self.SE_info_table['ellip'] if 'ellip' in self.SE_info_table.colnames \
                    else ( 1.0 - self.SE_info_table['B']/self.SE_info_table['A'] )
        # Apply filters
        if minval:
            res *= ellip >= float(minval)
        if maxval:
            res *= ellip <= float(maxval)
        
        return res
    
    @recordLog
    def starElgnCheck(self):
        '''
        Returns boolean array after elongation cuts
        '''
        assert self.SE_info_table is not None
        
        if not self.config.has_section('star_checks') : return None

        # Lower value of the check
        minval = self.config.get('star_checks','min_elgn')
        # Upper value of the check
        maxval = self.config.get('star_checks', 'max_elgn')
        # Resulting boolean array
        # star with all True
        res    = np.ones(len(self.SE_info_table)).astype(bool)
        # Compute ellipticity if not provided
        elgn   = self.SE_info_table['elgn'] if 'elgn' in self.SE_info_table.colnames \
                    else self.SE_info_table['A']/self.SE_info_table['B']
        # Apply filters
        if minval:
            res *= ( elgn >= float(minval) )
        if maxval:
            res *= ( elgn <= float(maxval) )
        
        return res
        
    @recordLog
    def starClassStarCheck(self):
        '''
        Returns boolean array with class star limit cut
        '''
        assert self.SE_info_table is not None
        
        if not self.config.has_section('star_checks') : return None
        
        # Lower value of the check
        minval = self.config.get('star_checks','min_class_star')
        # Upper value of the check
        maxval = self.config.get('star_checks', 'max_class_star')
        # Resulting boolean array
        res    = np.ones(len(self.SE_info_table)).astype(bool)
        # Apply filters
        if minval:
            res *= self.SE_info_table['class_st'] >= float(minval)
        if maxval:
            res *= self.SE_info_table['class_st'] <= float(maxval)
        
        return res
    
    @recordLog
    def starFWHMCheck(self):
        '''
        Returns boolean array with FWHM cut
        '''
        assert self.SE_info_table is not None
        
        if not self.config.has_section('star_checks') : return None
        
        # Lower value of the check
        minval = self.config.get('star_checks','min_FWHM')
        # Upper value of the check
        maxval = self.config.get('star_checks', 'max_FWHM')
        # Resulting boolean array
        res    = np.ones(len(self.SE_info_table)).astype(bool)
        # Apply filters
        if minval:
            res *= self.SE_info_table['fwhm'] >= float(minval)
        if maxval:
            res *= self.SE_info_table['fwhm'] <= float(maxval)
        
        return res
    
    ##########################################################
    ### Following checks perform checks which can be common to
    ### both stars and galaxies
    ##########################################################
    
    @recordLog
    def flagCheck(self):
        '''
        Return boolean array for object satisfying good flag
        '''
        assert self.SE_info_table is not None
        
        if not self.config.has_section('common_checks'): return None
        ### Resulting boolean array
        res = np.ones(len(self.SE_info_table)).astype(bool)
        allowed_flag = self.config.get('common_checks', 'allowed_flag')
        if allowed_flag:
            res *= ( self.SE_info_table['flags'] == int(allowed_flag) )
        
        return res
    
    @recordLog
    def edgeObjectCheck(self):
        '''
        Return boolean depending on objects NOT being at the edge of the field
        So Trues are objects which are not at the edge
        '''
        assert self.SE_info_table is not None
        
        if not self.config.has_section('common_checks'): return None
        ### Resulting boolean array
        res = np.ones(len(self.SE_info_table)).astype(bool)
        ### Read dimension of image
        x_dim = self.img_header['NAXIS1']
        y_dim = self.img_header['NAXIS2']
        
        edge_lim = self.config.get('common_checks', 'edge_limit')
        if edge_lim:
            x_vals   = self.SE_info_table['x']
            y_vals   = self.SE_info_table['y']
            edge_lim = int(edge_lim)
            res *= ( x_vals >= edge_lim ) * ( (x_dim - x_vals) >= edge_lim ) * \
                ( y_vals >= edge_lim ) * ( (y_dim - y_vals) >= edge_lim )
        
        return res
    
    @recordLog
    def enoughObjectsCheck(self):
        '''
        Return boolean VALUE depending on SE finding a minimum
        number of objects in the image
        '''
        assert self.SE_info_table is not None
        
        if not self.config.has_section('common_checks'): return None
        
        min_num_obj = self.config.get('common_checks', 'min_num_obj')
        if min_num_obj:
            min_num_obj = int(min_num_obj)
            return len(self.SE_info_table) >= min_num_obj
            ###FIXME: Could be modified to only flag == 0 objs
        return True
    
    def _performCommonChecks(self):
        '''
        Perform checks in _common_checks
        '''
        res = np.ones(len(self.SE_info_table)).astype('bool')
        for check in CalibratedImage._common_checks:
            ###FIXME: Could be fragile
            eval_str = 'self.'+check+'()'
            res *= eval(eval_str)
        return res
    
    def _performStarChecks(self):
        '''
        Perform checks in _star_checks
        '''
        res = np.ones(len(self.SE_info_table)).astype('bool')
        for check in CalibratedImage._common_checks:
            eval_str = 'self.'+check+'()'
            res *= eval(eval_str)
        for check in CalibratedImage._star_checks:
            ###FIXME: Could be fragile
            eval_str = 'self.'+check+'()'
            res *= eval(eval_str)
        return res
    
    def _performGalChecks(self):
        '''
        Perform checks in _gal_checks
        '''
        res = np.ones(len(self.SE_info_table)).astype('bool')
        for check in CalibratedImage._common_checks:
            eval_str = 'self.'+check+'()'
            res *= eval(eval_str)
        for check in CalibratedImage._gal_checks:
            ###FIXME: Could be fragile
            eval_str = 'self.'+check+'()'
            res *= eval(eval_str)
        return res
    
    def giveStarsTable(self):
        '''
        Return astropy table with potential star like objects
        '''
        assert self.SE_info_table is not None
        
        star_mask = self._performStarChecks()
        return self._giveTable(star_mask)
    
    def giveGalaxiesTable(self):
        '''
        Return astropy table with potential galaxy like objects
        '''
        assert self.SE_info_table is not None
        
        gal_mask = self._performGalChecks()
        return self._giveTable(gal_mask)
    
    def _giveTable(self, mask):
        '''
        Filters the astropy table based on the mask supplied
        '''
        assert mask.dtype == bool and len(mask) == len(self.SE_info_table)
        
        return self.SE_info_table[mask]
    
    def runSE(self, catalog_dir=None):
        '''
        Run source extractor on the image and store it to a catalog
        '''
        raise NotImplementedError           


class Injection(CalibratedImage):
    '''
    A fake image.
    '''
    name = "Injection"
        
    def __init__(self, fits_file, SE_cat,
                 num_inj=1,  
                 config=None, 
                 logger=None
                ):
        
        CalibratedImage.__init__(self,fits_file, SE_cat, config=config, logger=logger)
        self.num_inj = num_inj
        self._addFakeCard()
        
        
    def _sanity_checks(self):
        '''
        Try the sanity checks
        '''
        pass
    
    def _addFakeCard(self):
        '''
        Adds a header to the image to denote fake image
        '''
        if not self.config.has_section('fake_card_info'): return None
        
        name   = self.config.get('fake_card_info','card_name')
        value  = self.config.get('fake_card_info','card_value')
        comment=self.config.get('fake_card_info','card_comment')
        
        insert_after_card=self.config.get('fake_card_info','insert_after_card')
        ### Avoid multiple insertion of the fake card
        if name.upper() not in self.img_header.keys():
            card = fits.Card(name, value, comment)
            self.img_header.insert(insert_after_card, card)


class CloneStampedInjection(Injection):
    '''
    Performs injections according to clone stamping method.
    Brightest stars in the randomly chosen galaxies
    '''
    name = "CloneStampedInjection"
    
    def __init__(self, fits_file, SE_cat, 
                 num_inj,
                 config,
                 logger=None
                ):
        Injection.__init__(self, fits_file, SE_cat, num_inj, config, logger)
        
        ### presence of clone stamping section is required
        #assert self.config.has_section['clone_stamped_injection']
        
        ### read the dimensions of the stamping region
        self.delta_x = self.config.getint('clone_stamped_injection','delta_x')
        self.delta_y = self.config.getint('clone_stamped_injection','delta_y')
        ### make changes to the original image?
        self.stamp_original = self.config.getboolean('clone_stamped_injection','stamp_original')
        
    def getBrightStars(self):
        '''
        Returns astropy table with num_inj brightest stars, sorted based on flux
        '''
        star_table     = self.giveStarsTable()
        sorted_on_flux = star_table[np.argsort(-star_table['flux'])]
        return sorted_on_flux[:self.num_inj]
    
    def getBrightGalaxies(self):
        '''
        Returns astropy table with num_inj bright galaxies from the galaxy list
        '''
        gal_table      = self.giveGalaxiesTable()
        sorted_on_flux = gal_table[np.argsort(-gal_table['flux'])]
        return sorted_on_flux[:self.num_inj]

    def getDimStars(self):
        '''
        Returns astropy table with num_inj dimmest stars, sorted based on flux
        '''
        star_table     = self.giveStarsTable()
        sorted_on_flux = star_table[np.argsort(star_table['flux'])]
        return sorted_on_flux[:self.num_inj]
    
    def getDimGalaxies(self):
        '''
        Returns astropy table with num_inj dim galaxies from the galaxy list
        '''
        gal_table      = self.giveGalaxiesTable()
        sorted_on_flux = gal_table[np.argsort(gal_table['flux'])]
        return sorted_on_flux[:self.num_inj]
    
    def stampBrightsInDims(self):
        '''
        Stamp dim galaxies with bright stars
        '''
        stars    = self.getBrightStars()
        galaxies = self.getDimGalaxies()
        ### same number of stars and galaxies
        assert len(stars) == len(galaxies)
        ### create a copy of the image data
        duplicate_image = np.copy(self.img_data)
        ### randomize any one set, say stars, cast into Table object
        stars = Table(np.random.choice(stars, self.num_inj))
        ### make a list of offsets
        ###FIXME: offset choice currently random. Need rationale
        offset_list = np.arange(-5,6)
        ### Loop over and stamp
        for idx, star in enumerate(stars):
            galaxy = galaxies[idx]
            duplicate_image = \
            CloneStampedInjection._stamp(duplicate_image, star, galaxy,\
                                                           delta_x=self.delta_x,\
                                                           delta_y=self.delta_y,\
                                                           scale_fac=1.0,\
                                                           offset_x=np.random.choice(offset_list),\
                                                           offset_y=np.random.choice(offset_list)\
                                                          )
        
        if self.stamp_original:
            self.img_data = duplicate_image
        
        return duplicate_image
    
    @staticmethod
    def _stamp(image, source, target, \
                   delta_x=20,\
                   delta_y=20,\
                   scale_fac=1.0,\
                   offset_x=1,\
                   offset_y=1):
        '''
        FUNCTION :: Stamps source location into target location.
        INPUTS   :: image    - fits image data
                    source   - astropy row source to stamp
                    target   - astropy row where to stamp
                    scale_fac- factor by which to blow source up
                    offset_x - offset from targetx while stamping
                    offset_y - offset from targety while stamping
        OUTPUT   :: numpy 2D array of new image
        WARNING  :: Built ad-hoc, use with care.
        '''
        assert isinstance(source, astropy.table.row.Row) and \
                isinstance(target, astropy.table.row.Row)
        
        ### create a duplicate image
        #duplicate_image = np.copy(self.img_data)
        ### source patch
        source_x_lo = int( source['x'] - delta_x/2.0 )
        source_x_hi = int( source['x'] + delta_x/2.0 )
        source_y_lo = int( source['y'] - delta_y/2.0 )
        source_y_hi = int( source['y'] + delta_y/2.0 )
        ### target patch
        target_x_lo = int( target['x'] - delta_x/2.0 ) + offset_x
        target_x_hi = int( target['x'] + delta_x/2.0 ) + offset_x
        target_y_lo = int( target['y'] - delta_y/2.0 ) + offset_y
        target_y_hi = int( target['y'] + delta_y/2.0 ) + offset_y
        
        patch = ( image[source_y_lo:source_y_hi, source_x_lo:source_x_hi] \
            - source['bkgnd'] ) * scale_fac
        image[target_y_lo:target_y_hi, target_x_lo:target_x_hi] += patch
        
        return image
