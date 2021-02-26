import numpy as np
from astropy.io import fits

## Setup
num_data = fits.open("")[0].data


## Calculate the number of points needed for the turbulent power spectra
num_points = 4 * np.floor(num_data/2)) # Want even
print("Number of points: ", num_points)

b = 0

if not keyword_set(sdish):
    ra_tmp = 0
    dec_tmp = 0
    vis_tmp = 0

    # ===========================
    # CALCULATE BEAM PROFILE
    # (i.e. visibility)
    # ===========================

    # Dirty Beam File Name
    fits_B = '/Users/houde/Documents/IDL/Default/Dispersion/BLASTPol/'+fitsB
    print('Dirty beam input file = ', fitsB)


    # Load Dirty Beam File.
    fits_b_file = fits.open(fits_B)[0].data
    # fits_open, fitsB, fcb


    # What is this block doing?
    fits_read, fcb, imageB, hdrB, exten_no = 0
    epochI = sxpar(hdrB, 'EQUINOX')
    help, imageB
    fits_close, fcb

    # What is dimenB
    check_fits, imageB, hdrB, dimenB;,/update
    extast, hdrB, astrB

    corners = np.array([0, dimenB[0]-1, 0, dimenB[1]-1])
    print("corners: ", corners)

    dx = dimenB[0]
    dy = dimenB[1]

    data = replicate({x:0., y:0., ra:0.d, dec:0.d}, dx, dy)
    xr = findgen(dx)
    yc = findgen(dy)
    data.x = xr # (yc*0+1)
    data.y = (xr*0 + 1) # yc
    xy2ad, data.x, data.y, astrB, tmpra, tmpdec
    data.ra = tmpra-astrB.CRVAL[0]
    data.dec = tmpdec-astrB.CRVAL[1]

    old_pix_res = pix_res
    nbeam = astrB.NAXIS[0]
    pix_res = abs(astrB.CDELT[0]) * 3600. ; in arcsec
    if abs(pix_res-old_pix_res)/pix_res gt 1.e-3 then begin
      print
      print, '****** pix_res or beam not equal pix_res of data -- abort ******'
      print
      stop
    endif

    ;  imageB -= np.mean(imageB)

    n1 = 2*nbeam
    beam_vis = fltarr(n1,n1)

    ictr = long(astrB.CRPIX[0])
    jctr = long(astrB.CRPIX[1])
    for i = 0, nbeam-1 do begin
      for j = 0, nbeam-1 do begin
        beam_vis[i+n1/2-ictr,j+n1/2-jctr] = imageB[i,j] ;beam centred at n1/2
      endfor
    endfor
  endif else begin
    n1 = npts
    beam_vis = fltarr(n1,n1)
  endelse
