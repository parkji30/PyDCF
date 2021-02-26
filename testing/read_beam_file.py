
  ;; Read beamfile
  ;;===============

  openr, lun, beamfile, /get_lun

  readf, lun, sigx1, sigy1, theta1, coeff_sigma1    ; get the parameters of the synthetic Gaussian beam
  readf, lun, sigx2, sigy2, theta2, coeff_sigma2    ; get the parameters of the synthetic Gaussian beam
  readf, lun, fitsB                                 ; for interferrometer data get the fits file contaning the dirty beam

  sigx1=sigx1*coeff_sigma1
  sigy1=sigy1*coeff_sigma1
  sigx1 /= sqrt(8.*alog(2.))   ; beamwidth in arcsec
  sigy1 /= sqrt(8.*alog(2.))   ; beamwidth in arcsec
  W1 = sqrt(sigx1*sigy1)

  sigx2=sigx2*coeff_sigma2
  sigy2=sigy2*coeff_sigma2
  sigx2 /= sqrt(8.*alog(2.))   ; beamwidth in arcsec
  sigy2 /= sqrt(8.*alog(2.))   ; beamwidth in arcsec
  W2 = sqrt(sigx2*sigy2)

  free_lun, lun

  ; All in arcmin

  sigx1 /= 60.
  sigy1 /= 60.
  sigx2 /= 60.
  sigy2 /= 60.
  W1 /= 60.
  W2 /= 60.
  Dprime /= 60.

  print, 'W1, W2 = ', W1, W2
  print, 'Dprime = ', Dprime


  ;; Initialize parameters of "interferometer" fit
  ;;========================================

  name = ['delta','bratio2' ,'coeff2', 'coeff4', 'coeff6', 'coeff8', 'W1', 'W2', 'Dprime']

  np = n_elements(name)
  pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D], parname:[' ']}, np)
  ; parameter information
  ; Returns an vector of dimension np, filled with {fixed:0, limited:[0,0], limits:[0.D,0.D], parname:[' ']}.

  for i = 0, np-1 do begin
      pi[i].parname = name[i]
  endfor

  st =    [1., 0.1, 0.0, 0.0, 0.0, 0.0, W1,  W2, Dprime]
  fixed = [ 0, pf0, pf1, pf2, pf3, pf4,  1, pf5, 1]
  limited = [ $
            [1,0], $                      ; delta
            [1,0], $                      ; bratio2 = Bt^2/B0^2
            [1,0], $                      ; coeff2
            [0,0], $                      ; coeff4
            [0,0], $                      ; coeff6
            [0,0], $                      ; coeff8
            [1,0], $                      ; W1
            [1,0], $                      ; W2
            [1,0] ]                       ; Dprime
  limits =[ $
          [0.,10.], $                 ; delta
          [0.,10.], $                 ; bratio2
          [0.,0.], $                  ; coeff2
          [0.,0.], $                  ; coeff4
          [0.,0.], $                  ; coeff6
          [0.,0.], $                  ; coeff8
          [0.,0.], $                  ; W1
          [2.*W1,0.], $               ; W2
          [0.,0.] ]                   ; Dprime

  start = st
  pi.limits = limits
  pi.limited = limited
  pi.fixed = fixed
