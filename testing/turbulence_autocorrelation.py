 ;; Plot turbulence autocorrelation function and Gaussian fit
  ;;===========================================================

  if not keyword_set(nospec) then begin

  xmin = kpix[0]-0.2*dk
  xmax = srange[1]
  ymin = -0.1*max(kplot[1:krange[1]])
;  ymin = min(kplot[1:krange[1]]) - 0.1*max(kplot[1:krange[1]])
  ymax = 1.1*max(kplot[1:krange[1]])

  nz = indgen(n_elements(kpix))
  if keyword_set(g_wiener) then $
    nz = where(kpix ge srange[0] and kpix le srange[1] and kpix ne 0. and kplot ne 0. and abs(h2) ge b_cut and sig_noi ge snr_cut)

  ploterror, kpix[nz], kplot[nz], stdk_wiener[nz], psym=4, $    ;plot of deconvolved power spectrum tspec
        xrange = [xmin,xmax], /xsty,  xmarg=[10,8], $
        yrange = [ymin,ymax], ysty=9, ymarg=[3.5,2], $
         xtit='k/(2!4p!x) (arcmin!e-1!n)', ytit = '!13R!3!lt!n(k)/<!sB!r!a!3-!n!3!e2!n>'
  oplot, [-10,100], [0,0], /linesty
  oplot, [0,0], [-1,1000],/linesty

  ; The next few lines are for shading the low-drequency filtered part of the spectrum

  if not keyword_set(sdish) then begin
    y_shade[*] = ymax
    y_shade[0] = ymin
    y_shade[n_elements(kpix_shade)-1] = ymin
    polyfill, kpix_shade, y_shade, /line_fill, spacing=0.12, orientation=45.
  endif

; if keyword_set(late) then begin
;   if gplot eq 1 then oplot, smx, fitkol, linesty=3
; endif