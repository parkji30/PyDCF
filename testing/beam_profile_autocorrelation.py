import numpy as np

## Set data
nbeam = astrB.NAXIS[0]

n1 = 2*nbeam
beam_vis = np.zeros((n1,n1))
beam_tf = np.zeros((n1,n1))
gaus_vis = np.zeros((n1,n1))
gaus_tf = np.zeros((n1,n1))

pix_res = 0.30 #resoltuion in arcsec for BLAST-TNG

## Compute beam profile
b = 0
sigx1 =
sigy1 =
sigx2 =
sigy2 =
theta1 =
theta2 =

tot_gaus = 0

for i in range(n1):
    for j in range(n1):
        # Compute the synthetic Gaussian beam
        gaus_vis[i,j] = np.exp(-0.5*((j-n1/2)**2 * ((np.cos(theta1) / sigx1)**2 + (np.sin(theta1)/sigy1)**2) \
                        + (i - n1/2)**2 * (np.(cos(theta1) / sigy1)**2 + (np.sin(theta1) / sigx1)**2) \
                        -2.0 * (i - n1/2) * (j - n1/2) * np.sin(theta1) * np.cos(theta1) * (1.0 / sigx1**2 - 1.0 / sigy1**2))) \
                        / (2.0 * np.pi * sigx1 * sigy1)

        if not keyword_set(sdish):
        # Subtract a wider Gaussian beam to mimc the spatial filtering
            tot_gaus += gaus_vis[i,j]
            gaus_vis[i,j] -= np.exp(-0.5 * ((j - n1/2)**2 * ((np.cos(theta2) / sigx2)**2 + (np.sin(theta2) / sigy2)**2) \
                            + (i - n1/2)**2 * ((np.cos(theta2) / sigy2)**2 + (np.sin(theta2) / sigx2)**2) \
                            -2.0 * (i - n1/2) * (j - n1/2) * np.sin(theta2) * np.cos(theta2) * (1.0 / sigx2**2 - 1.0 / sigy2**2))) \
                            / (2.0 * np.pi * sigx2 * sigy2)
        else:
            beam_vis[i,j] = gaus_vis[i,j]     # in lieu of the dirty beam for an interferometer
        b++


if keyword_set(sdish):  # for single dish make sure beams are normalized
    beam_vis = beam_vis / np.sum(beam_vis)
    gaus_vis = gaus_vis / np.sum(gaus_vis)
else:
    print('tot_gaus, total(gaus_vis) = ', tot_gaus, np.sum(gaus_vis))
    gaus_vis = gaus_vis / tot_gaus


print('b, total(beam_vis), total(gaus_vis) :', b, np.sum(beam_vis), np.sum(gaus_vis))


if not keyword_set(sdish):
    print('b, total(imageB): ', b, np.sum(imageB))  # returns the sum of the elements in beam_vis array

if b != n1**2:
    print('Error!  - incompatible dimensions calculated. Quitting')
    exit #Stop Program.


## Calculate the wrapped version of the beam

for i in range(n1/2):
    for j in range(n1/2):
        beam_tf[i, j] = beam_vis[i + n1/2, j + n1/2]
        beam_tf[i + n1/2, j] = beam_vis[i, j + n1/2]
        beam_tf[i, j + n1/2] = beam_vis[i + n1/2, j]
        beam_tf[i + n1/2, j + n1/2] = beam_vis[i, j]

        gaus_tf[i, j] = gaus_vis[i + n1/2, j + n1/2]
        gaus_tf[i + n1/2, j] = gaus_vis[i, j + n1/2]
        gaus_tf[i, j + n1/2] = gaus_vis[i + n1/2, j]
        gaus_tf[i + n1/2, j + n1/2] = gaus_vis[i, j]


## Calculate the beam's autocorrelation
beam_fft = np.fft.fft(beam_vis) * n1 ** 2
beam_fft_tf = np.fft.fft(beam_tf) * n1**2 # wrapped version

correlation_beam = np.real(np.fft.fft(beam_fft * np.conjugate(beam_fft_tf),1))/n1**2  # autocorrelation
print('total(corr_beam), max(corr_beam): ', np.sum(corr_beam), np.max(corr_beam))


gaussian_fft = np.fft.fft(gaus_vis) * n1**2
gaussian_fft_tf = np.fft.fft(gaus_tf) * n1**2

correlation_gaussian = np.real(np.fft.fft(gaussian_fft * np.conjugate(gaussian_fft_tf)))/n1**2 # autocorrelation
print('total(corr_gaus), max(corr_gaus): ', np.sum(correlation_gaussian), max(correlation_gaussian))


## Construct the beam's mean Radial Profile

beam_x = 0
beam_y = 0
delta_r = np.zeros(n1**2)
beam = np.zeros(n1^2)
gbeam = np.zeros(n1^2)

c = 0
for i in range(n1):
    for j in range(n1):
        beam_x = (i - n1/2) * pix_res
        beam_y = (j - n1/2) * pix_res
        delta_r[c] = np.sqrt(beam_x**2 + beam_y**2)
        beam[c] = corr_beam[i, j]
        gbeam[c] = corr_gaus[i, j]

    if c != n1**2:
        print('Error!  - incompatible dimensions calculated. Quitting')
        exit #Stop Program.


print('Max delta_r: ', np.max(delta_r))
print('Min delta_r: ', np.min(delta_r))


bins = 3 * (n1/2) + 1
minb = -0.5 * pix_res
maxb = (2 * nb - 3) * np.abs(minb) # '-3' ensures that the binsize is pix_res

# CHECK THIS OVER
tmp = np.histogram(delta_r, max=maxb, min=minb, locations=ell, nbins=bins, reverse_ind = ind, omin=lmin, omax=lmax)
print('Bins, Minimum and Maximum ell bins = ', bins, lmin, lmax)


mean_beam = np.zeros(n_elements(ell))
mean_gaus = np.zeros(n_elements(ell))


for i in range(len(ell)):
    if ind[i] != ind[i+1]:
        indices = ind[ind[i] : ind[i+1] - 1]
        mean_beam[i] = np.mean(beam[indices])
        mean_gaus[i] = np.mean(gbeam[indices])


mean_beam /= np.max(mean_beam)
print('ell :', ell)
print('mean_beam :', mean_beam/mean_beam[0])
print('total mean_beam : ', total(mean_beam))
