import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from photutils.aperture import CircularAperture, aperture_photometry
from photutils.detection import DAOStarFinder
from astropy.stats import mad_std
from astropy.coordinates import SkyCoord, AltAz
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import EarthLocation
from scipy.stats import linregress
import io

# 관측소 위치 설정
location = EarthLocation(lat=31.9584 * u.deg, lon=-111.5967 * u.deg, height=2120 * u.m)

st.title("FITS 이미지 분석 웹앱")
st.markdown("FITS 이미지 업로드 후, 별 탐지 및 밝기 분석을 수행합니다.")

uploaded_file = st.file_uploader("FITS 파일 업로드", type=["fits", "fz"])
if uploaded_file is not None:
    hdul = fits.open(uploaded_file)
    image_data = hdul[0].data
    header = hdul[0].header
    hdul.close()

    st.subheader("1. 이미지 시각화")
    fig, ax = plt.subplots()
    ax.imshow(image_data, cmap='gray', origin='lower', 
              vmin=np.percentile(image_data, 5), 
              vmax=np.percentile(image_data, 99))
    st.pyplot(fig)

    st.subheader("2. 별 탐지 및 밝기 측정")
    sigma = mad_std(image_data)
    daofind = DAOStarFinder(fwhm=3.0, threshold=5.0 * sigma)
    sources = daofind(image_data)
    
    if sources is not None:
        positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
        apertures = CircularAperture(positions, r=5.)
        phot_table = aperture_photometry(image_data, apertures)

        for i in range(min(5, len(phot_table))):
            flux = phot_table[i]['aperture_sum']
            x = phot_table[i]['xcenter']
            y = phot_table[i]['ycenter']
            st.write(f"별 {i+1}: 위치=({x:.1f}, {y:.1f}), 밝기(Flux)={flux:.2f}")

    else:
        st.warning("별을 탐지하지 못했습니다.")

    st.subheader("3. Flux vs. Magnitude 회귀 분석")
    flux = np.array([2500, 4300, 1600, 3100, 5400])
    mag = np.array([14.7, 13.5, 15.3, 14.1, 13.2])
    log_flux = np.log10(flux)
    slope, intercept, r_value, _, _ = linregress(log_flux, mag)
    
    st.write(f"회귀식: mag = {slope:.2f} * log10(flux) + {intercept:.2f}")
    st.write(f"R² = {r_value**2:.4f}")
    
    fig2, ax2 = plt.subplots()
    ax2.scatter(log_flux, mag, label='Data')
    ax2.plot(log_flux, slope * log_flux + intercept, label='Fit', color='red')
    ax2.set_xlabel('log10(Flux)')
    ax2.set_ylabel('Magnitude')
    ax2.invert_yaxis()
    ax2.legend()
    ax2.grid(True)
    st.pyplot(fig2)

st.markdown("출처: [GPTOnline.ai](https://gptonline.ai/ko/) - AI와 함께 천문 데이터 분석을 시작하세요!")
