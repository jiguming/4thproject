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

location = EarthLocation(lat=31.9584 * u.deg, lon=-111.5967 * u.deg, height=2120 * u.m)

st.title("천문 이미지 분석 앱 - GPTOnline.ai")
st.markdown("FITS 파일 3개를 업로드하여 별 밝기 측정, 회귀 분석 및 대기 소광 계수를 계산합니다.")

# -------- 1. 첫 번째 파일 (단일 분석용) --------
st.header("1. 별 탐지 및 회귀 분석용 FITS 파일")
file1 = st.file_uploader("분석용 FITS 파일 업로드", type=["fits"], key="file1")

if file1:
    hdul = fits.open(file1)
    image_data = None
    for hdu in hdul:
        if hdu.data is not None:
            image_data = hdu.data
            break
    hdul.close()

    if image_data is None:
        st.error("이미지 데이터를 찾을 수 없습니다.")
        st.stop()

    st.subheader("• 이미지 보기")
    fig1, ax1 = plt.subplots()
    ax1.imshow(image_data, cmap='gray', origin='lower',
               vmin=np.percentile(image_data, 5),
               vmax=np.percentile(image_data, 99))
    st.pyplot(fig1)

    st.subheader("• 별 탐지 및 Flux 추출")
    sigma = mad_std(image_data)
    daofind = DAOStarFinder(fwhm=3.0, threshold=5.0 * sigma)
    sources = daofind(image_data)

    if sources is None or len(sources) == 0:
        st.warning("별을 탐지하지 못했습니다.")
    else:
        positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
        apertures = CircularAperture(positions, r=5.)
        phot_table = aperture_photometry(image_data, apertures)

        for i in range(min(5, len(phot_table))):
            flux = phot_table[i]['aperture_sum']
            x = phot_table[i]['xcenter']
            y = phot_table[i]['ycenter']
            st.write(f"별 {i+1}: 위치=({x:.1f}, {y:.1f}), Flux={flux:.2f}")

        st.subheader("• Flux vs. Magnitude 회귀 분석")
        flux = np.array([2500, 4300, 1600, 3100, 5400])
        mag = np.array([14.7, 13.5, 15.3, 14.1, 13.2])
        log_flux = np.log10(flux)
        slope, intercept, r_value, _, _ = linregress(log_flux, mag)

        st.write(f"회귀식: mag = {slope:.2f} * log10(flux) + {intercept:.2f}")
        st.write(f"R² = {r_value**2:.4f}")

        fig2, ax2 = plt.subplots()
        ax2.scatter(log_flux, mag, label='Data')
        ax2.plot(log_flux, slope * log_flux + intercept, color='red', label='Fit')
        ax2.set_xlabel("log10(Flux)")
        ax2.set_ylabel("Magnitude (u-band)")
        ax2.invert_yaxis()
        ax2.legend()
        ax2.grid(True)
        st.pyplot(fig2)

# -------- 2. 두 번째 & 세 번째 파일 (소광 계수용) --------
st.header("2. 소광 계수 분석용 FITS 파일 2개")
file2 = st.file_uploader("FITS 파일 1 (낮은 AIRMASS)", type=["fits", "fz"], key="file2")
file3 = st.file_uploader("FITS 파일 2 (높은 AIRMASS)", type=["fits", "fz"], key="file3")

def process_file(fits_file):
    hdul = fits.open(fits_file)
    data, header = None, None
    for hdu in hdul:
        if hdu.data is not None:
            data = hdu.data
            header = hdu.header
            break
    hdul.close()
    if data is None:
        raise ValueError("이미지 데이터를 찾을 수 없습니다.")
    
    ra = header.get("RA")
    dec = header.get("DEC")
    date_obs = header.get("DATE-OBS")
    if not (ra and dec and date_obs):
        raise ValueError("RA/DEC/DATE-OBS 정보가 부족합니다.")

    target = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    obs_time = Time(date_obs)
    altaz_frame = AltAz(obstime=obs_time, location=location)
    altaz = target.transform_to(altaz_frame)
    airmass = altaz.secz

    sigma = mad_std(data)
    daofind = DAOStarFinder(fwhm=3.0, threshold=5.0 * sigma)
    sources = daofind(data)
    if sources is None or len(sources) == 0:
        raise ValueError("별을 탐지하지 못했습니다.")
    brightest = sources[np.argmax(sources['flux'])]
    position = (brightest['xcentroid'], brightest['ycentroid'])
    aperture = CircularAperture(position, r=5.)
    phot = aperture_photometry(data, aperture)
    flux = phot[0]['aperture_sum']
    return flux, airmass

if file2 and file3:
    try:
        flux1, X1 = process_file(file2)
        flux2, X2 = process_file(file3)
        delta_mag = -2.5 * np.log10(flux2 / flux1)
        k = delta_mag / (X2 - X1)

        st.subheader("• 소광 계수 계산 결과")
        st.write(f"Flux1 = {flux1:.2f} @ AIRMASS = {X1:.4f}")
        st.write(f"Flux2 = {flux2:.2f} @ AIRMASS = {X2:.4f}")
        st.write(f"Δmag = {delta_mag:.4f}")
        st.success(f"소광 계수 (k) = {k:.4f}")
    except Exception as e:
        st.error(f"처리 중 오류 발생: {e}")

st.markdown("---")
st.markdown("🔗 분석 파트너: [GPTOnline.ai](https://gptonline.ai/ko/) | Streamlit + AI로 천문 분석을 자동화하세요.")

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import mad_std
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, aperture_photometry

# 기본 설정
st.title("별의 물리량 분석 앱 🌟 (u-band 기반)")
st.markdown("u-band FITS 이미지에서 밝기 데이터를 추출하여 별의 물리량을 추정합니다.")

# 별 물리량 추정 함수 (단일 u-band 기준)
def estimate_star_properties_u_band(mag_u, distance_pc, extinction=0.0):
    abs_mag = mag_u - extinction - 5 * np.log10(distance_pc / 10)
    luminosity = 10**(-0.4 * (abs_mag - 4.83))
    temperature = 5778 * (luminosity ** 0.25)  # 단순 가정: 태양 유사
    mass = luminosity ** (1/3.5)
    radius = np.sqrt(luminosity) * (5778 / temperature)**2
    return {
        "Absolute Magnitude (u-band)": abs_mag,
        "Luminosity (L\u2609)": luminosity,
        "Estimated Temperature (K)": temperature,
        "Mass (M\u2609)": mass,
        "Radius (R\u2609)": radius
    }

# 사용자 입력
file_u = st.file_uploader("u-band FITS 이미지 업로드", type=["fits"], key="u")
distance = st.number_input("별까지 거리 (parsec)", min_value=1.0, value=1000.0)
extinction = st.number_input("소광 계수 (A_u)", min_value=0.0, value=0.3)

# 이미지에서 flux 추출 함수
def extract_flux(fits_file):
    hdul = fits.open(fits_file)
    data = None
    for hdu in hdul:
        if hdu.data is not None:
            data = hdu.data
            break
    hdul.close()
    if data is None:
        raise ValueError("이미지 데이터를 찾을 수 없습니다.")

    sigma = mad_std(data)
    daofind = DAOStarFinder(fwhm=3.0, threshold=5. * sigma)
    sources = daofind(data)
    if sources is None or len(sources) == 0:
        raise ValueError("별을 탐지하지 못했습니다.")

    brightest = sources[np.argmax(sources['flux'])]
    position = (brightest['xcentroid'], brightest['ycentroid'])
    aperture = CircularAperture(position, r=5.)
    phot = aperture_photometry(data, aperture)
    return phot[0]['aperture_sum']

if file_u:
    try:
        flux_u = extract_flux(file_u)
        mag_u = -2.5 * np.log10(flux_u)

        result = estimate_star_properties_u_band(mag_u, distance, extinction)

        st.subheader("⭐ 분석 결과")
        for k, v in result.items():
            st.write(f"{k}: {v:.3f}")

        st.subheader("🔭 밝기 정보")
        st.write(f"u-band Flux: {flux_u:.2f} → Magnitude: {mag_u:.3f}")

    except Exception as e:
        st.error(f"처리 중 오류 발생: {e}")

def temperature_to_rgb(temp):
    if temp < 3700:
        return (1.0, 0.5, 0.5)  # Red
    elif temp < 5200:
        return (1.0, 0.7, 0.4)  # Orange
    elif temp < 6000:
        return (1.0, 1.0, 0.6)  # Yellow
    elif temp < 7500:
        return (1.0, 1.0, 1.0)  # White
    else:
        return (0.6, 0.6, 1.0)  # Blue

def show_colored_star(data, position, temperature):
    rgb_color = temperature_to_rgb(temperature)

    fig, ax = plt.subplots()
    ax.imshow(data, cmap='gray', origin='lower',
              vmin=np.percentile(data, 5), vmax=np.percentile(data, 99))
    ax.plot(position[0], position[1], marker='o', markersize=15,
            markerfacecolor=rgb_color, markeredgecolor='white', markeredgewidth=1.5)
    ax.set_title(f"별 위치 및 색상 (T={temperature:.0f}K)")
    st.pyplot(fig)
