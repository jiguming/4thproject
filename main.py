import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.stats import mad_std
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, aperture_photometry
from astropy.coordinates import SkyCoord, AltAz
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import EarthLocation
from scipy.stats import linregress

st.set_page_config(page_title="천문 이미지 분석", layout="wide")

# 위치 정보
location = EarthLocation(lat=31.9584 * u.deg, lon=-111.5967 * u.deg, height=2120 * u.m)

# ----------------- 사이드바 -----------------
st.sidebar.title("🧭 분석 설정")
tab_choice = st.sidebar.radio("기능 선택", ["별 밝기 분석", "소광 계수", "물리량 추정"])

# ----------------- 탭 구조 -----------------
st.title("🌌 천문 이미지 분석 앱")
st.markdown("FITS 파일을 분석하여 **별 밝기**, **소광 계수**, **물리량**을 추정합니다.")

# -------------------------------------
# 1. 별 밝기 분석 탭
# -------------------------------------
if tab_choice == "별 밝기 분석":
    file1 = st.sidebar.file_uploader("분석용 FITS 파일 업로드", type=["fits"], key="file1")
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

        st.subheader("🔭 이미지 미리보기")
        fig1, ax1 = plt.subplots()
        ax1.imshow(image_data, cmap='gray', origin='lower',
                   vmin=np.percentile(image_data, 5),
                   vmax=np.percentile(image_data, 99))
        st.pyplot(fig1)

        st.subheader("✨ 별 탐지 및 밝기 분석")
        sigma = mad_std(image_data)
        daofind = DAOStarFinder(fwhm=3.0, threshold=5.0 * sigma)
        sources = daofind(image_data)

        if sources is None or len(sources) == 0:
            st.warning("별을 탐지하지 못했습니다.")
        else:
            positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
            apertures = CircularAperture(positions, r=5.)
            phot_table = aperture_photometry(image_data, apertures)
            fluxes = []
            mags = []

            for i in range(min(5, len(phot_table))):
                flux = phot_table[i]['aperture_sum']
                x = phot_table[i]['xcenter']
                y = phot_table[i]['ycenter']
                st.write(f"⭐ 별 {i+1}: 위치=({x:.1f}, {y:.1f}), Flux={flux:.2f}")
                fluxes.append(flux)

            # 가상의 mag 값 사용
            mags = np.array([14.7, 13.5, 15.3, 14.1, 13.2])
            log_flux = np.log10(fluxes)
            slope, intercept, r_value, _, _ = linregress(log_flux, mags)

            st.success(f"회귀식: mag = {slope:.2f} * log10(flux) + {intercept:.2f}")
            st.write(f"R² = {r_value**2:.4f}")

            fig2, ax2 = plt.subplots()
            ax2.scatter(log_flux, mags, label='Data')
            ax2.plot(log_flux, slope * log_flux + intercept, color='red', label='Fit')
            ax2.set_xlabel("log10(Flux)")
            ax2.set_ylabel("Magnitude (u-band)")
            ax2.invert_yaxis()
            ax2.legend()
            ax2.grid(True)
            st.pyplot(fig2)

# -------------------------------------
# 2. 소광 계수 분석
# -------------------------------------
elif tab_choice == "소광 계수":
    file2 = st.sidebar.file_uploader("FITS 파일 1 (낮은 AIRMASS)", type=["fits"], key="file2")
    file3 = st.sidebar.file_uploader("FITS 파일 2 (높은 AIRMASS)", type=["fits"], key="file3")

    def process_file(fits_file):
        hdul = fits.open(fits_file)
        data, header = None, None
        for hdu in hdul:
            if hdu.data is not None:
                data = hdu.data
                header = hdu.header
                break
        hdul.close()

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

            st.subheader("🧮 소광 계수 계산 결과")
            st.write(f"Flux1 = {flux1:.2f} @ AIRMASS = {X1:.4f}")
            st.write(f"Flux2 = {flux2:.2f} @ AIRMASS = {X2:.4f}")
            st.write(f"Δmag = {delta_mag:.4f}")
            st.success(f"소광 계수 (k) = {k:.4f}")
        except Exception as e:
            st.error(f"오류 발생: {e}")

# -------------------------------------
# 3. 별 물리량 추정
# -------------------------------------
elif tab_choice == "물리량 추정":
    file_u = st.sidebar.file_uploader("u-band FITS 이미지 업로드", type=["fits"], key="u")
    distance = st.sidebar.number_input("별까지 거리 (pc)", min_value=1.0, value=1000.0)
    extinction = st.sidebar.number_input("소광 계수 (A_u)", min_value=0.0, value=0.3)

    def extract_flux(fits_file):
        hdul = fits.open(fits_file)
        data = None
        for hdu in hdul:
            if hdu.data is not None:
                data = hdu.data
                break
        hdul.close()

        sigma = mad_std(data)
        daofind = DAOStarFinder(fwhm=3.0, threshold=5. * sigma)
        sources = daofind(data)
        if sources is None or len(sources) == 0:
            raise ValueError("별을 탐지하지 못했습니다.")

        brightest = sources[np.argmax(sources['flux'])]
        position = (brightest['xcentroid'], brightest['ycentroid'])
        aperture = CircularAperture(position, r=5.)
        phot = aperture_photometry(data, aperture)
        return phot[0]['aperture_sum'], data, position

    def estimate_star_properties(mag_u, distance_pc, extinction=0.0):
        abs_mag = mag_u - extinction - 5 * np.log10(distance_pc / 10)
        luminosity = 10**(-0.4 * (abs_mag - 4.83))
        temperature = 5778 * (luminosity ** 0.25)
        mass = luminosity ** (1/3.5)
        radius = np.sqrt(luminosity) * (5778 / temperature)**2
        return {
            "Absolute Magnitude (u)": abs_mag,
            "Luminosity (L☉)": luminosity,
            "Temperature (K)": temperature,
            "Mass (M☉)": mass,
            "Radius (R☉)": radius
        }

    def temperature_to_rgb(temp):
        if temp < 3700: return (1.0, 0.5, 0.5)
        elif temp < 5200: return (1.0, 0.7, 0.4)
        elif temp < 6000: return (1.0, 1.0, 0.6)
        elif temp < 7500: return (1.0, 1.0, 1.0)
        else: return (0.6, 0.6, 1.0)

    def show_colored_star(data, position, temperature):
        rgb_color = temperature_to_rgb(temperature)
        fig, ax = plt.subplots()
        ax.imshow(data, cmap='gray', origin='lower',
                  vmin=np.percentile(data, 5), vmax=np.percentile(data, 99))
        ax.plot(position[0], position[1], marker='o', markersize=15,
                markerfacecolor=rgb_color, markeredgecolor='white', markeredgewidth=1.5)
        ax.set_title(f"별 위치 및 색상 (T={temperature:.0f}K)")
        st.pyplot(fig)

    if file_u:
        try:
            flux_u, image_data, star_pos = extract_flux(file_u)
            mag_u = -2.5 * np.log10(flux_u)
            result = estimate_star_properties(mag_u, distance, extinction)

            st.subheader("🌟 물리량 추정 결과")
            for k, v in result.items():
                st.write(f"{k}: {v:.3f}")
            st.subheader("🔭 u-band Flux 및 시각화")
            st.write(f"u-band Flux: {flux_u:.2f} → Magnitude: {mag_u:.3f}")
            show_colored_star(image_data, star_pos, result["Temperature (K)"])
        except Exception as e:
            st.error(f"오류 발생: {e}")
