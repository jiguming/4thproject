import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import mad_std
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, aperture_photometry
from astropy.coordinates import SkyCoord, AltAz
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import EarthLocation
from scipy.stats import linregress

st.set_page_config(page_title="우주 관측실 🌌", layout="wide")
st.title("🪐 천문 이미지 분석 앱")
st.caption("FITS 이미지를 활용한 별의 물리량 분석, 소광 계수 계산, 그리고 퀴즈까지!")

# 공통 위치 정보
location = EarthLocation(lat=31.9584 * u.deg, lon=-111.5967 * u.deg, height=2120 * u.m)

# 탭 구조 설정
tabs = st.tabs(["🎯 별 밝기 분석", "🌫️ 소광 계수 계산", "🔬 별의 물리량 추정", "🧠 퀴즈 풀기"])

# ------------------------------------------------------
# 🎯 별 밝기 분석
# ------------------------------------------------------
with tabs[0]:
    st.header("🎯 별 밝기 분석")
    file1 = st.file_uploader("FITS 파일 업로드", type=["fits", "fz"], key="brightness")
    if file1:
        with st.spinner("🔍 이미지를 분석하는 중..."):
            hdul = fits.open(file1)
            image_data = next((hdu.data for hdu in hdul if hdu.data is not None), None)
            hdul.close()

            if image_data is None:
                st.error("❌ 이미지 데이터가 없습니다.")
            else:
                fig, ax = plt.subplots()
                ax.imshow(image_data, cmap='gray', origin='lower',
                          vmin=np.percentile(image_data, 5), vmax=np.percentile(image_data, 99))
                st.pyplot(fig)

                sigma = mad_std(image_data)
                daofind = DAOStarFinder(fwhm=3.0, threshold=5. * sigma)
                sources = daofind(image_data)

                if sources is None or len(sources) == 0:
                    st.warning("🌌 별을 찾지 못했습니다.")
                else:
                    st.success("✅ 별 탐지 성공!")
                    phot_table = aperture_photometry(image_data, CircularAperture(
                        np.transpose((sources['xcentroid'], sources['ycentroid'])), r=5.))
                    fluxes = [phot_table[i]['aperture_sum'] for i in range(min(5, len(phot_table)))]
                    mags = np.array([14.7, 13.5, 15.3, 14.1, 13.2])  # 예시값
                    log_flux = np.log10(fluxes)

                    slope, intercept, r, _, _ = linregress(log_flux, mags)
                    st.code(f"mag = {slope:.2f} * log10(flux) + {intercept:.2f}   (R²={r**2:.4f})")

                    fig2, ax2 = plt.subplots()
                    ax2.scatter(log_flux, mags, label='측정 데이터')
                    ax2.plot(log_flux, slope * log_flux + intercept, 'r-', label='회귀선')
                    ax2.invert_yaxis()
                    ax2.set_xlabel("log10(Flux)")
                    ax2.set_ylabel("Magnitude")
                    ax2.legend()
                    st.pyplot(fig2)

# ------------------------------------------------------
# 🌫️ 소광 계수
# ------------------------------------------------------
with tabs[1]:
    st.header("🌫️ 대기 소광 계수 계산")
    file2 = st.file_uploader("🔻 낮은 AIRMASS 이미지", type=["fits", "fz"], key="airmass1")
    file3 = st.file_uploader("🔺 높은 AIRMASS 이미지", type=["fits", "fz"], key="airmass2")

    def process_fits(file):
        hdul = fits.open(file)
        data = next((hdu.data for hdu in hdul if hdu.data is not None), None)
        header = next((hdu.header for hdu in hdul if hdu.data is not None), None)
        hdul.close()

        ra = header.get("RA")
        dec = header.get("DEC")
        date_obs = header.get("DATE-OBS")
        if not (ra and dec and date_obs):
            raise ValueError("RA, DEC, DATE-OBS 정보가 누락됨.")

        coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
        time = Time(date_obs)
        airmass = coord.transform_to(AltAz(obstime=time, location=location)).secz

        sigma = mad_std(data)
        daofind = DAOStarFinder(fwhm=3.0, threshold=5. * sigma)
        sources = daofind(data)
        brightest = sources[np.argmax(sources['flux'])]
        flux = aperture_photometry(data, CircularAperture((brightest['xcentroid'], brightest['ycentroid']), r=5.))[0]['aperture_sum']
        return flux, airmass

    if file2 and file3:
        try:
            flux1, X1 = process_fits(file2)
            flux2, X2 = process_fits(file3)
            delta_mag = -2.5 * np.log10(flux2 / flux1)
            k = delta_mag / (X2 - X1)
            st.success(f"☑️ 소광 계수 (k): {k:.4f}")
            st.write(f"Flux1: {flux1:.2f} @ AIRMASS {X1:.3f}")
            st.write(f"Flux2: {flux2:.2f} @ AIRMASS {X2:.3f}")
        except Exception as e:
            st.error(f"🚫 오류 발생: {e}")

# ------------------------------------------------------
# 🔬 별의 물리량 추정
# ------------------------------------------------------
with tabs[2]:
    st.header("🔬 u-band 기반 별 물리량 추정")
    file_u = st.file_uploader("u-band 이미지 업로드", type=["fits"], key="u_band")
    distance = st.slider("⭐ 거리 (parsec)", 1, 10000, 1000)
    extinction = st.slider("🌫️ 소광 계수 (A_u)", 0.0, 2.0, 0.3)

    def temperature_to_rgb(temp):
        if temp < 3700: return (1.0, 0.5, 0.5)
        elif temp < 5200: return (1.0, 0.7, 0.4)
        elif temp < 6000: return (1.0, 1.0, 0.6)
        elif temp < 7500: return (1.0, 1.0, 1.0)
        else: return (0.6, 0.6, 1.0)

    def estimate_physics(mag, d, A=0.0):
        abs_mag = mag - A - 5 * np.log10(d / 10)
        L = 10 ** (-0.4 * (abs_mag - 4.83))
        T = 5778 * (L ** 0.25)
        M = L ** (1 / 3.5)
        R = np.sqrt(L) * (5778 / T) ** 2
        return abs_mag, L, T, M, R

    if file_u:
        hdul = fits.open(file_u)
        data = next((hdu.data for hdu in hdul if hdu.data is not None), None)
        hdul.close()

        sigma = mad_std(data)
        daofind = DAOStarFinder(fwhm=3.0, threshold=5. * sigma)
        sources = daofind(data)
        brightest = sources[np.argmax(sources['flux'])]
        flux = aperture_photometry(data, CircularAperture((brightest['xcentroid'], brightest['ycentroid']), r=5.))[0]['aperture_sum']
        mag = -2.5 * np.log10(flux)
        abs_mag, L, T, M, R = estimate_physics(mag, distance, extinction)

        st.success("⭐ 물리량 추정 결과")
        st.write(f"Absolute Magnitude: {abs_mag:.3f}")
        st.write(f"Luminosity (L☉): {L:.3f}")
        st.write(f"Temperature (K): {T:.1f}")
        st.write(f"Mass (M☉): {M:.3f}")
        st.write(f"Radius (R☉): {R:.3f}")

        fig, ax = plt.subplots()
        ax.imshow(data, cmap='gray', origin='lower',
                  vmin=np.percentile(data, 5), vmax=np.percentile(data, 99))
        ax.plot(brightest['xcentroid'], brightest['ycentroid'], 'o',
                markersize=15, markerfacecolor=temperature_to_rgb(T),
                markeredgecolor='white')
        ax.set_title(f"별 위치 및 색상 (T = {T:.0f}K)")
        st.pyplot(fig)

# ------------------------------------------------------
# 🧠 퀴즈
# ------------------------------------------------------
with tabs[3]:
    st.header("🧠 천문 퀴즈 타임!")
    st.markdown("별의 물리량과 관련된 퀴즈를 풀어보세요!")

    q1 = st.radio("1. 별의 겉보기 등급이 일정할 때, 거리가 멀수록 어떤가요?", 
                  ["더 밝아진다", "더 어두워진다", "변하지 않는다"], index=None)
    if q1:
        if q1 == "더 어두워진다":
            st.success("✅ 정답입니다!")
        else:
            st.error("❌ 오답입니다. 거리가 멀수록 빛이 퍼지기 때문에 어두워집니다.")

    q2 = st.radio("2. 별의 온도가 높을수록 어떤 색을 띨까요?", 
                  ["붉은색", "노란색", "푸른색"], index=None)
    if q2:
        if q2 == "푸른색":
            st.success("✅ 정확해요! 고온 별은 푸른 빛을 냅니다.")
        else:
            st.warning("💡 다시 생각해보세요! 고온 = 고에너지 파장입니다.")
