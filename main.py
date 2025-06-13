import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import mad_std
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, aperture_photometry

# 기본 설정
st.title("별의 물리량 분석 앱 🌟 (u-band 기반)")
st.markdown("u-band FITS 이미지에서 밝기 데이터를 추출하여 별의 물리량을 추정하고 예상 색상을 시각화합니다.")

# 별 물리량 추정 함수 (단일 u-band 기준)
def estimate_star_properties_u_band(mag_u, distance_pc, extinction=0.0):
    abs_mag = mag_u - extinction - 5 * np.log10(distance_pc / 10)
    luminosity = 10**(-0.4 * (abs_mag - 4.83))
    temperature = 5778 * (luminosity ** 0.25)  # 단순 가정: 태양 유사
    mass = luminosity ** (1/3.5)
    radius = np.sqrt(luminosity) * (5778 / temperature)**2
    return {
        "Absolute Magnitude (u-band)": abs_mag,
        "Luminosity (L☉)": luminosity,
        "Estimated Temperature (K)": temperature,
        "Mass (M☉)": mass,
        "Radius (R☉)": radius
    }

# 온도에 따른 색상 추정 함수
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

# 이미지에 색상 입히기

def show_colored_star(data, position, temperature):
    rgb_color = temperature_to_rgb(temperature)

    fig, ax = plt.subplots()
    ax.imshow(data, cmap='gray', origin='lower',
              vmin=np.percentile(data, 5), vmax=np.percentile(data, 99))
    ax.plot(position[0], position[1], marker='o', markersize=15,
            markerfacecolor=rgb_color, markeredgecolor='white', markeredgewidth=1.5)
    ax.set_title(f"별 위치 및 색상 (T={temperature:.0f}K)")
    st.pyplot(fig)

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
    return phot[0]['aperture_sum'], data, position

if file_u:
    try:
        flux_u, image_data, star_pos = extract_flux(file_u)
        mag_u = -2.5 * np.log10(flux_u)

        result = estimate_star_properties_u_band(mag_u, distance, extinction)

        st.subheader("⭐ 분석 결과")
        for k, v in result.items():
            st.write(f"{k}: {v:.3f}")

        st.subheader("🔭 밝기 정보")
        st.write(f"u-band Flux: {flux_u:.2f} → Magnitude: {mag_u:.3f}")

        # 색상 입힌 별 시각화
        show_colored_star(image_data, star_pos, result["Estimated Temperature (K)"])

    except Exception as e:
        st.error(f"처리 중 오류 발생: {e}")
