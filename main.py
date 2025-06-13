import streamlit as st
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import requests
from io import BytesIO

st.title("🔭 FITS 이미지 및 WCS 시각화 (HDU 1 기준)")

# GitHub raw URL
url = "https://raw.githubusercontent.com/jiguming/4thproject/main/k21i_100108_031209_ori.fits.fz"

try:
    response = requests.get(url)
    response.raise_for_status()

    with fits.open(BytesIO(response.content)) as hdul:
        hdu = hdul[1]  # ← 실제 이미지가 있는 HDU 번호로 바꿔야 함

        data = hdu.data
        header = hdu.header

        st.subheader("🧾 헤더")
        st.text(str(header))

        if data is not None and data.ndim == 2:
            st.subheader("🖼 이미지 시각화")
            # 이미지 정규화
            norm = (data - np.min(data)) / (np.max(data) - np.min(data))
            st.image(norm, use_column_width=True, clamp=True)

            # WCS 시도
            try:
                wcs = WCS(header)
                if wcs.has_celestial:
                    st.success("🌐 WCS (World Coordinate System) 정보가 포함되어 있음!")
                    st.code(wcs.to_header_string(), language="fits")
                else:
                    st.info("WCS 좌표 정보는 포함되어 있지 않음.")
            except Exception as wcs_err:
                st.warning(f"WCS 파싱 중 오류 발생: {wcs_err}")
        else:
            st.warning("⚠️ HDU에 2차원 데이터가 없습니다.")
except Exception as e:
    st.error(f"❌ 오류 발생: {e}")
