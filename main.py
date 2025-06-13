import streamlit as st
from astropy.io import fits
import numpy as np
import requests
from io import BytesIO

st.title("🌌 GitHub에서 FITS.FZ 파일 불러오기")

# URL 고정 (혹은 st.text_input으로 입력 가능)
url = "https://raw.githubusercontent.com/jiguming/4thproject/main/k21i_100108_031209_ori.fits.fz"

try:
    response = requests.get(url)
    response.raise_for_status()

    with fits.open(BytesIO(response.content)) as hdul:
        st.write("📁 HDU 구조:")
        st.text(hdul.info())

        hdu = hdul[0]
        header = hdu.header
        data = hdu.data

        st.subheader("🧾 헤더")
        st.text(str(header))

        if data is not None and data.ndim == 2:
            # 이미지 정규화
            data_norm = (data - np.min(data)) / (np.max(data) - np.min(data))
            st.image(data_norm, caption="FITS 이미지", use_column_width=True, clamp=True)
        else:
            st.warning("2차원 이미지 데이터가 아닙니다.")
except Exception as e:
    st.error(f"❌ 오류 발생: {e}")
