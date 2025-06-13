import streamlit as st
from astropy.io import fits
import numpy as np
import requests
from io import BytesIO

st.title("🌌 GitHub FITS.FZ 파일 - 모든 HDU 탐색")

# GitHub raw URL
url = "https://raw.githubusercontent.com/jiguming/4thproject/main/k21i_100108_031209_ori.fits.fz"

try:
    response = requests.get(url)
    response.raise_for_status()

    with fits.open(BytesIO(response.content)) as hdul:
        st.write("📁 HDU 구조:")
        st.text(hdul.info())

        found = False
        for i, hdu in enumerate(hdul):
            data = hdu.data
            if data is not None and data.ndim == 2:
                st.subheader(f"🖼 HDU {i} - 2D 이미지 데이터")
                # 데이터 정규화 후 시각화
                norm_data = (data - np.min(data)) / (np.max(data) - np.min(data))
                st.image(norm_data, caption=f"HDU {i}", use_column_width=True, clamp=True)

                st.subheader(f"🧾 HDU {i} 헤더")
                st.text(str(hdu.header))
                found = True

        if not found:
            st.warning("⚠️ 시각화 가능한 2차원 이미지가 포함된 HDU가 없습니다.")
except Exception as e:
    st.error(f"❌ 오류 발생: {e}")
