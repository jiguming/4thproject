import streamlit as st
from astropy.io import fits
import numpy as np
import requests
from io import BytesIO

st.title("🌌 GitHub에서 FITS.FZ 파일 불러오기 및 시각화")

# 고정된 GitHub raw URL
url = "https://raw.githubusercontent.com/jiguming/4thproject/main/kwb_190326_032609_ori.fits.fz"

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
                st.subheader(f"🖼 HDU {i} - 2차원 이미지")
                # 이미지 정규화
                norm_data = (data - np.min(data)) / (np.max(data) - np.min(data))
                st.image(norm_data, caption=f"HDU {i}", use_column_width=True, clamp=True)

                st.subheader(f"🧾 HDU {i} 헤더")
                st.text(str(hdu.header))
                found = True

        if not found:
            st.warning("⚠️ 2차원 이미지 데이터가 포함된 HDU가 없습니다.")
except Exception as e:
    st.error(f"❌ 오류 발생: {e}")
