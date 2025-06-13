import streamlit as st
from astropy.io import fits
import numpy as np

st.title("🛰️ 업로드된 FITS.FZ 파일 탐색 및 시각화")

# Streamlit에 업로드한 파일 받기
uploaded_file = st.file_uploader("FITS 파일 (.fits.fz) 업로드", type=["fits", "fz"])

if uploaded_file is not None:
    try:
        with fits.open(uploaded_file) as hdul:
            st.write("📁 HDU 구조:")
            st.text(hdul.info())

            found = False
            for i, hdu in enumerate(hdul):
                data = hdu.data
                if data is not None and data.ndim == 2:
                    st.subheader(f"🖼 HDU {i} - 2D 이미지")
                    # 정규화 및 시각화
                    norm = (data - np.min(data)) / (np.max(data) - np.min(data))
                    st.image(norm, caption=f"HDU {i}", use_column_width=True, clamp=True)

                    st.subheader(f"🧾 HDU {i} 헤더")
                    st.text(str(hdu.header))
                    found = True

            if not found:
                st.warning("⚠️ 2차원 이미지 데이터가 포함된 HDU가 없습니다.")
    except Exception as e:
        st.error(f"❌ 파일 처리 중 오류 발생: {e}")
