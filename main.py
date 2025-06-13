import streamlit as st
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

st.title("FITS 파일 분석 도구")

# 파일 업로드
uploaded_file = st.file_uploader("FITS 파일 (.fits 또는 .fits.fz)을 업로드하세요", type=["fits", "fz"])

if uploaded_file is not None:
    with fits.open(uploaded_file) as hdul:
        st.write("📁 파일 구조 (HDU 리스트):")
        hdul.info()

        # 첫 번째 HDU 선택
        hdu = hdul[0]
        header = hdu.header
        data = hdu.data

        st.subheader("🧾 헤더 정보")
        st.text(header)

        if data is not None:
            st.subheader("🖼 데이터 미리보기")
            if data.ndim == 2:
                fig, ax = plt.subplots()
                ax.imshow(data, cmap='gray', origin='lower')
                st.pyplot(fig)
            else:
                st.write(f"데이터 차원: {data.shape} (시각화 불가)")
        else:
            st.warning("데이터가 포함되지 않은 HDU입니다.")
