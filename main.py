import streamlit as st
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import gzip
import shutil
import os

st.title("FITS 파일 뷰어 🛰️")

# 파일 업로드
uploaded_file = st.file_uploader("FITS 파일(.fits, .fits.fz) 업로드", type=["fits", "fz"])

if uploaded_file is not None:
    # 저장 경로 설정
    temp_path = "temp_file.fits"

    # fz 확장자면 압축 해제
    if uploaded_file.name.endswith(".fz"):
        with open("temp_file.fits.fz", "wb") as f_out:
            f_out.write(uploaded_file.read())
        with gzip.open("temp_file.fits.fz", "rb") as f_in:
            with open(temp_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove("temp_file.fits.fz")
    else:
        with open(temp_path, "wb") as f:
            f.write(uploaded_file.read())

    # FITS 파일 열기
    try:
        with fits.open(temp_path) as hdul:
            st.subheader("📋 헤더 정보")
            st.text(hdul[0].header)

            # 이미지 데이터가 있으면 시각화
            if hdul[0].data is not None:
                st.subheader("🖼 이미지 데이터 시각화")
                data = hdul[0].data
                if data.ndim == 2:
                    fig, ax = plt.subplots()
                    ax.imshow(data, cmap='gray', origin='lower', vmin=np.percentile(data, 5), vmax=np.percentile(data, 95))
                    ax.set_title("FITS 이미지")
                    st.pyplot(fig)
                else:
                    st.warning(f"이미지 차원: {data.ndim}D. 2D 데이터만 시각화됩니다.")
            else:
                st.info("이미지 데이터 없음.")

    except Exception as e:
        st.error(f"FITS 파일을 여는 중 오류 발생: {e}")

    # 임시 파일 삭제
    if os.path.exists(temp_path):
        os.remove(temp_path)
