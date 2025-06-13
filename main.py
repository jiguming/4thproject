import streamlit as st
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import tempfile

st.title("FITS / FPACK (.fits.fz) 파일 뷰어 🛰️")

uploaded_file = st.file_uploader("FITS 파일(.fits 또는 .fits.fz) 업로드", type=["fits", "fz"])

if uploaded_file is not None:
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fits.fz") as tmp:
        tmp.write(uploaded_file.read())
        tmp_path = tmp.name

    try:
        with fits.open(tmp_path) as hdul:
            st.subheader("📋 헤더 정보")
            st.text(hdul[0].header)

            if hdul[0].data is not None:
                st.subheader("🖼 이미지 시각화")
                data = hdul[0].data
                if data.ndim == 2:
                    fig, ax = plt.subplots()
                    ax.imshow(data, cmap='gray', origin='lower', 
                              vmin=np.percentile(data, 5), vmax=np.percentile(data, 95))
                    ax.set_title("FITS 이미지")
                    st.pyplot(fig)
                else:
                    st.warning(f"{data.ndim}차원 데이터는 시각화되지 않습니다.")
            else:
                st.info("이미지 데이터 없음.")

    except Exception as e:
        st.error(f"FITS 파일 열기 오류: {e}")
