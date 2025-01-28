import streamlit as st
from pubchempy import get_compounds, Compound
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import streamlit.components.v1 as components
import py3Dmol
from stmol import showmol
from io import BytesIO
import sys
import pubchempy as pcp



# Configuración de la barra lateral
st.sidebar.image("ima.webp") #, caption="Autores:  -Guadalupe Enriquez -Dr.Jesus Alvarado"
st.sidebar.markdown("""
    <h1 style='color:green; font-size: 24px;'>COMPUESTOS SIMULADOS</h1>
""", unsafe_allow_html=True)
st.sidebar.markdown("Autor:\n- Dr. Jesus Alvarado\n- Bach. Guadalupe Enriquez")


# Streamlit Title
st.title("COMPUESTOS SIMULADOS")
def img_to_base64(image_path):
    """Convert image to base64."""
    try:
        with open(image_path, "rb") as img_file:
            return base64.b64encode(img_file.read()).decode()
    except Exception as e:
        logging.error(f"Error converting image to base64: {str(e)}")
        return None

@st.cache_data(show_spinner=False)
