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
