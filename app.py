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
st.sidebar.image("ima.webp") #, caption="Autor: -Dr.Jesus Alvarado 
-Guadalupe Enriquez"
st.sidebar.title("COMPUESTOS SIMULADOS")
st.sidebar.markdown("Autor: -Dr. Jesus Alvarado -Guadalupe Enriquez")


