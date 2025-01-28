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

#title probando

        
def load_data(file_path):
    
    try:
        return pd.read_csv(file_path)
    except Exception as e:
        st.error(f"Error al cargar el archivo CSV: {e}")
        return None

#práctica
def Home():
    st.header('De nombre común a 2D :cat:', divider='rainbow')
    st.sidebar.markdown("# Nombre clásico:")
    st.sidebar.markdown("Trivial name, non-systematic name for a chemical substance, son otras denominaciones en inglés")

    entrada = st.text_input("Escribe el nombre común en inglés:", "Aluminium hydroxide")

     st.markdown("### IUPAC")  
    nombreiupac = pcp.get_compounds(entrada,'name')
    st.text(nombreiupac[0].iupac_name)
    st.markdown("### SMILES Isomérico")
    smilesisomerico = get_compounds(entrada, 'name')
    st.text(smilesisomerico[0].isomeric_smiles)

    st.markdown("### Masa molecular (g/mol)")
    masamolecular = get_compounds(entrada, 'name')
    st.text(masamolecular[0].exact_mass)
  
    st.markdown("### Coeficiente de partición")
    coeficientedeparticion = get_compounds(entrada, 'name')
    st.text(coeficientedeparticion[0].xlogp)

    st.markdown("### PubChem ID")
    id_pubchem = pcp.get_compounds(entrada, 'name')
    st.text(id_pubchem)

    st.markdown("### Representación simplificada")
    m0 = Chem.MolFromSmiles(smilesisomerico[0].isomeric_smiles)    
    Draw.MolToFile(m0,'mol0.png')
    #st.pyplot()
    st.write('Molecule 2D :smiley:')
    st.image('mol0.png')


def page2():
    st.header('De SMILES a 2D :smiley:', divider='rainbow')
    st.sidebar.markdown("# Simplified Molecular Input Line Entry System")
    st.sidebar.markdown("Sistema de introducción molecular lineal simplificada")
    
    entrada = st.text_input("Escribe el nombre SMILES: ", "C1=CC2=C(C3=C(C=CC=N3)C=C2)N=C1")
    st.markdown("### PubChem ID:")
    identificador = pcp.get_compounds(entrada, 'smiles')
    st.text(identificador)

    st.markdown("### Nombre IUPAC")  
    nombreiupac = pcp.get_compounds(entrada,'smiles')
    st.text(nombreiupac[0].iupac_name)

    st.markdown("### Representación simplificada")
    m1 = Chem.MolFromSmiles(entrada)    
    Draw.MolToFile(m1,'mol1.png')
    st.write('Molecule 2D :smiley:')
    st.image('mol1.png')

# Generar archivo SDF
def generate_sdf(mol):
    """Genera un archivo SDF 3D de la molécula dada."""
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    
    sdf_data = BytesIO()
    writer = Chem.SDWriter(sdf_data)
    writer.write(mol)
    writer.close()
    sdf_data.seek(0)
    return sdf_data


csv_path = "base_datos_plant_smile.csv"  # Nombre del archivo
data = load_data(csv_path)

if data is not None:
    # Validar columnas necesarias
    if "SMILES" not in data.columns or "Name" not in data.columns:
        st.error("El archivo CSV debe contener las columnas 'SMILES' y 'Name'.")
    else:
        st.title("Predictor plant")
        
        # Mostrar un slider para seleccionar una molécula
        molecule_index = st.slider("Selecciona una molécula", 0, len(data) - 1, 0)
        selected_row = data.iloc[molecule_index]
        
        # Extraer datos de la molécula seleccionada
        smiles = selected_row["SMILES"]
        name = selected_row["Name"]

        #st.markdown("### Nombre IUPAC")  
        #nombreiupac = pcp.get_compounds(smiles,'smiles')
        #st.text(nombreiupac[0].iupac_name)

        st.markdown("### Coeficiente de partición")
        coeficientedeparticion = get_compounds(smiles, 'smiles')
        st.text(coeficientedeparticion[0].xlogp)
        
        #st.subheader(f"ID: {name}")
        st.text(f"Código SMILES: {smiles}")
        
        # Generar la representación 2D de la molécula
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            Draw.MolToFile(mol, "molecule.png")
            st.image("molecule.png", caption="Representación 2D de la molécula")
        else:
            st.error("No se pudo generar la representación molecular a partir del código SMILES.")
        
        # Visualización en 3D
        st.subheader("Visualización en 3D")
        def show_3d(smi):
            mol = Chem.MolFromSmiles(smi)
            mol = Chem.AddHs(mol)
            Chem.AllChem.EmbedMolecule(mol)
            Chem.AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            mblock = Chem.MolToMolBlock(mol)
            viewer = py3Dmol.view(width=500, height=400)
            viewer.addModel(mblock, "mol")
            viewer.setStyle({"stick": {}})
            viewer.zoomTo()
            showmol(viewer, height=400, width=500)

        show_3d(smiles)
        
    
else:
    st.warning("El archivo CSV no contiene datos válidos o no está cargado.")

