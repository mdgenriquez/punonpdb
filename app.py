import streamlit as st
import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import Draw

# Título de la aplicación
st.title("Visualización de Moléculas en 2D")

# URL del archivo
github_url = "https://raw.githubusercontent.com/mdgenriquez/punonpdb/main/MDG202_cl2.csv"

try:
    response = requests.get(github_url)
    response.raise_for_status()  # Verifica errores de descarga
    df = pd.read_csv(io.StringIO(response.text))
    
    if 'SMILES' not in df.columns:
        st.error("El archivo debe contener una columna llamada 'SMILES'.")
    else:
        st.success("Archivo cargado correctamente.")
        st.dataframe(df)

        st.subheader("Visualización de Moléculas en 2D")
        for index, row in df.iterrows():
            smiles = row['SMILES']
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                img = Draw.MolToImage(mol)
                st.image(img, caption=f"Molécula {index + 1}: {smiles}")
            else:
                st.warning(f"No se pudo procesar el SMILES: {smiles}")
except Exception as e:
    st.error(f"Error al descargar o procesar el archivo: {e}")

