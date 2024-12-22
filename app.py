import streamlit as st
import pandas as pd
import requests
from chemtools import ChemTools

# Título de la aplicación
st.title("Visualización de Moléculas en 2D con ChemTools")

# URL del archivo CSV en GitHub
github_url = "https://raw.githubusercontent.com/mdgenriquez/punonpdb/main/MDG202_cl2.csv"

try:
    # Descargar y leer el archivo CSV
    response = requests.get(github_url)
    response.raise_for_status()  # Verifica errores de descarga
    df = pd.read_csv(io.StringIO(response.text))

    # Verificar que la columna 'SMILES' esté presente
    if 'SMILES' not in df.columns:
        st.error("El archivo debe contener una columna llamada 'SMILES'.")
    else:
        st.success("Archivo cargado correctamente.")
        st.dataframe(df)

        # Visualización de las moléculas en 2D
        st.subheader("Visualización de Moléculas en 2D")
        for index, row in df.iterrows():
            smiles = row['SMILES']
            mol = ChemTools.from_smiles(smiles)
            if mol:
                # Generar imagen de la molécula en 2D
                img = mol.to_image(size=(300, 300))  # Genera la imagen en 2D
                st.image(img, caption=f"Molécula {index + 1}: {smiles}")
            else:
                st.warning(f"No se pudo procesar el SMILES: {smiles}")
except Exception as e:
    st.error(f"Error al descargar o procesar el archivo: {e}")
