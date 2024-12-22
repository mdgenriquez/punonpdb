import streamlit as st
import pandas as pd
import requests
import io  # Importar io para manejar StringIO
from rdkit import Chem
from rdkit.Chem import Draw

# Título de la aplicación
st.title("Visualización de Moléculas en 2D")

# URL del archivo
github_url = "https://raw.githubusercontent.com/mdgenriquez/punonpdb/main/MDG202_cl2.csv"

try:
    # Descargar el archivo desde GitHub
    response = requests.get(github_url)
    response.raise_for_status()  # Verifica errores de descarga
    
    # Leer el contenido como CSV
    df = pd.read_csv(io.StringIO(response.text))
    
    # Verificar si la columna SMILES está presente
    if 'SMILES' not in df.columns:
        st.error("El archivo debe contener una columna llamada 'SMILES'.")
    else:
        st.success("Archivo cargado correctamente.")
        
        # Mostrar la tabla de datos
        st.subheader("Tabla de Moléculas")
        st.dataframe(df)

        # Generar y mostrar imágenes de las moléculas
        st.subheader("Visualización de Moléculas en 2D")
        for index, row in df.iterrows():
            smiles = row['SMILES']
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    img = Draw.MolToImage(mol)
                    st.image(img, caption=f"Molécula {index + 1}: {smiles}")
                else:
                    st.warning(f"No se pudo procesar el SMILES: {smiles}")
            except Exception as e:
                st.error(f"Error al generar la imagen para el SMILES '{smiles}': {e}")
except Exception as e:
    st.error(f"Error al descargar o procesar el archivo: {e}")


