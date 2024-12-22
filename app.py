import streamlit as st
import pandas as pd
import requests

# Título de la aplicación
st.title("Visualización de Moléculas (Sin RDKit)")

# URL del archivo
github_url = "https://raw.githubusercontent.com/mdgenriquez/punonpdb/main/MDG202_cl2.csv"

try:
    response = requests.get(github_url)
    response.raise_for_status()
    df = pd.read_csv(io.StringIO(response.text))

    if 'SMILES' not in df.columns:
        st.error("El archivo debe contener una columna llamada 'SMILES'.")
    else:
        st.success("Archivo cargado correctamente.")
        st.dataframe(df)

except Exception as e:
    st.error(f"Error al procesar el archivo: {e}")
