import streamlit as st
import pandas as pd
import requests

# Título de la aplicación
st.title("Visualización de Moléculas en 2D con MolView")

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

        # Visualización de las moléculas en 2D con MolView
        st.subheader("Visualización de Moléculas en 2D")
        for index, row in df.iterrows():
            smiles = row['SMILES']
            molview_url = f"https://molview.org/?lang=es#mode=ballAndStick&url=smiles:{smiles}"
            st.markdown(f"### Molécula {index + 1}: {smiles}")
            st.markdown(f"Visualización en MolView: [Abrir MolView]( {molview_url} )")

except Exception as e:
    st.error(f"Error al descargar o procesar el archivo: {e}")
