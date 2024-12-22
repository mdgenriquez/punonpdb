import streamlit as st
import pandas as pd
import requests
from rdkit import Chem
#from rdkit.Chem import Draw

# Título de la aplicación
st.title("Visualización de Moléculas en 2D")
st.markdown("Visualiza moléculas desde el archivo `MDG202_cl2.csv` en tu repositorio de GitHub.")

# URL del archivo en GitHub (MDG202_cl2.csv)
github_url = "MDG202_cl2.csv"

# Descargar el archivo
try:
    response = requests.get(github_url)
    response.raise_for_status()  # Verifica que no haya errores en la descarga
    df = pd.read_csv(pd.compat.StringIO(response.text))
    
    if 'SMILES' not in df.columns:
        st.error("El archivo debe contener una columna llamada 'SMILES'.")
    else:
        st.success("Archivo cargado correctamente.")
        
        # Mostrar tabla con las moléculas
        st.subheader("Tabla de Moléculas")
        st.dataframe(df)
        
        # Generar imágenes de moléculas
        st.subheader("Visualización de Moléculas en 2D")
        
        for index, row in df.iterrows():
            smiles = row['SMILES']
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    img = Draw.MolToImage(mol, size=(300, 300))
                    st.image(img, caption=f"Molécula {index + 1}: {smiles}", use_column_width=False)
                else:
                    st.warning(f"No se pudo generar la estructura para el SMILES: {smiles}")
            except Exception as e:
                st.error(f"Error procesando la molécula en la fila {index + 1}: {e}")
except Exception as e:
    st.error(f"Error al descargar o procesar el archivo desde GitHub: {e}")
